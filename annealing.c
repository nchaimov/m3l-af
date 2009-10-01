#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"
#include "rates.h"
#include "annealing.h"
#include "numeric.h"
#include "spr.h"
#include "emig.h"

/**
* JSJ: Note that all external libraries, including math.h, time.h and gsl are included
* 		through utilities.h
*/

annealing anneal;

// "Set_Anneal" sets the default values of thermal annealing parameters. . .
void Set_Anneal(option *io){
	anneal.accept_ratio = 0.3; //this one gets set from io later... don't change this line.
	anneal.end_temp = io->temp_end;
	anneal.iters_per_temp = io->iters_per_stage;
	anneal.set_back = io->set_back;
	anneal.start_temp = io->temp_start;
	anneal.start_tau = io->tau_start;
	anneal.end_tau = io->tau_end;
	anneal.num_anneal_stages = io->num_anneal_stages;
	//printf("(annealing.c) anneal.num_anneal_stages = %d, %d", io->num_anneal_stages, anneal.num_anneal_stages);

	anneal.max_alpha = io->max_alpha;


	anneal.brlen_sigma = io->brlen_sigma;
	anneal.pinvar_sigma = io->pinvar_sigma;
	anneal.gamma_sigma = io->gamma_sigma;
	anneal.emig_sigma = io->emig_sigma;
	anneal.prob_emig = io->prob_emig;

	// prob_X = the probability of stepping parameter X during each SA iteration.
	anneal.prob_NNI = io->prob_NNI;
	anneal.prob_SPR = io->prob_SPR;
	//anneal.prob_TBR = 0.3;
	anneal.prob_brlen = io->prob_brlen;
	anneal.prob_gamma = io->prob_gamma;
	anneal.prob_kappa = io->prob_kappa;
	anneal.prob_lambda = io->prob_lambda;
	anneal.prob_rr = io->prob_rr;
	anneal.prob_pi = io->prob_pi;
	anneal.prob_rate_proportion = io->prob_rate_proportion;
	anneal.prob_topology = io->prob_topology;
	//	anneal.prob_trans_model = 0.1;
	anneal.prob_pinvar = io->prob_pinvar;

	anneal.random_seed = (int)time(NULL);
	anneal.rng = gsl_rng_alloc(gsl_rng_mt19937);
	anneal.no_change = 1;
	gsl_rng_set(anneal.rng, anneal.random_seed);
}

void Free_Anneal(){
	gsl_rng_free( anneal.rng );
}


void Swap_Tree_And_Mod(arbre *tree1, arbre *tree2){
	arbre *tmptree;
	model *tmpmod;
	tmptree = tree1;
	tmpmod = tree1->mod;

	tree1 = tree2;
	tree1->mod = tree2->mod;
	tree2 = tmptree;
	tree2->mod = tree1->mod;
}

// This method returns the scaled and estimated acceptance ratio for thermal annealing.
// The variable "tree" will be modified by this method, such that tree will be in a good
// starting position for thermal annealing after this method concludes.
m3ldbl Scale_Acceptance_Ratio(arbre *tree){
	PhyML_Printf(". Scaling acceptance ratio for simulated annealing...\n");

	/* psuedo-code:

	lnL_best = tree's lnL
	lnL_proposed = max_lnL

	fsum = 0.0 // function mean
	fsqsum  0.0 // function standard deviation
	n = 0 // number of points tried

	tree_proposed = copy of tree

	for iter in range(0, iters_per_temp):
		Get_Neighbor_Proposition(tree_proposed)
		lnL_proposed = tree_proposed's likelihood
		if (lnL_proposed > -DBL_MAX):
			n++
			fsum += working_lnL
			fsqsum += working_lnL * working_lnL
		if (lnL_proposed > lnL_best):
			lnL_best = lnL_proposed
			tree = copy of tree_proposed

	fsum = fsum / n
	fsqsum = fsqsum / n
	fsqsum = sqrt( fsqsum - fsum * fsum)
	aratio = accept_ratio * fsqsum / start_temp
	return aratio
	 */

	m3ldbl aratio;
	Lk(tree);
	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_proposed = tree->c_lnL;
	m3ldbl fsum = 0.0;
	m3ldbl fsqsum = 0.0;
	int n = 0;
	int i;
	//Optimiz_All_Free_Param(tree,0);
	arbre *best_tree = Make_Tree(tree->n_otu);
	Init_Tree(best_tree,tree->n_otu);
	Make_All_Tree_Nodes(best_tree, tree->mod->n_l);
	Make_All_Tree_Edges(best_tree, tree->mod->n_l);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree,best_tree);

	For(i,anneal.iters_per_temp){
		Get_Neighbor_Proposition(tree,anneal.start_temp);
		lnL_proposed = tree->c_lnL;
		if(lnL_proposed > UNLIKELY){
			n++;
			fsum += lnL_proposed;
			fsqsum += (lnL_proposed * lnL_proposed);
		}
		if(lnL_proposed > lnL_best){
			lnL_best = lnL_proposed;
			Copy_Tree(tree,best_tree);
			Record_Model(tree->mod,best_tree->mod);
		}else{
			//JSJ: restore tree and model
			Copy_Tree(best_tree,tree);
			Record_Model(best_tree->mod,tree->mod);
		}

	}
	fsum /= (double)n; //average proposed likelihood
	fsqsum /= (double)n; //average (squared) proposed likelihood
	fsqsum = sqrt(fsqsum - fsum * fsum);
	aratio = anneal.accept_ratio * fsqsum / anneal.start_temp;

	Free_Model(best_tree->mod);
	Free_Tree(best_tree);

	//now lets get a really good starting position...
	//Speed_Spr_Loop(tree);

	return aratio;
}

// This method modifies some or all parameters of "tree", including:
// * topology
// * branch lengths
// * branch length proportions
// * alpha (for gamma distributed ASRV)
// * gamma proportions
// * mu (evolutionary rate)
void Get_Neighbor_Proposition(arbre *tree,m3ldbl temp){
	// For now, we'll always perturb every parameter.
	// In the future, we'll do something more sophisticated, where the probability of perturbing
	// any particular parameter will be drawn from a probability distribution.
	tree->both_sides = 0; //set 1 if topology, brlen, or props
	anneal.no_change = 1;
	double x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_rate_proportion) Step_Brlen_Proportion(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_gamma)Step_Gamma(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_pinvar)Step_Pinvar(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_kappa)Step_Kappa(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_lambda)Step_Lambda(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_rr)Step_RR(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_pi)Step_Pi(tree);
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_topology){
		//PhyML_Printf("JSJ: Proposing a new topology\n");
		Step_Topology(tree);
	}
	x = gsl_rng_uniform(anneal.rng);
	if(x < anneal.prob_brlen)Step_Branch_Lengths(tree,temp);

	//	if(x == 1)Step_Pi(tree);
	//	if(x == 1)Step_Lambda(tree);
	//	if(x == 1)Step_Kappa(tree);
	//	if(x == 1)Step_RR(tree);
	//	x = gsl_rng_uniform(anneal.rng);
	//	if(x < anneal.prob_pinvar)Step_Pinvar(tree);

	// 2. set update_eigen to 0 (and in some cases, update_eigen will already == 0)


	// 1. Update the likelihood of tree
	if(anneal.no_change == 0)
	{
		Lk(tree);
	}
	tree->mod->update_eigen = 0;
	tree->both_sides = 1; // reset to 1
}

//
// Returns the absolute geometric distance between last_tree and tree, where
// distance accounts for changes in parameters: alpha, branch length proportions,
// prob. invar, kappa, lamda, RR, pi, topology, and branch lengths.
//
// This method is useful for traversing multi-parametric space.
//
m3ldbl Get_distance_between_trees(arbre *last_tree, arbre *tree){
	return 1.0;

	m3ldbl distance = 0.0; // distance is geometric.

	distance += (last_tree->mod->alpha - tree->mod->alpha)*(last_tree->mod->alpha - tree->mod->alpha);
	int i;
	for(i = 0; i < tree->mod->n_l; i++){
		distance += (last_tree->mod->bl_props[i] - tree->mod->bl_props[i])*(last_tree->mod->bl_props[i] - tree->mod->bl_props[i]);
	}
	distance += (last_tree->mod->pinvar - tree->mod->pinvar)*(last_tree->mod->pinvar - tree->mod->pinvar);
	distance += (last_tree->mod->kappa - tree->mod->kappa)*(last_tree->mod->kappa - tree->mod->kappa);
	distance += (last_tree->mod->lambda - tree->mod->lambda)*(last_tree->mod->lambda - tree->mod->lambda);
	for(i = 0; i < tree->mod->n_diff_rr; i++){
		distance += (last_tree->mod->rr_val[i] - tree->mod->rr_val[i])*(last_tree->mod->rr_val[i] - tree->mod->rr_val[i]);
	}
	for(i = 0; i < tree->mod->ns; i++){
		distance += (last_tree->mod->pi_unscaled[i] - tree->mod->pi_unscaled[i])*(last_tree->mod->pi_unscaled[i] - tree->mod->pi_unscaled[i]);
	}
	// What about topology distance? and branch length distance?
	m3ldbl d = Get_distance_between_topologies(last_tree->noeud[0], tree->noeud[0]);
	distance += d*d;

	return distance;
}

//
// This is a helper method for Get_distance_between_trees
// nodea and nodeb
//
//
m3ldbl Get_distance_between_topologies(node *nodea, node *nodeb){
	return 0.0;
}

// helper for "Get_Neighbor_Proposition"
void Step_Brlen_Proportion(arbre *tree){
	// For now, do something simple and stupid, like stepping branch lengths by +/- 0.001.
	// We'll eventually do something more sophisticated, like selected the new branch
	// lengths from a Dirichlet distribution.
	// when we update proportions we need to recalculate the likelihoods on the whole tree...
	if(tree->mod->s_opt->opt_props == 1){
		tree->both_sides = 1;
		anneal.no_change = 0;
		//		int i,j;
		//		double r = (((double)rand() + 1.0) / ((double)(RAND_MAX)+ 1.0));
		//		int prange = (int)(Rand_Int(1,(tree->mod->n_l)) * r);
		//		For(i,prange){
		//			j = Rand_Int(0,(tree->mod->n_l - 1));
		//			r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0; //r is in the range -0.04999 to 0.049999
		//			tree->props[j] += r;
		//			if(tree->props[j] > 1.0) tree->props[j] = 1.0;
		//			else if(tree->props[j] < 0.0) tree->props[j] = 0.0;
		//		}
		double *alpha;
		alpha = (m3ldbl *)mCalloc(tree->mod->n_l, sizeof(m3ldbl));
		int i;
		For(i,tree->mod->n_l){
			alpha[i] = anneal.max_alpha;
			PhyML_Printf("(annealing.c 201) alpha[%d] = %f", i, alpha[i]);
		}
		gsl_ran_dirichlet(anneal.rng,tree->mod->n_l,alpha,tree->mod->bl_props);
		Normalize_Props(tree->mod);

		Free(alpha);
	}
}

void Step_Pinvar(arbre *tree){
	if(tree->mod->s_opt->opt_pinvar == 1){
		double r;
		anneal.no_change = 0;
		r = gsl_ran_gaussian(anneal.rng,anneal.pinvar_sigma);
		tree->mod->pinvar += r;
		if(tree->mod->pinvar < 0.0001) tree->mod->pinvar = 0.0001;
		else if(tree->mod->pinvar > 0.9999) tree->mod->pinvar = 0.9999;
	}
}

void Step_RR(arbre * tree){
	if(((tree->mod->whichmodel == GTR) && tree->mod->n_diff_rr > 1) ||
			((tree->mod->whichmodel == CUSTOM) && (tree->mod->s_opt->opt_rr == 1) && (tree->mod->n_diff_rr > 1)))
	{
		int i;
		m3ldbl tmp = 0;
		anneal.no_change = 0;
		double *alpha = calloc(tree->mod->n_diff_rr,sizeof(double));
		if(tree->mod->n_diff_rr > 5){
			tmp = tree->mod->rr_val[5]; //save this number to restore later
		}
		for(i = 0; i < tree->mod->n_diff_rr; i++){
			alpha[i] = anneal.max_alpha;
		}
		gsl_ran_dirichlet(anneal.rng,tree->mod->n_diff_rr,alpha,tree->mod->rr_val);
		if(tree->mod->n_diff_rr > 5){
			tree->mod->rr_val[5] = tmp; //restore this number
			//JSJ: not exactly sure why but this is what is done in Optimize_All_Free_Params
		}
		free(alpha);
//		int i,j;
//		int range = Rand_Int(1,(tree->mod->n_diff_rr));
//		For(i,range){
//			j = Rand_Int(0,(tree->mod->n_diff_rr - 1));
//			if(j == 5) j = (j+1)%(tree->mod->n_diff_rr);
//			r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0;
//			tree->mod->rr_val[j] += r;
//			if(tree->mod->rr_val[j] < 0.01) tree->mod->rr_val[j] = 0.01;
//			else if(tree->mod->rr_val[j] > 100.0) tree->mod->rr_val[j] = 99.999;
//		}
	}
}

//.1 to 100.0
void Step_Kappa(arbre * tree){
	if(tree->mod->s_opt->opt_kappa == 1){
		double r;
		tree->mod->update_eigen = 1;
		anneal.no_change = 0;
		r = gsl_ran_gaussian(anneal.rng,tree->mod->kappa);
		tree->mod->kappa += r;
		if(tree->mod->kappa < 0.1) tree->mod->kappa = 0.1;
		else if(tree->mod->kappa > 100.0) tree->mod->kappa = 99.9999;
//		double r;
//		tree->mod->update_eigen = 1;
//		r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0;
//		tree->mod->kappa += r;
//		if(tree->mod->kappa < 0.1) tree->mod->kappa = 0.1;
//		else if(tree->mod->kappa > 100.0) tree->mod->kappa = 99.9999;
	}
}

//.001 to 100
void Step_Lambda(arbre * tree){
	if(tree->mod->s_opt->opt_lambda == 1){
		double r;
		anneal.no_change = 0;
		r = gsl_ran_gaussian(anneal.rng,tree->mod->lambda);
		tree->mod->lambda += r;
		if(tree->mod->lambda < 0.001) tree->mod->lambda = 0.001;
		else if(tree->mod->lambda > 100.0) tree->mod->lambda = 99.9999;
	}
}

// helper for "Get_Neighbor_Proposition"
void Step_Gamma(arbre *tree){
	// For now, do something stupid like incrementing gamma +/- 0.001.
	// In reality, we'll want to perturb the gamma parameter based on a Gaussian distribution
	// with mean equal to the current value of gamma and standard deviation equal to some
	// user-specified parameter sigma.
	if(tree->mod->s_opt->opt_alpha == 1){
		anneal.no_change = 0;
		double r = gsl_ran_gaussian(anneal.rng,anneal.gamma_sigma);
		tree->mod->alpha += r;
		if(tree->mod->alpha < 0.01) tree->mod->alpha = 0.01;
		else if(tree->mod->alpha > 100.0) tree->mod->alpha = 100.0;
	}
}

void Step_Pi(arbre * tree){
	if((tree->mod->s_opt->opt_state_freq) && (tree->mod->datatype == NT)){
		tree->mod->update_eigen = 1;
		anneal.no_change = 0;
		int i;
		double *alpha = calloc(tree->mod->ns, sizeof(double));
		for(i = 0; i < tree->mod->ns; i++){
			alpha[i] = anneal.max_alpha;
		}
		gsl_ran_dirichlet(anneal.rng,tree->mod->ns,alpha,tree->mod->pi_unscaled);
		for(i = 0; i < tree->mod->ns; i++){
			if(tree->mod->pi_unscaled[i] < -1000.0) tree->mod->pi_unscaled[i] = -999.9999;
			else if(tree->mod->pi_unscaled[i] > 1000.0) tree->mod->pi_unscaled[i] = 999.9999;
		}
	}
}

// helper for "Get_Neighbor_Proposition"
void Step_Branch_Lengths(arbre *tree, m3ldbl temp){
	if(tree->mod->s_opt->opt_bl == 1){
		tree->both_sides = 1;
		anneal.no_change = 0;
		int i,j,edge_range,set_range,m,n,max;
		int n_edges = (tree->n_otu * 2) - 3;
		double rand_gauss;
		/* r is a random floating point value in the range (0,1) {not including 0,
		 * or 1}. Note we must convert rand() and/or RAND_MAX+1 to
		 * floating point values to avoid integer division. In addition, Sean
		 * Scanlon pointed out the possibility that RAND_MAX may be the largest
		 * positive integer the architecture can represent, so (RAND_MAX+1)
		 * may result in an overflow, or more likely the value will end up being
		 * the largest negative integer the architecture can represent, so
		 * to avoid this we convert RAND_MAX and 1 to doubles before adding. */
		max = (n_edges + 1) * temp;
		if (max < 1) max = 1;
		edge_range = gsl_rng_uniform_int(anneal.rng,max);
		if(edge_range < 1) edge_range = 1;

		max = (tree->mod->n_l + 1) * temp;
		if(max < 1) max = 1;
		set_range = gsl_rng_uniform_int(anneal.rng,max);
		if(set_range < 1) set_range = 1;
		For(i,edge_range){
			j = gsl_rng_uniform_int(anneal.rng,n_edges);
			For(m,set_range){
				rand_gauss = gsl_ran_gaussian(anneal.rng,anneal.brlen_sigma);
				n = gsl_rng_uniform_int(anneal.rng,tree->mod->n_l);
				tree->t_edges[j]->l[n] += rand_gauss;
				if(tree->t_edges[j]->l[n] < BL_MIN) tree->t_edges[j]->l[n] = BL_MIN;
				else if(tree->t_edges[j]->l[n] > BL_MAX) tree->t_edges[j]->l[n] = BL_MAX;


				//Update_Lk_At_Given_Edge(tree->t_edges[j],tree); //calls update_p_lk on appropriate nodes and this edge

				//JSJ: just for fun...
				//				Br_Len_Brent_Iter(10.*tree->t_edges[j]->l[n],tree->t_edges[j]->l[n],BL_MIN,
				//									0.00000001,
				//									tree->t_edges[j],tree,
				//									10,
				//									0,n);
			}
//			if(anneal.full_lk == 0){
//				Update_PMat_At_Given_Edge(tree->t_edges[j],tree);
//				if(!tree->t_edges[j]->left->tax) Update_P_Lk(tree,tree->t_edges[j],tree->t_edges[j]->left);
//				if(!tree->t_edges[j]->rght->tax) Update_P_Lk(tree,tree->t_edges[j],tree->t_edges[j]->rght);
//				Update_Lk_At_Given_Edge(tree->t_edges[j],tree);
//			}
		}


	}
}
//
//void Random_NNI(arbre *tree){
//
//}


// helper for "Get_Neighbor_Proposition"
void Step_Topology(arbre *tree){
	// invoke NNI or SPR to create a new topology.
	// JSJ: pick a random edge, from that edge find a path connecting 4 nodes
	//       swap those four nodes.
	/**
	* Find an edge where both nodes are internal (n_edges > 3 or opt_topo should be set false)
	* When both nodes are internal
	*
	*/
	if(tree->mod->s_opt->opt_topo){
		tree->mod->update_eigen = 1;
		tree->both_sides = 1;
		anneal.no_change = 0;
		double p = gsl_rng_uniform(anneal.rng);
		if(p <= anneal.prob_NNI){
			node *b,*c;
			node *a = NULL;
			node *d = NULL;
			int i,j;
			int n_edges = (tree->n_otu * 2) - 3;
			int edge = gsl_rng_uniform_int(anneal.rng,n_edges);
			//PhyML_Printf("JSJ: Edge: %i, Number of Edges: %i\n",edge,n_edges);
			double r = gsl_rng_uniform(anneal.rng);
			while(tree->t_edges[edge]->left->tax == 1 || tree->t_edges[edge]->rght->tax == 1){
				if(r < 0.5) edge = ((edge + 1) % (n_edges)); //starting from here, look through edge list
				else{
					if(edge == 0) edge = n_edges - 1;
					else edge = ((edge - 1) % (n_edges));
				}
			}// now we have an internal edge...
			b = tree->t_edges[edge]->left;
			c = tree->t_edges[edge]->rght;

			/**
			* now we need to find nodes a (from b) and d (from c) such that
			* a != c, d != a, and either a,d == terminal or a,d == internal
			*/
			For(i,3){
				if(b->v[i] != c){
					a = b->v[i];
					For(j,3){
						if(c->v[j] != b){
							d = c->v[j];
							if(a->tax == d->tax){
								break;
							}
						}
					}
					if(a->tax == d->tax){
						break;
					}
				}
			}

			Swap(a,b,c,d,tree);
//			if(anneal.full_lk == 0){
//				tree->both_sides = 1;
//				Update_PMat_At_Given_Edge(tree->t_edges[edge],tree);
//				if(!tree->t_edges[edge]->left->tax){
//					For(i,3)
//					if(tree->t_edges[edge]->left->v[i] != tree->t_edges[edge]->rght)
//						Update_P_Lk(tree,tree->t_edges[edge]->left->b[i],tree->t_edges[edge]->left);
//				}
//				if(!tree->t_edges[edge]->rght->tax){
//					For(i,3)
//					if(tree->t_edges[edge]->rght->v[i] != tree->t_edges[edge]->left)
//						Update_P_Lk(tree,tree->t_edges[edge]->rght->b[i],tree->t_edges[edge]->rght);
//				}
//				Update_Lk_At_Given_Edge(tree->t_edges[edge],tree);
//			}
		}else if((p - anneal.prob_NNI) <= anneal.prob_SPR){ //otherwise do SPR
			//int num_spr = gsl_rng_uniform_int(anneal.rng,3);
			Random_Spr(1,tree);
		}else{
			Migrate_One_Edge(tree,&(anneal));
		}
	}
}


void Get_QA_Neighbor_Proposition(arbre *tree){

}

// Boltzmann probability sampler, for thermal annealing
//
// INPUT: Given the likelihood (lnl_curr) of our current state, the likelihood (lnl_new) of our new proposed state,
// and the temperature (temperature) of the system, calculate the probability of the new state,
// where the probability is drawn from a Boltzmann distribution.
//
// OUTPUT: a floating point decimal, between 0.0 and 1.0
m3ldbl Boltzmann_P_TA(m3ldbl lnl_curr, m3ldbl lnl_new, m3ldbl temperature){
	m3ldbl result = 1.0;
	if(lnl_new < lnl_curr){
		result = exp((lnl_new - lnl_curr)/(anneal.accept_ratio * temperature));
	}
	//PhyML_Printf("In Boltzmann_P: accept_ratio= %lf, lnl_curr = %lf, lnl_new = %lf, temp= %lf, result = %lf\n",anneal.accept_ratio,lnl_curr,lnl_new,temperature,result);
	return result;
}

//Boltzmann probability sampler, for quantum annealing
//
// distance = the distance between the proposed state and the current state
// tfield = the transverse field value (tau), which is QA's analog to the temperature in TA.
//
m3ldbl Boltzmann_P_QA(m3ldbl lnl_curr, m3ldbl lnl_new, m3ldbl distance, m3ldbl tfield){
	m3ldbl result = 1.0;
	if(lnl_new < lnl_curr){
		result = exp(distance / (anneal.accept_ratio * tfield));
	}
	return result;
}

// INPUT: a tree structure, where tree->mod->s_opt contains the
//      parameters which are free to be optimized.
// OUTPUT: the likelihood of the best found tree
//
// At the end of this method, the object "tree" will contain the best-found topology,
//      branch lengths, and model parameters.  In other words, tree is mutable!
m3ldbl Thermal_Anneal_All_Free_Params(arbre *tree, int verbose){

	Set_Anneal(tree->io);
	//m3ldbl result = 1.0;
	tree->mod->update_eigen = 1;
	tree->both_sides = 1; //search both pre and post order on all subtrees
	int n_edges = (tree->n_otu * 2) - 3;
	if (n_edges <= 3){
		tree->mod->s_opt->opt_topo = 0; //make sure that opt_topo is false if there are no meaningful branch swaps.
	}
	m3ldbl temp = anneal.start_temp;
	m3ldbl tempmult = exp(log(anneal.end_temp/anneal.start_temp)/(((double)anneal.num_anneal_stages) - 1.0));
	if(tree->io->acc_ratio < 0.0) anneal.accept_ratio = Scale_Acceptance_Ratio(tree);
	else{
		anneal.accept_ratio = tree->io->acc_ratio; //if positive, user input a value, don't estimate.
		Lk(tree);
#ifdef MEASURE
		Count_Mean_Compressability(tree);
#endif
	}

	// VHS: do we need this? I think it gets us to a good starting location.
	//Optimiz_All_Free_Param(tree,1);

	arbre *best_tree = Make_Tree(tree->n_otu);
	Init_Tree(best_tree,tree->n_otu);
	Make_All_Tree_Nodes(best_tree, tree->mod->n_l);
	Make_All_Tree_Edges(best_tree, tree->mod->n_l);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree,best_tree);

	//Speed_Spr_Loop(tree);
	arbre *last_tree = Make_Tree(tree->n_otu);
	Init_Tree(last_tree,tree->n_otu);
	Make_All_Tree_Nodes(last_tree, tree->mod->n_l);
	Make_All_Tree_Edges(last_tree, tree->mod->n_l);
	last_tree->mod = Copy_Model(best_tree->mod);

	time_t start = time(NULL);
	time_t now;

	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_current = tree->c_lnL;
	m3ldbl lnL_proposed = tree->c_lnL;
	m3ldbl acc_prob;
	double r;
	int itemp,iter;
	int steps_tried = 0;
	int steps_accepted = 0;

	m3ldbl temp_of_best = 0; // VHS: for optimization purposes, let's record the temperature at which we discovered the best tree.
	m3ldbl iter_of_best = 0;

	PhyML_Printf("\n\n Starting simulated thermal annealing with the following parameters:\n");
	PhyML_Printf("\t. acceptance ratio = %f\n", anneal.accept_ratio);
	PhyML_Printf("\t. start temperature = %f\n", anneal.start_temp);
	PhyML_Printf("\t. end temperature = %f\n", anneal.end_temp);
	PhyML_Printf("\t. temperature count = %d\n", anneal.num_anneal_stages);
	PhyML_Printf("\t. temperature multiplier = %f\n", tempmult);

	PhyML_Printf("\t. P of adjusting alpha = %f\n", anneal.prob_gamma );
	PhyML_Printf("\t. maximum alpha  = %f\n", anneal.max_alpha);
	PhyML_Printf("\t. a.s.r.v. gamma sigma = %f\n", anneal.gamma_sigma);

	PhyML_Printf("\t. P of adjusting branch lengths = %f\n", anneal.prob_brlen );
	PhyML_Printf("\t. branch length sigma = %f\n", anneal.brlen_sigma);

	PhyML_Printf("\t. P of adjusting pinvar = %f\n", anneal.prob_pinvar );
	PhyML_Printf("\t. pinvar. sigma = %f\n", anneal.pinvar_sigma);

	PhyML_Printf("\t. P of changing topology = %f\n", anneal.prob_topology );
	PhyML_Printf("\t. P of using TBR = %f\n", anneal.prob_TBR );
	PhyML_Printf("\t. P of using NNI = %f\n", anneal.prob_NNI );
	PhyML_Printf("\t. P of using SPR = %f\n", anneal.prob_SPR );

	PhyML_Printf("\t. P of adjusting branch length rate proportions = %f\n", anneal.prob_rate_proportion );


	for(itemp = 0; itemp < anneal.num_anneal_stages; itemp++){
		//recenter the search at each temperature.
		//make the tree our best tree so far
		//
		// VHS: This re-centering step might be meretricious.
		// Let's examine the run without this step.
		Copy_Tree(best_tree,tree);
		Record_Model(best_tree->mod,tree->mod);
		Copy_Tree(tree,last_tree);
		Record_Model(tree->mod,last_tree->mod);
		lnL_current = lnL_best;

		//steps_tried = 0;
		//steps_accepted = 0;

		//Optimiz_All_Free_Param(tree,1);

		for(iter = 0; iter < anneal.iters_per_temp; iter++){
			steps_tried++;
			//get our next proposition.
			Get_Neighbor_Proposition(tree,temp);
			lnL_proposed = tree->c_lnL;
#ifdef MEASURE
			Count_Mean_Compressability(tree);
#endif

			now = time(NULL);
			// Some useful debugging statements:
			//Print_Tree_Screen(tree);
			PhyML_Printf("T: %f iter: %d current lnL = %f proposed lnL = %f\n", temp, iter, lnL_current, lnL_proposed);
			//PhyML_Printf("plot1 %d %lf\n", (steps_tried * (itemp + 1)), lnL_current);
			//PhyML_Printf("plot2 %ld %d\n", (long)now, (steps_tried * (itemp + 1)));
			PhyML_Printf("plot1 %d %lf\n", steps_tried, lnL_current);
			PhyML_Printf("plot2 %ld %d\n", (long)difftime(now,start), steps_tried);
			//PhyML_Printf("proposed tree = %s\n", Write_Tree(tree) );

			if(lnL_proposed > lnL_best){
				//save this tree into best_tree
				//PhyML_Printf("(annealing.c, 465): Found a new best tree with lnL = %f\n", lnL_best);
				Copy_Tree(tree,best_tree);
				Record_Model(tree->mod, best_tree->mod);
				lnL_best = lnL_proposed;
				temp_of_best = temp;
				iter_of_best = iter;

				// This temperature is yielding good results: let's stay here
				// longer, thus increasing our chances of reaching the best ground
				// state at this temperature.
				iter -= anneal.set_back;
				if(iter < 0) iter = 0; //make sure not set back into negative...
				acc_prob = 1.0;
				r = 0.0;
			}
			else if(lnL_proposed > lnL_current){
				acc_prob = 1.0;
				r = 0.0;
			}else{
				acc_prob = Boltzmann_P_TA(lnL_current, lnL_proposed, temp);
				r = gsl_rng_uniform(anneal.rng);
			}

			if(acc_prob >= r){
				//save the current tree
				if(anneal.no_change == 0){
					Copy_Tree(tree,last_tree);
					Record_Model(tree->mod,last_tree->mod);
				}
				steps_accepted++;
				lnL_current = tree->c_lnL;

				//PhyML_Printf("JSJ: Proposed Likelihood Accepted, acc_prob = %lf\n",acc_prob);

			}else{
				//restore the current tree to be our last tree
				if(anneal.no_change == 0){
					Copy_Tree(last_tree,tree);
					Record_Model(last_tree->mod,tree->mod);
				}
			}
		}//end inner for loop.
		temp *= tempmult;

	}
	if(verbose){
		PhyML_Printf("Annealing finished.\n");
		PhyML_Printf("In temperature range [%f, %f], the best tree was found at %f, iteration %f\n", temp, anneal.start_temp, temp_of_best, iter_of_best);
	}
	Copy_Tree(best_tree,tree);
	Record_Model(best_tree->mod,tree->mod);

	Free_Model(best_tree->mod);
	Free_Tree(best_tree);
	Free_Model(last_tree->mod);
	Free_Tree(last_tree);
	Free_Anneal(); //local function that frees the gsl_rng;
	Lk(tree);

	return tree->c_lnL;
}


// INPUT: a tree structure, where tree->mod->s_opt contains the
//      parameters which are free to be optimized.
// OUTPUT: the likelihood of the best found tree
//
// At the end of this method, the object "tree" will contain the best-found topology,
//      branch lengths, and model parameters.  In other words, tree is mutable!
m3ldbl Quantum_Anneal_All_Free_Params(arbre *tree, int verbose){

	Set_Anneal(tree->io);
	//m3ldbl result = 1.0;
	tree->mod->update_eigen = 1;
	tree->both_sides = 1; //search both pre and post order on all subtrees
	int n_edges = (tree->n_otu * 2) - 3;
	if (n_edges <= 3){
		tree->mod->s_opt->opt_topo = 0; //make sure that opt_topo is false if there are no meaningful branch swaps.
	}
	m3ldbl tau = anneal.start_tau;
	m3ldbl taumult = exp(log(anneal.end_tau/anneal.start_tau)/(((double)anneal.num_anneal_stages) - 1.0));
	if(tree->io->acc_ratio < 0.0) anneal.accept_ratio = Scale_Acceptance_Ratio(tree);
	else{
		anneal.accept_ratio = tree->io->acc_ratio; //if positive, user input a value, don't estimate.
		Lk(tree);
	}

	Optimiz_All_Free_Param(tree,1);
	//printf("(annealing.c 742) anneal.num_anneal_stages = %d", anneal.num_anneal_stages);

	arbre *best_tree = Make_Tree(tree->n_otu);
	Init_Tree(best_tree,tree->n_otu);
	Make_All_Tree_Nodes(best_tree, tree->mod->n_l);
	Make_All_Tree_Edges(best_tree, tree->mod->n_l);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree,best_tree);

	//Speed_Spr_Loop(tree);
	arbre *last_tree = Make_Tree(tree->n_otu);
	Init_Tree(last_tree,tree->n_otu);
	Make_All_Tree_Nodes(last_tree, tree->mod->n_l);
	Make_All_Tree_Edges(last_tree, tree->mod->n_l);
	last_tree->mod = Copy_Model(best_tree->mod);
	time_t start = time(NULL);
	time_t now;

	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_current = tree->c_lnL;
	m3ldbl lnL_proposed = tree->c_lnL;
	m3ldbl acc_prob;
	double r;
	int itau,iter;
	int steps_tried = 0;
	int steps_accepted = 0;

	m3ldbl tau_of_best = 0; // VHS: for optimization purposes, let's record the tau at which we discovered the best tree.
	m3ldbl iter_of_best = 0;

	PhyML_Printf("\n\n Starting simulated quantum annealing with the following parameters:\n");
	PhyML_Printf("\t. acceptance ratio = %f\n", anneal.accept_ratio);
	PhyML_Printf("\t. starting tau value = %f\n", anneal.start_tau);
	PhyML_Printf("\t. ending tau value = %f\n", anneal.end_tau);
	PhyML_Printf("\t. tau count = %d\n", anneal.num_anneal_stages);
	PhyML_Printf("\t. tau multiplier = %f\n", taumult);

	PhyML_Printf("\t. P of adjusting alpha = %f\n", anneal.prob_gamma );
	PhyML_Printf("\t. maximum alpha  = %f\n", anneal.max_alpha);
	PhyML_Printf("\t. a.s.r.v. gamma sigma = %f\n", anneal.gamma_sigma);

	PhyML_Printf("\t. P of adjusting branch lengths = %f\n", anneal.prob_brlen );
	PhyML_Printf("\t. branch length sigma = %f\n", anneal.brlen_sigma);

	PhyML_Printf("\t. P of adjusting pinvar = %f\n", anneal.prob_pinvar );
	PhyML_Printf("\t. pinvar. sigma = %f\n", anneal.pinvar_sigma);

	PhyML_Printf("\t. P of changing topology = %f\n", anneal.prob_topology );
	PhyML_Printf("\t. P of using TBR = %f\n", anneal.prob_TBR );
	PhyML_Printf("\t. P of using NNI = %f\n", anneal.prob_NNI );
	PhyML_Printf("\t. P of using SPR = %f\n", anneal.prob_SPR );

	PhyML_Printf("\t. P of adjusting branch length rate proportions = %f\n", anneal.prob_rate_proportion );

	//m3ldbl  prob_pi;
	//m3ldbl  prob_kappa;
	//m3ldbl  prob_lambda;
	//m3ldbl  prob_rr;
	//m3ldbl  prob_emig;


	for(itau = 0; itau < anneal.num_anneal_stages; itau++){
		//recenter the search at each tau.
		//make the tree our best tree so far
		//
		// VHS: This recentering step might be meretricious.
		// Let's examine the run without this step.
		Copy_Tree(best_tree,tree);
		Record_Model(best_tree->mod,tree->mod);
		Copy_Tree(tree,last_tree);
		Record_Model(tree->mod,last_tree->mod);
		lnL_current = lnL_best;

		//steps_tried = 0;
		//steps_accepted = 0;

		//Optimiz_All_Free_Param(tree,1);

		for(iter = 0; iter < anneal.iters_per_temp; iter++){
			steps_tried++;
			//get our next proposition.
			Get_Neighbor_Proposition(tree,tau);
			lnL_proposed = tree->c_lnL;

			now = time(NULL);
			// Some useful debugging statements:
			//Print_Tree_Screen(tree);
			//if(verbose){
			PhyML_Printf("T: %f iter: %d current lnL = %f proposed lnL = %f\n", tau, iter, lnL_current, lnL_proposed);
			//PhyML_Printf("plot1 %d %lf\n", (steps_tried * (itau + 1)), lnL_current);
			//PhyML_Printf("plot2 %ld %d\n", (long)now, (steps_tried * (itau + 1)));
			PhyML_Printf("plot1 %d %lf\n", steps_tried, lnL_current);
			PhyML_Printf("plot2 %ld %d\n", (long)difftime(now,start), steps_tried);
			//}

			//      PhyML_Printf("JSJ: Proposed Likelihood at iter %i: %lf\n",iter,lnL_current);

			if(lnL_proposed > lnL_best){
				//save this tree into best_tree

				//PhyML_Printf("(annealing.c, 465): Found a new best tree with lnL = %f\n", lnL_best);
				Copy_Tree(tree,best_tree);
				Record_Model(tree->mod, best_tree->mod);

				lnL_best = lnL_proposed;

				tau_of_best = tau;
				iter_of_best = iter;

				// This tau is yielding good results: let's stay here
				// longer, thus increasing our chances of reaching the best ground
				// state at this tau.
				iter -= anneal.set_back;
				if(iter < 0) iter = 0; //make sure not set back into negative...

				acc_prob = 1.0;
				r = 0.0;
			}
			else if(lnL_proposed > lnL_current){
				acc_prob = 1.0;
				r = 0.0;
			}else{
				m3ldbl distance = Get_distance_between_trees(last_tree, tree);
				acc_prob = Boltzmann_P_QA(lnL_current, lnL_proposed, distance, tau); // VHS: change this to use QA-based Boltzmann
				r = gsl_rng_uniform(anneal.rng);
			}


			if(acc_prob >= r){
				//save the current tree
				if(anneal.no_change == 0){
					Copy_Tree(tree,last_tree);
					Record_Model(tree->mod,last_tree->mod);
				}
				steps_accepted++;
				lnL_current = tree->c_lnL;

				//PhyML_Printf("JSJ: Proposed Likelihood Accepted, acc_prob = %lf\n",acc_prob);

			}else{
				//restore the current tree to be our last tree
				if(anneal.no_change == 0){
					Copy_Tree(last_tree,tree);
					Record_Model(last_tree->mod,tree->mod);
				}
			}
		}//end inner for loop.
		tau *= taumult;

	}
	if(verbose){
		PhyML_Printf("Annealing finished.\n");
		PhyML_Printf("In tau range [%f, %f], the best tree was found at %f, iteration %f\n", tau, anneal.start_tau, tau_of_best, iter_of_best);
	}
	Copy_Tree(best_tree,tree);
	Record_Model(best_tree->mod,tree->mod);

	Free_Model(best_tree->mod);
	Free_Tree(best_tree);
	Free_Model(last_tree->mod);
	Free_Tree(last_tree);
	Free_Anneal(); //local function that frees the gsl_rng;
	Lk(tree);
	/* Here is psuedocode for the algorithm:
	 *
	 * m3ldbl temp = start_temp // start_temp is a global
	 * m3ldbl tempmult = exp(log(end_temp / start_temp) / double(num_anneal_stages - 1) )
	 *
	 * aratio = Scale_Acceptance_Ratio(tree) // accept_ratio is a global value
	 *
	 * lnL_best = tree's lnL
	 * tree_current = copy of tree
	 * lnL_current = tree_current's lnL
	 * tree_proposed = copy of tree
	 * lnL_proposed = tree_proposed's lnL
	 *
	 * for (int itemp = 0; itemp < num_anneal_stages; itemp++): // num_anneal_stages is a global
	 *
	 *              // Here we recenter the search at each tau, using the best values found thus far
	 *              tree_current = copy of tree
	 *              lnL_current = lnL_best
	 *
	 *              steps_tries = 0
	 *              steps_accepted = 0
	 *
	 *              for (int iter = 0; iter < iters_per_temp; iter++):
	 *                      steps_tried++;
	 *                      tree_proposed = copy of tree_current
	 *                      Get_Neighbor_Proposition(tree_proposed)
	 *                      lnL_proposed = tree_proposed's likelihood
	 *                      if lnL_proposed > lnL_best:
	 *                              lnL_best = lnL_proposed
	 *                              tree = copy of tree_proposed
	 *                              iter -= set_back // keep going at this tau if we are improving the likelihood
	 *                              if (iter < 0):
	 *                                      iter = 0
	 *                      acc_prob = Boltzmann_P(lnL_proposed, lnL_current, itemp)
	 *                      if acc_prob > random(0,1):
	 *                              steps_accepted++
	 *                              lnL_current = lnL_proposed
	 *                              tree_current = copy of tree_proposed
	 *
	 *              temp *= tempmult // reduce the tau with geometric descent
	 */
	return tree->c_lnL;

}

