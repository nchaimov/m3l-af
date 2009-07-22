#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"
#include "rates.h"
#include "annealing.h"
#include "numeric.h"

annealing anneal;

void Set_Anneal(){
	anneal.accept_ratio = 0.5;
	anneal.end_temp = 0.001;
	anneal.iters_per_temp = 100;
	anneal.set_back = 5;
	anneal.start_temp = 1.0;
	anneal.temp_count = 100;
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
	/* psuedo-code:

	lnL_best = tree's lnL
	lnL_proposed = max_lnL

	fsum = 0.0 // function mean
	fsqsum  0.0 // function standard deviation
	n = 0 // number of points tried

	tree_proposed = copy of tree

	for iter in range(0, iters_per_temp):
		Get_TA_Neighbor_Proposition(tree_proposed)
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
	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_proposed = tree->c_lnL;
	m3ldbl fsum = 0.0;
	m3ldbl fsqsum = 0.0;
	int n = 0;
	int i;
	arbre *best_tree = Make_Tree(tree->n_otu,tree->n_l);
	Init_Tree(best_tree,tree->n_otu, tree->n_l);
	Make_All_Tree_Nodes(best_tree);
	Make_All_Tree_Edges(best_tree);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree,best_tree);
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);

	For(i,anneal.iters_per_temp){
		Get_TA_Neighbor_Proposition(tree);
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
	fsum /= (double)n;
	fsqsum /= (double)n;
	fsqsum = sqrt(fsqsum - (fsum * fsum));
	aratio = anneal.accept_ratio * fsqsum / anneal.start_temp;

	Free_Model(best_tree->mod);
	Free_Tree(best_tree);

	return aratio;
}

// This method modifies some or all parameters of "tree", including:
// * topology
// * branch lengths
// * branch length proportions
// * alpha (for gamma distributed ASRV)
// * gamma proportions
// * mu (evolutionary rate)
void Get_TA_Neighbor_Proposition(arbre *tree){
	// For now, we'll always perturb every parameter.
	// In the future, we'll do something more sophisticated, where the probability of perturbing
	// any particular parameter will be drawn from a probability distribution.
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	Step_Brlen_Proportion(tree);
	Step_Branch_Lengths(tree);
	Step_Gamma(tree);
	Step_Topology(tree);
}

// helper for "Get_TA_Neighbor_Proposition"
void Step_Brlen_Proportion(arbre *tree){
	// For now, do something simple and stupid, like stepping branch lengths by +/- 0.001.
	// We'll eventually do something more sophisticated, like selected the new branch
	// lengths from a Dirichlet distribution.
	// when we update proportions we need to recalculate the likelihoods on the whole tree...
	if(tree->mod->s_opt->opt_props == 1){
		int i,j;
		double r = (((double)rand() + 1.0) / ((double)(RAND_MAX)+ 1.0));
		int prange = (int)(Rand_Int(1,(tree->n_l - 1)) * r);
		For(i,prange){
			j = Rand_Int(0,(tree->n_l - 1));
			r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0; //r is in the range -0.04999 to 0.049999
			tree->props[j] += r;
		}
		Normalize_Props(tree);
		//JSJ: recalculate the likelihood of the tree
		tree->both_sides = 1;
		Lk(tree);
	}
}

// helper for "Get_TA_Neighbor_Proposition"
void Step_Gamma(arbre *tree){
	// For now, do something stupid like incrementing gamma +/- 0.001.
	// In reality, we'll want to perturb the gamma parameter based on a normal distribution
	// with mean equal to the current value of gamma and standard deviation equal to some
	// user-specified parameter sigma.
	if(tree->mod->s_opt->opt_alpha == 1){
		m3ldbl r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0;
		tree->mod->alpha += r;
	}
	tree->both_sides = 1;
	Lk(tree);
}

// helper for "Get_TA_Neighbor_Proposition"
void Step_Branch_Lengths(arbre *tree){
	// For now, do something stupid like incrementing lengths by +/- 0.001.
	if(tree->mod->s_opt->opt_bl == 1){
		int i,j,edge_range,set_range,m,n;
		int n_edges = (tree->n_otu * 2) - 3;
		double r = (((double)rand() + 1.0) / ((double)(RAND_MAX)+ 1.0));
		/* r is a random floating point value in the range (0,1) {not including 0,
		 * or 1}. Note we must convert rand() and/or RAND_MAX+1 to
		 * floating point values to avoid integer division. In addition, Sean
		 * Scanlon pointed out the possibility that RAND_MAX may be the largest
		 * positive integer the architecture can represent, so (RAND_MAX+1)
		 * may result in an overflow, or more likely the value will end up being
		 * the largest negative integer the architecture can represent, so
		 * to avoid this we convert RAND_MAX and 1 to doubles before adding. */
		edge_range = (int)(Rand_Int(1,(n_edges - 1)) * r); //choose a random range shifted toward 1
		r = (((double)rand() + 1.0) / ((double)(RAND_MAX)+ 1.0));
		set_range = (int)(Rand_Int(1,(tree->n_l - 1)) * r);
		For(i,edge_range){
			j = Rand_Int(0,(n_edges - 1));
			For(m,set_range){
				r = ((((double)rand() + 1.0) / ((double)(RAND_MAX)+1.0)) - 0.5)/10.0; //r is in the range -0.04999 to 0.049999
				n = Rand_Int(0,(tree->n_l - 1));
				tree->t_edges[j]->l[n] += r;
				Update_Lk_At_Given_Edge(tree->t_edges[j],tree); //calls update_p_lk on appropriate nodes and this edge
			}
		}
	}
}

// helper for "Get_TA_Neighbor_Proposition"
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
		node *a,*b,*c,*d;
		int i,j;
		int n_edges = (tree->n_otu * 2) - 3;
		int edge = Rand_Int(0,(n_edges - 1));
		while(tree->t_edges[edge]->left->tax == 1 || tree->t_edges[edge]->rght->tax == 1){
			edge = ((edge + 1) % (n_edges - 1)); //starting from here, look through edge list
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
		tree->both_sides = 1;
		Lk(tree);


	}
}


void Get_QA_Neighbor_Proposition(arbre *tree){

}

// INPUT: Given the likelihood (lnl_curr) of our current state, the likelihood (lnl_new) of our new proposed state,
// and the temperature (temperature) of the system, calculate the probability of the new state,
// where the probability is drawn from a Boltzmann distribution.
//
// OUTPUT: a floating point decimal, between 0.0 and 1.0
m3ldbl Boltzmann_P(m3ldbl lnl_curr, m3ldbl lnl_new, m3ldbl temperature){
	m3ldbl result = 1.0;
	/* psuedocode:
	if (lnl_new < lnl_curr)
	{
		result = exp( (lnl_new - lnl_curr)/(accept_ratio * temperature) );
	}
	 */
	if(lnl_new < lnl_curr){
		result = exp((lnl_new - lnl_curr)/(anneal.accept_ratio * temperature));
	}
	return result;
}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
//
// At the end of this method, the object "tree" will contain the best-found topology, branch lengths,
// and model parameters.
m3ldbl Thermal_Anneal_All_Free_Params(arbre *tree, int verbose){
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	Set_Anneal();
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	m3ldbl result = 1.0;
	int n_edges = (tree->n_otu * 2) - 3;
	if (n_edges <= 3){
		tree->mod->s_opt->opt_topo = 0; //make sure that opt_topo is false if there are no meaningful branch swaps.
	}
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	m3ldbl temp = anneal.start_temp;
	m3ldbl tempmult = exp(log(anneal.end_temp/anneal.start_temp)/(((double)anneal.temp_count) - 1.0));
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	m3ldbl aratio = Scale_Acceptance_Ratio(tree);
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
	arbre *best_tree = Make_Tree(tree->n_otu,tree->n_l);
	Init_Tree(best_tree,tree->n_otu, tree->n_l);
	Make_All_Tree_Nodes(best_tree);
	Make_All_Tree_Edges(best_tree);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree,best_tree);

	arbre *last_tree = Make_Tree(tree->n_otu,tree->n_l);
	Init_Tree(last_tree,tree->n_otu, tree->n_l);
	Make_All_Tree_Nodes(last_tree);
	Make_All_Tree_Edges(last_tree);

	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_current = tree->c_lnL;
	m3ldbl lnL_proposed = tree->c_lnL;
	m3ldbl acc_prob;
	double r;
	int itemp,iter;
	int steps_tried;
	int steps_accepted;
	PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);

	For(itemp,anneal.temp_count){
		//recenter the search at each temperature.
		last_tree->mod = Copy_Model(best_tree->mod);
		Copy_Tree(best_tree,tree);
		Record_Model(best_tree->mod,tree->mod);
		lnL_current = lnL_best;

		steps_tried = 0;
		steps_accepted = 0;

		For(iter,anneal.iters_per_temp){
			steps_tried++;
			//record current tree and model
			Copy_Tree(tree,last_tree);
			Record_Model(tree->mod,last_tree->mod);
			//get our next proposition.
			Get_TA_Neighbor_Proposition(tree);

			lnL_proposed = tree->c_lnL;
			if(lnL_proposed > lnL_best){
				//save this tree into best_tree
				PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
				Copy_Tree(tree,best_tree);
				Record_Model(tree->mod, best_tree->mod);
				lnL_best = lnL_proposed;

				iter -= anneal.set_back;
				if(iter < 0) iter = 0; //make sure not set back into negative...
			}
			r = ((double)rand() / ((double)(RAND_MAX)+ 1.0));
			acc_prob = Boltzmann_P(lnL_current, lnL_proposed, temp);

			if(acc_prob > r){
				steps_accepted++;
				lnL_current = tree->c_lnL;
				PhyML_Printf("JSJ: Made it to line %d, in file %s\n",__LINE__,__FILE__);
			}else{
				//restore the current tree
				Copy_Tree(last_tree,tree);
				Record_Model(last_tree->mod,tree->mod);
			}
		}//end inner for loop.
		temp *= tempmult;

	}
	/* Here is psuedocode for the algorithm:
	 *
	 * m3ldbl temp = start_temp // start_temp is a global
	 * m3ldbl tempmult = exp(log(end_temp / start_temp) / double(temp_count - 1) )
	 *
	 * aratio = Scale_Acceptance_Ratio(tree) // accept_ratio is a global value
	 *
	 * lnL_best = tree's lnL
	 * tree_current = copy of tree
	 * lnL_current = tree_current's lnL
	 * tree_proposed = copy of tree
	 * lnL_proposed = tree_proposed's lnL
	 *
	 * for (int itemp = 0; itemp < temp_count; itemp++): // temp_count is a global
	 *
	 *		// Here we recenter the search at each temperature, using the best values found thus far
	 * 		tree_current = copy of tree
	 * 		lnL_current = lnL_best
	 *
	 *		steps_tries = 0
	 * 		steps_accepted = 0
	 *
	 *		for (int iter = 0; iter < iters_per_temp; iter++):
	 *			steps_tried++;
	 *			tree_proposed = copy of tree_current
	 *			Get_TA_Neighbor_Proposition(tree_proposed)
	 *			lnL_proposed = tree_proposed's likelihood
	 *			if lnL_proposed > lnL_best:
	 *				lnL_best = lnL_proposed
	 *				tree = copy of tree_proposed
	 *				iter -= set_back // keep going at this temperature if we are improving the likelihood
	 *				if (iter < 0):
	 *					iter = 0
	 *			acc_prob = Boltzmann_P(lnL_proposed, lnL_current, itemp)
	 *			if acc_prob > random(0,1):
	 *				steps_accepted++
	 *				lnL_current = lnL_proposed
	 *				tree_current = copy of tree_proposed
	 *
	 *		temp *= tempmult // reduce the temperature with geometric descent
	 */
	return result;
}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
m3ldbl Quantum_Anneal_All_Free_Params(arbre *tree, int verbose){
	Set_Anneal();
	m3ldbl result = 1.0;
	int n_edges = (tree->n_otu * 2) - 3;
	if (n_edges <= 3){
		tree->mod->s_opt->opt_topo = 0; //make sure that opt_topo is false if there are no meaningful branch swaps.
	}

	return result;

}

