/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

 */

/*

Implementation of aLRT branch tests, and 5-branchs NNI research optimization.

Authors : Jean-Francois Dufayard & Stephane Guindon.

 */


#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "simu.h"
#include "alrt.h"

/*********************************************************/

/*
 * Check every testable branch of the tree,
 * check for NNIs, optimizing 5 branches,
 * param tree : the tree to check
 */

int Check_NNI_Five_Branches(arbre *tree)
{
	int best_edge;
	m3ldbl best_gain;
	int best_config;
	int i;
	int better_found; /* = 1 if a phylogeny with greater likelihood than current one was found */
	int result;
	m3ldbl init_lnL;

	init_lnL     = UNLIKELY;
	better_found = 1;

	//While there is at least one NNI to do...
	while(better_found)
	{
		Update_Dirs(tree);

		//Interface output
		if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n\n. Checking for NNIs, optimizing five branches...\n");

		better_found  =  0;
		result        = -1;
		best_edge     = -1;
		best_gain     = 1.E+10;
		best_config   = -1;

		tree->both_sides = 1;
		init_lnL = Return_Lk(tree);

		//For every branch
		For(i,2*tree->n_otu-3)
		{
			//if this branch is not terminal
			if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax))
			{
				result = -1;

				//Try every NNIs for the tested branch
				result = NNI_Neigh_BL(tree->t_edges[i],tree);

				//Look for possible NNI to do, and check if it is the best one
				switch(result)
				{
				case 1 : /* lk1 > lk0 > lk2 */
				{
					if((tree->t_edges[i]->nni->lk0 - tree->t_edges[i]->nni->lk1) < best_gain)
					{
						better_found = 1;
						best_edge    = i;
						best_gain    = tree->t_edges[i]->nni->lk0-tree->t_edges[i]->nni->lk1;
						best_config  = 1;
					}
					break;
				}
				case 2 : /* lk2 > lk0 > lk1 */
				{
					if((tree->t_edges[i]->nni->lk0 - tree->t_edges[i]->nni->lk2) < best_gain)
					{
						better_found = 1;
						best_edge    = i;
						best_gain    = tree->t_edges[i]->nni->lk0-tree->t_edges[i]->nni->lk2;
						best_config  = 2;
					}
					break;
				}
				case 3 : /* lk1 > lk2 > lk0 */
				{
					if((tree->t_edges[i]->nni->lk2 - tree->t_edges[i]->nni->lk1) < best_gain)
					{
						better_found = 1;
						best_edge    = i;
						best_gain    = tree->t_edges[i]->nni->lk2-tree->t_edges[i]->nni->lk1;
						best_config  = 1;
					}
					break;
				}
				case 4 : /* lk2 > lk1 > lk0 */
				{
					if((tree->t_edges[i]->nni->lk1 - tree->t_edges[i]->nni->lk2) < best_gain)
					{
						better_found = 1;
						best_edge    = i;
						best_gain    = tree->t_edges[i]->nni->lk1-tree->t_edges[i]->nni->lk2;
						best_config  = 2;
					}
					break;
				}
				default : /* lk2 = lk1 = lk0 */
				{
					if(best_gain > .0) best_gain = .0;
					break;
				}
				}
			}
		}

		if((tree->c_lnL < init_lnL - tree->mod->s_opt->min_diff_lk_local) || (tree->c_lnL > init_lnL + tree->mod->s_opt->min_diff_lk_local))
		{
			PhyML_Printf("\n\n. tree->c_lnL = %f init_lnL = %f.",tree->c_lnL,init_lnL);
			PhyML_Printf("\n. Err in file %s at line %d.\n",__FILE__,__LINE__);
			Warn_And_Exit("\n");
		}

		//Don't do any NNI if the user doesn't want to optimize topology
		if(!tree->mod->s_opt->opt_topo) better_found = 0;
		/*       if(fabs(best_gain) <= tree->mod->s_opt->min_diff_lk_move) better_found = 0; */

		//If there is one swap to do, make the best one.
		if(better_found)
		{
			Make_Target_Swap(tree,tree->t_edges[best_edge],best_config);
			Lk(tree);

			if(tree->c_lnL < init_lnL)
			{
				PhyML_Printf("\n\n. tree->c_lnL = %f init_lnL = %f.",tree->c_lnL,init_lnL);
				PhyML_Printf("\n. Err in file %s at line %d.\n",__FILE__,__LINE__);
				Warn_And_Exit("\n");
			}

			if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Topology           ]");
			return 1;
		}
	}
	return 0;
}

/*********************************************************/

/* Compute aLRT supports */
void aLRT(arbre *tree)
{
	int i;


	/* aLRT support will label each internal branch */
	tree->print_alrt_val = 1;

	/* The topology will not be modified when assessing the branch support. We make sure that it will
     not be modified afterwards by locking the topology */

	tree->both_sides = 1;
	Lk(tree);


	For(i,2*tree->n_otu-3)
	if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax))
	{
		/* Compute likelihoods for each of the three configuration */
		NNI_Neigh_BL(tree->t_edges[i],tree);
		/* Compute the corresponding statistical support */
		Compute_Likelihood_Ratio_Test(tree->t_edges[i],tree);
	}

	tree->lock_topo = 1;
}

/*********************************************************/
/*
 * Launch one branch testing,
 * analyse the result
 * and convert supports as wished by the user.
 * param tree : the tree to check
 * param tested_edge : the tested edge of the tree
 * param old_loglk : the initial likelihood, before any aLRT analysis
 * param isBoot : 1 if used from the Bootstrap procedure, 0 if not
 * return an integer, informative to analyse the results and potential NNIs to do
 */
int Compute_Likelihood_Ratio_Test(edge *tested_edge, arbre *tree)
{
	int result=0;

	tested_edge->ratio_test     =  0.0;
	tested_edge->alrt_statistic = -1.0;

	if((tested_edge->nni->lk0 > tested_edge->nni->lk1) && (tested_edge->nni->lk0 > tested_edge->nni->lk2))
	{
		if(tested_edge->nni->lk1 < tested_edge->nni->lk2)
		{
			//lk0 > lk2 > lk1
			tested_edge->alrt_statistic = 2*(tested_edge->nni->lk0 - tested_edge->nni->lk2);
		}
		else
		{
			//lk0 > lk1 >= lk2
			tested_edge->alrt_statistic = 2*(tested_edge->nni->lk0 - tested_edge->nni->lk1);
		}

		if (tested_edge->alrt_statistic < 0.0)
		{
			tested_edge->alrt_statistic = 0.0;
			tested_edge->ratio_test = 0.0;
		}
		else
		{
			//aLRT statistic is valid, compute the wished support
			if (tree->io->ratio_test == 2)
				tested_edge->ratio_test = Statistics_To_Probabilities(tested_edge->alrt_statistic);

			else if(tree->io->ratio_test == 3)
			{
				m3ldbl sh_support;
				m3ldbl param_support;

				sh_support    = Statistics_To_SH(tree);
				param_support = Statistics_To_Probabilities(tested_edge->alrt_statistic);

				if(sh_support < param_support) tested_edge->ratio_test = sh_support;
				else                           tested_edge->ratio_test = param_support;
			}

			else if(tree->io->ratio_test == 1)
			{
				tested_edge->ratio_test=tested_edge->alrt_statistic;
			}
			else
			{
				tested_edge->ratio_test = Statistics_To_SH(tree);
			}
		}
	}
	//statistic is not valid, give the negative statistics if wished, or a zero support.
	else if((tested_edge->nni->lk1 > tested_edge->nni->lk0) && (tested_edge->nni->lk1 > tested_edge->nni->lk2))
	{
		/*       tested_edge->ratio_test = 2*(tested_edge->nni->lk0-tested_edge->nni->lk1); */
		tested_edge->ratio_test = 0.0;
		if(tree->io->ratio_test > 1) tested_edge->alrt_statistic = 0.0;
	}
	else if((tested_edge->nni->lk2 > tested_edge->nni->lk0) && (tested_edge->nni->lk2 > tested_edge->nni->lk1))
	{
		/*       tested_edge->ratio_test = 2*(tested_edge->nni->lk0-tested_edge->nni->lk2); */
		tested_edge->ratio_test = 0.0;
		if(tree->io->ratio_test > 1) tested_edge->alrt_statistic = 0.0;
	}
	else // lk0 ~ lk1 ~ lk2
	{
		tested_edge->ratio_test = 0.0;
		if(tree->io->ratio_test > 1) tested_edge->ratio_test = 0.0;
	}
	return result;
}

/*********************************************************/
/*
 * Test the 3 NNI positions for one branch.
 * param tree : the tree to check
 * param tested_edge : the tested edge of the tree
 * param old_loglk : the initial likelihood, before any aLRT analysis
 * param isBoot : 1 if used from the Bootstrap procedure, 0 if not
 */
int NNI_Neigh_BL(edge *b_fcus, arbre *tree)
{
	int n_patterns;
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4, site;
	node *v1,*v2,*v3,*v4;
	edge *e1,*e2,*e3,*e4;
	m3ldbl lk0, lk1, lk2;
	m3ldbl lk_init, lk_temp;
	int i,j;
	int result,counter,wei;

	m3ldbl len_e1[MAX_BL_SET];
	m3ldbl len_e2[MAX_BL_SET];
	m3ldbl len_e3[MAX_BL_SET];
	m3ldbl len_e4[MAX_BL_SET];
	m3ldbl bl_init[MAX_BL_SET];
	m3ldbl l1[MAX_BL_SET];
	m3ldbl l2[MAX_BL_SET];
	m3ldbl l3[MAX_BL_SET];
//	m3ldbl l_infa[MAX_BL_SET];
//	m3ldbl l_infb[MAX_BL_SET];
//	m3ldbl l_max[MAX_BL_SET];



	result = 0;

	//Initialization
	//JSJ: temp fix to l
	For(j,tree->n_l) bl_init[j] = b_fcus->l[j];
	lk_init = tree->c_lnL;
	lk_temp = UNLIKELY;
	n_patterns = tree->n_pattern;

	b_fcus->nni->score = .0;

	lk0 = lk1 = lk2    = UNLIKELY;
	v1 = v2 = v3 = v4  = NULL;

	l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

	l_r = b_fcus->l_r;
	r_l = b_fcus->r_l;

	//vertices
	v1 = b_fcus->left->v[b_fcus->l_v1];
	v2 = b_fcus->left->v[b_fcus->l_v2];
	v3 = b_fcus->rght->v[b_fcus->r_v1];
	v4 = b_fcus->rght->v[b_fcus->r_v2];

	if(v1->num < v2->num)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	if(v3->num < v4->num)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	For(j,tree->n_l){
		l1[j] = l2[j] = l3[j] = -1.;
		l1[j] = b_fcus->l[j];
	}

	//edges
	e1 = b_fcus->left->b[b_fcus->l_v1];
	e2 = b_fcus->left->b[b_fcus->l_v2];
	e3 = b_fcus->rght->b[b_fcus->r_v1];
	e4 = b_fcus->rght->b[b_fcus->r_v2];

	//JSJ: fixed l
	//record initial branch lengths
	For(j,tree->n_l){
		len_e1[j] = e1->l[j];
		len_e2[j] = e2->l[j];
		len_e3[j] = e3->l[j];
		len_e4[j] = e4->l[j];
	}

	//Optimize branch lengths and update likelihoods for
	//the original configuration.

	/**
	 * JSJ: for testing purposes
	 */
	PhyML_Printf("\n. Printing tree from file %s at line %d\n",__FILE__,__LINE__);
	Print_Tree_Screen(tree);
	/**
	 * JSJ: end tree test print section
	 */

	do
	{
		lk0 = lk_temp;

		For(i,3)
		if(b_fcus->left->v[i] != b_fcus->rght)
		{ //JSJ: fixed l
			Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->left->b[i]->l[j];
//				l_max[j]   = b_fcus->left->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				/**
				 * JSJ: for testing purposes
				 */
				PhyML_Printf("\n. Printing tree from file %s at line %d\n",__FILE__,__LINE__);
				Print_Tree_Screen(tree);
				/**
				 * JSJ: end tree test print section
				 */
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->left->b[i]->l[j],b_fcus->left->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->left->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->left->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}

		Update_P_Lk(tree,b_fcus,b_fcus->left);
		//JSJ: temp fix to l
		For(j,tree->n_l){
//			l_infa[j]  = 10.*b_fcus->l[j];
//			l_max[j]   = b_fcus->l[j];
//			l_infb[j]  = BL_MIN;
			/**
			 * JSJ: for testing purposes
			 */
			PhyML_Printf("\n. Printing tree from file %s at line %d\n",__FILE__,__LINE__);
			Print_Tree_Screen(tree);
			/**
			 * JSJ: end tree test print section
			 */
			lk_temp = Br_Len_Brent_Iter(10.*b_fcus->l[j],b_fcus->l[j],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,0,j);
		}

//		lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//				tree->mod->s_opt->min_diff_lk_local,
//				b_fcus,tree,
//				tree->mod->s_opt->brent_it_max,0);


		For(i,3)
		if(b_fcus->rght->v[i] != b_fcus->left)
		{
			Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);
			//JSJ: temp fix to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->rght->b[i]->l[j];
//				l_max[j]   = b_fcus->rght->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				/**
				 * JSJ: for testing purposes
				 */
				PhyML_Printf("\n. Printing tree from file %s at line %d\n",__FILE__,__LINE__);
				Print_Tree_Screen(tree);
				/**
				 * JSJ: end tree test print section
				 */
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->rght->b[i]->l[j],b_fcus->rght->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->rght->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);

			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->rght->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}

		Update_P_Lk(tree,b_fcus,b_fcus->rght);

		if(lk_temp < lk0 - (tree->mod->s_opt->min_diff_lk_local))
		{
			/**
			 * JSJ: for testing purposes
			 */
			PhyML_Printf("\n. Printing tree from file %s at line %d\n",__FILE__,__LINE__);
			Print_Tree_Screen(tree);
			/**
			 * JSJ: end tree test print section
			 */
			PhyML_Printf("\n. lk_temp = %f lk0 = %f\n",lk_temp,lk0);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	//		Warn_And_Exit(""); //JSJ: temporarily comment out this line...
		}
	}
	while(fabs(lk_temp-lk0) > tree->mod->s_opt->min_diff_lk_global);


	lk0 = tree->c_lnL;
	counter=0;
	For(site,n_patterns)
	{
		wei=0;
		For(wei,tree->data->wght[site])
		{
			tree->log_lks_aLRT[0][counter]= tree->c_lnL_sorted[site] / tree->data->wght[site];
			counter++;
		}
	}
	//JSJ: fixed l
	/* Go back to initial branch lengths */
	For(j,tree->n_l){
		e1->l[j]     = len_e1[j];
		e2->l[j]     = len_e2[j];
		e3->l[j]     = len_e3[j];
		e4->l[j]     = len_e4[j];
		b_fcus->l[j] = bl_init[j];
	}
	Update_PMat_At_Given_Edge(e1,tree);
	Update_PMat_At_Given_Edge(e2,tree);
	Update_PMat_At_Given_Edge(e3,tree);
	Update_PMat_At_Given_Edge(e4,tree);
	Update_PMat_At_Given_Edge(b_fcus,tree);


	//Do first possible swap
	Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);


	tree->both_sides = 1;
	lk1 = Update_Lk_At_Given_Edge(b_fcus,tree);

	For(i,3)
	if(b_fcus->left->v[i] != b_fcus->rght)
		Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);

	For(i,3)
	if(b_fcus->rght->v[i] != b_fcus->left)
		Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);


	//Optimize branch lengths and update likelihoods
	lk_temp = UNLIKELY;
	do
	{
		lk1 = lk_temp;

		For(i,3)
		if(b_fcus->left->v[i] != b_fcus->rght)
		{
			Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);
			//JSJ: temporary fixes to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->left->b[i]->l[j];
//				l_max[j]   = b_fcus->left->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->left->b[i]->l[j],b_fcus->left->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->left->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->left->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}


		Update_P_Lk(tree,b_fcus,b_fcus->left);
		//JSJ: temp fixes to l
		For(j,tree->n_l){
//			l_infa[j]  = 10.*b_fcus->l[j];
//			l_max[j]   = b_fcus->l[j];
//			l_infb[j]  = BL_MIN;
			lk_temp = Br_Len_Brent_Iter(10.*b_fcus->l[j],b_fcus->l[j],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,0,j);
		}
//		lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//				tree->mod->s_opt->min_diff_lk_local,
//				b_fcus,tree,
//				tree->mod->s_opt->brent_it_max,0);



		For(i,3)
		if(b_fcus->rght->v[i] != b_fcus->left)
		{
			Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);
			//JSJ: temp fixes to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->rght->b[i]->l[j];
//				l_max[j]   = b_fcus->rght->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->rght->b[i]->l[j],b_fcus->rght->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->rght->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->rght->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}

		Update_P_Lk(tree,b_fcus,b_fcus->rght);


		if(lk_temp < lk1 - (tree->mod->s_opt->min_diff_lk_local))
		{
			PhyML_Printf("\n. lk_temp = %f lk1 = %f\n",lk_temp,lk1);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		//	Warn_And_Exit(""); //JSJ: commented out temporarily
		}
	}
	while(fabs(lk_temp-lk1) > tree->mod->s_opt->min_diff_lk_global);
	//until no significative improvement is detected

	lk1 = tree->c_lnL;

	//save likelihood of each sites, in order to compute branch supports
	counter=0;
	For(site,n_patterns)
	{
		wei=0;
		For(wei,tree->data->wght[site])
		{
			tree->log_lks_aLRT[1][counter]= tree->c_lnL_sorted[site] / tree->data->wght[site];
			counter++;
		}
	}


	//undo the swap
	Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);

	/* Go back to initial branch lengths */
	For(j,tree->n_l){
		l2[j]  = b_fcus->l[j];
		e1->l[j]     = len_e1[j];
		e2->l[j]     = len_e2[j];
		e3->l[j]     = len_e3[j];
		e4->l[j]     = len_e4[j];
		b_fcus->l[j] = bl_init[j];
	}

	Update_PMat_At_Given_Edge(e1,tree);
	Update_PMat_At_Given_Edge(e2,tree);
	Update_PMat_At_Given_Edge(e3,tree);
	Update_PMat_At_Given_Edge(e4,tree);
	Update_PMat_At_Given_Edge(b_fcus,tree);
	//JSJ: temp fixes to l
	/***********/
	//do the second possible swap
	Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	For(j,tree->n_l) b_fcus->l[j] = bl_init[j];
	tree->both_sides = 1;

	lk2 = Update_Lk_At_Given_Edge(b_fcus,tree);

	For(i,3)
	if(b_fcus->left->v[i] != b_fcus->rght)
		Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);

	For(i,3)
	if(b_fcus->rght->v[i] != b_fcus->left)
		Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);


	//Optimize branch lengths and update likelihoods
	lk_temp = UNLIKELY;
	do
	{
		lk2 = lk_temp;

		For(i,3)
		if(b_fcus->left->v[i] != b_fcus->rght)
		{
			Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);
			//JSJ: temp fixes to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->left->b[i]->l[j];
//				l_max[j]   = b_fcus->left->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->left->b[i]->l[j],b_fcus->left->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->left->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->left->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}


		Update_P_Lk(tree,b_fcus,b_fcus->left);
		//JSJ: temp fixes to l
		For(j,tree->n_l){
//			l_infa[j]  = 10.*b_fcus->l[j];
//			l_max[j]   = b_fcus->l[j];
//			l_infb[j]  = BL_MIN;
			lk_temp = Br_Len_Brent_Iter(10.*b_fcus->l[j],b_fcus->l[j],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,0,j);
		}
//		lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//				tree->mod->s_opt->min_diff_lk_local,
//				b_fcus,tree,
//				tree->mod->s_opt->brent_it_max,0);

		For(i,3)
		if(b_fcus->rght->v[i] != b_fcus->left)
		{
			Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);
			//JSJ: temp fixes to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->rght->b[i]->l[j];
//				l_max[j]   = b_fcus->rght->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->rght->b[i]->l[j],b_fcus->rght->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->rght->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->rght->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}

		Update_P_Lk(tree,b_fcus,b_fcus->rght);

		if(lk_temp < lk2 - (tree->mod->s_opt->min_diff_lk_local))
		{
			PhyML_Printf("\n. lk_temp = %f lk2 = %f\n",lk_temp,lk2);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	//		Warn_And_Exit(""); //JSJ: temporarily comment out
		}
	}
	while(fabs(lk_temp-lk2) > tree->mod->s_opt->min_diff_lk_global);
	//until no significative improvement is detected

	lk2 = tree->c_lnL;

	//save likelihood of each sites, in order to compute branch supports
	counter=0;
	For(site,n_patterns)
	{
		wei=0;
		For(wei,tree->data->wght[site])
		{
			tree->log_lks_aLRT[2][counter]= tree->c_lnL_sorted[site] / tree->data->wght[site];
			counter++;
		}
	}
	//undo the swap
	Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
	/***********/

	For(j,tree->n_l){
		//JSJ: temp fixes to l
		//save current length
		l3[j]  = b_fcus->l[j];
		//restore the initial branch length values
		e1->l[j]     = len_e1[j];
		e2->l[j]     = len_e2[j];
		e3->l[j]     = len_e3[j];
		e4->l[j]     = len_e4[j];
		b_fcus->l[j] = bl_init[j];
	}

	//recompute likelihoods
	Update_PMat_At_Given_Edge(e1,tree);
	Update_PMat_At_Given_Edge(e2,tree);
	Update_PMat_At_Given_Edge(e3,tree);
	Update_PMat_At_Given_Edge(e4,tree);
	Update_PMat_At_Given_Edge(b_fcus,tree);

	Update_P_Lk(tree,b_fcus,b_fcus->rght);
	Update_P_Lk(tree,b_fcus,b_fcus->left);

	For(i,3)
	if(b_fcus->left->v[i] != b_fcus->rght)
		Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);

	For(i,3)
	if(b_fcus->rght->v[i] != b_fcus->left)
		Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);

	lk_temp = Lk_At_Given_Edge(b_fcus,tree);


	if((lk_temp > lk_init + tree->mod->s_opt->min_diff_lk_local) || (lk_temp < lk_init - tree->mod->s_opt->min_diff_lk_local))
	{
		PhyML_Printf("\n. lk_temp = %f lk_init = %f\n",lk_temp,lk_init);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//		Warn_And_Exit("");
	}


	//save likelihoods in NNI structures
	b_fcus->nni->lk0 = lk0;
	b_fcus->nni->lk1 = lk1;
	b_fcus->nni->lk2 = lk2;

	b_fcus->nni->score = lk0 - MAX(lk1,lk2);

	if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
			(b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
	{
		b_fcus->nni->score = .0;
		b_fcus->nni->lk1 = b_fcus->nni->lk2 = b_fcus->nni->lk0;
	}


	if((b_fcus->nni->lk1 > b_fcus->nni->lk0) && (b_fcus->nni->lk1 > b_fcus->nni->lk2))
	{
		if(b_fcus->nni->lk0 > b_fcus->nni->lk2) result = 1; //lk1 > lk0 > lk2
		else                                    result = 3; //lk1 > lk2 > lk0
	}
	else if((b_fcus->nni->lk2 > b_fcus->nni->lk0) && (b_fcus->nni->lk2 > b_fcus->nni->lk1))
	{
		if(b_fcus->nni->lk0 > b_fcus->nni->lk1) result = 2; //lk2 > lk0 > lk1
		else                                    result = 4; //lk2 > lk1 > lk0
	}

	/*   int counter=0; */
	/*   For(site,n_patterns) */
	/*     { */
	/*       PhyML_Printf("%3d %3d %15f ",b_fcus->num,tree->data->wght[site],b_fcus->nni->score); */
	/*       PhyML_Printf("%15f %15f\n", */
	/* 	     tree->log_lks_aLRT[0][counter], */
	/* 	     (lk1 > lk2)?(tree->log_lks_aLRT[1][counter]):(tree->log_lks_aLRT[2][counter])); */
	/*       counter+=tree->data->wght[site]; */
	/*     } */
	return result;
}

/*********************************************************/
/*
 * Make one target swap, optimizing five branches.
 * param tree : the tree to check
 * param tested_edge : the swaping edge of the tree
 * param swapToDo : 1 or 2, to select the first or the second swap to do
 */

void Make_Target_Swap(arbre *tree, edge *b_fcus, int swaptodo)
{
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
	node *v1,*v2,*v3,*v4;
	edge *e1,*e2,*e3,*e4;
	m3ldbl lktodo;
	m3ldbl bl_init;
	m3ldbl lk_init, lk_temp;
	int i,j;

	m3ldbl l_infa[MAX_BL_SET];
	m3ldbl l_infb[MAX_BL_SET];
	m3ldbl l_max[MAX_BL_SET];

	if(swaptodo < 0)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	//JSJ: temp fix to l
	//Initialization
	bl_init = b_fcus->l[0];
	lk_init = tree->c_lnL;

	b_fcus->nni->score = .0;

	lktodo = UNLIKELY;
	v1 = v2 = v3 = v4 = NULL;

	l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

	l_r = b_fcus->l_r;
	r_l = b_fcus->r_l;

	v1 = b_fcus->left->v[b_fcus->l_v1];
	v2 = b_fcus->left->v[b_fcus->l_v2];
	v3 = b_fcus->rght->v[b_fcus->r_v1];
	v4 = b_fcus->rght->v[b_fcus->r_v2];

	if(v1->num < v2->num)
	{
		Free(l_infa);
		Free(l_max);
		Free(l_infb);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	if(v3->num < v4->num)
	{
		Free(l_infa);
		Free(l_max);
		Free(l_infb);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	e1 = b_fcus->left->b[b_fcus->l_v1];
	e2 = b_fcus->left->b[b_fcus->l_v2];
	e3 = b_fcus->rght->b[b_fcus->r_v1];
	e4 = b_fcus->rght->b[b_fcus->r_v2];


	/***********/
	//make the selected swap
	if(swaptodo==1)
	{
		Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
	}
	else
	{
		Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	}

	tree->both_sides = 1;
	lktodo = Update_Lk_At_Given_Edge(b_fcus,tree);

	For(i,3)
	if(b_fcus->left->v[i] != b_fcus->rght)
		Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);

	For(i,3)
	if(b_fcus->rght->v[i] != b_fcus->left)
		Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);


	//Optimize branch lengths and update likelihoods
	lk_temp = UNLIKELY;
	do
	{
		lktodo = lk_temp;

		For(i,3)
		if(b_fcus->left->v[i] != b_fcus->rght)
		{
			Update_P_Lk(tree,b_fcus->left->b[i],b_fcus->left);
			//JSJ: temp fixes to l
			For(j,tree->n_l){
//				l_infa[j] = 10.*b_fcus->left->b[i]->l[j];
//				l_max[j]  = b_fcus->left->b[i]->l[j];
//				l_infb[j] = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->left->b[i]->l[j],b_fcus->left->b[i]->l[j],BL_MIN,
									tree->mod->s_opt->min_diff_lk_local,
									b_fcus->left->b[i],tree,
									tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->left->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}


		Update_P_Lk(tree,b_fcus,b_fcus->left);
		//JSJ: temp fixes to l
		For(j,tree->n_l){
//			l_infa[j]  = 10.*b_fcus->l[j];
//			l_max[j]   = b_fcus->l[j];
//			l_infb[j]  = BL_MIN;
			lk_temp = Br_Len_Brent_Iter(10.*b_fcus->l[j],b_fcus->l[j],BL_MIN,
							tree->mod->s_opt->min_diff_lk_local,
							b_fcus,tree,
							tree->mod->s_opt->brent_it_max,0,j);
		}
//		lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//				tree->mod->s_opt->min_diff_lk_local,
//				b_fcus,tree,
//				tree->mod->s_opt->brent_it_max,0);



		For(i,3)
		if(b_fcus->rght->v[i] != b_fcus->left)
		{
			Update_P_Lk(tree,b_fcus->rght->b[i],b_fcus->rght);
			//JSJ: temp fixes to l
			For(j,tree->n_l){
//				l_infa[j]  = 10.*b_fcus->rght->b[i]->l[j];
//				l_max[j]   = b_fcus->rght->b[i]->l[j];
//				l_infb[j]  = BL_MIN;
				lk_temp = Br_Len_Brent_Iter(10.*b_fcus->rght->b[i]->l[j],b_fcus->rght->b[i]->l[j],BL_MIN,
						tree->mod->s_opt->min_diff_lk_local,
						b_fcus->rght->b[i],tree,
						tree->mod->s_opt->brent_it_max,0,j);
			}
//			lk_temp = Br_Len_Brent(l_infa,l_max,l_infb,
//					tree->mod->s_opt->min_diff_lk_local,
//					b_fcus->rght->b[i],tree,
//					tree->mod->s_opt->brent_it_max,0);
		}

		Update_P_Lk(tree,b_fcus,b_fcus->rght);


		if(lk_temp < lktodo - tree->mod->s_opt->min_diff_lk_local)
		{
			PhyML_Printf("\n. Edge %3d lk_temp = %f lktodo = %f\n",b_fcus->num,lk_temp,lktodo);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}
	}
	while(fabs(lk_temp-lktodo) > tree->mod->s_opt->min_diff_lk_global);
	//until no significative improvement is detected


	/*   PhyML_Printf("\n.<< [%3d] l=%f lk_init=%f tree->c_lnL=%f score=%12f v1=%3d v2=%3d v3=%3d v4=%3d >>", */
	/* 	 b_fcus->num, */
	/* 	 b_fcus->l, */
	/* 	 lk_init, */
	/* 	 tree->c_lnL, */
	/* 	 lk_init-tree->c_lnL, */
	/* 	 v1->num,v2->num,v3->num,v4->num);       */


	if(tree->c_lnL < lk_init)
	{
		PhyML_Printf("\n. [%3d] v1=%d v2=%d v3=%d v4=%d",
				b_fcus->num,v1->num,v2->num,v3->num,v4->num);
		PhyML_Printf("\n. tree->c_lnL = %f lk_init = %f\n",tree->c_lnL,lk_init);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

}

/*********************************************************/
/**
* Convert an aLRT statistic to a none parametric support
* param in: the statistic
*/

m3ldbl Statistics_To_Probabilities(m3ldbl in)
{
	m3ldbl rough_value=0.0;
	m3ldbl a=0.0;
	m3ldbl b=0.0;
	m3ldbl fa=0.0;
	m3ldbl fb=0.0;

	if(in>=0.000000393 && in<0.00000157)
	{
		a=0.000000393;
		b=0.00000157;
		fa=0.0005;
		fb=0.001;
	}
	else if(in>=0.00000157 && in<0.0000393)
	{
		a=0.00000157;
		b=0.0000393;
		fa=0.001;
		fb=0.005;
	}
	else if(in>=0.0000393 && in<0.000157)
	{
		a=0.0000393;
		b=0.000157;
		fa=0.005;
		fb=0.01;
	}
	else if(in>=0.000157 && in<0.000982)
	{
		a=0.000157;
		b=0.000982;
		fa=0.01;
		fb=0.025;
	}
	else if(in>0.000982 && in<0.00393)
	{
		a=0.000982;
		b=0.00393;
		fa=0.025;
		fb=0.05;
	}
	else if(in>=0.00393 && in<0.0158)
	{
		a=0.00393;
		b=0.0158;
		fa=0.05;
		fb=0.1;
	}
	else if(in>=0.0158 && in<0.0642)
	{
		a=0.0158;
		b=0.0642;
		fa=0.1;
		fb=0.2;
	}
	else if(in>=0.0642 && in<0.148)
	{
		a=0.0642;
		b=0.148;
		fa=0.2;
		fb=0.3;
	}
	else if(in>=0.148 && in<0.275)
	{
		a=0.148;
		b=0.275;
		fa=0.3;
		fb=0.4;
	}
	else if(in>=0.275 && in<0.455)
	{
		a=0.275;
		b=0.455;
		fa=0.4;
		fb=0.5;
	}
	else if(in>=0.455 && in<0.708)
	{
		a=0.455;
		b=0.708;
		fa=0.5;
		fb=0.6;
	}
	else if(in>=0.708 && in<1.074)
	{
		a=0.708;
		b=1.074;
		fa=0.6;
		fb=0.7;
	}
	else if(in>=1.074 && in<1.642)
	{
		a=1.074;
		b=1.642;
		fa=0.7;
		fb=0.8;
	}
	else if(in>=1.642 && in<2.706)
	{
		a=1.642;
		b=2.706;
		fa=0.8;
		fb=0.9;
	}
	else if(in>=2.706 && in<3.841)
	{
		a=2.706;
		b=3.841;
		fa=0.9;
		fb=0.95;
	}
	else if(in>=3.841 && in<5.024)
	{
		a=3.841;
		b=5.024;
		fa=0.95;
		fb=0.975;
	}
	else if(in>=5.024 && in<6.635)
	{
		a=5.024;
		b=6.635;
		fa=0.975;
		fb=0.99;
	}
	else if(in>=6.635 && in<7.879)
	{
		a=6.635;
		b=7.879;
		fa=0.99;
		fb=0.995;
	}
	else if(in>=7.879 && in<10.828)
	{
		a=7.879;
		b=10.828;
		fa=0.995;
		fb=0.999;
	}
	else if(in>=10.828 && in<12.116)
	{
		a=10.828;
		b=12.116;
		fa=0.999;
		fb=0.9995;
	}
	if (in>=12.116)
	{
		rough_value=0.9999;
	}
	else if(in<0.000000393)
	{
		rough_value=0.0001;
	}
	else
	{
		rough_value=(b-in)/(b-a)*fa + (in - a)/(b-a)*fb;
	}
	rough_value=rough_value+(1.0-rough_value)/2.0;
	rough_value=rough_value*rough_value*rough_value;
	return rough_value;
}

/*********************************************************/
/**
* deprecated
* Compute a RELL support, using the latest tested branch
* param tree: the tested tree
*/
m3ldbl Statistics_to_RELL(arbre *tree)
{
	int i;
	int occurence=1000;
	m3ldbl nb=0.0;
	m3ldbl res;
	int site;
	m3ldbl lk0=0.0;
	m3ldbl lk1=0.0;
	m3ldbl lk2=0.0;
	m3ldbl buff = -1.;
	int position = -1;

	//1000 times
	For(i,occurence)
	{
		lk0=0.0;
		lk1=0.0;
		lk2=0.0;
		//Shuffle the data and increment the support, if needed
		For(site, tree->data->init_len)
		{
			buff  = rand();
			buff /= (RAND_MAX+1.);
			/* 	  buff *= tree->data->crunch_len; */
			buff *= tree->data->init_len;
			position = (int)floor(buff);

			lk0+=tree->log_lks_aLRT[0][position];
			lk1+=tree->log_lks_aLRT[1][position];
			lk2+=tree->log_lks_aLRT[2][position];
		}
		if (lk0>=lk1 && lk0>=lk2) nb++;
	}

	res= nb/(m3ldbl)occurence;

	return res;
}
/*********************************************************/
/**
* deprecated
* Compute a SH-like support, using the latest tested branch
* param tree: the tested tree
*/
m3ldbl Statistics_To_SH(arbre *tree)
{
	int i;
	int occurence=1000;
	m3ldbl nb=0.0;
	m3ldbl res;
	int site;
	m3ldbl lk0=0.0;
	m3ldbl lk1=0.0;
	m3ldbl lk2=0.0;
	m3ldbl c0=0.0;
	m3ldbl c1=0.0;
	m3ldbl c2=0.0;
	m3ldbl buff=-1.;
	int position=-1;
	m3ldbl delta_local=-1.;


	//Compute the total log-lk of each NNI position
	For(site, tree->data->init_len)
	{
		c0+=tree->log_lks_aLRT[0][site];
		c1+=tree->log_lks_aLRT[1][site];
		c2+=tree->log_lks_aLRT[2][site];
	}
	m3ldbl delta=0.0;
	if (c0>=c1 && c0>=c2)
	{
		if (c1>=c2)
		{
			delta=c0-c1;
		}
		else
		{
			delta=c0-c2;
		}
	}
	else if(c1>=c0 && c1>=c2)
	{
		if (c0>=c2)
		{
			delta=c1-c0;
		}
		else
		{
			delta=c1-c2;
		}
	}
	else
	{
		if (c1>=c0)

		{
			delta=c2-c1;
		}
		else
		{
			delta=c2-c0;
		}
	}
	//1000 times
	For(i,occurence)
	{
		lk0=0.0;
		lk1=0.0;
		lk2=0.0;

		//Shuffle the data
		For(site, tree->data->init_len)
		{
			buff  = rand();
			buff /= (RAND_MAX+1.);
			buff *= tree->data->init_len;
			position = (int)floor(buff);

			lk0+=tree->log_lks_aLRT[0][position];
			lk1+=tree->log_lks_aLRT[1][position];
			lk2+=tree->log_lks_aLRT[2][position];
		}

		//return to null hypothesis
		lk0=lk0-c0;
		lk1=lk1-c1;
		lk2=lk2-c2;

		//compute results and increment if needed
		delta_local=0.0;
		if (lk0>=lk1 && lk0>=lk2)
		{
			if (lk1>=lk2)
			{
				delta_local=lk0-lk1;
			}
			else
			{
				delta_local=lk0-lk2;
			}
		}
		else if(lk1>=lk0 && lk1>=lk2)
		{
			if (lk0>=lk2)
			{
				delta_local=lk1-lk0;
			}
			else
			{
				delta_local=lk1-lk2;
			}
		}
		else
		{
			if (lk1>=lk0)

			{
				delta_local=lk2-lk1;
			}
			else
			{
				delta_local=lk2-lk0;
			}
		}

		if (delta>(delta_local+0.1)) {
			nb++;
		}
	}

	res= nb/occurence;

	return res;
}

/*********************************************************/
/**
* deprecated
* Compute one side likelihood
* param b_fcus : concerned edge
* param tree : b_fcus tree
* param exclude :  side to exclude for computation
*/
m3ldbl Update_Lk_At_Given_Edge_Excluding(edge *b_fcus, arbre *tree, node *exclude)
{
	if((!b_fcus->left->tax) && (exclude == NULL || exclude != b_fcus->left))
		Update_P_Lk(tree,b_fcus,b_fcus->left);
	if((!b_fcus->rght->tax) && (exclude == NULL || exclude != b_fcus->rght))
		Update_P_Lk(tree,b_fcus,b_fcus->rght);

	tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);

	return tree->c_lnL;
}


/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
