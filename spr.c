/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

 */

/*
 ** spr.c: Routines for performing SPR moves on the tree.
 **
 ** Wim Hordijk   Last modified: 28 August 2006
 ** Stephane Guindon 2007
 */

#include "utilities.h"
#include "models.h"
#include "lk.h"
#include "free.h"
#include "optimiz.h"
#include "spr.h"
#include "alrt.h"
#include "pars.h"
#include "simu.h"


/*
 ** BIG: Some big number.
 */

#define BIG  1e05

/*
 ** Global vars.
 **
 **   - cur_lk:        The current likelihood of the tree.
 **   - subtree_dist:  The average subtree distances matrix.
 **   - seq_dist:      The sequence distance matrix.
 **   - optim_cand:    Array for holding candidate moves for local and global branch
 **                    length optimization.
 **   - rgrft_cand:    Array for holding candidate regraft positions.
 **   - v_tmp:         The central node of the temporary regraft structure for
 **                    estimating changes in likelihood.
 **   - path:          The path through the tree during the recursive tree length
 **                    calculation.
 **   - sum_scale_tmp  Array for temporarily storing scaling factors.
 **   - p_lk_tmp:      Temporary partial likelihood storage.
 **   - e_brent:       A temporary edge to use for estimating distances using Brent.

 **   - tree->mod->s_opt->wim_n_rgrft:       Number of promising regraft positions to consider when
                      performing all improving SPR moves.
 **   - tree->mod->s_opt->wim_n_optim:       Number of candidate moves on which to perform local branch
                      length optimization.
 **   - tree->mod->s_opt->wim_max_dist:      Maximum regraft distance to consider.
 **   - tree->mod->s_opt->wim_n_globl:       Number of candidates moves on which to perform global branch
                      length optimization.
 **   - tree->mod->s_opt->wim_n_best:        Number of promising regraft positions to consider when
                      performing only the best SPR move.

 **   - nr_d_l:        Total number of change in tree length calculations done.
 **   - nr_d_lk:       Total number of change in likelihood calculations done.
 **   - nr_loc:        Total number of local branch length optimizations done.
 **   - nr_glb:        Total number of global branch length optimizations done.
 */

m3ldbl   cur_lk, **subtree_dist, *sum_scale_tmp, *p_lk_tmp;
matrix  *seq_dist;
_move_ **optim_cand, **rgrft_cand;
node    *v_tmp=NULL, **path;
edge    *e_brent=NULL;
int      nr_d_L, nr_d_lk, nr_loc, nr_glb;

/*
 ** Init_SPR: Initialize the SPR algorithm: allocate memory and set variables.
 **
 ** Parameters:
 **   - tree: The current tree to use for initialization.
 */

//JSJ: Modified Make_Node_Light() and Make_Edge_Light() so they take the n_l as an argument
void Init_SPR (arbre *tree)
{
	int   i, nr_nodes, nr_edges;
	node *u_0, *u_1, *u_2;

	/*
	 ** Get the SPR parameter values.
	 */
	nr_edges = 2*tree->n_otu-3;

	if(tree->mod->s_opt->wim_n_rgrft  < 0) tree->mod->s_opt->wim_n_rgrft = 1 + nr_edges / 5;
	if(tree->mod->s_opt->wim_n_globl  < 0) tree->mod->s_opt->wim_n_globl = 1 + nr_edges / 10;
	if(tree->mod->s_opt->wim_max_dist < 0) tree->mod->s_opt->wim_max_dist = 1 + nr_edges / 10;
	if(tree->mod->s_opt->wim_n_optim  < 0) tree->mod->s_opt->wim_n_optim = 100;
	if(tree->mod->s_opt->wim_n_best   < 0) tree->mod->s_opt->wim_n_best = tree->mod->s_opt->wim_n_rgrft; /* can't
	 *  be
	 *  anything else
	 */


	/*
	 ** If it doesn't exist yet, create the temporary regraft structure:
	 ** a central node with three edges and tip nodes adjacent to it.
	 */
	if (v_tmp == NULL)
	{

		v_tmp=Make_Node_Light(0, tree->n_l);
		v_tmp->tax = 0;
		u_0=Make_Node_Light(1, tree->n_l);
		u_0->tax = 1;
		u_1=Make_Node_Light(2, tree->n_l);
		u_1->tax = 1;
		u_2=Make_Node_Light(3, tree->n_l);
		u_2->tax = 1;

		v_tmp->v[0] = u_0;
		v_tmp->v[1] = u_1;
		v_tmp->v[2] = u_2;
		u_0->v[0]   = v_tmp;
		u_1->v[0]   = v_tmp;
		u_2->v[0]   = v_tmp;

		edge *edge_0 = Make_Edge_Light (v_tmp, u_0, 0, tree->n_l);
		Make_Edge_Lk (edge_0, tree);
		edge *edge_1 = Make_Edge_Light (v_tmp, u_1,1, tree->n_l);
		Make_Edge_Lk (edge_1, tree);
		edge *edge_2 = Make_Edge_Light (v_tmp, u_2,2, tree->n_l);
		Make_Edge_Lk (edge_2, tree);


		/*     For(i,tree->data->crunch_len) */
		/*       { */
		/* 	For(j,tree->mod->n_catg) */
		/* 	  { */
		/* 	    Free(edge_0->p_lk_rght[i][j]); */
		/* 	  } */
		/* 	Free(edge_0->p_lk_rght[i]); */
		/*       } */
		Free(edge_0->p_lk_rght);

		if(!edge_0->rght->tax) Free(edge_0->sum_scale_f_rght);


		/*     For(i,tree->data->crunch_len) */
		/*       { */
		/* 	For(j,tree->mod->n_catg) */
		/* 	  { */
		/* 	    Free(edge_1->p_lk_rght[i][j]); */
		/* 	  } */
		/* 	Free(edge_1->p_lk_rght[i]); */
		/*       } */
		Free(edge_1->p_lk_rght);

		if(!edge_1->rght->tax) Free(edge_1->sum_scale_f_rght);


		/*     For(i,tree->data->crunch_len) */
		/*       { */
		/* 	For(j,tree->mod->n_catg) */
		/* 	  { */
		/* 	    Free(edge_2->p_lk_rght[i][j]); */
		/* 	  } */
		/* 	Free(edge_2->p_lk_rght[i]); */
		/*       } */
		Free(edge_2->p_lk_rght);

		if(!edge_2->rght->tax) Free(edge_2->sum_scale_f_rght);
	}

	/*
	 ** If it doesn't exist yet, create the temporary edge.
	 */
	if (e_brent == NULL)
	{

		u_1=Make_Node_Light(4, tree->n_l);
		u_1->tax = 1;
		u_2=Make_Node_Light(5, tree->n_l);
		u_2->tax = 1;
		u_1->v[0] = u_2;
		u_2->v[0] = u_1;
		edge *edge_4 = Make_Edge_Light (u_1, u_2, 3, tree->n_l);
		Make_Edge_Lk (edge_4, tree);
		e_brent = u_1->b[0];


		/*     For(i,tree->data->crunch_len) */
		/*       { */
		/* 	For(j,tree->mod->n_catg) */
		/* 	  { */
		/* 	    Free(edge_4->p_lk_rght[i][j]); */
		/* 	  } */
		/* 	Free(edge_4->p_lk_rght[i]); */
		/*       } */
		Free(edge_4->p_lk_rght);

		if(!edge_4->rght->tax) Free(edge_4->sum_scale_f_rght);


		/*     For(i,tree->data->crunch_len) */
		/*       { */
		/* 	For(j,tree->mod->n_catg) */
		/* 	  { */
		/* 	    Free(edge_4->p_lk_left[i][j]); */
		/* 	  } */
		/* 	Free(edge_4->p_lk_left[i]); */
		/*       } */
		Free(edge_4->p_lk_left);

		if(!edge_4->left->tax) Free(edge_4->sum_scale_f_left);
	}

	/*
	 ** Allocate memory for temporarily storing partial likelihoods and
	 ** scaling factors.
	 */
	p_lk_tmp = (m3ldbl *)mCalloc (tree->n_pattern*tree->mod->n_catg*tree->mod->ns, sizeof (m3ldbl));
	/*   p_lk_tmp = (m3ldbl ***)mCalloc (tree->n_pattern, sizeof (m3ldbl **)); */
	/*   for (i = 0; i < tree->n_pattern; i++) */
	/*   { */
	/*     p_lk_tmp[i] = (m3ldbl **)mCalloc (tree->mod->n_catg, sizeof (m3ldbl *)); */
	/*     for (j = 0; j < tree->mod->n_catg; j++) */
	/*     { */
	/*       p_lk_tmp[i][j] = (m3ldbl *)mCalloc (tree->mod->ns, sizeof (m3ldbl)); */
	/*     } */
	/*   } */
	sum_scale_tmp = (m3ldbl *)mCalloc (tree->n_pattern, sizeof (m3ldbl));

	/*
	 ** Allocate memory for storing the average subtree distances.
	 */
	nr_nodes = 2*tree->n_otu-2;
	subtree_dist = (m3ldbl **)malloc (nr_nodes * sizeof (m3ldbl *));
	for (i = 0; i < nr_nodes; i++)
	{
		subtree_dist[i] = (m3ldbl *)malloc (nr_nodes * sizeof (m3ldbl));
	}

	/*
	 ** Allocate memory for storing the candidate regraft positions and
	 ** edge length optimization moves.
	 */
	rgrft_cand = (_move_ **)malloc (MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best) * sizeof (_move_ *));
	for (i = 0; i < MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best); i++)
	{
		rgrft_cand[i] = (_move_ *)malloc (sizeof (_move_));
		rgrft_cand[i]->path = (node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (node *));
	}
	optim_cand = (_move_ **)malloc (tree->mod->s_opt->wim_n_optim * sizeof (_move_ *));
	for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
	{
		optim_cand[i] = (_move_ *)malloc (sizeof (_move_));
		optim_cand[i]->path = (node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (node *));
	}
	path = (node **)malloc ((tree->mod->s_opt->wim_max_dist+2) * sizeof (node *));

	if(!tree->mat)
	{
		seq_dist = ML_Dist (tree->data, tree->mod);
		tree->mat = seq_dist;
	}
	else
		seq_dist = tree->mat;

	/*
	 ** Set variables.
	 */
	nr_d_L = 0;
	nr_d_lk = 0;
	nr_loc = 0;
	nr_glb = 0;
}


/*
 ** Clean_SPR: Free up the used memory.
 **
 ** Parameters:
 **   - tree: The current tree.
 */

void Clean_SPR (arbre *tree)
{
	int i;

	/*
	 ** Clean up the temporary regraft structure.
	 */
	Free_Node (v_tmp->v[0]);
	Free_Node (v_tmp->v[1]);
	Free_Node (v_tmp->v[2]);
	v_tmp->b[0]->p_lk_rght = NULL;
	Free_Edge_Lk (tree, v_tmp->b[0]);
	Free_Edge (v_tmp->b[0]);
	v_tmp->b[1]->p_lk_rght = NULL;
	Free_Edge_Lk(tree, v_tmp->b[1]);
	Free_Edge(v_tmp->b[1]);
	v_tmp->b[2]->p_lk_rght = NULL;
	Free_Edge_Lk (tree, v_tmp->b[2]);
	Free_Edge (v_tmp->b[2]);
	Free_Node (v_tmp);
	v_tmp = NULL;

	/*
	 ** Clean up the temporary edge.
	 */
	Free_Node (e_brent->left);
	Free_Node (e_brent->rght);
	e_brent->p_lk_left = NULL;
	e_brent->p_lk_rght = NULL;
	Free_Edge_Lk (tree, e_brent);
	Free_Edge (e_brent);
	e_brent = NULL;

	/*
	 ** Free the temporary partial likelihood and scaling memory.
	 */
	/*   for (i = 0; i < tree->n_pattern; i++) */
	/*   { */
	/*     for (j = 0; j < tree->mod->n_catg; j++) */
	/*     { */
	/*       free (p_lk_tmp[i][j]); */
	/*     } */
	/*     free (p_lk_tmp[i]); */
	/*   } */
	free (p_lk_tmp);
	free (sum_scale_tmp);

	/*
	 ** Free the subtree distance matrix.
	 */
	for (i = 0; i < 2*tree->n_otu - 2; i++)
	{
		free (subtree_dist[i]);
	}
	free (subtree_dist);

	/*
	 ** Free the arrays for storing the candidate regrafting positions and
	 ** edge length optimization moves.
	 */
	for (i = 0; i < MAX(tree->mod->s_opt->wim_n_rgrft,tree->mod->s_opt->wim_n_best); i++)
	{
		free (rgrft_cand[i]->path);
		free (rgrft_cand[i]);
	}
	free (rgrft_cand);
	for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
	{
		free (optim_cand[i]->path);
		free (optim_cand[i]);
	}
	free (optim_cand);
	free (path);


	/*
	 ** Print some statistics (for "research" purposes only).
	 */
	/*   PhyML_Printf ("nr_d_L:  %d\n", nr_d_L); */
	/*   PhyML_Printf ("nr_d_lk: %d\n", nr_d_lk); */
	/*   PhyML_Printf ("nr_loc:  %d\n", nr_loc); */
	/*   PhyML_Printf ("nr_glb:  %d\n", nr_glb); */
}


/*
 ** Optim_SPR: Optimize the tree using SPR moves.
 **
 ** Parameters:
 **   - tree:     The tree to optimize.
 **   - max_size: The maximum size (= number of taxa) of the subtrees to be
 **               pruned. If m=0 or m>ntax, all possible prunings will be
 **               considered.
 **   - method:   The optimization method to use ("ALL" or "BEST").
 */

void Optim_SPR (arbre *tree, int max_size, int method)
{
	int   nr_moves, improvement;
	node *root;

	if(tree->mod->s_opt->print) PhyML_Printf("\n\n. Starting SPR moves...\n");

	/*
	 ** Calculate the current likelihood value.
	 */
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	time(&(tree->t_current));
	if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

	/*
	 ** Optimize all edge lengths and calculate the new likelihood value.
	 */
	/*   PhyML_Printf("\n. Optimizing edge lengths."); */
	root = tree->noeud[0];
	Optimize_Br_Len_Serie (root, root->v[0], root->b[0], tree, tree->data);
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	time(&(tree->t_current));
	if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

	/*
	 ** While improvements were found, perform another round of SPR moves.
	 */
	nr_moves = 0;
	improvement = 1;
	while (improvement)
	{
		/*
		 ** Perform one round of SPR moves.
		 */
		if (method == ALL)
		{
			improvement = Perform_SPR_Moves (tree, max_size);
		}
		else if (method == BEST)
		{
			improvement = Perform_Best_SPR (tree, max_size);
		}
		else if (method == ONE)
		{
			improvement = Perform_One_SPR (tree, max_size);
		}
		else
		{
			PhyML_Printf ("\n. Unknown SPR optimization method, bailing out...\n");
			exit (1);
		}

		/*     If an improvement was found, update statistics. */
		if(improvement)
		{
			nr_moves++;
			if((nr_moves == 1) || (nr_moves % 4 == 1))
			{
				/*
				 ** Optimize model parameters.
				 */
				Optimiz_All_Free_Param (tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
				tree->both_sides = 1;
				Lk(tree);
			}
		}

		/* Beg SG 28 May 2007 */
		if(method == BEST || method == ONE) break;
		/* Beg SG 28 May 2007 */
	}

	if(tree->mod->s_opt->print) PhyML_Printf ("\n\n. Number of SPR moves: %d\n", nr_moves);

	/*
	 ** Perform a last round of optimization steps (for edge lengths).
	 */
	Round_Optimize(tree,tree->data,ROUND_MAX);
	Check_NNI_Five_Branches(tree);
}


/*
 ** Perform_SPR_Moves: Perform a round of SPR moves on the tree. Prune each subtree in
 **                    turn and calculate the change in tree length for each candidate
 **                    regraft position. Estimate change in likelihood for the most
 **                    promising moves, and perform all moves that result in an
 **                    improvement. If no improvements were found at all, try local edge
 **                    length optimization. If still no improvement, try global edge
 **                    length optimization.
 **
 ** Parameters:
 **   - tree:     The tree to perform the SPR moves on.
 **   - max_size: The maximum size (= number of taxa) of the subtrees to be
 **               pruned. If m=0 or m>ntax, all possible prunings will be
 **               considered.
 **
 ** Returns:
 **   If the current tree could be improved: 1.
 **   Otherwise:                             0.
 */

int Perform_SPR_Moves (arbre *tree, int max_size)
{
	int   nr_edges, i, j, candidate, improvement;
	node *root, *v_prune;
	edge *e_prune;

	/*
	 ** Calculate the average subtree distances.
	 */
	root = tree->noeud[0];
	PostOrder_v (tree, root->v[0], root->b[0]);

	/*
	 ** Initialize the array of optimization candidates.
	 */
	for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
	{
		optim_cand[i]->delta_lk = -1.0*BIG;
		optim_cand[i]->d_L = -1.0*BIG;
	}

	/*
	 ** Try all possible SPR moves and perform the ones that give an improvement.
	 */
	nr_edges = 2*tree->n_otu - 3;
	cur_lk = tree->c_lnL;
	improvement = 0;


	/*   PhyML_Printf("\n >>>>>>>>>>>>>>>>>>"); */
	/*   PhyML_Printf("\n. cur_lk = %f %f",cur_lk,Return_Lk(tree)); */

	/*   PhyML_Printf ("\n. Trying SPR moves"); */
	/*   PhyML_Printf ("\n.  - calculating tree distances and estimating likelihoods"); */

	for(i = 0; i < nr_edges; i++)
	{
		/*
		 ** Get the next prune edge.
		 */
		e_prune = tree->t_edges[i];
		/*
		 ** Try right subtree if appropriate.
		 */
		if (!e_prune->left->tax)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->left;
			if ((max_size == 0) || (e_prune->num_tax_rght <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Perform moves that give an improvement.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
				{
					improvement = 1;
					Make_Move (rgrft_cand[candidate],0,tree);
					/* 		  PhyML_Printf("\n. Make simple move"); */
					/* 		  PhyML_Printf("\n. lk after simple move = %f",Return_Lk(tree)); */
				}
			}
		}
		/*
		 ** Try left subtree if appropriate.
		 */
		if (!e_prune->rght->tax)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->rght;
			if ((max_size == 0) || (e_prune->num_tax_left <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Perform moves that give an improvement.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
				{
					improvement = 1;
					Make_Move (rgrft_cand[candidate],0,tree);
					/* 		  PhyML_Printf("\n. Make simple move"); */
					/* 		  PhyML_Printf("\n. lk after simple move = %f",Return_Lk(tree)); */
				}
			}
		}
	}

	/*
	 ** If there was no improvement at all, try local edge length optimization at the
	 ** regraft position.
	 */

	/*   PhyML_Printf("\n. before local = %f %f",tree->c_lnL,Return_Lk(tree)); */

	if (!improvement)
	{
		/*       PhyML_Printf ("\n.  - performing local edge length optimizations"); */
		if ((candidate = Find_Optim_Local (tree)) >= 0)
		{
			/*  	  PhyML_Printf("\n. make local move"); */
			improvement = 1;
			Make_Move (optim_cand[candidate],1,tree);
			/* 	  PhyML_Printf("\n. lk after local move = %f",Return_Lk(tree)); */
		}
	}


	/*
	 ** If there was still no improvement, try global edge length optimization.
	 */
	/*   PhyML_Printf("\n. before global = %f %f",tree->c_lnL,Return_Lk(tree)); */

	if (!improvement)
	{
		/*       PhyML_Printf ("\n.  - performing global edge length optimization"); */
		if ((candidate = Find_Optim_Globl (tree)) >= 0)
		{
			/*  	  PhyML_Printf("\n. make global move"); */
			improvement = 1;
			Make_Move (optim_cand[candidate],2,tree);
			/* 	  PhyML_Printf("\n. lk after global move = %f",Return_Lk(tree)); */
		}
	}

	/*   PhyML_Printf("\n. after all = %f %f",tree->c_lnL,Return_Lk(tree)); */

	/*
	 ** Optimize all edge lengths again to make sure we got an updated
	 ** likelihood value.
	 */
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	root = tree->noeud[0];
	Optimize_Br_Len_Serie (root, root->v[0], root->b[0], tree, tree->data);
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	time(&(tree->t_current));
	if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

	/*
	 ** Return the result.
	 */
	return (improvement);
}


/*
 ** Perform_Best_SPR: Perform the best SPR move on the tree. Prune each subtree in
 **                   turn and calculate the change in tree length for each candidate
 **                   regraft position. Estimate change in likelihood for the most
 **                   promising regraft positions, and store the best one. Then choose
 **                   the best candidate over all moves. If no improving move can be
 **                   found, try local edge length optimization, and if necessary
 **                   global edge length optimization.
 **
 ** Parameters:
 **   - tree:     The tree to perform the SPR moves on.
 **   - max_size: The maximum size (= number of taxa) of the subtrees to be
 **               pruned. If m=0 or m>ntax, all possible prunings will be
 **               considered.
 **
 ** Returns:
 **   If an improving move could be performed: 1.
 **   Otherwise:                               0.
 */

int Perform_Best_SPR (arbre *tree, int max_size)
{
	int   nr_edges, i, j, candidate, improvement;
	node *root, *v_prune;
	edge *e_prune;

	/*
	 ** Calculate the average subtree distances.
	 */
	root = tree->noeud[0];
	PostOrder_v (tree, root->v[0], root->b[0]);

	/*
	 ** Initialize the array of optimization candidates.
	 */
	for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
	{
		optim_cand[i]->delta_lk = -1.0*BIG;
		optim_cand[i]->d_L = -1.0*BIG;
	}

	/*
	 ** Try all possible SPR moves and perform the best one.
	 */
	nr_edges = 2*tree->n_otu - 3;
	cur_lk = tree->c_lnL;
	improvement = 0;
	/*   PhyML_Printf ("\n. Trying SPR moves"); */
	/*   PhyML_Printf ("\n.  -calculating tree distances and estimating likelihoods"); */
	for (i = 0; i < nr_edges; i++)
	{
		/*
		 ** Get the next prune edge.
		 */
		e_prune = tree->t_edges[i];
		/*
		 ** Try right subtree if appropriate.
		 */
		if (!e_prune->left->tax)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_best; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->left;
			if ((max_size == 0) || (e_prune->num_tax_rght <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Store the best one in the optimization list.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				candidate = Best_Lk_Change (e_prune, v_prune, tree);
			}
		}
		/*
		 ** Try left subtree if appropriate.
		 */
		if (!e_prune->rght->tax)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->rght;
			if ((max_size == 0) || (e_prune->num_tax_left <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Perform moves that give an improvement.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				candidate = Best_Lk_Change (e_prune, v_prune, tree);
			}
		}
	}

	/* If the best candidate has a positive estimated change in
	 ** likelihood, perform that move.
	 */
	if (optim_cand[0]->delta_lk > 1.0/BIG)
	{
		improvement = 1;
		Make_Move (optim_cand[0],0,tree);
	}

	/*
	 ** If there was no improvement at all, try local edge length optimization at the
	 ** regraft position.
	 */
	if (!improvement)
	{
		/*     PhyML_Printf ("\n.  - performing local edge length optimizations"); */
		if ((candidate = Find_Optim_Local (tree)) >= 0)
		{
			improvement = 1;
			Make_Move (optim_cand[candidate],1,tree);
		}
	}

	/*
	 ** If there was still no improvement, try global edge length optimization.
	 */
	if (!improvement)
	{
		/*     PhyML_Printf ("\n.  - performing global edge length optimization"); */
		if ((candidate = Find_Optim_Globl (tree)) >= 0)
		{
			improvement = 1;
			Make_Move (optim_cand[candidate],2,tree);
		}
	}

	/*
	 ** Optimize all edge lengths again to make sure we got an updated
	 ** likelihood value.
	 */
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	root = tree->noeud[0];
	Optimize_Br_Len_Serie (root, root->v[0], root->b[0], tree, tree->data);
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	time(&(tree->t_current));
	if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

	/*
	 ** Return the result.
	 */
	return (improvement);
}


/*
 ** Perform_One_Moves: Perform a round of SPR moves on the tree. Prune each subtree in
 **                    turn and calculate the change in tree length for each candidate
 **                    regraft position. Estimate change in likelihood for the most
 **                    promising moves, and perform the first move that results in an
 **                    improvement. If no improvements were found at all, try local edge
 **                    length optimization. If still no improvement, try global edge
 **                    length optimization.
 **
 ** Parameters:
 **   - tree:     The tree to perform the SPR moves on.
 **   - max_size: The maximum size (= number of taxa) of the subtrees to be
 **               pruned. If m=0 or m>ntax, all possible prunings will be
 **               considered.
 **
 ** Returns:
 **   If the current tree could be improved: 1.
 **   Otherwise:                             0.
 */

int Perform_One_SPR(arbre *tree, int max_size)
{
	int   nr_edges, i, j, candidate, improvement;
	node *root, *v_prune;
	edge *e_prune;

	/*
	 ** Calculate the average subtree distances.
	 */

	root = tree->noeud[0];
	PostOrder_v (tree, root->v[0], root->b[0]);

	/*
	 ** Initialize the array of optimization candidates.
	 */
	for (i = 0; i < tree->mod->s_opt->wim_n_optim; i++)
	{
		optim_cand[i]->delta_lk = -1.0*BIG;
		optim_cand[i]->d_L = -1.0*BIG;
	}

	/*
	 ** Try all possible SPR moves and perform the ones that give an improvement.
	 */
	nr_edges = 2*tree->n_otu - 3;
	cur_lk = tree->c_lnL;
	improvement = 0;

	/*   PhyML_Printf("\n >>>>>>>>>>>>>>>>>>"); */
	/*   PhyML_Printf("\n. cur_lk = %f %f",cur_lk,Return_Lk(tree)); */

	/*   PhyML_Printf ("\n. Trying SPR moves"); */
	/*   PhyML_Printf ("\n.  - calculating tree distances and estimating likelihoods"); */

	for(i = 0; i < nr_edges; i++)
	{
		/*
		 ** Get the next prune edge.
		 */
		e_prune = tree->t_edges[i];
		/*
		 ** Try right subtree if appropriate.
		 */
		if (!e_prune->left->tax)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->left;
			if ((max_size == 0) || (e_prune->num_tax_rght <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Perform moves that give an improvement.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
				{
					improvement = 1;
					Make_Move (rgrft_cand[candidate],0,tree);
				}
			}
		}
		/*
		 ** Try left subtree if appropriate.
		 */
		if (!e_prune->rght->tax && !improvement)
		{
			/*
			 ** Clear the regraft candidate list.
			 */
			for (j = 0; j < tree->mod->s_opt->wim_n_rgrft; j++)
			{
				rgrft_cand[j]->d_L = -1.0*BIG;
			}
			v_prune = e_prune->rght;
			if ((max_size == 0) || (e_prune->num_tax_left <= max_size))
			{
				/*
				 ** Calculate changes in tree length, and estimate changes in likelihood for
				 ** the most promising candidates. Perform moves that give an improvement.
				 */
				Calc_Tree_Length (e_prune, v_prune, tree);
				if ((candidate = Est_Lk_Change (e_prune, v_prune, tree)) >= 0)
				{
					improvement = 1;
					Make_Move (rgrft_cand[candidate],0,tree);
					/* 		  PhyML_Printf("\n. Make simple move"); */
					/* 		  PhyML_Printf("\n. lk after simple move = %f",Return_Lk(tree)); */
				}
			}
		}
		if(improvement) break;
	}

	/*
	 ** If there was no improvement at all, try local edge length optimization at the
	 ** regraft position.
	 */

	/*   PhyML_Printf("\n. before local = %f %f",tree->c_lnL,Return_Lk(tree)); */

	if (!improvement)
	{
		/*       PhyML_Printf ("\n.  - performing local edge length optimizations"); */
		if ((candidate = Find_Optim_Local (tree)) >= 0)
		{
			/*  	  PhyML_Printf("\n. make local move"); */
			improvement = 1;
			Make_Move (optim_cand[candidate],1,tree);
			/* 	  PhyML_Printf("\n. lk after local move = %f",Return_Lk(tree)); */
		}
	}


	/*
	 ** If there was still no improvement, try global edge length optimization.
	 */
	/*   PhyML_Printf("\n. before global = %f %f",tree->c_lnL,Return_Lk(tree)); */

	if (!improvement)
	{
		/*       PhyML_Printf ("\n.  - performing global edge length optimization"); */
		if ((candidate = Find_Optim_Globl (tree)) >= 0)
		{
			/*  	  PhyML_Printf("\n. make global move"); */
			improvement = 1;
			Make_Move (optim_cand[candidate],2,tree);
			/* 	  PhyML_Printf("\n. lk after global move = %f",Return_Lk(tree)); */
		}
	}

	/*   PhyML_Printf("\n. after all = %f %f",tree->c_lnL,Return_Lk(tree)); */

	/*
	 ** Optimize all edge lengths again to make sure we got an updated
	 ** likelihood value.
	 */
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	root = tree->noeud[0];
	Optimize_Br_Len_Serie (root, root->v[0], root->b[0], tree, tree->data);
	tree->both_sides = 1;
	cur_lk = Return_Lk (tree);
	time(&(tree->t_current));
	if(tree->mod->s_opt->print) Print_Lk(tree,"topology");

	/*
	 ** Return the result.
	 */
	return (improvement);
}


/*
 **  Calc_Tree_Length: Calculate the change in tree length, given a pruned subtree,
 **                    for each possible regraft position.
 **
 ** Parameters:
 **   - e_prune: The edge at which the subtree is pruned.
 **   - v_prune: The root of the pruned subtree.
 */

void Calc_Tree_Length (edge *e_prune, node *v_prune, arbre *tree)
{
	int     i, d0, d1, d2;
	m3ldbl  d_uu;
	node   *u_prune, *u1, *u2;

	/*
	 ** Get the directions from node v_prune.
	 */
	d0 = -1;
	u_prune = NULL;
	for (i = 0; i < 3; i++)
	{
		if (v_prune->b[i] == e_prune)
		{
			d0 = i;
			u_prune = v_prune->v[i];
			break;
		}
	}
	d1 = (d0 + 1) % 3;
	d2 = 3 - d0 - d1;

	/*
	 ** Get the relevant average subtree distance within the pruned subtree.
	 */
	if (!u_prune->tax)
	{
		u1 = u2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (u_prune->b[i] != e_prune)
			{
				if (u1 == NULL)
				{
					u1 = u_prune->v[i];
				}
				else
				{
					u2 = u_prune->v[i];
				}
			}
		}
		d_uu = subtree_dist[u1->num][u2->num];
	}
	else
	{
		d_uu = 0.0;
	}

	/*
	 ** Recursively calculate the change in tree length for each
	 ** possible regraft position.
	 **
	 ** First recurse into direction d1.
	 */
	if (!v_prune->v[d1]->tax)
	{
		u1 = u2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (v_prune->v[d1]->b[i] != v_prune->b[d1])
			{
				if (u1 == NULL)
				{
					u1 = v_prune->v[d1]->v[i];
				}
				else
				{
					u2 = v_prune->v[d1]->v[i];
				}
			}
		}
		Tree_Length(v_prune, u_prune, v_prune->v[d1], v_prune->v[d2], u1, v_prune->v[d2],
				u2, subtree_dist[u_prune->num][v_prune->v[d2]->num], d_uu, 0.0, 1, tree);
		Tree_Length(v_prune, u_prune, v_prune->v[d1], v_prune->v[d2], u2, v_prune->v[d2],
				u1, subtree_dist[u_prune->num][v_prune->v[d2]->num], d_uu, 0.0, 1, tree);
	}
	/*
	 ** Next recurse into direction d2.
	 */
	if (!v_prune->v[d2]->tax)
	{
		u1 = u2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (v_prune->v[d2]->b[i] != v_prune->b[d2])
			{
				if (u1 == NULL)
				{
					u1 = v_prune->v[d2]->v[i];
				}
				else
				{
					u2 = v_prune->v[d2]->v[i];
				}
			}
		}
		Tree_Length(v_prune, u_prune, v_prune->v[d2], v_prune->v[d1], u1, v_prune->v[d1],
				u2, subtree_dist[u_prune->num][v_prune->v[d1]->num], d_uu, 0.0, 1, tree);
		Tree_Length(v_prune, u_prune, v_prune->v[d2], v_prune->v[d1], u2, v_prune->v[d1],
				u1, subtree_dist[u_prune->num][v_prune->v[d1]->num], d_uu, 0.0, 1, tree);
	}
}


/*
 ** Tree_Length: Recursively calculate the change in tree length for a given pruned
 **              subtree and regraft position.
 **
 ** Parameters:
 **   - v_prune: The root of the pruned subtree.
 **   - u_prune: The node adjacent to v_p along the pruned edge.
 **   - v_n:     The node adjacent to the regraft edge in the "backward" direction.
 **   - v_n_1:   The previous v_n.
 **   - v_nx1:   The node adjacent to the regrafting edge in the "forward" direction.
 **   - v_0:     The other node originally adjacent to v_p;
 **   - u_n:     The other node adjecent to v_n (besides v_n_1 and v_nx1);
 **   - d_uv_1:  The distance between u_p and v_n_1;
 **   - d_uu:    The subtree distance between descendants of u_prune.
 **   - d_L_1:   The previous change in tree length.
 **   - n:       The current distance from the prune position.
 */

void Tree_Length (node *v_prune, node *u_prune, node *v_n, node *v_n_1, node *v_nx1,
		node *v_0, node *u_n, m3ldbl d_up_v_1, m3ldbl d_L_1, m3ldbl d_uu,
		int n, arbre *tree)
{
	int     i, j;
	m3ldbl  d_un_v, d_up_v, d_L;
	node   *u1, *u2;
	edge   *e_prune, *e_regraft;
	_move_ *tmp_cand;

	/*
	 ** Update the path and number of calculations.
	 */
	path[n] = v_n;
	nr_d_L++;
	e_prune = NULL;
	e_regraft = NULL;

	/*
	 ** Calculate the change in tree length for the current pruned subtree and regraft
	 ** position.
	 */
	if (n == 1)
	{
		d_un_v = subtree_dist[u_n->num][v_0->num];
	}
	else
	{
		d_un_v = subtree_dist[u_n->num][v_n_1->num] -
				(pow (0.5, n) * subtree_dist[u_n->num][u_prune->num]) +
				(pow (0.5, n) * subtree_dist[u_n->num][v_0->num]);
	}
	d_up_v = 0.5 * (d_up_v_1 + subtree_dist[u_prune->num][u_n->num]);
	/*
	 ** Alternative method for calculating d_up_v. Just kept it around for reference...
	 **
  d_up_v = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
  if (!u_n->tax)
  {
    u1 = u2 = NULL;
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] != v_n)
      {
	if (u1 == NULL)
	{
	  u1 = u_n->v[i];
	}
	else
	{
	  u2 = u_n->v[i];
	}
      }
    }
    d_up_v -= 0.5 * subtree_dist[u1->num][u2->num];
  }
  for (i = 0; i < 3; i++)
  {
    if (u_n->v[i] == v_n)
    {
      d_up_v -= u_n->b[i]->l;
      break;
    }
  }
	 */
	d_L = d_L_1 + 0.25*((d_up_v_1 + subtree_dist[u_n->num][v_nx1->num]) -
			(d_un_v + subtree_dist[u_prune->num][v_nx1->num]));

	/*
	 ** If the change is within the tree->mod->s_opt->wim_n_rgrft best ones so far, save it.
	 */
	if (d_L > rgrft_cand[tree->mod->s_opt->wim_n_rgrft-1]->d_L)
	{
		for (i = 0; i < 3; i++)
		{
			if (v_prune->v[i] == u_prune)
			{
				e_prune = v_prune->b[i];
			}
			if (v_n->v[i] == v_nx1)
			{
				e_regraft = v_n->b[i];
			}
		}
		i = tree->mod->s_opt->wim_n_rgrft-1;
		rgrft_cand[i]->v_prune = v_prune;
		rgrft_cand[i]->u_prune = u_prune;
		rgrft_cand[i]->v_n = v_n;
		rgrft_cand[i]->v_nx1 = v_nx1;
		rgrft_cand[i]->u_n = u_n;
		rgrft_cand[i]->e_prune = e_prune;
		rgrft_cand[i]->e_regraft = e_regraft;
		rgrft_cand[i]->d_L = d_L;
		rgrft_cand[i]->d_up_v = d_up_v;
		rgrft_cand[i]->d_un_v = d_un_v;
		rgrft_cand[i]->dist = n;
		for (j = 1; j <= n; j++)
		{
			rgrft_cand[i]->path[j] = path[j];
		}

		rgrft_cand[i]->path[n+1] = v_nx1;
		/*
		 ** Move the candidate to the appropriate position in the list, so the list
		 ** remains sorted in decreasing d_L value.
		 */
		while ((i > 0) && (rgrft_cand[i]->d_L > rgrft_cand[i-1]->d_L))
		{
			tmp_cand = rgrft_cand[i];
			rgrft_cand[i] = rgrft_cand[i-1];
			rgrft_cand[i-1] = tmp_cand;
			i--;
		}
	}

	/*
	 ** Recurse.
	 */
	if (n < tree->mod->s_opt->wim_max_dist)
	{
		if (!v_nx1->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_nx1->v[i] != v_n)
				{
					if (u1 == NULL)
					{
						u1 = v_nx1->v[i];
					}
					else
					{
						u2 = v_nx1->v[i];
					}
				}
			}
			Tree_Length (v_prune, u_prune, v_nx1, v_n, u1, v_0, u2, d_up_v, d_uu, d_L, n+1, tree);
			Tree_Length (v_prune, u_prune, v_nx1, v_n, u2, v_0, u1, d_up_v, d_uu, d_L, n+1, tree);
		}
	}
}


/*
 ** Est_Lk_Change: Estimate the changes in likelihood for the most promising candidate
 **                regraft positions given a pruned subtree.
 **
 ** Parameters:
 **   - e_prune: The edge at which the subtree was pruned.
 **   - v_prune: The root of the pruned subtree.
 **   - tree:    The tree on which to do the calculations.
 **
 ** Returns:
 **   If an improvement as found: The candidate which gives the improvement (which
 **                               will be the first one found).
 **   Otherwise:                  -1.
 */

int Est_Lk_Change (edge *e_prune, node *v_prune, arbre *tree)
{
	int     i, j, k, cand, best_cand, d0, d1, d2, n, pat, cat, ste;
	m3ldbl  d_uu[MAX_BL_SET], best_d_lk, l_connect[MAX_BL_SET], l_01[MAX_BL_SET], l_02[MAX_BL_SET], l_12[MAX_BL_SET], l_est[MAX_BL_SET][3], new_lk,
	l_simple[MAX_BL_SET][3], l_dist[MAX_BL_SET][3];
	plkflt *p_lk1_tmp, *p_lk2_tmp, *p_lk, *p_sum;
	node   *u_prune, *v_n, *v_nx1, *u_n, *u1, *u2;
	edge   *e_regraft, *e_tmp;
	_move_ *tmp_cand;
	int dim1, dim2;


	dim1 = tree->mod->ns * tree->mod->n_catg;
	dim2 = tree->mod->n_catg;

	/*
	 ** Get the directions from node v_prune.
	 */
	d0 = -1;
	u_prune = NULL;
	for (i = 0; i < 3; i++)
	{
		if (v_prune->b[i] == e_prune)
		{
			d0 = i;
			u_prune = v_prune->v[i];
			break;
		}
	}
	d1 = (d0 + 1) % 3;
	d2 = 3 - d0 - d1;

	/*
	 ** Copy the relevant partial likelihoods to the temporary regraft structure.
	 ** We can point to the original matrices, cos they won't be changed anyway.
	 */
	if (v_prune == e_prune->left)
	{
		v_tmp->b[0]->p_lk_rght = e_prune->p_lk_rght;
		v_tmp->b[0]->sum_scale_f_rght = e_prune->sum_scale_f_rght;
	}
	else
	{
		v_tmp->b[0]->p_lk_rght = e_prune->p_lk_left;
		v_tmp->b[0]->sum_scale_f_rght = e_prune->sum_scale_f_left;
	}
	v_tmp->num = v_prune->num;
	v_tmp->v[0]->num = u_prune->num;
	v_tmp->b[0]->num = e_prune->num;

	/*
	 ** Estimate the length of the edge that will connect the two "detached" nodes
	 ** after pruning. (The average of the sum of the lengths of the original two
	 ** edges and the average subtree distance based estimate.)
	 */
	For(k,tree->n_l){
		l_connect[k] = subtree_dist[v_prune->v[d1]->num][v_prune->v[d2]->num];
		if (!v_prune->v[d1]->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_prune->v[d1]->b[i] != v_prune->b[d1])
				{
					if (u1 == NULL)
					{
						u1 = v_prune->v[d1]->v[i];
					}
					else
					{
						u2 = v_prune->v[d1]->v[i];
					}
				}
			}
			l_connect[k] -= 0.5 * subtree_dist[u1->num][u2->num];
		}
		if (!v_prune->v[d2]->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_prune->v[d2]->b[i] != v_prune->b[d2])
				{
					if (u1 == NULL)
					{
						u1 = v_prune->v[d2]->v[i];
					}
					else
					{
						u2 = v_prune->v[d2]->v[i];
					}
				}
			}
			l_connect[k] -= 0.5 * subtree_dist[u1->num][u2->num];
		} //JSJ: temp fixes to l
		l_connect[k] += (v_prune->b[d1]->l[k] + v_prune->b[d2]->l[k]);
		l_connect[k] /= 2.0;
		if (l_connect[k] < BL_MIN)
		{
			l_connect[k] = BL_MIN;
		}
	}

	/*
	 ** Temporarily swap the relevant partial likelihoods at the prune site.
	 **
	 ** Direction d1.
	 */
	if (v_prune == v_prune->b[d1]->left)
	{
		p_lk1_tmp = v_prune->b[d1]->p_lk_left;
		if (v_prune == v_prune->b[d2]->left)
		{
			v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_rght;
		}
		else
		{
			v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_left;
		}
	}
	else
	{
		p_lk1_tmp = v_prune->b[d1]->p_lk_rght;
		if (v_prune == v_prune->b[d2]->left)
		{
			v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_rght;
		}
		else
		{
			v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_left;
		}
	}
	/*
	 ** Direction d2.
	 */
	if (v_prune == v_prune->b[d2]->left)
	{
		p_lk2_tmp = v_prune->b[d2]->p_lk_left;
		if (v_prune == v_prune->b[d1]->left)
		{
			v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_rght;
		}
		else
		{
			v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_left;
		}
	}
	else
	{
		p_lk2_tmp = v_prune->b[d2]->p_lk_rght;
		if (v_prune == v_prune->b[d1]->left)
		{
			v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_rght;
		}
		else
		{
			v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_left;
		}
	}

	/*
	 ** Temporarily set the edge lengths and update transition prob's at the
	 ** prune site.
	 */ //JSJ: temp fixes to l
	For(k,tree->n_l){
		v_prune->b[d1]->l_old[k] = v_prune->b[d1]->l[k];
		v_prune->b[d2]->l_old[k] = v_prune->b[d2]->l[k];
		v_prune->b[d1]->l[k] = l_connect[k];
		v_prune->b[d2]->l[k] = l_connect[k];
	}
	Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
	Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

	/*
	 ** Get the relevant average subtree distance within the pruned subtree.
	 */
	if (!u_prune->tax)
	{
		u1 = u2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (u_prune->b[i] != e_prune)
			{
				if (u1 == NULL)
				{
					u1 = u_prune->v[i];
				}
				else
				{
					u2 = u_prune->v[i];
				}
			}
		}
		For(k,tree->n_l)d_uu[k] = subtree_dist[u1->num][u2->num];
	}
	else
	{
		For(k,tree->n_l) d_uu[k] = 0.0;
	}

	/*
	 ** Try each candidate SPR and estimate the change in likelihood.
	 */
	best_d_lk = 1.0/BIG;
	best_cand = -1;
	for (cand = 0; cand < tree->mod->s_opt->wim_n_rgrft; cand++)
	{
		/*
		 ** If there are no more candidates, bail out...
		 */
		if (rgrft_cand[cand]->d_L == -1.0*BIG)
		{
			break;
		}
		else
		{
			nr_d_lk++;
		}

		/*
		 ** Get the relevant nodes and edges.
		 */
		v_n = rgrft_cand[cand]->v_n;
		v_nx1 = rgrft_cand[cand]->v_nx1;
		u_n = rgrft_cand[cand]->u_n;
		e_regraft = rgrft_cand[cand]->e_regraft;

		/*
		 ** Update the relevant partial likelihoods along the path between the prune
		 ** and regraft positions (temporarily save the first one).
		 */
		n = rgrft_cand[cand]->dist;
		e_tmp = NULL;
		p_lk = NULL;
		p_sum = NULL;
		for (i = 1; i <= n; i++)
		{
			/*
			 ** Get the next edge along the path.
			 */
			for (j = 0; j < 3; j++)
			{
				if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
				{
					e_tmp = rgrft_cand[cand]->path[i]->b[j];
					break;
				}
			}
			if (i == 1)
			{
				/*
				 ** Save the first partial likelihood along the path.
				 */
				if (rgrft_cand[cand]->path[i] == e_tmp->left)
				{
					p_lk  = e_tmp->p_lk_left;
					p_sum = e_tmp->sum_scale_f_left;
				}
				else
				{
					p_lk  = e_tmp->p_lk_rght;
					p_sum = e_tmp->sum_scale_f_rght;
				}

				for (pat = 0; pat < tree->n_pattern; pat++)
				{
					sum_scale_tmp[pat] = p_sum[pat];
					for (cat = 0; cat < tree->mod->n_catg; cat++)
					{
						for (ste = 0; ste < tree->mod->ns; ste++)
						{
							p_lk_tmp[pat*dim1+cat*dim2+ste] = p_lk[pat*dim1+cat*dim2+ste];
						}
					}
				}
			}
			Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
		}
		if (v_n == e_regraft->left)
		{
			v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_left;
			v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_rght;
			v_tmp->b[1]->sum_scale_f_rght = e_regraft->sum_scale_f_left;
			v_tmp->b[2]->sum_scale_f_rght = e_regraft->sum_scale_f_rght;
		}
		else
		{
			v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_rght;
			v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_left;
			v_tmp->b[1]->sum_scale_f_rght = e_regraft->sum_scale_f_rght;
			v_tmp->b[2]->sum_scale_f_rght = e_regraft->sum_scale_f_left;
		}

		/*
		 ** Estimate edge lengths of the three relevant regraft edges based on
		 ** average subtree distances.
		 **
		 ** l_01
		 */
		/*
		 ** Alternative method of estimating l_01. Kept it around for reference...
		 **
    l_01 = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
    if (!u_n->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
	if (u_n->v[i] != v_n)
	{
	  if (u1 == NULL)
	  {
	    u1 = u_n->v[i];
	  }
	  else
	  {
	    u2 = u_n->v[i];
	  }
	}
      }
      l_01 -= 0.5 * subtree_dist[u1->num][u2->num];
    }
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] == v_n)
      {
	l_01 -= u_n->b[i]->l;
	break;
      }
    }
		 */
		For(k,tree->n_l) l_01[k] = rgrft_cand[cand]->d_up_v - (0.5 * rgrft_cand[cand]->d_un_v) - (0.5 * d_uu[k]);
		/*
		 ** l_02
		 */
		For(k,tree->n_l) l_02[k] = subtree_dist[u_prune->num][v_nx1->num] - (0.5 * d_uu[k]);
		if (!v_nx1->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_nx1->v[i] != v_n)
				{
					if (u1 == NULL)
					{
						u1 = v_nx1->v[i];
					}
					else
					{
						u2 = v_nx1->v[i];
					}
				}
			}
			For(k,tree->n_l) l_02[k] -= (0.5 * subtree_dist[u1->num][u2->num]);
		}
		/*
		 ** l_12
		 */ //JSJ: temp fixes to l
		For(k,tree->n_l){
			l_12[k] = e_regraft->l[k];

			/*
			 ** Simple estimates.
			 */
			l_simple[k][0] = l_02[k] - (0.5*e_regraft->l[k]);
			l_simple[k][1] = 0.5 * e_regraft->l[k];
			l_simple[k][2] = 0.5 * e_regraft->l[k];
			for (i = 0; i < 3; i++)
			{
				if (l_simple[k][i] < BL_MIN)
				{
					l_simple[k][i] = BL_MIN;
				}
			}
			/*
			 ** Average subtree distance based estimates.
			 */
			l_dist[k][0] = 0.5 * ( l_01[k] + l_02[k] - l_12[k]);
			l_dist[k][1] = 0.5 * ( l_01[k] - l_02[k] + l_12[k]);
			l_dist[k][2] = 0.5 * (-l_01[k] + l_02[k] + l_12[k]);
			for (i = 0; i < 3; i++)
			{
				if (l_dist[k][i] < BL_MIN)
				{
					l_dist[k][i] = BL_MIN;
				}
			}
			/*
			 ** Take the average of the two estimates.
			 */
			l_est[k][0] = (l_simple[k][0] + l_dist[k][0]) / 2.0;
			l_est[k][1] = (l_simple[k][1] + l_dist[k][1]) / 2.0;
			l_est[k][2] = (l_simple[k][2] + l_dist[k][2]) / 2.0;
		}

		/*
		 ** Set the edge lengths and update the relevant transition prob's and
		 ** partial likelihoods in the temporary regraft structure.
		 */
		//JSJ: temp fixes to l
		for (i = 0; i < 3; i++)
		{
			For(k,tree->n_l){
				v_tmp->b[i]->l[k] = l_est[k][i]; /* TO DO */
				if (v_tmp->b[i]->l[k] < BL_MIN)
				{
					v_tmp->b[i]->l[k] = BL_MIN;
				}
			}
			Update_PMat_At_Given_Edge (v_tmp->b[i], tree);
		}

		/* Beg SG 18 May 2007 */
		if(tree->mod->s_opt->wim_inside_opt)
		{ //JSJ: temp fixes to l
			Triple_Dist(v_tmp,tree,0);
			For(k,tree->n_l){
				For(i,3) l_est[k][i] = v_tmp->b[i]->l[k];
			}
		}
		/* End SG 18 May 2007 */


		/*
		 ** Calculate the change in likelihood locally. Save it and the estimated edge
		 ** lengths in the current candidate in the list.
		 */
		Update_P_Lk (tree, v_tmp->b[0], v_tmp);
		new_lk = Lk_At_Given_Edge (v_tmp->b[0],tree);
		/*     PhyML_Printf("\n. new_lk = %f",new_lk); */





		rgrft_cand[cand]->delta_lk = new_lk - cur_lk;
		rgrft_cand[cand]->rgrft_rank = cand;
		rgrft_cand[cand]->optim_rank = -1;
		rgrft_cand[cand]->globl_rank = -1;
		For(k,tree->n_l) rgrft_cand[cand]->l_connect[k] = l_connect[k];
		for (i = 0; i < 3; i++)
		{//JSJ: temp fix to l
			For(k,tree->n_l) rgrft_cand[cand]->l_est[k][i] = v_tmp->b[i]->l[k];
		}
		if (rgrft_cand[cand]->delta_lk > best_d_lk)
		{
			best_d_lk = rgrft_cand[cand]->delta_lk;
			best_cand = cand;
		}

		/*
		 ** If the change is within the tree->mod->s_opt->wim_n_optim best ones, save it in the list of
		 ** optimization candidates.
		 */
		if (rgrft_cand[cand]->delta_lk > optim_cand[tree->mod->s_opt->wim_n_optim-1]->delta_lk)
		{
			i = tree->mod->s_opt->wim_n_optim-1;
			optim_cand[i]->v_prune = rgrft_cand[cand]->v_prune;
			optim_cand[i]->u_prune = rgrft_cand[cand]->u_prune;
			optim_cand[i]->v_n = rgrft_cand[cand]->v_n;
			optim_cand[i]->v_nx1 = rgrft_cand[cand]->v_nx1;
			optim_cand[i]->u_n = rgrft_cand[cand]->u_n;
			optim_cand[i]->e_prune = rgrft_cand[cand]->e_prune;
			optim_cand[i]->e_regraft = rgrft_cand[cand]->e_regraft;
			optim_cand[i]->d_L = rgrft_cand[cand]->d_L;
			optim_cand[i]->dist = rgrft_cand[cand]->dist;
			optim_cand[i]->rgrft_rank = rgrft_cand[cand]->rgrft_rank;
			optim_cand[i]->optim_rank = rgrft_cand[cand]->optim_rank;
			optim_cand[i]->globl_rank = rgrft_cand[cand]->globl_rank;
			For(k,tree->n_l) optim_cand[i]->l_connect[k] = rgrft_cand[cand]->l_connect[k];
			for (j = 0; j < 3; j++)
			{
				For(k,tree->n_l) optim_cand[i]->l_est[k][j] = rgrft_cand[cand]->l_est[k][j];
			}
			optim_cand[i]->delta_lk = rgrft_cand[cand]->delta_lk;
			/*
			 ** Move the candidate to the appropriate position in the list, so the list
			 ** remains sorted in decreasing delta_Lk value.
			 */
			while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
			{
				tmp_cand = optim_cand[i];
				optim_cand[i] = optim_cand[i-1];
				optim_cand[i-1] = tmp_cand;
				i--;
			}
		}

		/*
		 ** Reset the partial likelihoods along the path.
		 */
		for (pat = 0; pat < tree->n_pattern; pat++)
		{
			p_sum[pat] = sum_scale_tmp[pat];
			for (cat = 0; cat < tree->mod->n_catg; cat++)
			{
				for (ste = 0; ste < tree->mod->ns; ste++)
				{
					p_lk[pat*dim1+cat*dim2+ste] = p_lk_tmp[pat*dim1+cat*dim2+ste];
				}
			}
		}
		n = rgrft_cand[cand]->dist;
		for (i = 2; i <= n; i++)
		{
			for (j = 0; j < 3; j++)
			{
				if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
				{
					e_tmp = rgrft_cand[cand]->path[i]->b[j];
					break;
				}
			}
			Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
		}

		/*
		 ** If an improvement was found, forget the other candidates...
		 */
		if (best_cand >= 0)
		{
			break;
		}
	}

	/*
	 ** Swap back the relevant partial likelihoods at the prune site.
	 */
	if (v_prune == v_prune->b[d1]->left)
	{
		v_prune->b[d1]->p_lk_left = p_lk1_tmp;
	}
	else
	{
		v_prune->b[d1]->p_lk_rght = p_lk1_tmp;
	}
	if (v_prune == v_prune->b[d2]->left)
	{
		v_prune->b[d2]->p_lk_left = p_lk2_tmp;
	}
	else
	{
		v_prune->b[d2]->p_lk_rght = p_lk2_tmp;
	}

	/*
	 ** Reset the relevant edge lengths and transition prob's at the prune site.
	 */
	For(k,tree->n_l){
		v_prune->b[d1]->l[k] = v_prune->b[d1]->l_old[k];
		v_prune->b[d2]->l[k] = v_prune->b[d2]->l_old[k];
	}
	Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
	Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

	/*
	 ** Return the best candidate.
	 */
	return (best_cand);
}


/*
 ** Best_Lk_Change: Estimate the changes in likelihood for the most promising candidate
 **                 regraft positions given a pruned subtree and save the best one.
 **
 ** Parameters:
 **   - e_prune: The edge at which the subtree was pruned.
 **   - v_prune: The root of the pruned subtree.
 **   - tree:    The tree on which to do the calculations.
 **
 ** Returns:
 **   The candidate which gives the best (possibly negative) improvement.
 */

int Best_Lk_Change (edge *e_prune, node *v_prune, arbre *tree)
{
	int     i, j, k, cand, best_cand, d0, d1, d2, n, pat, cat, ste;
	m3ldbl  d_uu[MAX_BL_SET], best_d_lk, l_connect[MAX_BL_SET], l_01[MAX_BL_SET], l_02[MAX_BL_SET], l_12[MAX_BL_SET], l_est[MAX_BL_SET][3], new_lk,
	l_simple[MAX_BL_SET][3], l_dist[MAX_BL_SET][3];
	plkflt *p_lk1_tmp, *p_lk2_tmp, *p_lk, *p_sum;
	node   *u_prune, *v_n, *v_nx1, *u_n, *u1, *u2;
	edge   *e_regraft, *e_tmp;
	_move_ *tmp_cand;
	int dim1, dim2;

	dim1 = tree->mod->ns * tree->mod->n_catg;
	dim2 = tree->mod->ns ;

	/*
	 ** Get the directions from node v_prune.
	 */
	d0 = -1;
	u_prune = NULL;
	for (i = 0; i < 3; i++)
	{
		if (v_prune->b[i] == e_prune)
		{
			d0 = i;
			u_prune = v_prune->v[i];
			break;
		}
	}
	d1 = (d0 + 1) % 3;
	d2 = 3 - d0 - d1;

	/*
	 ** Copy the relevant partial likelihoods to the temporary regraft structure.
	 ** We can point to the original matrices, cos they won't be changed anyway.
	 */
	if (v_prune == e_prune->left)
	{
		v_tmp->b[0]->p_lk_rght = e_prune->p_lk_rght;
		v_tmp->b[0]->sum_scale_f_rght = e_prune->sum_scale_f_rght;
	}
	else
	{
		v_tmp->b[0]->p_lk_rght = e_prune->p_lk_left;
		v_tmp->b[0]->sum_scale_f_rght = e_prune->sum_scale_f_left;
	}
	v_tmp->num = v_prune->num;
	v_tmp->v[0]->num = u_prune->num;
	v_tmp->b[0]->num = e_prune->num;

	/*
	 ** Estimate the length of the edge that will connect the two "detached" nodes
	 ** after pruning. (The average of the sum of the lengths of the original two
	 ** edges and the average subtree distance based estimate.)
	 */
	For(k,tree->n_l){
		l_connect[k] = subtree_dist[v_prune->v[d1]->num][v_prune->v[d2]->num];
		if (!v_prune->v[d1]->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_prune->v[d1]->b[i] != v_prune->b[d1])
				{
					if (u1 == NULL)
					{
						u1 = v_prune->v[d1]->v[i];
					}
					else
					{
						u2 = v_prune->v[d1]->v[i];
					}
				}
			}
			l_connect[k] -= 0.5 * subtree_dist[u1->num][u2->num];
		}
		if (!v_prune->v[d2]->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_prune->v[d2]->b[i] != v_prune->b[d2])
				{
					if (u1 == NULL)
					{
						u1 = v_prune->v[d2]->v[i];
					}
					else
					{
						u2 = v_prune->v[d2]->v[i];
					}
				}
			}
			l_connect[k] -= 0.5 * subtree_dist[u1->num][u2->num];
		} //JSJ: temp fix to l
		l_connect[k] += (v_prune->b[d1]->l[k] + v_prune->b[d2]->l[k]);
		l_connect[k] /= 2.0;
		if (l_connect[k] < BL_MIN)
		{
			l_connect[k] = BL_MIN;
		}
	}

	/*
	 ** Temporarily swap the relevant partial likelihoods at the prune site.
	 **
	 ** Direction d1.
	 */
	if (v_prune == v_prune->b[d1]->left)
	{
		p_lk1_tmp = v_prune->b[d1]->p_lk_left;
		if (v_prune == v_prune->b[d2]->left)
		{
			v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_rght;
		}
		else
		{
			v_prune->b[d1]->p_lk_left = v_prune->b[d2]->p_lk_left;
		}
	}
	else
	{
		p_lk1_tmp = v_prune->b[d1]->p_lk_rght;
		if (v_prune == v_prune->b[d2]->left)
		{
			v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_rght;
		}
		else
		{
			v_prune->b[d1]->p_lk_rght = v_prune->b[d2]->p_lk_left;
		}
	}
	/*
	 ** Direction d2.
	 */
	if (v_prune == v_prune->b[d2]->left)
	{
		p_lk2_tmp = v_prune->b[d2]->p_lk_left;
		if (v_prune == v_prune->b[d1]->left)
		{
			v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_rght;
		}
		else
		{
			v_prune->b[d2]->p_lk_left = v_prune->b[d1]->p_lk_left;
		}
	}
	else
	{
		p_lk2_tmp = v_prune->b[d2]->p_lk_rght;
		if (v_prune == v_prune->b[d1]->left)
		{
			v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_rght;
		}
		else
		{
			v_prune->b[d2]->p_lk_rght = v_prune->b[d1]->p_lk_left;
		}
	}

	/*
	 ** Temporarily set the edge lengths and update transition prob's at the
	 ** prune site.
	 */ //JSJ: temp fixes to l
	For(k,tree->n_l){
		v_prune->b[d1]->l_old[k] = v_prune->b[d1]->l[k];
		v_prune->b[d2]->l_old[k] = v_prune->b[d2]->l[k];
		v_prune->b[d1]->l[k] = l_connect[k];
		v_prune->b[d2]->l[k] = l_connect[k];
	}
	Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
	Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

	/*
	 ** Get the relevant average subtree distance within the pruned subtree.
	 */
	if (!u_prune->tax)
	{
		u1 = u2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (u_prune->b[i] != e_prune)
			{
				if (u1 == NULL)
				{
					u1 = u_prune->v[i];
				}
				else
				{
					u2 = u_prune->v[i];
				}
			}
		}
		For(k,tree->n_l) d_uu[k] = subtree_dist[u1->num][u2->num];
	}
	else
	{
		For(k,tree->n_l) d_uu[k] = 0.0;
	}

	/*
	 ** Try the best candidate SPRs and estimate the change in likelihood.
	 */
	best_d_lk = -1.0*BIG;
	best_cand = 0;
	for (cand = 0; cand < tree->mod->s_opt->wim_n_best; cand++)
	{
		/*
		 ** If there are no more candidates, bail out...
		 */
		if (rgrft_cand[cand]->d_L == -1.0*BIG)
		{
			break;
		}
		else
		{
			nr_d_lk++;
		}

		/*
		 ** Get the relevant nodes and edges.
		 */
		v_n = rgrft_cand[cand]->v_n;
		v_nx1 = rgrft_cand[cand]->v_nx1;
		u_n = rgrft_cand[cand]->u_n;
		e_regraft = rgrft_cand[cand]->e_regraft;

		/*
		 ** Update the relevant partial likelihoods along the path between the prune
		 ** and regraft positions (temporarily save the first one).
		 */
		n = rgrft_cand[cand]->dist;
		e_tmp = NULL;
		p_lk = NULL;
		p_sum = NULL;
		for (i = 1; i <= n; i++)
		{
			/*
			 ** Get the next edge along the path.
			 */
			for (j = 0; j < 3; j++)
			{
				if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
				{
					e_tmp = rgrft_cand[cand]->path[i]->b[j];
					break;
				}
			}
			if (i == 1)
			{
				/*
				 ** Save the first partial likelihood along the path.
				 */
				if (rgrft_cand[cand]->path[i] == e_tmp->left)
				{
					p_lk = e_tmp->p_lk_left;
					p_sum = e_tmp->sum_scale_f_left;
				}
				else
				{
					p_lk = e_tmp->p_lk_rght;
					p_sum = e_tmp->sum_scale_f_rght;
				}
				for (pat = 0; pat < tree->n_pattern; pat++)
				{
					sum_scale_tmp[pat] = p_sum[pat];
					for (cat = 0; cat < tree->mod->n_catg; cat++)
					{
						for (ste = 0; ste < tree->mod->ns; ste++)
						{
							p_lk_tmp[pat*dim1+cat*dim2+ste] = p_lk[pat*dim1+cat*dim2+ste];
						}
					}
				}
			}
			Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
		}
		if (v_n == e_regraft->left)
		{
			v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_left;
			v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_rght;
			v_tmp->b[1]->sum_scale_f_rght = e_regraft->sum_scale_f_left;
			v_tmp->b[2]->sum_scale_f_rght = e_regraft->sum_scale_f_rght;
		}
		else
		{
			v_tmp->b[1]->p_lk_rght = e_regraft->p_lk_rght;
			v_tmp->b[2]->p_lk_rght = e_regraft->p_lk_left;
			v_tmp->b[1]->sum_scale_f_rght = e_regraft->sum_scale_f_rght;
			v_tmp->b[2]->sum_scale_f_rght = e_regraft->sum_scale_f_left;
		}

		/*
		 ** Estimate edge lengths of the three relevant regraft edges based on
		 ** average subtree distances.
		 **
		 ** l_01
		 */
		/*
		 ** Alternative method of estimating l_01. Kept it around for reference...
		 **
    l_01 = subtree_dist[u_prune->num][u_n->num] - (0.5 * d_uu);
    if (!u_n->tax)
    {
      u1 = u2 = NULL;
      for (i = 0; i < 3; i++)
      {
	if (u_n->v[i] != v_n)
	{
	  if (u1 == NULL)
	  {
	    u1 = u_n->v[i];
	  }
	  else
	  {
	    u2 = u_n->v[i];
	  }
	}
      }
      l_01 -= 0.5 * subtree_dist[u1->num][u2->num];
    }
    for (i = 0; i < 3; i++)
    {
      if (u_n->v[i] == v_n)
      {
	l_01 -= u_n->b[i]->l;
	break;
      }
    }
		 */
		For(k,tree->n_l)l_01[k] = rgrft_cand[cand]->d_up_v - (0.5 * rgrft_cand[cand]->d_un_v) - (0.5 * d_uu[k]);
		/*
		 ** l_02
		 */
		For(k,tree->n_l) l_02[k] = subtree_dist[u_prune->num][v_nx1->num] - (0.5 * d_uu[k]);
		if (!v_nx1->tax)
		{
			u1 = u2 = NULL;
			for (i = 0; i < 3; i++)
			{
				if (v_nx1->v[i] != v_n)
				{
					if (u1 == NULL)
					{
						u1 = v_nx1->v[i];
					}
					else
					{
						u2 = v_nx1->v[i];
					}
				}
			}
			For(k,tree->n_l) l_02[k] -= (0.5 * subtree_dist[u1->num][u2->num]);
		}

		For(k,tree->n_l){
			/*
			 ** l_12
			 */ //JSJ: temp fixes to l
			l_12[k] = e_regraft->l[k];
			/*
			 ** Simple estimates.
			 */
			l_simple[k][0] = l_02[k] - (0.5*e_regraft->l[k]);
			l_simple[k][1] = 0.5 * e_regraft->l[k];
			l_simple[k][2] = 0.5 * e_regraft->l[k];
			for (i = 0; i < 3; i++)
			{
				if (l_simple[k][i] < BL_MIN)
				{
					l_simple[k][i] = BL_MIN;
				}
			}
			/*
			 ** Average subtree distance based estimates.
			 */
			l_dist[k][0] = 0.5 * ( l_01[k] + l_02[k] - l_12[k]);
			l_dist[k][1] = 0.5 * ( l_01[k] - l_02[k] + l_12[k]);
			l_dist[k][2] = 0.5 * (-l_01[k] + l_02[k] + l_12[k]);
			for (i = 0; i < 3; i++)
			{
				if (l_dist[k][i] < BL_MIN)
				{
					l_dist[k][i] = BL_MIN;
				}
			}
			/*
			 ** Take the average of the two estimates.
			 */
			l_est[k][0] = (l_simple[k][0] + l_dist[k][0]) / 2.0;
			l_est[k][1] = (l_simple[k][1] + l_dist[k][1]) / 2.0;
			l_est[k][2] = (l_simple[k][2] + l_dist[k][2]) / 2.0;

		}

		/*
		 ** Set the edge lengths and update the relevant transition prob's and
		 ** partial likelihoods in the temporary regraft structure.
		 */
		for (i = 0; i < 3; i++)
		{ //JSJ: temp fixes to l
			For(k,tree->n_l){
				v_tmp->b[i]->l[k] = l_est[k][i];
				if (v_tmp->b[i]->l[k] < BL_MIN)
				{
					v_tmp->b[i]->l[k] = BL_MIN;
				}
			}
			Update_PMat_At_Given_Edge (v_tmp->b[i], tree);
		}
		Update_P_Lk (tree, v_tmp->b[0], v_tmp);

		/*
		 ** Calculate the change in likelihood locally. Save it and the estimated edge
		 ** lengths in the current candidate in the list.
		 */
		new_lk = Lk_At_Given_Edge (v_tmp->b[0],tree);
		rgrft_cand[cand]->delta_lk = new_lk - cur_lk;
		rgrft_cand[cand]->rgrft_rank = cand;
		rgrft_cand[cand]->optim_rank = -1;
		rgrft_cand[cand]->globl_rank = -1;
		For(k,tree->n_l) rgrft_cand[cand]->l_connect[k] = l_connect[k];
		for (i = 0; i < 3; i++)
		{//JSJ: temp fixes to l
			For(k,tree->n_l) rgrft_cand[cand]->l_est[k][i] = v_tmp->b[i]->l[k];
		}
		if (rgrft_cand[cand]->delta_lk > best_d_lk)
		{
			best_d_lk = rgrft_cand[cand]->delta_lk;
			best_cand = cand;
		}

		/*
		 ** Reset the partial likelihoods along the path.
		 */
		for (pat = 0; pat < tree->n_pattern; pat++)
		{
			p_sum[pat] = sum_scale_tmp[pat];
			for (cat = 0; cat < tree->mod->n_catg; cat++)
			{
				for (ste = 0; ste < tree->mod->ns; ste++)
				{
					p_lk[pat*dim1+cat*dim2+ste] = p_lk_tmp[pat*dim1+cat*dim2+ste];
				}
			}
		}
		n = rgrft_cand[cand]->dist;
		for (i = 2; i <= n; i++)
		{
			for (j = 0; j < 3; j++)
			{
				if (rgrft_cand[cand]->path[i]->v[j] == rgrft_cand[cand]->path[i+1])
				{
					e_tmp = rgrft_cand[cand]->path[i]->b[j];
					break;
				}
			}
			Update_P_Lk (tree, e_tmp, rgrft_cand[cand]->path[i]);
		}
	}

	/*
	 ** If the best candidate is within the tree->mod->s_opt->wim_n_optim best ones, save it in the list of
	 ** optimization candidates.
	 */
	if (rgrft_cand[best_cand]->delta_lk > optim_cand[tree->mod->s_opt->wim_n_optim-1]->delta_lk)
	{
		i = tree->mod->s_opt->wim_n_optim-1;
		optim_cand[i]->v_prune = rgrft_cand[best_cand]->v_prune;
		optim_cand[i]->u_prune = rgrft_cand[best_cand]->u_prune;
		optim_cand[i]->v_n = rgrft_cand[best_cand]->v_n;
		optim_cand[i]->v_nx1 = rgrft_cand[best_cand]->v_nx1;
		optim_cand[i]->u_n = rgrft_cand[best_cand]->u_n;
		optim_cand[i]->e_prune = rgrft_cand[best_cand]->e_prune;
		optim_cand[i]->e_regraft = rgrft_cand[best_cand]->e_regraft;
		optim_cand[i]->d_L = rgrft_cand[best_cand]->d_L;
		optim_cand[i]->dist = rgrft_cand[best_cand]->dist;
		optim_cand[i]->rgrft_rank = rgrft_cand[best_cand]->rgrft_rank;
		optim_cand[i]->optim_rank = rgrft_cand[best_cand]->optim_rank;
		optim_cand[i]->globl_rank = rgrft_cand[best_cand]->globl_rank;
		For(k,tree->n_l) optim_cand[i]->l_connect[k] = rgrft_cand[best_cand]->l_connect[k];

		for (j = 0; j < 3; j++)
		{
			For(k,tree->n_l) optim_cand[i]->l_est[k][j] = rgrft_cand[best_cand]->l_est[k][j];
		}
		optim_cand[i]->delta_lk = rgrft_cand[best_cand]->delta_lk;
		/*
		 ** Move the candidate to the appropriate position in the list, so the list
		 ** remains sorted in decreasing delta_Lk value.
		 */
		while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
		{
			tmp_cand = optim_cand[i];
			optim_cand[i] = optim_cand[i-1];
			optim_cand[i-1] = tmp_cand;
			i--;
		}
	}

	/*
	 ** Swap back the relevant partial likelihoods at the prune site.
	 */
	if (v_prune == v_prune->b[d1]->left)
	{
		v_prune->b[d1]->p_lk_left = p_lk1_tmp;
	}
	else
	{
		v_prune->b[d1]->p_lk_rght = p_lk1_tmp;
	}
	if (v_prune == v_prune->b[d2]->left)
	{
		v_prune->b[d2]->p_lk_left = p_lk2_tmp;
	}
	else
	{
		v_prune->b[d2]->p_lk_rght = p_lk2_tmp;
	}

	/*
	 ** Reset the relevant edge lengths and transition prob's at the prune site.
	 */
	For(k,tree->n_l){
		v_prune->b[d1]->l[k] = v_prune->b[d1]->l_old[k];
		v_prune->b[d2]->l[k] = v_prune->b[d2]->l_old[k];
	}
	Update_PMat_At_Given_Edge (v_prune->b[d1], tree);
	Update_PMat_At_Given_Edge (v_prune->b[d2], tree);

	/*
	 ** Return the best candidate.
	 */
	return (best_cand);
}


/*
 ** Make_Move: Perform an actual SPR move and calculate the new likelihood.
 **
 ** Parameters:
 **   - candidate: The candidate move to perform.
 **   - tree:      The tree on which to perform the move.
 **
 */

void Make_Move (_move_ *move, int type, arbre *tree)
{
	int     i,j;
	node   *v_prune, *u_prune, *v_n, *root;
	edge   *e_prune, *e_regraft, *e_connect, *e_avail;
	m3ldbl  new_lk;

	/*
	 ** Get the relevant nodes and edges.
	 */
	v_prune = move->v_prune;
	u_prune = move->u_prune;
	v_n = move->v_n;
	e_prune = move->e_prune;
	e_regraft = move->e_regraft;
	/*   PhyML_Printf ("  making move: %d -> %d (%f)\n", e_prune->num, e_regraft->num, move->delta_lk); */
	/*
	 ** Perform the move and set edge lengths.
	 */
	Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
	Regraft (e_regraft, v_prune, e_avail, tree);
	For(j,tree->n_l) e_connect->l[j] = move->l_connect[j]; //JSJ: temp fix to l

	for (i = 0; i < 3; i++)
	{
		if (v_prune->v[i] == u_prune)
		{//JSJ: temp fixes to l
			For(j,tree->n_l) v_prune->b[i]->l[j] = move->l_est[j][0];
		}
		else if (v_prune->v[i] == v_n)
		{
			For(j,tree->n_l)v_prune->b[i]->l[j] = move->l_est[j][1];
		}
		else
		{
			For(j,tree->n_l)v_prune->b[i]->l[j] = move->l_est[j][2];
		}
	}


	if(type > 0) /* local or global move */
	{
		Restore_Br_Len(NULL,tree);
	}

	/*
	 ** Calculate the new likelihood.
	 */
	tree->both_sides = 1;
	new_lk = Return_Lk (tree);

	if(tree->c_lnL < cur_lk-tree->mod->s_opt->min_diff_lk_local)
	{
		PhyML_Printf("\n. tree->c_lnL = %f cur_lk = %f",tree->c_lnL,cur_lk);
		PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	cur_lk = new_lk;

	/*
	 ** Recalculate the average distances between all (non-overlapping) subtrees.
	 */
	root = tree->noeud[0];
	PostOrder_v (tree, root->v[0], root->b[0]);
}


/*
 ** Find_Optim_Local: Perform local edge length optimization on the candidates in the
 **                   optimization list, and return the first one that gives an
 **                   improvement in likelihood.
 **
 ** Parameters:
 **   - tree: The tree on which to check the moves.
 **
 ** Returns:
 **   If an improvement was found: The candidate that gives the improvement.
 **   Otherwise:                   -1.
 */

int Find_Optim_Local (arbre *tree)
{
	int     best_cand, cand, i, j;
	node   *v_prune, *u_prune, *v_n;
	edge   *e_prune, *e_regraft, *e_connect, *e_avail;
	m3ldbl  max_change, new_lk;
	_move_ *move, *tmp_cand;

	/*
	 ** Try all candidate moves starting from the first one.
	 */
	best_cand = -1;
	max_change = 1.0/BIG;
	for(cand = 0; cand < tree->mod->s_opt->wim_n_optim; cand++)
	{
		move = optim_cand[cand];
		if(move->delta_lk > -1.0*BIG)
		{
			/*
			 ** Get the relevant nodes and edges.
			 */
			v_prune   = move->v_prune;
			u_prune   = move->u_prune;
			v_n       = move->v_n;
			e_prune   = move->e_prune;
			e_regraft = move->e_regraft;

			/*
			 ** Perform the move and set edge lengths.
			 */ //JSJ: temp fixes of l
			Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
			Regraft (e_regraft, v_prune, e_avail, tree);
			For(j,tree->n_l){
				e_connect->l_old[j] = e_connect->l[j];
				e_connect->l[j] = move->l_connect[j];
			}

			for (i = 0; i < 3; i++)
			{
				For(j,tree->n_l){
					v_prune->b[i]->l_old[j] = v_prune->b[i]->l[j];

					if (v_prune->v[i] == u_prune)
					{
						v_prune->b[i]->l[j] = move->l_est[j][0];//JSJ: havn't yet modified move
					}
					else if (v_prune->v[i] == v_n)
					{
						v_prune->b[i]->l[j] = move->l_est[j][1];
					}
					else
					{
						v_prune->b[i]->l[j] = move->l_est[j][2];
					}
				}
			}

			tree->both_sides = 1;
			Lk(tree);  // Not sure anymore whether this is required...

			/*
			 ** Use Brent optimization on the relevant edges at the regraft position
			 ** and calculate the new likelihood value.
			 */
			int tmp;
//			m3ldbl max[3][MAX_BL_SET];
//			m3ldbl min[3][MAX_BL_SET];
//			For(tmp,3){
//				For(j,tree->n_l){
//					max[tmp][j] = v_prune->b[tmp]->l[j];
//					max[tmp][j] *= 10.0;
//					min[tmp][j] =BL_MIN;
//
//				}
//			}
			For(tmp,3){
				For(j,tree->n_l){
					Br_Len_Brent_Iter(v_prune->b[tmp]->l[j] * 10.0, v_prune->b[tmp]->l[j], BL_MIN, 1.e-10, v_prune->b[tmp], tree, 250, 0,j);
				}
			}

			//			Br_Len_Brent (10.*(v_prune->b[0]->l[0]), v_prune->b[0]->l[0], BL_MIN, 1.e-10,
			//					v_prune->b[0], tree, 250, 0);
			//			Br_Len_Brent (10.*(v_prune->b[1]->l[0]), v_prune->b[1]->l[0], BL_MIN, 1.e-10,
			//					v_prune->b[1], tree, 250, 0);
			//			Br_Len_Brent (10.*(v_prune->b[2]->l[0]), v_prune->b[2]->l[0], BL_MIN, 1.e-10,
			//					v_prune->b[2], tree, 250, 0);

			/* 	  Update_PMat_At_Given_Edge (v_prune->b[0], tree); */
			/* 	  Update_PMat_At_Given_Edge (v_prune->b[1], tree); */
			/* 	  Update_PMat_At_Given_Edge (v_prune->b[2], tree); */

			Update_P_Lk (tree, v_prune->b[0], v_prune);
			new_lk = Lk_At_Given_Edge (v_prune->b[0],tree);

			/* 	  PhyML_Printf("\n. local new_lk = %f",new_lk); */

			/*
			 ** Save the change in likelihood and move the current candidate to the
			 ** appropriate place in the list.
			 */
			move->delta_lk = new_lk - cur_lk;
			move->optim_rank = cand;
			i = cand;
			while ((i > 0) && (optim_cand[i]->delta_lk > optim_cand[i-1]->delta_lk))
			{
				tmp_cand = optim_cand[i];
				optim_cand[i] = optim_cand[i-1];
				optim_cand[i-1] = tmp_cand;
				i--;
			}
			if (move->delta_lk > max_change)
			{
				best_cand = i;
				max_change = move->delta_lk;
				Record_Br_Len(NULL,tree);
			}


			/*
			 ** Undo the move again.
			 */
			Prune (e_prune, v_prune, &(e_regraft), &(e_avail), tree);
			Regraft (e_connect, v_prune, e_avail, tree);
			For(j,tree->n_l) e_regraft->l[j] = e_regraft->l_old[j];
			for (i = 0; i < 3; i++)
			{
				For(j,tree->n_l) v_prune->b[i]->l[j] = v_prune->b[i]->l_old[j];
			}
			tree->both_sides = 1;
			Lk(tree);
			nr_loc++;
			/* 	  PhyML_Printf("\n. local back to = %f",tree->c_lnL); */
		}
		else
		{
			break;
		}

		/*
		 ** If an improvement was found, forget the other candidates...
		 */
		if (best_cand >= 0)
		{
			break;
		}
	}

	/*
	 ** Return the best candidate.
	 */
	return (best_cand);
}


/*
 ** Find_Optim_Globl: Perform global edge length optimization on the candidates in the
 **                   optimization list, and return the first one that gives an
 **                   improvement in likelihood.
 **
 ** Parameters:
 **   - tree: The tree on which to check the moves.
 **
 ** Returns:
 **   If an improvement is found: The candidate that gives the improvement.
 **   Otherwise:                  -1.
 */

int Find_Optim_Globl (arbre *tree)
{
	int     best_cand, cand, i, j;
	node   *v_prune, *u_prune, *v_n, *root;
	edge   *e_prune, *e_regraft, *e_connect, *e_avail;
	m3ldbl  max_change, new_lk;
	_move_ *move;

	/*
	 ** Try all moves.
	 */
	best_cand = -1;
	max_change = 1.0/BIG;
	for (cand = 0; cand < tree->mod->s_opt->wim_n_globl; cand++)
	{
		move = optim_cand[cand];
		if (move->delta_lk > -1.0*BIG)
		{
			Record_Br_Len(NULL,tree);

			/*
			 ** Get the relevant nodes and edges.
			 */
			v_prune = move->v_prune;
			u_prune = move->u_prune;
			v_n = move->v_n;
			e_prune = move->e_prune;
			e_regraft = move->e_regraft;

			/*
			 ** Perform the move and optimize all edge lengths.
			 */
			Prune (e_prune, v_prune, &(e_connect), &(e_avail), tree);
			Regraft (e_regraft, v_prune, e_avail, tree);
			For(j,tree->n_l){
				e_connect->l_old[j] = e_connect->l[j]; //JSJ: temp fixes to l
				e_connect->l[j] = move->l_connect[j];
			}

			for (i = 0; i < 3; i++)
			{
				For(j,tree->n_l) v_prune->b[i]->l_old[j] = v_prune->b[i]->l[j];
				if (v_prune->v[i] == u_prune)
				{
					For(j,tree->n_l) v_prune->b[i]->l[j] = move->l_est[j][0];
				}
				else if (v_prune->v[i] == v_n)
				{
					For(j,tree->n_l) v_prune->b[i]->l[j] = move->l_est[j][1];
				}
				else
				{
					For(j,tree->n_l) v_prune->b[i]->l[j] = move->l_est[j][2];
				}
			}

			tree->both_sides = 1;
			Lk(tree);
			root = tree->noeud[0];
			Optimize_Br_Len_Serie (root, root->v[0], root->b[0], tree, tree->data);
			tree->both_sides = 1;
			new_lk = Return_Lk (tree);

			/*       PhyML_Printf("\n. global new_lk = %f\n",tree->c_lnL); */

			/*
			 ** Save the change in likelihood and undo the move.
			 */

			move->delta_lk = new_lk - cur_lk;
			move->globl_rank = cand;
			if (move->delta_lk > max_change)
			{
				best_cand = cand;
				max_change = move->delta_lk;
				Record_Br_Len(NULL,tree);
			}
			//JSJ: temp fixes to l
			Prune (e_prune, v_prune, &(e_regraft), &(e_avail), tree);
			Regraft (e_connect, v_prune, e_avail, tree);
			For(j,tree->n_l){
				e_regraft->l[j] = e_regraft->l_old[j];
				for (i = 0; i < 3; i++)
				{ //JSJ: temp fixes to l
					v_prune->b[i]->l[j] = v_prune->b[i]->l_old[j];
				}
			}
			tree->both_sides = 1;
			Restore_Br_Len(NULL,tree);
			Lk(tree);
			nr_glb++;
			/*       PhyML_Printf("\n. global back to = %f",tree->c_lnL); */

		}
		else break;

		/*
		 ** If an improvement was found, forget the other candidates...
		 */
		if (best_cand >= 0) break;
	}

	/*
	 ** Return the best candidate.
	 */
	return (best_cand);
}


/*
 ** Prune: Prune the subtree at a certain edge and node. Note that edge
 **        lengths are not set and partial likelihoods are not updated.
 **
 ** Parameters:
 **   - e:         The edge at which to prune the subtree.
 **   - v:         The node adjacent to edge e which forms the root of the subtree.
 **   - e_connect: An edge pointer which will point to the edge that was left
 **                after pruning.
 **   - e_avail:   The edge that is "left over" and should be used in the
 **                regrafting step.
 **   - tree:      The tree on which the pruning is done.
 **
 **
 **
 **          \ /
 **           o
 **           |
 **           | e    --> subtree to be pruned
 **           |
 **           o v
 **          / \
 **       e1/   \e2
 **        /     \
 **       o       o
 **      u1       u2   --> such that u1->num < u2->num
 */

void Prune (edge *e, node *v, edge **e_connect, edge **e_avail, arbre *tree)
{
	int     dir0, dir1, dir2, v0, v1, v2, tmp_dir, i, j, k, m;
	node   *u1, *u2, *tmp_node;
	edge   *e1, *e2;
	plkflt *sum_scale_f, *p_lk;
	int dim1, dim2;


	dim1 = tree->mod->ns * tree->mod->n_catg;
	dim2 = tree->mod->ns;

	/*
	 ** Get the relevant directions, nodes and edges.
	 ** Make sure that node u1 is the node with the smaller number.
	 */
	dir0 = -1;
	for (i = 0; i < 3; i++)
	{
		if (v->b[i] == e)
		{
			dir0 = i;
			break;
		}
	}
	dir1 = (dir0 + 1) % 3;
	dir2 = 3 - dir0 - dir1;
	u1 = v->v[dir1];
	u2 = v->v[dir2];
	if (u1->num > u2->num)
	{
		tmp_node = u1;
		u1 = u2;
		u2 = tmp_node;
		tmp_dir = dir1;
		dir1 = dir2;
		dir2 = tmp_dir;
	}
	e1 = v->b[dir1];
	e2 = v->b[dir2];

	/*
	 ** Detach node v from the tree.
	 */
	v->v[dir1] = NULL;
	v->v[dir2] = NULL;
	v->b[dir1] = NULL;
	v->b[dir2] = NULL;

	/*
	 ** Connect nodes u1 and u2 via edge e1 and copy relevant partial likelihoods.
	 */
	if (u2 == e2->left)
	{
		v0 = e2->l_r;
		v1 = e2->l_v1;
		v2 = e2->l_v2;
		sum_scale_f = e2->sum_scale_f_left;
		p_lk = e2->p_lk_left;
	}
	else
	{
		v0 = e2->r_l;
		v1 = e2->r_v1;
		v2 = e2->r_v2;
		sum_scale_f = e2->sum_scale_f_rght;
		p_lk = e2->p_lk_rght;
	}
	if (u1 == e1->left)
	{
		e1->rght = u2;
		e1->r_l = v0;
		e1->r_v1 = v1;
		e1->r_v2 = v2;
		for (i = 0; i < tree->n_pattern; i++)
		{
			e1->sum_scale_f_rght[i] = sum_scale_f[i];
			for (j = 0; j < tree->mod->n_catg; j++)
			{
				for (k = 0; k < tree->mod->ns; k++)
				{
					e1->p_lk_rght[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
				}
			}
		}
	}
	else
	{
		e1->left = u2;
		e1->l_r = v0;
		e1->l_v1 = v1;
		e1->l_v2 = v2;
		for (i = 0; i < tree->n_pattern; i++)
		{
			e1->sum_scale_f_left[i] = sum_scale_f[i];
			for (j = 0; j < tree->mod->n_catg; j++)
			{
				for (k = 0; k < tree->mod->ns; k++)
				{
					e1->p_lk_left[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
				}
			}
		}
	}
	for (i = 0; i < 3; i++)
	{
		if (u1->v[i] == v)
		{
			u1->v[i] = u2;
		}
		if (u2->v[i] == v)
		{
			u2->v[i] = u1;
			u2->b[i] = e1;
			For(m,tree->n_l) u2->l[m][i] = e1->l[m];
		}
	}

	/*
	 ** Make sure that a possible tip node is still on the right side.
	 */
	if (e1->left->tax)
	{
		/*
		 ** Swap left and right.
		 */
		tmp_node = e1->left;
		e1->left = e1->rght;
		e1->rght = tmp_node;
		tmp_dir = e1->l_r;
		e1->l_r = e1->r_l;
		e1->r_l = tmp_dir;
		tmp_dir = e1->l_v1;
		e1->l_v1 = e1->r_v1;
		e1->r_v1 = tmp_dir;
		tmp_dir = e1->l_v2;
		e1->l_v2 = e1->r_v2;
		e1->r_v2 = tmp_dir;
		p_lk = e1->p_lk_left;
		e1->p_lk_left = e1->p_lk_rght;
		e1->p_lk_rght = p_lk;
		sum_scale_f = e1->sum_scale_f_left;
		e1->sum_scale_f_left = e1->sum_scale_f_rght;
		e1->sum_scale_f_rght = sum_scale_f;
	}

	/*
	 ** Set the connecting and available edges.
	 */
	*(e_connect) = e1;
	*(e_avail) = e2;
}


/*
 ** Regraft: Regraft a subtree at a certain edge. Note that edge lengths
 **          are not set and partial likelihoods are not updated.
 **
 ** Parameters:
 **   - e:     The edge to regraft the subtree on.
 **   - v:     The root of the subtree to regraft.
 **   - avail: A previously deleted edge now available for insertion again.
 **   - tree:  The tree on which the regrafting is done.
 **
 **
 **      \ /
 **       o
 **       |
 **       |   --> subtree to regraft
 **       |
 **       o v
 **
 **   o--------o
 **  u1   e    u2   --> such that u1->num < u2->num
 */

void Regraft (edge *e, node *v, edge *avail, arbre *tree)
{
	int     dir0, dir1, dir2, i, j, k;
	plkflt *sum_scale_f, *p_lk;
	node   *u1, *u2;
	int dim1, dim2;

	dim1 = tree->mod->ns * tree->mod->n_catg;
	dim2 = tree->mod->ns ;

	/*
	 ** Get the relevant directions and nodes.
	 */
	dir0 = -1;
	for (i = 0; i < 3; i++)
	{
		if (v->b[i] != NULL)
		{
			dir0 = i;
			break;
		}
	}
	dir1 = (dir0 + 1) % 3;
	dir2 = 3 - dir0 - dir1;
	if (e->left->num < e->rght->num)
	{
		u1 = e->left;
		u2 = e->rght;
		sum_scale_f = e->sum_scale_f_rght;
		p_lk = e->p_lk_rght;
	}
	else
	{
		u1 = e->rght;
		u2 = e->left;
		sum_scale_f = e->sum_scale_f_left;
		p_lk = e->p_lk_left;
	}

	/*
	 ** Connect nodes v and u2 via the available edge 'avail' and copy the
	 ** relevant partial likelihood.
	 ** (We want to do this first, cos we need some of the values of edge
	 **  e before changing them below).
	 */
	avail->left = v;
	avail->rght = u2;
	avail->l_r = dir2;
	avail->l_v1 = dir0;
	avail->l_v2 = dir1;
	v->v[dir2] = u2;
	v->b[dir2] = avail;
	for (i = 0; i < 3; i++)
	{
		if (e == u2->b[i])
		{
			u2->v[i] = v;
			u2->b[i] = avail;
			avail->r_l = i;
			avail->r_v1 = (i + 1) % 3;
			avail->r_v2 = 3 - i - avail->r_v1;
			break;
		}
	}
	for (i = 0; i < tree->n_pattern; i++)
	{
		avail->sum_scale_f_rght[i] = sum_scale_f[i];
		for (j = 0; j < tree->mod->n_catg; j++)
		{
			for (k = 0; k < tree->mod->ns; k++)
			{
				avail->p_lk_rght[i*dim1+j*dim2+k] = p_lk[i*dim1+j*dim2+k];
			}
		}
	}

	/*
	 ** Connect nodes v and u1 via edge e.
	 */
	if (u1 == e->left)
	{
		e->rght = v;
		e->r_l = dir1;
		e->r_v1 = dir0;
		e->r_v2 = dir2;
	}
	else
	{
		e->left = v;
		e->l_r = dir1;
		e->l_v1 = dir0;
		e->l_v2 = dir2;
	}
	v->v[dir1] = u1;
	v->b[dir1] = e;
	for (i = 0; i < 3; i++)
	{
		if (u1->v[i] == u2)
		{
			u1->v[i] = v;
			break;
		}
	}
}


/*
 ** PostOrder_v: Recursively visit all nodes v in postorder to calculate
 **              the average distance between subtrees.
 **
 ** Parameters:
 **   - tree: The tree for which to calculate the average distances.
 **   - v:    The current node.
 **   - e:    The edge we came from.
 */

void PostOrder_v (arbre *tree, node *v, edge *e)
{
	int   i;
	node *w;

	/*
	 ** If v is not a taxon, recurse.
	 */
	if (!v->tax)
	{
		for (i = 0; i < 3; i++)
		{
			if (v->b[i] != e)
			{
				PostOrder_v (tree, v->v[i], v->b[i]);
			}
		}
	}

	/*
	 ** Recurse over all nodes w not in the current subtree and calculate
	 ** the average distance between the current subtree and all others.
	 */
	if (v == e->left)
	{
		w = e->rght;
	}
	else
	{
		w = e->left;
	}
	PostOrder_w (tree, v, e, w, e);
}


/*
 ** PostOrder_w: Recursively visit all nodes w not in the subtree of v in
 **              postorder and calculate the average distance between the
 **              subtree of v and all others.
 **
 ** Parameters:
 **   - tree: The tree for which to calculate the average distances.
 **   - v:    The root node of the first subtree.
 **   - v_e:  The edge adjacent to the first subtree.
 **   - w:    The current node.
 **   - e:    The edge we came from.
 */

void PostOrder_w (arbre *tree, node *v, edge *v_e, node *w, edge *e)
{
	int   i;
	node *w1, *w2, *v1, *v2;

	/*
	 ** If w is not a taxon, recurse.
	 */
	if (!w->tax)
	{
		for (i = 0; i < 3; i++)
		{
			if (w->b[i] != e)
			{
				PostOrder_w (tree, v, v_e, w->v[i], w->b[i]);
			}
		}
	}

	/*
	 ** Calculate the average distance between the subtrees defined by
	 ** nodes v and w.
	 */
	if (v->tax && w->tax)
	{
		subtree_dist[v->num][w->num] = seq_dist->dist[v->num][w->num];
	}
	else if (!v->tax)
	{
		v1 = v2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (v->b[i] != v_e)
			{
				if (v1 == NULL)
				{
					v1 = v->v[i];
				}
				else
				{
					v2 = v->v[i];
				}
			}
		}
		subtree_dist[v->num][w->num] = 0.5*(subtree_dist[v1->num][w->num] +
				subtree_dist[v2->num][w->num]);
	}
	else
	{
		w1 = w2 = NULL;
		for (i = 0; i < 3; i++)
		{
			if (w->b[i] != e)
			{
				if (w1 == NULL)
				{
					w1 = w->v[i];
				}
				else
				{
					w2 = w->v[i];
				}
			}
		}
		subtree_dist[v->num][w->num] = 0.5*(subtree_dist[v->num][w1->num] +
				subtree_dist[v->num][w2->num]);
	}
	subtree_dist[w->num][v->num] = subtree_dist[v->num][w->num];
}


/*********************************************************/
/*********************************************************/
/*********************************************************/
/* Below are my functions for SPR search (Stephane Guindon, 2007) */


void Randomize_Spr_List(arbre *tree)
{
	int i,j;
	spr *buff;

	For(i,tree->size_spr_list)
	{
		j = (int)floor(rand()/(RAND_MAX+1.)*tree->size_spr_list);
		buff              = tree->spr_list[i];
		tree->spr_list[i] = tree->spr_list[j];
		tree->spr_list[j] = buff;
	}
}

/*********************************************************/
//JSJ: currently a bug somewhere in this function
int Spr(m3ldbl init_lnL, arbre *tree)
{
	int br;
	int pars_diff, max_pars_diff, new_pars, old_pars;
	edge *b;

	tree->both_sides = 1;
	pars_diff        = -1;
	max_pars_diff    = -1;
	new_pars         = -1;
	old_pars         = -1;

	Reset_Spr_List(tree);
	//JSJ: OK to here... The bug is in the following for loop!!!
	For(br,2*tree->n_otu-3)
	{
		b = tree->t_edges[br];

		old_pars = tree->c_pars;
		Spr_Subtree(b,b->left,tree);
		new_pars = tree->c_pars;
		//JSJ: gets weird after second time through above...
		pars_diff =  new_pars - old_pars;
		if(pars_diff > max_pars_diff) max_pars_diff = pars_diff;

		old_pars = tree->c_pars;
		Spr_Subtree(b,b->rght,tree);
		new_pars = tree->c_pars;
		//JSJ: not altered further by this point at second itteration...
		pars_diff = new_pars - old_pars;
		if(pars_diff > max_pars_diff) max_pars_diff = pars_diff;
	}

	/*   tree->mod->s_opt->pars_thresh = MAX(5,max_pars_diff); */

	return 1;
}

/*********************************************************/

void Spr_Subtree(edge *b, node *link, arbre *tree)
{
	int i;
	int n_moves_pars, n_moves, curr_pars, min_pars, best_move;
	spr *best_pars_move;
	edge *target, *residual;
	int ret_val;

	best_move     = -1;
	tree->n_moves = 0;
	curr_pars     = tree->c_pars;
	ret_val       = 0;

	if((link != b->left) && (link != b->rght))
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	else
	{
		if(!link->tax) Test_All_Spr_Targets(b,link,tree);

		if(tree->n_moves)
		{
			n_moves_pars = 0;
			n_moves      = 0;

			For(i,tree->n_moves)
			if(curr_pars - tree->spr_list[i]->pars >= -tree->mod->s_opt->pars_thresh)
				n_moves_pars++;
			n_moves_pars = MAX(n_moves_pars,1);

			if(tree->mod->s_opt->spr_lnL) n_moves = 15;
			else                          n_moves = n_moves_pars;

			n_moves = MIN(n_moves,2*tree->n_otu-3);

			if(tree->mod->s_opt->spr_pars)
			{
				min_pars = 1E+8;
				best_pars_move = NULL;

				For(i,n_moves)
				if(tree->spr_list[i]->pars < min_pars)
				{
					best_pars_move = tree->spr_list[i];
					min_pars = tree->spr_list[i]->pars;
				}

				if(best_pars_move->pars < tree->best_pars)
				{
					Prune_Subtree(best_pars_move->n_link,best_pars_move->n_opp_to_link,&target,&residual,tree);
					Graft_Subtree(best_pars_move->b_target,best_pars_move->n_link,residual,tree);
					tree->both_sides = 1;
					Pars(tree);
					tree->best_pars = tree->c_pars;
					if(tree->best_pars != best_pars_move->pars)
					{
						printf("\n. best_pars = %d move_pars = %d",tree->best_pars,best_pars_move->pars);
						PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
						Warn_And_Exit("");
					}
					tree->n_improvements++;
				}
				else
					Pars(tree);
			}
			else
			{
				best_move = Evaluate_List_Of_Regraft_Pos_Triple(tree->spr_list,n_moves,tree);

				if(tree->spr_list[best_move]->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
				{
					Try_One_Spr_Move_Triple(tree->spr_list[best_move],tree);
					ret_val = 1;
				}
				else
				{
					Pars(tree);
				}
			}
		}
		Reset_Spr_List(tree);
	}
}


/*********************************************************/

int Test_All_Spr_Targets(edge *b_pulled, node *n_link, arbre *tree)
{
	node *n_opp_to_link,*n_v1,*n_v2,*n_up;
	edge *b_target,*b_residual;
	int i,k,dir1,dir2;
	m3ldbl init_len_v1[MAX_BL_SET], init_len_v2[MAX_BL_SET], init_len_pulled[MAX_BL_SET];
	int best_found,approx;
	m3ldbl init_lnL;


	init_lnL = tree->c_lnL;
	n_up = NULL;
	b_target = b_residual = NULL;
	n_opp_to_link  = (n_link == b_pulled->rght)?(b_pulled->left):(b_pulled->rght);
	approx = 1;
	//JSJ: temp fix to l
	For(k,tree->n_l) init_len_pulled[k] = b_pulled->l[k];
	dir1 = dir2 = -1;
	For(i,3)
	if(n_link->v[i] != n_opp_to_link)
	{
		if(dir1<0) dir1 = i;
		else       dir2 = i;
	}

	if(n_link->v[dir1]->num < n_link->v[dir2]->num)
	{
		n_v1        = n_link->v[dir1];
		n_v2        = n_link->v[dir2]; //JSJ: temp fixes to l
		For(k,tree->n_l){
			init_len_v1[k] = n_link->b[dir1]->l[k];
			init_len_v2[k] = n_link->b[dir2]->l[k];
		}
	}
	else
	{
		n_v1        = n_link->v[dir2];
		n_v2        = n_link->v[dir1];
		For(k,tree->n_l){
			init_len_v1[k] = n_link->b[dir2]->l[k];
			init_len_v2[k] = n_link->b[dir1]->l[k];
		}
	}

	if(!(n_v1->tax && n_v2->tax)) /* Pruning is meaningless otherwise */
	{
		Prune_Subtree(n_link,n_opp_to_link,&b_target,&b_residual,tree);

		if(tree->mod->s_opt->spr_lnL)
		{
			Fast_Br_Len(b_target,tree,0);
			/* 	  Update_PMat_At_Given_Edge(b_target,tree); */
		}

		best_found = 0;
		tree->depth_curr_path = 0;
		tree->curr_path[0] = b_target->left;
		Test_One_Spr_Target_Recur(b_target->rght,
				b_target->left,
				b_pulled,n_link,b_residual,&best_found,tree);

		tree->depth_curr_path = 0;
		tree->curr_path[0] = b_target->rght;
		Test_One_Spr_Target_Recur(b_target->left,
				b_target->rght,
				b_pulled,n_link,b_residual,&best_found,tree);

		Graft_Subtree(b_target,n_link,b_residual,tree);

		if((n_link->v[dir1] != n_v1) || (n_link->v[dir2] != n_v2))
			PhyML_Printf("\n. Warning : -- SWITCH NEEDED -- ! \n");
		//JSJ: temp fixes to l
		For(k,tree->n_l){
			n_link->b[dir1]->l[k] = init_len_v1[k];
			n_link->b[dir2]->l[k] = init_len_v2[k];
			b_pulled->l[k] = init_len_pulled[k];
		}


		Update_PMat_At_Given_Edge(n_link->b[dir1],tree);
		Update_PMat_At_Given_Edge(n_link->b[dir2],tree);
		Update_PMat_At_Given_Edge(b_pulled,tree);


		if(tree->mod->s_opt->spr_lnL)
		{
			Update_P_Lk(tree,b_pulled,  n_link);
			Update_P_Lk(tree,b_target,  n_link);
			Update_P_Lk(tree,b_residual,n_link);
		}
		else
		{
			Update_P_Pars(tree,b_pulled,  n_link);
			Update_P_Pars(tree,b_target,  n_link);
			Update_P_Pars(tree,b_residual,n_link);
		}

		For(i,3)
		if(n_link->v[i] != n_opp_to_link)
		{
			if(tree->mod->s_opt->spr_lnL) Pre_Order_Lk(n_link,n_link->v[i],tree);
			else                          Pre_Order_Pars(n_link,n_link->v[i],tree);
		}
	}

	tree->c_lnL = init_lnL;

	return 0;

}

/*********************************************************/

void Test_One_Spr_Target_Recur(node *a, node *d, edge *pulled, node *link, edge *residual, int *best_found, arbre *tree)
{
	int i;

	if(*best_found) return;

	if(d->tax) return;
	else
	{
		m3ldbl move_lnL;

		For(i,3)
		{
			if(d->v[i] != a)
			{
				if(tree->mod->s_opt->spr_lnL) Update_P_Lk(tree,d->b[i],d);
				else                          Update_P_Pars(tree,d->b[i],d);

				tree->depth_curr_path++;
				tree->curr_path[tree->depth_curr_path] = d->v[i];

				if((tree->depth_curr_path <= tree->mod->s_opt->max_depth_path) &&
						(tree->depth_curr_path >= tree->mod->s_opt->min_depth_path))
				{
					move_lnL = Test_One_Spr_Target(d->b[i],pulled,link,residual,tree);
					if(move_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) *best_found = 1;
				}

				if(tree->depth_curr_path < tree->mod->s_opt->max_depth_path)
					Test_One_Spr_Target_Recur(d,d->v[i],pulled,link,residual,best_found,tree);

				tree->depth_curr_path--;
			}
		}
	}
}

/*********************************************************/

m3ldbl Test_One_Spr_Target(edge *b_target, edge *b_arrow, node *n_link, edge *b_residual, arbre *tree)
{
	m3ldbl init_target_len[MAX_BL_SET], init_arrow_len[MAX_BL_SET], init_residual_len[MAX_BL_SET];
	int i,k,dir_v0,dir_v1,dir_v2;
	m3ldbl l0[MAX_BL_SET],l1[MAX_BL_SET],l2[MAX_BL_SET];
	node *v1, *v2;
	m3ldbl init_lnL, move_lnL;
	int init_pars,move_pars;
	int approx;

	tree->n_moves++;

	approx    = 1;
	move_lnL  = UNLIKELY;
	init_lnL  = tree->c_lnL;
	init_pars = tree->c_pars;

	Graft_Subtree(b_target,n_link,b_residual,tree);
	//JSJ: temp fixes to l
	For(k,tree->n_l){
		init_target_len[k]   = b_target->l[k];
		init_arrow_len[k]    = b_arrow->l[k];
		init_residual_len[k] = b_residual->l[k];
	}

	if(tree->mod->s_opt->spr_lnL)
	{
		/*       move_lnL = Triple_Dist(n_link,tree,1); */
		Update_PMat_At_Given_Edge(b_target,tree);
		Update_PMat_At_Given_Edge(b_arrow,tree);
		Update_P_Lk(tree,b_residual,n_link);
		move_lnL = Lk_At_Given_Edge(b_residual,tree);
	}
	else
	{
		Update_P_Pars(tree,b_residual,n_link);
		move_pars = Pars_At_Given_Edge(b_residual,tree);
	}

	v1 = (b_residual->left == n_link)?(b_residual->rght):(b_residual->left);
	v2 = (b_target->left   == n_link)?(b_target->rght):(b_target->left);
	dir_v1 = dir_v2 = dir_v0 = -1;
	For(i,3)
	{
		if(n_link->v[i]      == v1) dir_v1 = i;
		else if(n_link->v[i] == v2) dir_v2 = i;
		else                        dir_v0 = i;
	} //JSJ: temp fixes to l

	For(k,tree->n_l) l0[k] = n_link->b[dir_v0]->l[k];
	if(n_link->v[dir_v1]->num > n_link->v[dir_v2]->num)
	{
		For(k,tree->n_l){
			l1[k] = n_link->b[dir_v2]->l[k];
			l2[k] = n_link->b[dir_v1]->l[k];
		}
	}
	else
	{
		For(k,tree->n_l){
			l1[k] = n_link->b[dir_v1]->l[k];
			l2[k] = n_link->b[dir_v2]->l[k];
		}
	}


	For(i,tree->depth_curr_path+1) tree->spr_list[tree->size_spr_list]->path[i] = tree->curr_path[i];
	tree->spr_list[tree->size_spr_list]->depth_path    = tree->depth_curr_path;
	tree->spr_list[tree->size_spr_list]->pars          = tree->c_pars;
	tree->spr_list[tree->size_spr_list]->lnL           = tree->c_lnL;
	tree->spr_list[tree->size_spr_list]->b_target      = b_target;
	tree->spr_list[tree->size_spr_list]->n_link        = n_link;
	tree->spr_list[tree->size_spr_list]->n_opp_to_link = (n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left);
	tree->spr_list[tree->size_spr_list]->b_opp_to_link = b_arrow;
	For(k,tree->n_l){
		tree->spr_list[tree->size_spr_list]->l0[k]            = l0[k];
		tree->spr_list[tree->size_spr_list]->l1[k]            = l1[k];
		tree->spr_list[tree->size_spr_list]->l2[k]            = l2[k];
	}
	tree->spr_list[tree->size_spr_list]->dist          = b_target->topo_dist_btw_edges;

	Include_One_Spr_To_List_Of_Spr(tree->spr_list[tree->size_spr_list],tree);
	//JSJ: temp fixes to l
	For(k,tree->n_l){
		b_target->l[k]   = init_target_len[k];
		b_arrow->l[k]    = init_arrow_len[k];
		b_residual->l[k] = init_residual_len[k];
	}

	Prune_Subtree(n_link,
			(n_link==b_arrow->left)?(b_arrow->rght):(b_arrow->left),
					&b_target,
					&b_residual,
					tree);

	if(tree->mod->s_opt->spr_lnL) Update_PMat_At_Given_Edge(b_target,tree);

	tree->c_lnL   = init_lnL;
	tree->c_pars  = init_pars;

	return move_lnL;
}

/*********************************************************/

void Speed_Spr_Loop(arbre *tree)
{
	m3ldbl lk_old;
	int init_thresh;

	init_thresh                      = tree->mod->s_opt->pars_thresh;
	tree->best_pars                  = 1E+8;
	tree->mod->s_opt->spr_lnL        = 0;
	tree->mod->s_opt->spr_pars       = 0;
	tree->mod->s_opt->quickdirty     = 0;

	if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n. Maximizing likelihood (using SPR moves)...\n");


	Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	tree->best_lnL = tree->c_lnL;

	/*****************************/
	lk_old = UNLIKELY;
	tree->mod->s_opt->max_depth_path = 2*tree->n_otu-3;
	tree->mod->s_opt->spr_lnL        = 0;
	do
	{
		lk_old = tree->c_lnL;
		Speed_Spr(tree,1);
		if(tree->n_improvements) Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
		if((!tree->n_improvements) || (fabs(lk_old-tree->c_lnL) < 1.)) break;
	}
	while(1);
	/*****************************/



	/*****************************/
	if(tree->mod->datatype == NT)
	{
		lk_old = UNLIKELY;
		tree->mod->s_opt->max_depth_path = 10;
		tree->mod->s_opt->spr_lnL        = 1;
		do
		{
			lk_old = tree->c_lnL;
			Speed_Spr(tree,1);
			if(tree->n_improvements) Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
			if((!tree->n_improvements) || (fabs(lk_old-tree->c_lnL) < 1.)) break;
		}
		while(1);
	}
	/*****************************/



	/*****************************/
	lk_old = UNLIKELY;
	do
	{
		lk_old = tree->c_lnL;
		Simu(tree,10);
		Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	}
	while(fabs(lk_old - tree->c_lnL) > tree->mod->s_opt->min_diff_lk_global);
	/*****************************/

	/*****************************/
	do
	{
		if(!Check_NNI_Five_Branches(tree)) break;
	}while(1);
	/*****************************/

	if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n");

}
/*********************************************************/

void Speed_Spr(arbre *tree, int max_cycles)
{
	int step,old_pars;
	m3ldbl old_lnL;

	if(tree->lock_topo)
	{
		PhyML_Printf("\n. The tree topology is locked.");
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}


	tree->both_sides = 1;
	Pars(tree);
	Lk(tree);
	Record_Br_Len(NULL,tree);

	tree->mod->s_opt->deepest_path  = 0;
	tree->best_pars                 = tree->c_pars;
	tree->best_lnL                  = tree->c_lnL;
	old_lnL                         = tree->c_lnL;
	old_pars                        = tree->c_pars;
	step                            = 0;

	//JSJ: there is a bug somewhere between here...
	do
	{
		++step;

		old_lnL  = tree->c_lnL;
		old_pars = tree->c_pars;

		tree->n_improvements         = 0;
		tree->perform_spr_right_away = 1;
		Spr(UNLIKELY,tree);
		//JSJ: And here..... so it is probably in SPR...
		if(!tree->mod->s_opt->spr_pars)
		{

			/* Optimise branch lengths */
			Optimize_Br_Len_Serie(tree->noeud[0],
					tree->noeud[0]->v[0],
					tree->noeud[0]->b[0],
					tree,
					tree->data);

			/* Update partial likelihoods */
			tree->both_sides = 1;
			Lk(tree);

			/* Print log-likelihood and parsimony scores */
			if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Branch lengths     ]");
		}
		else
		{
			if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Pars(tree);
		}

		Pars(tree);

		if(tree->io->print_trace)
		{
			PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,Write_Tree(tree)); fflush(tree->io->fp_out_trace);
			if((tree->io->print_site_lnl) && (!tree->mod->s_opt->spr_pars)) Print_Site_Lk(tree,tree->io->fp_out_lk); fflush(tree->io->fp_out_lk);
		}


		/* Record the current best log-likelihood and parsimony */
		tree->best_lnL  = tree->c_lnL;
		tree->best_pars = tree->c_pars;

		if(!tree->mod->s_opt->spr_pars)
		{
			if(tree->c_lnL < old_lnL-tree->mod->s_opt->min_diff_lk_local)
			{
				PhyML_Printf("\n. old_lnL = %f c_lnL = %f",old_lnL,tree->c_lnL);
				PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
				Warn_And_Exit("");
			}
		}
		else
		{
			if(tree->c_pars > old_pars)
			{
				PhyML_Printf("\n. old_pars = %d c_pars = %d",old_pars,tree->c_pars);
				PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
				Warn_And_Exit("");
			}
		}

		/* Record the current best branch lengths  */
		Record_Br_Len(NULL,tree);

		/* Exit if no improvements after complete optimization */
		if(step+1 > max_cycles) break;
		if((!tree->mod->s_opt->spr_pars) && (fabs(old_lnL-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global)) break;
		/*       if(( tree->mod->s_opt->spr_pars) && (fabs(old_pars-tree->c_pars) < 1.)) break; */
		if(!tree->n_improvements) break;

	}while(1);
}

/*********************************************************/

int Evaluate_List_Of_Regraft_Pos_Triple(spr **spr_list, int list_size, arbre *tree)
{
	spr *move;
	edge *init_target, *b_residual;
	int i,j,k,best_move;
	int dir_v0, dir_v1, dir_v2;
	m3ldbl recorded_l[MAX_BL_SET];
	m3ldbl best_lnL,init_lnL,delta_lnL;
	m3ldbl max_improv;

	best_lnL = UNLIKELY;
	delta_lnL = 0.0;
	max_improv = 0.0;
	init_target = b_residual = NULL;
	best_move = -1;
	init_lnL = tree->c_lnL;

	if(!list_size)
	{
		printf("\n. List size is 0 !");
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	For(k,tree->n_l) recorded_l[k] = -1.0;
	For(i,list_size)
	{
		move = spr_list[i];

		if(!move)
		{
			printf("\n. move is NULL\n");
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}

		if(move->b_target)
		{

			/* Record edge lengths */
			Record_Br_Len(NULL,tree);

			/* Prune subtree */
			Prune_Subtree(move->n_link,move->n_opp_to_link,&init_target,&b_residual,tree);
			For(k,tree->n_l){
				if(recorded_l[k] < 0.0)
				{
					/* Rough optimization of the branch length at prune site
					 * We only need to perform this optimization for the first
					 * element of spr_list because the pruned subtree is the
					 * same across all the elements of spr_list. It would not
					 * be true in the general case
					 */
					if(k == 0)Fast_Br_Len(init_target,tree,0);

					/* Record branch length at prune site */ //JSJ: temp fixes to l
					move->init_target_l[k] = init_target->l[k];
					recorded_l[k]          = init_target->l[k];
				}
				else
				{
					init_target->l[k]      = recorded_l[k];
					move->init_target_l[k] = recorded_l[k];
				}
			}

			/* Update the change proba matrix at prune position */
			Update_PMat_At_Given_Edge(init_target,tree);

			/* Update conditional likelihoods along the path from the prune to
	     the regraft position */
			Update_P_Lk_Along_A_Path(move->path,move->depth_path+1,tree);

			/* Regraft subtree */
			Graft_Subtree(move->b_target,move->n_link,b_residual,tree);


			dir_v1 = dir_v2 = dir_v0 = -1;
			For(j,3)
			{
				if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
				else if(dir_v1 < 0)                           dir_v1 = j;
				else                                          dir_v2 = j;
			}
			//JSJ: temp fixes to l
			For(k,tree->n_l){
				move->n_link->b[dir_v0]->l[k] = move->l0[k];

				if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
				{
					move->n_link->b[dir_v2]->l[k] = move->l1[k];
					move->n_link->b[dir_v1]->l[k] = move->l2[k];
				}
				else
				{
					move->n_link->b[dir_v1]->l[k] = move->l1[k];
					move->n_link->b[dir_v2]->l[k] = move->l2[k];
				}
			}

			/* 	  if(!tree->mod->s_opt->spr_lnL) */
			/* 	    { */
			/* 	      Update_PMat_At_Given_Edge(move->b_target,tree); */
			/* 	      Update_PMat_At_Given_Edge(b_residual,tree); */
			/* 	      Update_P_Lk(tree,move->b_opp_to_link,move->n_link); */
			/* 	      move->lnL = Lk_At_Given_Edge(move->b_opp_to_link,tree); */
			/* 	    } */
			/* 	  else */
			/* 	    { */
			move->lnL = Triple_Dist(move->n_link,tree,-1);
			/* 	    } */

			delta_lnL = 0.0;
			if((move->lnL < best_lnL) && (move->lnL > best_lnL - tree->mod->s_opt->max_delta_lnL_spr))
			{
				/* Estimate the three edge lengths at the regraft site */
				delta_lnL = Triple_Dist(move->n_link,tree,0) - move->lnL;
				move->lnL += delta_lnL;
			}
			//JSJ: temp fixes to l
			/* Record updated branch lengths for this move */
			For(k,tree->n_l){
				move->l0[k] = move->n_link->b[dir_v0]->l[k];

				if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
				{
					move->l1[k] = move->n_link->b[dir_v2]->l[k];
					move->l2[k] = move->n_link->b[dir_v1]->l[k];
				}
				else
				{
					move->l1[k] = move->n_link->b[dir_v1]->l[k];
					move->l2[k] = move->n_link->b[dir_v2]->l[k];
				}
			}

			if(move->lnL > best_lnL + tree->mod->s_opt->min_diff_lk_move)
			{
				best_lnL  = move->lnL;
				best_move = i;
				if(delta_lnL > max_improv) max_improv = delta_lnL;
			}


			/* Regraft the subtree at its original position */
			Prune_Subtree(move->n_link,
					move->n_opp_to_link,
					&move->b_target,
					&b_residual,
					tree);

			Graft_Subtree(init_target,
					move->n_link,
					b_residual,
					tree);

			/* Restore branch lengths */
			Restore_Br_Len(NULL,tree);

			/* Update relevant change proba matrices */
			/* 	  Update_PMat_At_Given_Edge(move->n_link->b[0],tree); */
			/* 	  Update_PMat_At_Given_Edge(move->n_link->b[1],tree); */
			/* 	  Update_PMat_At_Given_Edge(move->n_link->b[2],tree); */
			Update_PMat_At_Given_Edge(move->b_target,tree);

			/* 	  Update_P_Lk(tree,move->n_link->b[0],move->n_link); */
			/* 	  Update_P_Lk(tree,move->n_link->b[1],move->n_link); */
			/* 	  Update_P_Lk(tree,move->n_link->b[2],move->n_link); */

			/* Update conditional likelihoods along the path from the prune to
	     the regraft position */
			/* 	  Update_P_Lk_Along_A_Path(move->path,move->depth_path+1,tree); */

			tree->c_lnL = init_lnL;
		}

		/* Bail out as soon as you've found a true improvement */
		if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) break;
	}

	/*   printf("\n. [ %4d/%4d ]",i,list_size); */

	/*   printf("\n. max_improv = %f",max_improv); */

	For(i,list_size)
	{
		move = spr_list[i];
		if(move->b_target)
		{
			For(j,3) Update_PMat_At_Given_Edge(move->n_link->b[j],tree);
			For(j,3) Update_P_Lk(tree,move->n_link->b[j],move->n_link);

			/* TO DO : we don't need to update all these partial likelihoods here.
	     Would need to record only those that were along the paths examined
	     above */

			For(j,3)
			if(move->n_link->v[j] != move->n_opp_to_link)
				Pre_Order_Lk(move->n_link,move->n_link->v[j],tree);
			break;
		}
	}


#ifdef DEBUG
	if(best_move < 0)
	{
		PhyML_Printf("\n\n. Best_move < 0 !");

		PhyML_Printf("\n. List size = %d",list_size);
		For(i,list_size)
		{
			move = spr_list[i];
			PhyML_Printf("\n. %p %p",move,move->b_target);
		}

		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
#endif


	return best_move;
}

/*********************************************************/

int Try_One_Spr_Move_Triple(spr *move, arbre *tree)
{
	edge *init_target, *b_residual;
	int j,k;
	int dir_v0, dir_v1, dir_v2;


	Record_Br_Len(NULL,tree);

	Prune_Subtree(move->n_link,
			move->n_opp_to_link,
			&init_target,
			&b_residual,
			tree);
	//JSJ: temp fix to l
	For(k,tree->n_l){
		init_target->l[k] = move->init_target_l[k];
	}

	Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

	dir_v1 = dir_v2 = dir_v0 = -1;
	For(j,3)
	{
		if(move->n_link->v[j] == move->n_opp_to_link) dir_v0 = j;
		else if(dir_v1 < 0)                           dir_v1 = j;
		else                                          dir_v2 = j;
	}
	//JSJ temp fixes to l
	For(k,tree->n_l){
		move->n_link->b[dir_v0]->l[k] = move->l0[k];

		if(move->n_link->v[dir_v1]->num > move->n_link->v[dir_v2]->num)
		{
			move->n_link->b[dir_v2]->l[k] = move->l1[k];
			move->n_link->b[dir_v1]->l[k] = move->l2[k];
		}
		else
		{
			move->n_link->b[dir_v1]->l[k] = move->l1[k];
			move->n_link->b[dir_v2]->l[k] = move->l2[k];
		}
	}

	if(move->lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move) /* Apply the move */
	{
		time(&(tree->t_current));
		tree->both_sides = 1;
		Lk(tree);
		Pars(tree);

		if(fabs(tree->c_lnL - move->lnL) > 1.E-3)
		{
			if(tree->mod->s_opt->print) PhyML_Printf("\n. c_lnL = %f move_lnL = %f",
					tree->c_lnL,move->lnL);
			PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
			Warn_And_Exit("");
		}

		if((tree->mod->s_opt->print) && (!tree->io->quiet))
		{
			Print_Lk_And_Pars(tree);
			printf(" [depth=%5d]",move->depth_path); fflush(NULL);
		}

		tree->n_improvements++;
		tree->best_lnL = tree->c_lnL;
		Record_Br_Len(NULL,tree);

		if(move->depth_path > tree->mod->s_opt->deepest_path)
			tree->mod->s_opt->deepest_path = move->depth_path;

		return 1;
	}

	Prune_Subtree(move->n_link,
			move->n_opp_to_link,
			&move->b_target,
			&b_residual,
			tree);

	Graft_Subtree(init_target,
			move->n_link,
			b_residual,
			tree);

	Restore_Br_Len(NULL,tree);
	tree->both_sides = 1;
	Lk(tree);
	Pars(tree);
	return 0;
}

/*********************************************************/

int Try_One_Spr_Move_Full(spr *move, arbre *tree)
{
	edge *init_target, *b_residual;

	Record_Br_Len(NULL,tree);

	Prune_Subtree(move->n_link,
			move->n_opp_to_link,
			&init_target,
			&b_residual,
			tree);

	Graft_Subtree(move->b_target,move->n_link,b_residual,tree);

	tree->both_sides = 1;
	Lk(tree);

	Optimize_Br_Len_Serie(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],
			tree,
			tree->data);

	tree->both_sides = 1;
	Lk(tree);

	if(tree->c_lnL > tree->best_lnL + tree->mod->s_opt->min_diff_lk_move)
	{
		Pars(tree);
		if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Topology           ]");
		tree->n_improvements++;
		tree->best_lnL = tree->c_lnL;
		Record_Br_Len(NULL,tree);
		return 1;
	}
	else
	{
		Prune_Subtree(move->n_link,
				move->n_opp_to_link,
				&move->b_target,
				&b_residual,
				tree);

		Graft_Subtree(init_target,
				move->n_link,
				b_residual,
				tree);

		Restore_Br_Len(NULL,tree);
		tree->both_sides = 1;
		Lk(tree);
		Pars(tree);
		return 0;
	}

	return -1;
}

/*********************************************************/

void Include_One_Spr_To_List_Of_Spr(spr *move, arbre *tree)
{
	int i,j;
	spr *buff_spr;

	if((( tree->mod->s_opt->spr_lnL) && (move->lnL  > tree->spr_list[tree->size_spr_list-1]->lnL)) ||
			((!tree->mod->s_opt->spr_lnL) && (move->pars <= tree->spr_list[tree->size_spr_list-1]->pars)))
	{
		tree->spr_list[tree->size_spr_list-1]->depth_path    = move->depth_path;
		tree->spr_list[tree->size_spr_list-1]->pars          = move->pars;
		tree->spr_list[tree->size_spr_list-1]->lnL           = move->lnL;
		tree->spr_list[tree->size_spr_list-1]->b_target      = move->b_target;
		tree->spr_list[tree->size_spr_list-1]->n_link        = move->n_link;
		tree->spr_list[tree->size_spr_list-1]->n_opp_to_link = move->n_opp_to_link;
		tree->spr_list[tree->size_spr_list-1]->b_opp_to_link = move->b_opp_to_link;
		For(j,tree->n_l){
			tree->spr_list[tree->size_spr_list-1]->l0[j]            = move->l0[j]; //JSJ: temp fixes to l
			tree->spr_list[tree->size_spr_list-1]->l1[j]            = move->l1[j];
			tree->spr_list[tree->size_spr_list-1]->l2[j]            = move->l2[j];
		}
		tree->spr_list[tree->size_spr_list-1]->dist          = move->dist;

		For(i,tree->spr_list[tree->size_spr_list-1]->depth_path+1)
		tree->spr_list[tree->size_spr_list-1]->path[i] = move->path[i];

		for(i=tree->size_spr_list-1;i>0;i--)
		{
			if((( tree->mod->s_opt->spr_lnL) && (tree->spr_list[i]->lnL > tree->spr_list[i-1]->lnL)) ||
					((!tree->mod->s_opt->spr_lnL) && (tree->spr_list[i]->pars <=  tree->spr_list[i-1]->pars)))
			{
				buff_spr            = tree->spr_list[i-1];
				tree->spr_list[i-1] = tree->spr_list[i];
				tree->spr_list[i]   = buff_spr;
			}
			else  break;
		}
	}
}

/*********************************************************/

void Random_Spr(int n_moves, arbre *tree)
{
	int i;
	int br_pulled, br_target;
	spr *spr_struct;
	edge *target, *residual;

	spr_struct = Make_One_Spr(tree);
	Init_One_Spr(spr_struct,tree->n_l);
	target = residual = NULL;
	For(i,n_moves)
	{
		br_pulled = (int)((m3ldbl)rand()/RAND_MAX * (2*tree->n_otu-3-1));
		do
		{
			br_target = (int)((m3ldbl)rand()/RAND_MAX * (2*tree->n_otu-3-1));
		}while(br_target == br_pulled);

		spr_struct->n_link        = tree->t_edges[br_pulled]->left;
		spr_struct->n_opp_to_link = tree->t_edges[br_pulled]->rght;
		spr_struct->b_opp_to_link = tree->t_edges[br_pulled];
		spr_struct->b_target      = tree->t_edges[br_target];
		spr_struct->b_init_target = NULL;

		if(!Check_Spr_Move_Validity(spr_struct,tree))
		{
			spr_struct->n_link        = tree->t_edges[br_pulled]->rght;
			spr_struct->n_opp_to_link = tree->t_edges[br_pulled]->left;
		}

#ifdef DEBUG
		if(!Check_Spr_Move_Validity(spr_struct,tree))
		{
			Warn_And_Exit("\n. Could not find a valid move...\n");
		}
#endif

		Prune_Subtree(spr_struct->n_link,
				spr_struct->n_opp_to_link,
				&target,
				&residual,
				tree);

		Graft_Subtree(spr_struct->b_target,
				spr_struct->n_link,
				residual,tree);
	}
	Free(spr_struct);
}

/*********************************************************/

void Reset_Spr_List(arbre *tree)
{
	int i;

	For(i,tree->size_spr_list)
	{
		tree->spr_list[i]->depth_path     = 0;
		tree->spr_list[i]->pars           = MAX_PARS;
		tree->spr_list[i]->lnL            = UNLIKELY;
		tree->spr_list[i]->n_link         = NULL;
		tree->spr_list[i]->n_opp_to_link  = NULL;
		tree->spr_list[i]->b_target       = NULL;
	}
}

/*********************************************************/

int Check_Spr_Move_Validity(spr *this_spr_move, arbre *tree)
{
	int match;

	match = 0;
	Found_In_Subtree(this_spr_move->n_link,
			this_spr_move->n_opp_to_link,
			this_spr_move->b_target->left,
			&match,
			tree);

	if(match) return 0;
	else      return 1;
}

/*********************************************************/

void Make_Best_Spr(arbre *tree)
{
	tree->best_spr = Make_One_Spr(tree);
	Init_One_Spr(tree->best_spr,tree->n_l);
}

/*********************************************************/

void Make_Spr_List(arbre *tree)
{
	int i;

	tree->size_spr_list = 2*tree->n_otu-3;
	tree->spr_list = (spr **)mCalloc(2*tree->n_otu-2,sizeof(spr *));

	For(i,2*tree->n_otu-2)
	{
		tree->spr_list[i] = Make_One_Spr(tree);
		Init_One_Spr(tree->spr_list[i],tree->n_l);
	}
	tree->perform_spr_right_away = 0;
}

/*********************************************************/

void Init_One_Spr(spr *a_spr, int n_l)
{
	a_spr->lnL             = UNLIKELY;
	a_spr->pars            = 1E+5;
	a_spr->depth_path      = 0;
	a_spr->dist            = 0;
	a_spr->n_l			   = n_l; //JSJ: store the number of rate categories
	a_spr->n_link          = NULL;
	a_spr->n_opp_to_link   = NULL;
	a_spr->b_opp_to_link   = NULL;
	a_spr->b_target        = NULL;
	a_spr->b_init_target   = NULL;
	/**
	* JSJ: Initialize our arrays of branch lengths
	*/
	int i;
	For(i,n_l){
		a_spr->init_target_l[i]   = -1.;
		a_spr->l0[i]              = -1.;
		a_spr->l1[i]              = -1.;
		a_spr->l2[i]              = -1.;
	}
}

/*********************************************************/

spr *Make_One_Spr(arbre *tree)
{
	spr *a_spr;
	a_spr       = (spr *)mCalloc(1,sizeof(spr));
	a_spr->path = (node **)mCalloc(tree->n_otu,sizeof(node *));
	return a_spr;
}

/*********************************************************/

void Spr_Pars(arbre *tree)
{

	PhyML_Printf("\n. Minimizing parsimony...\n");

	tree->best_pars = 1E+8;
	tree->best_lnL  = UNLIKELY;
	tree->mod->s_opt->spr_lnL = 0;

	tree->mod->s_opt->spr_pars = 1;
	do
	{
		Speed_Spr(tree,1);
	}while(tree->n_improvements);
	tree->mod->s_opt->spr_pars = 0;

	PhyML_Printf("\n");
}

/*********************************************************/














/*
 ** EOF: spr.c
 */
