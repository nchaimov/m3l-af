/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

 */

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "m4.h"
#include "mc.h"

#ifdef COMPRESS_SUBALIGNMENTS
#include "compress.h"
#endif

#ifdef MC
#include "rates.h"
#endif



/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, plkflt *p_lk)
{
	switch(state)
	{
	case 'A' : p_lk[pos+0]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=.0;
	break;
	case 'C' : p_lk[pos+1]=1.; p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=.0;
	break;
	case 'G' : p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+0]=p_lk[pos+3]=.0;
	break;
	case 'T' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
	break;
	case 'U' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
	break;
	case 'M' : p_lk[pos+0]=p_lk[pos+1]=1.; p_lk[pos+2]=p_lk[pos+3]=.0;
	break;
	case 'R' : p_lk[pos+0]=p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+3]=.0;
	break;
	case 'W' : p_lk[pos+0]=p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=.0;
	break;
	case 'S' : p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+0]=p_lk[pos+3]=.0;
	break;
	case 'Y' : p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+2]=.0;
	break;
	case 'K' : p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+1]=.0;
	break;
	case 'B' : p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=.0;
	break;
	case 'D' : p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+1]=.0;
	break;
	case 'H' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+2]=.0;
	break;
	case 'V' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+3]=.0;
	break;
	case 'N' : case 'X' : case '?' : case 'O' : case '-' :
		p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.;break;
	default :
	{
		PhyML_Printf("\n. Unknown character state : %c\n",state);
		Exit("\n. Init failed (check the data type)\n");
		break;
	}
	}
}

/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars)
{
	switch(state)
	{
	case 'A' : p_pars[pos+0]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=0;
	break;
	case 'C' : p_pars[pos+1]=1; p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=0;
	break;
	case 'G' : p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+0]=p_pars[pos+3]=0;
	break;
	case 'T' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
	break;
	case 'U' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
	break;
	case 'M' : p_pars[pos+0]=p_pars[pos+1]=1; p_pars[pos+2]=p_pars[pos+3]=0;
	break;
	case 'R' : p_pars[pos+0]=p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+3]=0;
	break;
	case 'W' : p_pars[pos+0]=p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=0;
	break;
	case 'S' : p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+0]=p_pars[pos+3]=0;
	break;
	case 'Y' : p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+2]=0;
	break;
	case 'K' : p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+1]=0;
	break;
	case 'B' : p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=0;
	break;
	case 'D' : p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+1]=0;
	break;
	case 'H' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+2]=0;
	break;
	case 'V' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+3]=0;
	break;
	case 'N' : case 'X' : case '?' : case 'O' : case '-' :
		p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1;break;
	default :
	{
		PhyML_Printf("\n. Unknown character state : %c\n",state);
		Exit("\n. Init failed (check the data type)\n");
		break;
	}
	}
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Float(char aa, int pos, plkflt *p_lk)
{
	int i;

	For(i,20) p_lk[pos+i] = .0;

	switch(aa){
	case 'A' : p_lk[pos+0]= 1.; break;/* Alanine */
	case 'R' : p_lk[pos+1]= 1.; break;/* Arginine */
	case 'N' : p_lk[pos+2]= 1.; break;/* Asparagine */
	case 'D' : p_lk[pos+3]= 1.; break;/* Aspartic acid */
	case 'C' : p_lk[pos+4]= 1.; break;/* Cysteine */
	case 'Q' : p_lk[pos+5]= 1.; break;/* Glutamine */
	case 'E' : p_lk[pos+6]= 1.; break;/* Glutamic acid */
	case 'G' : p_lk[pos+7]= 1.; break;/* Glycine */
	case 'H' : p_lk[pos+8]= 1.; break;/* Histidine */
	case 'I' : p_lk[pos+9]= 1.; break;/* Isoleucine */
	case 'L' : p_lk[pos+10]=1.; break;/* Leucine */
	case 'K' : p_lk[pos+11]=1.; break;/* Lysine */
	case 'M' : p_lk[pos+12]=1.; break;/* Methionine */
	case 'F' : p_lk[pos+13]=1.; break;/* Phenylalanin */
	case 'P' : p_lk[pos+14]=1.; break;/* Proline */
	case 'S' : p_lk[pos+15]=1.; break;/* Serine */
	case 'T' : p_lk[pos+16]=1.; break;/* Threonine */
	case 'W' : p_lk[pos+17]=1.; break;/* Tryptophan */
	case 'Y' : p_lk[pos+18]=1.; break;/* Tyrosine */
	case 'V' : p_lk[pos+19]=1.; break;/* Valine */

	case 'B' : p_lk[pos+2]= 1.; break;/* Asparagine */
	case 'Z' : p_lk[pos+5]= 1.; break;/* Glutamine */

	case 'X' : case '?' : case '-' : For(i,20) p_lk[pos+i] = 1.; break;
	default :
	{
		PhyML_Printf("\n. Unknown character state : %c\n",aa);
		Exit("\n. Init failed (check the data type)\n");
		break;
	}
	}
}

/*********************************************************/

void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars)
{

	int i;

	For(i,20) p_pars[pos+i] = .0;

	switch(aa){
	case 'A' : p_pars[pos+0]  = 1; break;/* Alanine */
	case 'R' : p_pars[pos+1]  = 1; break;/* Arginine */
	case 'N' : p_pars[pos+2]  = 1; break;/* Asparagine */
	case 'D' : p_pars[pos+3]  = 1; break;/* Aspartic acid */
	case 'C' : p_pars[pos+4]  = 1; break;/* Cysteine */
	case 'Q' : p_pars[pos+5]  = 1; break;/* Glutamine */
	case 'E' : p_pars[pos+6]  = 1; break;/* Glutamic acid */
	case 'G' : p_pars[pos+7]  = 1; break;/* Glycine */
	case 'H' : p_pars[pos+8]  = 1; break;/* Histidine */
	case 'I' : p_pars[pos+9]  = 1; break;/* Isoleucine */
	case 'L' : p_pars[pos+10] = 1; break;/* Leucine */
	case 'K' : p_pars[pos+11] = 1; break;/* Lysine */
	case 'M' : p_pars[pos+12] = 1; break;/* Methionine */
	case 'F' : p_pars[pos+13] = 1; break;/* Phenylalanin */
	case 'P' : p_pars[pos+14] = 1; break;/* Proline */
	case 'S' : p_pars[pos+15] = 1; break;/* Serine */
	case 'T' : p_pars[pos+16] = 1; break;/* Threonine */
	case 'W' : p_pars[pos+17] = 1; break;/* Tryptophan */
	case 'Y' : p_pars[pos+18] = 1; break;/* Tyrosine */
	case 'V' : p_pars[pos+19] = 1; break;/* Valine */

	case 'B' : p_pars[pos+2]  = 1; break;/* Asparagine */
	case 'Z' : p_pars[pos+5]  = 1; break;/* Glutamine */

	case 'X' : case '?' : case '-' : For(i,20) p_pars[pos+i] = 1; break;
	default :
	{
		PhyML_Printf("\n. Unknown character state : %c\n",aa);
		Exit("\n. Init failed (check the data type)\n");
		break;
	}
	}
}

/*********************************************************/
void Get_All_Partial_Lk_Scale(arbre *tree, edge *b_fcus, node *d)
{
	if(d->tax) return;
	else Update_P_Lk(tree,b_fcus,d);
}


/*********************************************************/
void Post_Order_Lk(node *a, node *d, arbre *tree)
{
	//PhyML_Printf("Post_Order_Lk(a = %d, d = %d)\n", a->num, d->num);
	int i,dir;

	dir = -1;

	if(d->tax) return;
	else
	{
		For(i,3)
		{
			if(d->v[i] != a)
				Post_Order_Lk(d,d->v[i],tree);
			else dir = i;
		}
		Get_All_Partial_Lk_Scale(tree,d->b[dir],d);
	}
}

/*********************************************************/
void Pre_Order_Lk(node *a, node *d, arbre *tree)
{
	//PhyML_Printf("Pre_Order_Lk(a = %d, d = %d)\n", a->num, d->num);

	int i;

	if(d->tax) return;
	else
	{
		For(i,3)
		{
			if(d->v[i] != a)
			{
				Get_All_Partial_Lk_Scale(tree,d->b[i],d);
				Pre_Order_Lk(d,d->v[i],tree);
			}
		}
	}
}

//
// THis is basically the Lk method without the COMPRESS_SUBALIGNMENT pragmas.
//
// It is for debugging tests, and should be removed from the production code.
//
void debug_Lk_nocompress(arbre *tree)
{
	int br,site;
	int n_patterns;

	//PhyML_Printf(" . debug: entered debug_Lk_nocompress\n");
	//Print_Tree_Screen(tree);

	n_patterns = tree->n_pattern;
	tree->number_of_lk_calls++;
	Set_Model_Parameters(tree->mod);

#ifdef MC
	if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Get_Br_Len(tree);
	if(tree->bl_from_node_stamps) MC_Bl_From_T(tree);
#endif

	//	chunk = (2*tree->n_otu-3)/omp_get_num_procs();
	//	printf("chunk: %i total: %i\n",chunk,(2*tree->n_otu-3));

	// #pragma omp parallel
	//for shared(tree,n_patterns,chunk) schedule(static,chunk)
	for(br=0; br < 2*tree->n_otu-3; br++)
	{
		if(!tree->t_edges[br]->rght->tax)
			For(site,n_patterns) tree->t_edges[br]->sum_scale_f_rght[site] = .0;

		if(!tree->t_edges[br]->left->tax)
			For(site,n_patterns) tree->t_edges[br]->sum_scale_f_left[site] = .0;

		Update_PMat_At_Given_Edge(tree->t_edges[br],tree);
	}

	// VHS: the post- and pre-order traversals begin at a random node (the 0th node, to be specific).
	// This random strategy is okay, because the Felsenstein pruning algorithm operates on unrooted trees.
	Post_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
	if(tree->both_sides)
	{
		Pre_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
	}

	// 10.11.2009: here.
	return;

	tree->c_lnL     = .0;
	tree->curr_catg =  0;
	tree->curr_site =  0;
#ifdef USE_OPENMP
	int chunk = n_patterns/omp_get_num_procs();
#pragma omp parallel for \
		shared(tree,n_patterns,chunk) private(site)\
		schedule(static,chunk)
#endif
	for(site = 0; site < n_patterns; site++)
	{
		tree->c_lnL_sorted[site] = .0;
		tree->site_lk[site]      = .0;
		tree->curr_site          = site;
		Site_Lk(tree,site);
	}

	/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

	tree->c_lnL = .0;

	//#pragma omp parallel
	//for default(shared)
	//private(site) schedule(static,chunk)
	For(site,n_patterns)
	{
		if(tree->c_lnL_sorted[site] < .0) /* WARNING : change cautiously */
		{
			//PhyML_Printf(" . debug: lk.c line 372: tree->c_lnL_sorted[%d] = %f\n", site, tree->c_lnL_sorted[site]);
			tree->c_lnL += tree->c_lnL_sorted[site];
		}
	}

	// debug:
	//PhyML_Printf(" . debug: Lk is returning %f\n", tree->c_lnL);
}

/*********************************************************/

void Lk(arbre *tree)
{
	int br,site;
	int n_patterns;

	n_patterns = tree->n_pattern;
	tree->number_of_lk_calls++;
	Set_Model_Parameters(tree->mod);

#ifdef MC
	if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Get_Br_Len(tree);
	if(tree->bl_from_node_stamps) MC_Bl_From_T(tree);
#endif

	int chunk = (2*tree->n_otu-3)/omp_get_num_procs();
	//	printf("chunk: %i total: %i\n",chunk,(2*tree->n_otu-3));

#pragma omp parallel for shared(tree,n_patterns,chunk) schedule(static,chunk)
	for(br=0; br < 2*tree->n_otu-3; br++)
	{
		if(!tree->t_edges[br]->rght->tax)
			For(site,n_patterns) tree->t_edges[br]->sum_scale_f_rght[site] = .0;

		if(!tree->t_edges[br]->left->tax)
			For(site,n_patterns) tree->t_edges[br]->sum_scale_f_left[site] = .0;

		Update_PMat_At_Given_Edge(tree->t_edges[br],tree);

		// debug:
		//PhyML_Printf(" . debug: P matrix at edge %d:\n", tree->t_edges[br]->num);
		//Print_Pij( tree->t_edges[br]->Pij_rr, tree->mod);
	}


#ifdef COMPRESS_SUBALIGNMENTS
	//PhyML_Printf("\n. Compressing sub-alignments based on new phylogeny...\n");
	Init_All_Nodes_Red(tree);
	Compute_Red_Arrays(tree->noeud[0], tree->noeud[0]->v[0], tree);
	tree->red_arrays_invalid = 0; // now the red arrays are valid again.

	for(site = 0; site < n_patterns; site++)
	{
		Post_Order_Lk_Red(tree->noeud[0],tree->noeud[0]->v[0],tree,site);
	}
	if(tree->both_sides)
	{
		for(site = 0; site < n_patterns; site++)
		{
			Pre_Order_Lk_Red(tree->noeud[0],tree->noeud[0]->v[0],tree,site);
		}
	}
#endif

#ifndef COMPRESS_SUBALIGNMENTS
	// VHS: the post- and pre-order traversals begin at a random node (the 0th node, to be specific).
	// This random strategy is okay, because the Felsenstein pruning algorithm operates on unrooted trees.
	Post_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
	if(tree->both_sides)
	{
		Pre_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
	}
#endif

	  // debugging:

	  //PhyML_Printf(" . debug: returned from Post_Order_Lk and Pre_Order_Lk, likelihoods = \n");
	  /*
	  for(br = 1; br < 2*tree->n_otu-3; br++)
	  {
		PhyML_Printf("edge %d:\n", br);
		int is_tax = 1;
		if(!tree->t_edges[br]->left->tax)
			is_tax = 0;
		Print_Plk(tree->t_edges[br]->p_lk_left, tree, is_tax);
	  }
	  */


	tree->c_lnL     = .0;
	tree->curr_catg =  0;
	tree->curr_site =  0;
#ifdef USE_OPENMP
	int chunk = n_patterns/omp_get_num_procs();
#pragma omp parallel for \
		shared(tree,n_patterns,chunk) private(site)\
		schedule(static,chunk)
#endif
	for(site = 0; site < n_patterns; site++)
	{
		// Here we reset these three variables, so that we can update them in Site_Lk
		tree->c_lnL_sorted[site] = .0;
		tree->site_lk[site]      = .0;
		tree->curr_site          = site;
		Site_Lk(tree,site);
	}

	/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

	tree->c_lnL = .0;

	//#pragma omp parallel
	//for default(shared)
	//private(site) schedule(static,chunk)
	For(site,n_patterns)
	{
		//PhyML_Printf(" . debug: in lk.c 258: c_lnL_sorted[%d] = %f\n", site, tree->c_lnL_sorted[site]);
		if(tree->c_lnL_sorted[site] < .0) /* WARNING : change cautiously */
		{
			//PhyML_Printf(" . debug: lk.c line 372: tree->c_lnL_sorted[%d] = %f\n", site, tree->c_lnL_sorted[site]);
			tree->c_lnL += tree->c_lnL_sorted[site];
		}
	}

	// debug:
	//PhyML_Printf(" . debug: Lk is returning %f\n", tree->c_lnL);
}

/*********************************************************/

// VHS: This method calculates the likelihood for the entire tree, given the site
// specified by *tree->curr_site
void Site_Lk(arbre *tree, int site)
{
	edge *eroot;

	eroot = tree->noeud[0]->b[0];

	if(!eroot->rght->tax)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	if(tree->data->wght[tree->curr_site] > MDBL_MIN)
	{	Lk_Core(eroot,tree,site);
	}
	else
	{	tree->c_lnL_sorted[tree->curr_site] = 1.; /* WARNING : change cautiously */
	}
}

/*********************************************************/
m3ldbl Lk_At_Given_Edge(edge *b_fcus, arbre *tree)
{
	int n_patterns,site;
	tree->number_of_branch_lk_calls++;
	n_patterns = tree->n_pattern;

#ifdef MC
	if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Get_Br_Len(tree);
	if(tree->bl_from_node_stamps) MC_Bl_From_T(tree);
#endif

	Update_PMat_At_Given_Edge(b_fcus,tree);

	if(b_fcus->left->tax)
	{
		PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	tree->c_lnL = .0;
#ifdef USE_OPENMP
	int chunk = n_patterns/omp_get_num_procs();
#pragma omp parallel for \
		shared(tree,n_patterns,chunk,b_fcus) private(site) \
		schedule(static,chunk)
#endif
	for(site = 0; site < n_patterns; site++)
	{
		//printf("JSJ: Iterating over state pattern %i in Lk_At_Given_Edge\n",tree->curr_site);

		// debug: x-inf

		if(tree->data->wght[site] > MDBL_MIN)
		{
			//PhyML_Printf(" . debug: lk.c 533: calling Lk_Core for edge %d, site %d\n", b_fcus->num, site);

			// debug stuff:
			/*
			int dima = tree->mod->ns;
			int dimb = dima * tree->mod->n_l;
			int dimc = dimb * tree->mod->n_catg;
	*/

			//debug:
//			if ( (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] < 0.0 || (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] > 1.0 )
//			{
//				PhyML_Printf(" . debug: lk.c 543: the p_lk_left is invalid for edge %d at site %d\n", b_fcus->num, site);
//			}


			Lk_Core(b_fcus,tree,site);
		}
		else tree->c_lnL_sorted[site] = 1.; /* WARNING : change cautiously */
		//PhyML_Printf(" . debug: lk.c 537: returned from Lk_Core: tree->c_lnL_sorted[%d] = %f\n", site, tree->c_lnL_sorted[site]);
	}

	/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

	tree->c_lnL = .0;
	For(tree->curr_site,n_patterns)
	if(tree->c_lnL_sorted[tree->curr_site] < .0) /* WARNING : change cautiously */
		tree->c_lnL += tree->c_lnL_sorted[tree->curr_site];

	return tree->c_lnL;
}

/*********************************************************/
// VHS: This method calculates the likelihood for the entire tree (we presume) rooted at edge *b,
// for the site indicated by site
m3ldbl Lk_Core(edge *b, arbre *tree, int site)
{
	/**
	* Game plan: iterate through the entirety of this function
	* on each branch length set. Don't return the lk val in each
	* iteration though. Instead store an array of those values.
	* At the end of the loop deal with reducing the array and return
	* the likelihood given the branch length set.
	*/

	//PhyML_Printf(" . debug lk.c 563: entered Lk_Core with edge %d\n", b->num);
	//PhyML_Printf("   connecting nodes %d on left and %d on right\n", b->left->num, b->rght->num);

	/**
	* Global variables (in the scope of this loop)
	*/
	m3ldbl site_lk, site_lk_cat, site_lk_set, sum_site_lk, log_site_lk;
	int i,j; //loop counters
	m3ldbl sum;
	int ambiguity_check,state;
	int catg,ns,k,l;
	float scale_left, scale_rght;

	site_lk = site_lk_cat = site_lk_set = sum_site_lk = log_site_lk = 0.0;

	int dima, dimb, dimc, dimaa, dimd;
	dima = tree->mod->ns;
	dimb = dima * tree->mod->n_l;
	dimc = dimb * tree->mod->n_catg;
	dimaa = dima * dima; // = ns^2
	dimd = dimb * dima; // = number of BL sets * ns^2

	ambiguity_check = state = -1;
	ns = tree->mod->ns;
	scale_left = (b->sum_scale_f_left)? (b->sum_scale_f_left[site]): (0.0);
	scale_rght = (b->sum_scale_f_rght)? (b->sum_scale_f_rght[site]): (0.0);

	if((b->rght->tax) && (!tree->mod->s_opt->greedy))
	{
		ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];
		if(!ambiguity_check)
		{
			state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dima,tree);
		}
	}

	if(tree->mod->use_m4mod) ambiguity_check = 1;

	For(catg,tree->mod->n_catg)
	{
		site_lk_cat = .0;

		For(i,tree->mod->n_l)
		{
			site_lk_set = .0;

			if((b->rght->tax) && (!tree->mod->s_opt->greedy))
			{
				if(!ambiguity_check)
				{
					sum = .0;
					For(l,ns)
					{
						sum +=	b->Pij_rr[catg*dimd + i*dimaa + state*dima + l] *
								(m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l];
					}
					site_lk_set += sum * tree->mod->pi[state];
					//debug:
					//PhyML_Printf(" . debug: lk 519: site_lk_set = %f\n", site_lk_set);
				}
				else
				{
					For(k,ns)
					{
						sum = .0;
						if(b->p_lk_tip_r[site*dima+k] > .0)
						{
							For(l,ns)
							{
								sum +=	b->Pij_rr[catg*dimd + i*dimaa + k*dima + l] *
										(m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l];
								// debug:
								if ( (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l] < 0.0 || (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l] > 1.0)
								{
									PhyML_Printf(" . debug: lk.c 637: in Lk_Core: problem with edge %d, site %d, catg %d, bl %d, k %d, l%d\n", b->num, site, catg, i, k, l);
								}
							}
							site_lk_set +=
									sum *
									tree->mod->pi[k] *
									(m3ldbl)b->p_lk_tip_r[site*dima+k];
							//debug:
							//PhyML_Printf(" . debug: lk 536: site_lk_set = %e\n", site_lk_set);
						}
					}
				}
			}
			else
			{
				For(k,ns)
				{
					sum = .0;
					if(b->p_lk_rght[site*dimc + catg*dimb + i*dima + k] > .0)
					{
						For(l,ns)
						{
							sum += b->Pij_rr[catg*dimd + i*dimaa + k*dima + l] *
									(m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l];
							// note on 10.11.2009: VHS: the value in (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l];
							// equals nan for edge 3, site 3 in the apgm83 test suite.
							// Why?
							//PhyML_Printf("sum = %e += %e * %e\n", sum, b->Pij_rr[catg*dimd + i*dimaa + k*dima + l], (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l]);
							// debug:
							//if ( (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l] < 0.0 || (m3ldbl)b->p_lk_left[site*dimc + catg*dimb + i*dima + l] > 1.0)
							//{
							//	PhyML_Printf(" . debug: lk.c 668: in Lk_Core: problem with edge %d, site %d, catg %d, bl %d, k %d, l%d\n", b->num, site, catg, i, k, l);
							//}
						}
						site_lk_set +=
								sum *
								tree->mod->pi[k] *
								(m3ldbl)b->p_lk_rght[site*dimc + catg*dimb + i*dima + k];
						//debug:
						//PhyML_Printf("tree->mod->pi[%d] = %d\n",k, tree->mod->pi[k]);
						//PhyML_Printf(" . debug: lk 556: site_lk_set = %f\n", site_lk_set);
					}
				}
			}
			//debug:
			//PhyML_Printf(" . debug: tree->log_site_lk_set[BL set %d][site %d] = %f\n", i, site, site_lk_set);

			site_lk_cat += site_lk_set * tree->mod->bl_props[i];
			//PhyML_Printf( " . debug: lk 672: %lf * %lf = %lf\n", site_lk_set, tree->mod->bl_props[i], site_lk_set * tree->mod->bl_props[i]);
		} // end for BL

		// debug:
		//PhyML_Printf(" . debug: tree->log_site_lk_cat[catg %d][site %d] = %f\n", catg, site, site_lk_cat);

		// we do the log-ing further down in this method.
		tree->log_site_lk_cat[catg][site] = site_lk_cat;

		site_lk += site_lk_cat * tree->mod->gamma_r_proba[catg];
	} // end for gamma

	if(site_lk < 1.E-300)
	{
		site_lk = 1.E-300;
		PhyML_Printf(" . debug: WARNING, site_lk is too small! new site_lk = %e\n", site_lk);
	}


	//debug test
	//PhyML_Printf(" lk.c 685 log(site_lk) = %e\n", log(site_lk) );


	if(!tree->mod->invar)
	{
		log_site_lk = log(site_lk) + (m3ldbl)scale_left + (m3ldbl)scale_rght;

		//PhyML_Printf(" . debug: lk.c 689: log_site_lk = %f + %f + %f  = %f\n", log(site_lk), (m3ldbl)scale_left + (m3ldbl)scale_rght, log_site_lk);
	}
	else
	{
		if((m3ldbl)tree->data->invar[site] > -0.5)
		{
			if((scale_left + scale_rght > 0.0) || (scale_left + scale_rght < 0.0))
					site_lk *= (m3ldbl)exp(scale_left + scale_rght);
			log_site_lk = (m3ldbl)log(site_lk*(1.0-tree->mod->pinvar) +
							tree->mod->pinvar*tree->mod->pi[tree->data->invar[site]]);
		}
		else
		{
			log_site_lk = (m3ldbl)log(site_lk*(1.0-tree->mod->pinvar)) + (m3ldbl)scale_left + (m3ldbl)scale_rght;
		}
	}

/*
	if(log_site_lk < -MDBL_MAX)
	{
		PhyML_Printf("\n . debug: lk.c near 740: log_site_lk = %f\n", log_site_lk);
		Warn_And_Exit("\nlog_site_lk < -MDBL_MAX\n");
	}
*/

	//PhyML_Printf(" . debug: lk.c 712: tree->c_lnL_sorted[%d] = %d * %f = %f\n", site, tree->data->wght[site], log_site_lk, tree->c_lnL_sorted[site]);
	//PhyML_Printf(" . debug: line 713: scale_left = %f, scale_right = %f\n", scale_left, scale_rght);


	// to-do: VHS: we need to add an extra dimension to log_site_lk_cat to account for mixed BL.
	// Here we're just taking the log of non-logged likelihoods
	For(j, tree->mod->n_catg){
		tree->log_site_lk_cat[j][site] =
			(m3ldbl)log(tree->log_site_lk_cat[j][site]) +
			(m3ldbl)scale_left +
			(m3ldbl)scale_rght;
	}

	tree->site_lk[site] = log_site_lk;
	tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;

	//PhyML_Printf(" . debug: line 724: tree->c_lnL_sorted[%d] = %d * %f = %f\n", site, tree->data->wght[site], log_site_lk, tree->c_lnL_sorted[site]);


	return log_site_lk;
}

/*********************************************************/

m3ldbl Return_Lk(arbre *tree)
{
	Lk(tree);
	return tree->c_lnL;
}

/*********************************************************/

m3ldbl Return_Abs_Lk(arbre *tree)
{
	Lk(tree);
	return fabs(tree->c_lnL);
}

/*********************************************************/

matrix *ML_Dist(allseq *data, model *mod)
{
	int i,j,k,l;
	m3ldbl init;
	int n_catg;
	m3ldbl d_max,sum;
	matrix *mat;
	allseq *twodata,*tmpdata;
	int state0, state1,len;
	m3ldbl *F;
	eigen *eigen_struct;

	tmpdata             = (allseq *)mCalloc(1,sizeof(allseq));
	tmpdata->c_seq      = (seq **)mCalloc(2,sizeof(seq *));
	tmpdata->b_frq      = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl));
	tmpdata->ambigu     = (short int *)mCalloc(data->crunch_len,sizeof(short int));
	F                   = (m3ldbl *)mCalloc(mod->ns*mod->ns,sizeof(m3ldbl ));
	eigen_struct        = (eigen *)Make_Eigen_Struct(mod);

	tmpdata->n_otu      = 2;

	tmpdata->crunch_len = data->crunch_len;
	tmpdata->init_len   = data->init_len;

	mat = (mod->datatype == NT) ? ((mod->whichmodel < 10)?(K80_dist(data,2000)):(JC69_Dist(data,mod))):(JC69_Dist(data,mod));

	/*   Print_Mat(mat); */
	/*   Exit("\n"); */

	For(i,mod->n_catg) /* Don't use the discrete gamma distribution */
	{
		mod->gamma_rr[i]      = 1.0;
		mod->gamma_r_proba[i] = 1.0;
	}

	n_catg = mod->n_catg;
	mod->n_catg = 1;

	For(j,data->n_otu-1)
	{
		tmpdata->c_seq[0]       = data->c_seq[j];
		tmpdata->c_seq[0]->name = data->c_seq[j]->name;
		tmpdata->wght           = data->wght;

		for(k=j+1;k<data->n_otu;k++)
		{
			tmpdata->c_seq[1]       = data->c_seq[k];
			tmpdata->c_seq[1]->name = data->c_seq[k]->name;

			twodata = Compact_CSeq(tmpdata,mod);

			For(l,mod->ns) twodata->b_frq[l] = data->b_frq[l];
			Check_Ambiguities(twodata,mod->datatype,1);
			Hide_Ambiguities(twodata);

			init = mat->dist[j][k];
			if((init == DIST_MAX) || (init < .0)) init = 0.1;

			d_max = init;

			For(i,mod->ns*mod->ns) F[i]=.0;
			len = 0;
			For(l,twodata->c_seq[0]->len)
			{
				state0 = Assign_State(twodata->c_seq[0]->state+l,mod->datatype,mod->stepsize);
				state1 = Assign_State(twodata->c_seq[1]->state+l,mod->datatype,mod->stepsize);
				if((state0 > -1) && (state1 > -1))
				{
					F[mod->ns*state0+state1] += twodata->wght[l];
					len += (int)twodata->wght[l];
				}
			}
			if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (m3ldbl)len;}

			sum = 0.;
			For(i,mod->ns*mod->ns) sum += F[i];

			if(sum < .001) d_max = -1.;
			else if((sum > 1. - .001) && (sum < 1. + .001)){
				//				m3ldbl * tmp = mCalloc(tree->mod->n_l,sizeof(m3ldbl));
				//				/**
				//				 * JSJ: this is now in need of a code review
				//				 */
				//				int tmp0;
				//				For(tmp0,tree->mod->n_l) tmp[tmp0] = d_max;
				Opt_Dist_F_No_Bl(&(d_max),F,mod);
				//				Free(tmp);
			}
			else
			{
				PhyML_Printf("\n. sum = %f\n",sum);
				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
				Exit("");
			}


			/* /\* 	  BRENT *\/ */
			/* 	  d_max = Optimize_Dist(mod,init,twodata); */

			/* 	  PhyML_Printf("\n. Warning : not using the ML pairwise distances..."); */
			/* 	  d_max = init; */

			if(d_max >= DIST_MAX)
			{
				/* 	      PhyML_Printf("\n. Large distance encountered between %s and %s sequences.", */
				/* 		     tmpdata->c_seq[1]->name, */
				/* 		     tmpdata->c_seq[0]->name); */
				d_max = DIST_MAX;
			}

			/* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
			 *  will not be called
			 */

			mat->dist[j][k] = d_max;
			mat->dist[k][j] = mat->dist[j][k];
			Free_Cseq(twodata);
		}
	}

	mod->n_catg = n_catg;

	Free(tmpdata->ambigu);
	Free(tmpdata->b_frq);
	Free(tmpdata->c_seq);
	free(tmpdata);
	Free_Eigen(eigen_struct);
	Free(F);

	return mat;
}

/*********************************************************/


//m3ldbl Lk_Given_Two_Seq(allseq *data, int numseq1, int numseq2, m3ldbl dist, model *mod, m3ldbl *loglk)
//{
//
//	seq *seq1,*seq2;
//	m3ldbl site_lk,log_site_lk;
//	int i,j,k,l;
//	/*   plkflt **p_lk_l,**p_lk_r; */
//	plkflt *p_lk_l,*p_lk_r;
//	m3ldbl len;
//	int dim1,dim2;
//	dim1 = mod->ns;
//	dim2 = mod->ns * mod->ns;
//
//
//	DiscreteGamma(mod->gamma_r_proba, mod->gamma_rr, mod->alpha,
//			mod->alpha,mod->n_catg,mod->gamma_median);
//
//	seq1 = data->c_seq[numseq1];
//	seq2 = data->c_seq[numseq2];
//
//
//	p_lk_l = (plkflt *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(plkflt));
//	p_lk_r = (plkflt *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(plkflt));
//
//	if(dist < BL_MIN) dist = BL_START;
//	else if(dist > BL_MAX) dist = BL_START;
//
//	For(i,mod->n_catg)
//	{
//		len = dist*mod->gamma_rr[i];
//		if(len < BL_MIN) len = BL_MIN;
//		else if(len > BL_MAX) len = BL_MAX;
//		PMat(len,mod,dim2*i,mod->Pij_rr);
//	}
//
//	if(mod->datatype == NT)
//	{
//		For(i,data->c_seq[0]->len)
//		{
//			Init_Tips_At_One_Site_Nucleotides_Float(seq1->state[i],i*mod->ns,p_lk_l);
//			Init_Tips_At_One_Site_Nucleotides_Float(seq2->state[i],i*mod->ns,p_lk_r);
//		}
//	}
//	else
//	{
//		For(i,data->c_seq[0]->len)
//		{
//			Init_Tips_At_One_Site_AA_Float(seq1->state[i],i*mod->ns,p_lk_l);
//			Init_Tips_At_One_Site_AA_Float(seq2->state[i],i*mod->ns,p_lk_r);
//		}
//	}
//
//
//	site_lk = .0;
//	*loglk = 0;
//
//	For(i,data->c_seq[0]->len)
//	{
//		if(data->wght[i])
//		{
//			site_lk = log_site_lk = .0;
//			if(!data->ambigu[i])
//			{
//				For(k,mod->ns) {if(p_lk_l[i*mod->ns+k] > .0001) break;}
//				For(l,mod->ns) {if(p_lk_r[i*mod->ns+l] > .0001) break;}
//				For(j,mod->n_catg)
//				{
//					site_lk +=
//							mod->gamma_r_proba[j] *
//							mod->pi[k] *
//							//p_lk_l[i*dim1+k] *
//							p_lk_l[i*dimc+]
//							mod->Pij_rr[j*dim2+k*dim1+l] *
//							p_lk_r[i*dim1+l];
//				}
//			}
//			else
//			{
//				For(j,mod->n_catg)
//				{
//					For(k,mod->ns) /*sort sum terms ? No global effect*/
//					{
//						For(l,mod->ns)
//						{
//							site_lk +=
//									mod->gamma_r_proba[j] *
//									mod->pi[k] *
//									p_lk_l[i*dim1+k] *
//									mod->Pij_rr[j*dim2+k*dim1+l] *
//									p_lk_r[i*dim1+l];
//						}
//					}
//				}
//			}
//
//			if(site_lk <= .0)
//			{
//				PhyML_Printf("'%c' '%c'\n",seq1->state[i],seq2->state[i]);
//				Exit("\n. Err: site lk <= 0\n");
//			}
//
//			log_site_lk += (m3ldbl)log(site_lk);
//
//			*loglk += data->wght[i] * log_site_lk;/* sort sum terms ? No global effect*/
//		}
//	}
//
//	/*   For(i,data->c_seq[0]->len) */
//	/*     { */
//	/*       Free(p_lk_l[i]); */
//	/*       Free(p_lk_r[i]); */
//	/*     } */
//
//	Free(p_lk_l); Free(p_lk_r);
//	return *loglk;
//}

/*********************************************************/

void Unconstraint_Lk(arbre *tree)
{
	int i;

	tree->unconstraint_lk = .0;

	For(i,tree->data->crunch_len)
	{
		tree->unconstraint_lk +=
				tree->data->wght[i]*(m3ldbl)log(tree->data->wght[i]);
	}
	tree->unconstraint_lk -=
			tree->data->init_len*(m3ldbl)log(tree->data->init_len);
}


/*********************************************************/
void Update_P_Lk(arbre *tree, edge *b, node *d)
{

	/*
				   |
				   |<- b
				   |
				   d
				  / \
		 dir1 -> /   \ <- dir2
				/     \
			  n_v1    n_v2

	 */

	//PhyML_Printf("entered Update_P_Lk(..., edge %d, node %d)\n", b->num, d->num);
	//PhyML_Printf("left = %d, right = %d\n", b->left->num, b->rght->num);


	node *n_v1, *n_v2;
	m3ldbl p1_lk1,p2_lk2;
	plkflt *p_lk,*p_lk_v1,*p_lk_v2;
	double *Pij1, *Pij2;
	plkflt max_p_lk;
	plkflt *sum_scale, *sum_scale_v1, *sum_scale_v2;
	plkflt scale_v1, scale_v2;
	int i,j,k;
	int catg, site;
	int dir1,dir2;
	int n_patterns;
	/**
	* ambiguity_check_v1 is used to store the site specific ambiguity of the neighbor
	* node in direction dir1. Ambiguity check v2 does the same for dir2
	*/
	int ambiguity_check_v1,ambiguity_check_v2;
	int state_v1,state_v2;

	int dima, dimb, dimc, dimaa, dimd;
	dima = tree->mod->ns;
	dimb = dima * tree->mod->n_l;
	dimc = dimb * tree->mod->n_catg;
	dimaa = dima * dima; // = ns^2
	dimd = dimb * dima; // = number of BL sets * ns^2

	state_v1 = state_v2 = -1;
	ambiguity_check_v1 = ambiguity_check_v2 = -1;
	scale_v1 = scale_v2 = 0.0;
	p1_lk1 = p2_lk2 = .0;

	if(d->tax)
	{
		PhyML_Printf("\n. node %d is a leaf...",d->num);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("\n");
	}

	n_patterns = tree->n_pattern;

	dir1=dir2=-1;
	For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i); //VHS: here we set dir1 and dir2 to point to the two
															// edges (of three possible edges) which aren't
															// the edge pointed to by edge *b.

	if((dir1 == -1) || (dir2 == -1))
	{
		PhyML_Printf("\n. d = %d",d->num);
		PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num);
		PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num);
		PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	}

	n_v1 = d->v[dir1];
	n_v2 = d->v[dir2];

	// VHS: In the remainder of this method, we work on computing p_lk as the product of p_lk_v1 and p_lk_v2.
	// However, first we must establish where p_lk_v1 and p_lk_v2 point within our tree.  .

	// p_lk becomes a shortcut pointer to b->p_lk_left OR b->p_lk_right.
	if(d == b->left)
	{
		p_lk = b->p_lk_left;
		sum_scale = b->sum_scale_f_left;
	}
	else
	{
		p_lk = b->p_lk_rght;
		sum_scale = b->sum_scale_f_rght;
	}

	// VHS: p_lk_v1 is the partial likelihood for the subtree attached to branch dir1
	if(d == d->b[dir1]->left)
	{
		p_lk_v1 = d->b[dir1]->p_lk_rght;
		sum_scale_v1 = d->b[dir1]->sum_scale_f_rght;
	}
	else
	{
		p_lk_v1 = d->b[dir1]->p_lk_left;
		sum_scale_v1 = d->b[dir1]->sum_scale_f_left;
	}

	// VHS: p_lk_v2 is the partial likelihood for the subtree attached to branch dir2
	if(d == d->b[dir2]->left)
	{
		p_lk_v2 = d->b[dir2]->p_lk_rght;
		sum_scale_v2 = d->b[dir2]->sum_scale_f_rght;
	}
	else
	{
		p_lk_v2 = d->b[dir2]->p_lk_left;
		sum_scale_v2 = d->b[dir2]->sum_scale_f_left;
	}

	Pij1 = d->b[dir1]->Pij_rr;
	Pij2 = d->b[dir2]->Pij_rr;

	/**
	* For each site given all possible site patterns
	* (deal with site patterns rather than sites due
	* 	to phyml's sequence compression...)
	*/

#ifdef USE_OPENMP
	int chunk = n_patterns/omp_get_num_procs();
	//	int chunk = n_patterns/2;
	//printf("Chunk size: %i\n",chunk);
#pragma omp parallel for\
		default(shared) private(k,catg,i,j,site,scale_v1,scale_v2,\
				max_p_lk,state_v1,state_v2,ambiguity_check_v1,\
				ambiguity_check_v2,p1_lk1,p2_lk2)\
				schedule(static,chunk)
#endif
	for(site = 0; site < n_patterns; site++)
	{

		//printf("JSJ: In Update_P_Lk iterating over state pattern %i\n",site);

		/**
		* JSJ: If sum_scale_v1 was assigned a non-null value in the above if/else cascade,
		* the scale_v1 value is assigned as the site specific sum_scale_v1 from above.
		*/
		scale_v1 = (sum_scale_v1)?(sum_scale_v1[site]):(0.0);
		scale_v2 = (sum_scale_v2)?(sum_scale_v2[site]):(0.0);

		// VHS: As of October 2009, I don't think sum_scale[i] is ever set to a value
		// other than 0.0.  I'm not sure why sum_scale is included in this code at all.
		sum_scale[site] = scale_v1 + scale_v2;

		max_p_lk = -MDBL_MAX;
		state_v1 = state_v2 = -1; //just ints
		ambiguity_check_v1 = ambiguity_check_v2 = -1;

		/**
		* JSJ: If the model is set to greedy, and the node is terminal,
		*  then we check if the site is ambiguous
		* 	otherwise if it is not ambiguous we default to true
		*/
		if(!tree->mod->s_opt->greedy)
		{
			if(n_v1->tax)
			{
				ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
				if(!ambiguity_check_v1) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dima,tree);
			}

			if(n_v2->tax)
			{
				ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
				if(!ambiguity_check_v2) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dima,tree);
			}
		}
		//JSJ: if using Markov modulated Markov Model then the ambiguity checks are true
		if(tree->mod->use_m4mod)
		{
			ambiguity_check_v1 = 1;
			ambiguity_check_v2 = 1;
		}

		For(catg,tree->mod->n_catg)
		{
			For(k, tree->mod->n_l)
			{
				For(i,tree->mod->ns)
				{
					p1_lk1 = .0;

					// if n_v1 is terminal, then use it's known state...
					if((n_v1->tax) && (!tree->mod->s_opt->greedy))
					{
						if(!ambiguity_check_v1)
						{
							p1_lk1 = Pij1[catg*dimd + k*dimaa + i*dima + state_v1];  //Pij_rr[k][catg*dim3+i*dim2+state_v1];
						}
						else
						{
							For(j,tree->mod->ns)
							{
								p1_lk1 += Pij1[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)n_v1->b[0]->p_lk_tip_r[site*dima+j];
							}
						}
					}
					else // otherwise, n_v1 is non-terminal, so we need to consider all possible states...
					{
						For(j,tree->mod->ns)
						{
							p1_lk1 += Pij1[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)p_lk_v1[site*dimc + catg*dimb + k*dima + j];   //[site*dim1+catg*dim2+j];
						}
					}

					p2_lk2 = .0;

					if((n_v2->tax) && (!tree->mod->s_opt->greedy))
					{
						if(!ambiguity_check_v2)
						{
							p2_lk2 = Pij2[catg*dimd + k*dimaa + i*dima + state_v2];
						}
						else
						{
							For(j,tree->mod->ns)
							{
								p2_lk2 += Pij2[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)n_v2->b[0]->p_lk_tip_r[site*dima+j];
							}
						}
					}
					else
					{
						For(j,tree->mod->ns)
						{
							p2_lk2 += Pij2[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)p_lk_v2[site*dimc + catg*dimb + k*dima + j];//[site*dim1+catg*dim2+j];
						}
					}

					p_lk[site*dimc + catg*dimb + k*dima + i] = (plkflt)(p1_lk1 * p2_lk2);

					if( p_lk[site*dimc + catg*dimb + k*dima + i] > max_p_lk )
					{
						max_p_lk = p_lk[site*dimc + catg*dimb + k*dima + i];
					}
				}// end for ns
			}// end for n_l
		}// end for catg


		if((max_p_lk < LIM_SCALE_VAL) || (max_p_lk > (1./LIM_SCALE_VAL)))
		{
			For(catg,tree->mod->n_catg)
			{
				for(k = 0; k < tree->mod->n_l; k++)
				{
					For(i,tree->mod->ns)
					{
						p_lk[site*dimc + catg*dimb + k*dima + i] /= max_p_lk;
					}
				}
			}
			// sum_scale[site] gets uses later in Lk_Core to calculate the final likelihood of the tree.
			sum_scale[site] += (plkflt)log(max_p_lk);
			//PhyML_Printf(" . debug: sum_scale[%d] = %f\n", site, sum_scale[site]);
		}

		//PhyML_Printf(" . debug lk.c 1180: edge->num=%d, sum_scale[%d]=%f\n", b->num, site, sum_scale[site]);

	} // end For(site,patterns)
}

/*********************************************************/


void Make_Tree_4_Lk(arbre *tree, allseq *alldata, int n_site)
{
	int i;

	tree->c_lnL_sorted = (m3ldbl *)mCalloc(tree->n_pattern, sizeof(m3ldbl));
	tree->site_lk      = (m3ldbl *)mCalloc(alldata->crunch_len,sizeof(m3ldbl));

	tree->log_site_lk_cat      = (m3ldbl **)mCalloc(tree->mod->n_catg,sizeof(m3ldbl *));
	For(i,tree->mod->n_catg)
	tree->log_site_lk_cat[i] = (m3ldbl *)mCalloc(alldata->crunch_len,sizeof(m3ldbl));

	tree->log_lks_aLRT = (m3ldbl **)mCalloc(3,sizeof(m3ldbl *));
	For(i,3) tree->log_lks_aLRT[i] = (m3ldbl *)mCalloc(tree->data->init_len,sizeof(m3ldbl));

	For(i,2*tree->n_otu-3)
	{
		Make_Edge_Lk(tree->t_edges[i],tree);
		Make_Edge_NNI(tree->t_edges[i]);
	}

	For(i,2*tree->n_otu-2) Make_Node_Lk(tree->noeud[i]);

	if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
	else Init_P_Lk_Tips_Int(tree);
}

/*********************************************************/

void Init_P_Lk_Tips_Double(arbre *tree)
{
	int curr_site,i,j,k,l;

	//int dim1,dim2;
	//dim1 = tree->mod->n_catg * tree->mod->ns;
	//dim2 = tree->mod->ns;

	int dima, dimb, dimc;
	dima = tree->mod->ns;
	dimb = dima * tree->mod->n_l;
	dimc = dimb * tree->mod->n_catg;

	Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
	{
		For(i,tree->n_otu)
		{
			if (tree->mod->datatype == NT)
				Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site],
						curr_site*dimc+0*dimb+0*dima,
						//curr_site*dim1+0*dim2,
						tree->noeud[i]->b[0]->p_lk_rght);
			else
				Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site],
						curr_site*dimc+0*dimb+0*dima,
						//curr_site*dim1+0*dim2,
						tree->noeud[i]->b[0]->p_lk_rght);

			for(j=0;j<tree->mod->n_catg;j++)
			{
				for(l=0; l<tree->mod->n_l; l++)
				{
					For(k,tree->mod->ns)
					{
						tree->noeud[i]->b[0]->p_lk_rght[curr_site*dimc+j*dimb+l*dima+k] =
							tree->noeud[i]->b[0]->p_lk_rght[curr_site*dimc+0*dimb+l*dima+k];
					}
				}
			}
		}
	}

#ifdef M4
	if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Double(tree);
#endif
}

/*********************************************************/

void Init_P_Lk_Tips_Int(arbre *tree)
{
	int curr_site,i,dim1;

	dim1 = tree->mod->ns;

	Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
	{
		For(i,tree->n_otu)
		{
			if(tree->mod->datatype == NT)
			{
				Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						curr_site*dim1,
						tree->noeud[i]->b[0]->p_lk_tip_r);
			}
			else
			{
				Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
						curr_site*dim1,
						tree->noeud[i]->b[0]->p_lk_tip_r);
			}
		}
	}

#ifdef M4
	if(tree->mod->m4mod) M4_Init_P_Lk_Tips_Int(tree);
#endif

}

/*********************************************************/

void Update_PMat_At_Given_Edge(edge *b_fcus, arbre *tree)
{
	int i,j;
	m3ldbl len;

	len = -1.0;
	//JSJ: fixes to b_fcus_l, Pir_rr and has_zero_br_len
	For(i,tree->mod->n_l){
		if(b_fcus->l[i] < BL_MIN)      b_fcus->l[i] = BL_MIN;
		else if(b_fcus->l[i] > BL_MAX) b_fcus->l[i] = BL_MAX;
	}

	For(i,tree->mod->n_catg) // foreach gamma category
	{
		For(j,tree->mod->n_l) // foreach branch length category
		{
			//PhyML_Printf(" . debug: Update_PMat_At_Given_Edge: edge %d, catg %d, BL %d, edge->l[j] = %f\n", b_fcus->num, i, j, b_fcus->l[j]);

			if(b_fcus->has_zero_br_len[j])
			{
				//PhyML_Printf( " . debug: this branch/catg/BL set has zero length\n");
				len = -1.0;
			}
			else
			{
				len = b_fcus->l[j]*tree->mod->gamma_rr[i];
				//PhyML_Printf( " . debug: lk.c 1420: ncatg=%d, bl=%d, len=%f\n", i, j, len);
				if(len < BL_MIN)      len = BL_MIN;
				else if(len > BL_MAX) len = BL_MAX;
			}

			PMat(len,
					tree->mod,
					tree->mod->ns*tree->mod->ns*tree->mod->n_l*i + tree->mod->ns*tree->mod->ns*j,
					b_fcus->Pij_rr);
		}
	}
}

/*********************************************************/

/* void Update_P_Lk_On_A_Path(node *a, node *d, edge *b, node *target_one_side, node *target_other_side, arbre *tree) */
/* { */


/*   /\* */
/*                 \               / */
/* 	   target\___________b_/ */
/* 		 /  \	\  \   \ */
/* 		/    \	 \  \	\ */

/*     target is the root of the subtree at which we want */
/*     the likelihood to be updated */
/*   *\/ */



/* /\*   PhyML_Printf("Update_p_lk on (%d %d) at %d (target=%d %d)\n", *\/ */
/* /\* 	 b->left->num, *\/ */
/* /\* 	 b->rght->num, *\/ */
/* /\* 	 a->num, *\/ */
/* /\* 	 target_one_side->num, *\/ */
/* /\* 	 target_other_side->num); *\/ */

/*   Update_P_Lk(tree,b,a); */
/*   if((a == target_one_side) && (d == target_other_side))  */
/*     return; */
/*   else */
/*     { */
/*       Update_P_Lk_On_A_Path(d, */
/* 			    d->v[tree->t_dir[d->num][target_other_side->num]], */
/* 			    d->b[tree->t_dir[d->num][target_other_side->num]], */
/* 			    target_one_side, */
/* 			    target_other_side, */
/* 			    tree);  */
/*     } */
/* } */

void Update_P_Lk_Along_A_Path(node **path, int path_length, arbre *tree)
{
	int i,j;

	For(i,path_length-1)
	{
		For(j,3)
		if(path[i]->v[j] == path[i+1])
		{
			if(path[i] == path[i]->b[j]->left)
			{
				Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);
			}

			else if(path[i] == path[i]->b[j]->rght)
			{
				Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
			}
			else
			{
				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
				Exit("");
			}
			break;
		}
#ifdef DEBUG
		if(j == 3)
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Exit("");
		}
#endif
	}
}

// VHS: I think F is an array of length ns*ns
m3ldbl Lk_Dist_No_Bl(m3ldbl *F, m3ldbl dist, model *mod)
{
	int i,j,k,l;
	m3ldbl lnL,len;
	int dim1,dim2,dim3;

	For(k,mod->n_catg)
	{
		For (l, mod->n_l)
		{
			len = dist*mod->gamma_rr[k];
			if(len < BL_MIN)      len = BL_MIN;
			else if(len > BL_MAX) len = BL_MAX;
			PMat(len,mod, mod->ns*mod->ns*mod->n_l*k + mod->ns*mod->ns*l, mod->Pij_rr);
		}
	}

	dim1 = mod->ns*mod->ns;
	dim2 = mod->ns;
	dim3 = mod->n_l * dim1;
	lnL = .0;
	For(k,mod->n_catg)
	{
		For (l, mod->n_l)
		{
			For(i,mod->ns)
			{
				For(j,mod->ns)
				{
					lnL += F[k*dim3 + l*dim1 + dim2*i + j] * log(mod->Pij_rr[k*dim3 + l*dim1 + dim2*i + j]);
				}
			}
		}
	}

	return lnL;
}
/*********************************************************/

m3ldbl Lk_Dist(m3ldbl *F, m3ldbl *dist, model *mod, arbre *tree)
{
	int i,j,k,l;
	m3ldbl lnL,len;
	int dim1,dim2,dim3;
	lnL = .0;
	dim1 = mod->ns*mod->ns;
	dim2 = mod->ns;
	dim3 = dim1 * mod->n_l;

	For(k,mod->n_catg)
	{
		For(l, mod->n_l)
		{
			len = dist[k]*mod->gamma_rr[k];
			if(len < BL_MIN)      len = BL_MIN;
			else if(len > BL_MAX) len = BL_MAX;
			PMat(len, mod, k*mod->ns*mod->ns*mod->n_l + l*mod->ns*mod->ns, mod->Pij_rr);
		}
	}

	For(k,mod->n_catg)
	{
		For(l, mod->n_l)
		{
			For(i,mod->ns)
			{
				For(j,mod->ns)
				{
					lnL += F[k*dim3 + l*dim1 + i*dim2 +j] * log(mod->Pij_rr[k*dim3 + l*dim1 + i*dim2 +j]);
				}
			}
		}

	}

	return lnL;
}

/*********************************************************/

m3ldbl Update_Lk_At_Given_Edge(edge *b_fcus, arbre *tree)
{
	//PhyML_Printf(" . debug: lk.c 1607: tree->c_lnL = %f\n", tree->c_lnL);

	// debug stuff:
	int dima = tree->mod->ns;
	int dimb = dima * tree->mod->n_l;
	int dimc = dimb * tree->mod->n_catg;
	int site = 3;
	if ( (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] < 0.0 || (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] > 1.0 )
	{
		PhyML_Printf(" . debug: lk.c 1648: the p_lk_left is invalid for edge %d at site %d\n", b_fcus->num, site);
	}


	if(!b_fcus->left->tax) Update_P_Lk(tree,b_fcus,b_fcus->left);
	if(!b_fcus->rght->tax) Update_P_Lk(tree,b_fcus,b_fcus->rght);

	// debug stuff:
	/*
	if ( (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] < 0.0 || (m3ldbl)b_fcus->p_lk_left[site*dimc + 0*dimb + 0*dima + 0] > 1.0 )
	{
		PhyML_Printf(" . debug: lk.c 1650: the p_lk_left is invalid for edge %d at site %d\n", b_fcus->num, site);
	}
	*/

	tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
	//PhyML_Printf(" . debug: lk.c 1615: tree->c_lnL = %f\n", tree->c_lnL);

	return tree->c_lnL;
}

/*********************************************************/

/*
 * JSJ: this function is currently not used, so don't bother
 * changing it more than needed to simply make it compile...
 *      root
           \
           /
          a
	  |
	  |
	  d
	 / \
        /   \
       w     x

       d->t has changed and we need to compute
       the likelihood.
       (1) update the three branch lengths l(ad), l(dw) and l(dx)
       (2) update the change proba matrices along these branches
       (3) update the likelihood of subtree (w,x) (WARNING: (a,x) and (a,w) are not updated)
 */
//m3ldbl Lk_Triplet(node *a, node *d, arbre *tree)
//{
//	int i;
//	m3ldbl max_height;
//	m3ldbl up_bound, low_bound;
//
//	if(d->tax)
//	{
//		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
//		Warn_And_Exit("");
//	}
//
//	up_bound = low_bound = -1.0;
//	max_height = -1.0;
//	For(i,3)
//	{
//		if((d->v[i] != a) && (d->b[i] != tree->e_root))
//		{
//			if(tree->rates->nd_t[d->v[i]->num] > max_height)
//			{
//				max_height = tree->rates->nd_t[d->v[i]->num];
//			}
//		}
//		else
//		{
//			up_bound =
//					(a == tree->n_root)?
//							(tree->rates->nd_t[a->num]):
//								(tree->rates->nd_t[d->v[i]->num]);
//		}
//	}
//
//	low_bound = max_height;
//
//	if(up_bound < low_bound - 1.E-10)
//	{
//		PhyML_Printf("\n. a->num=%d d->num=%d",a->num,d->num);
//		PhyML_Printf("\n. up_bound = %f, low_bound = %f",up_bound,low_bound);
//		Warn_And_Exit("\n");
//	}
//
//	if(tree->rates->nd_t[d->num] < low_bound) tree->rates->nd_t[d->num] = low_bound;
//	else if(tree->rates->nd_t[d->num] > up_bound) tree->rates->nd_t[d->num] = up_bound;
//
//	/* Step (1) */
//	For(i,3)
//	{ //JSJ: temporarily fixed b[i]->l
//		if((d->v[i] != a) && (d->b[i] != tree->e_root))
//		{
//			d->b[i]->l[0] =
//					(tree->rates->nd_t[d->num] - tree->rates->nd_t[d->v[i]->num]) *
//					tree->rates->clock_r *
//					tree->rates->br_r[d->b[i]->num];
//		}
//		else
//		{
//			if(a == tree->n_root)
//			{
//				d->b[i]->l[0] =
//						(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num] +
//								tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]) * tree->rates->clock_r;
//			}
//			else
//			{
//				d->b[i]->l[0] = (tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->clock_r * tree->rates->br_r[d->b[i]->num];
//			}
//		}
//	}
//
//	/* Step (2) */
//	For(i,3) Update_PMat_At_Given_Edge(d->b[i],tree);
//
//	For(i,3)
//	if((d->v[i] == a) || (d->b[i] == tree->e_root))
//	{
//		Update_P_Lk(tree,d->b[i],d);
//		Lk_At_Given_Edge(d->b[i],tree);
//		break;
//	}
//
//	return tree->c_lnL;
//}

/*********************************************************/

void Print_Lk_Given_Edge_Recurr(node *a, node *d, edge *b, arbre *tree)
{
	PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
			b->num,
			b->left->num,
			b->rght->num,
			Lk_At_Given_Edge(b,tree));

	if(d->tax) return;
	else
	{
		int i;
		For(i,3)
		if(d->v[i] != a)
			Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
	}
}

/*********************************************************/

/* Returns a vector containing the posterior probabilities of
   the different branch rate classes
 */
m3ldbl *Post_Prob_Rates_At_Given_Edge(edge *b, m3ldbl *post_prob, arbre *tree)
{
	m3ldbl norm_factor;
	int rcat, scale_int;
	m3ldbl sum,log2,lnL,worst_lnL,best_lnL,mid_lnL;


	For(rcat,tree->mod->n_rr_branch) post_prob[rcat] = .0;

	log2 = 0.6931472;

	best_lnL  = UNLIKELY;
	worst_lnL = .0;
	For(rcat,tree->mod->n_rr_branch)
	{
		tree->rates->br_r[b->num] = tree->mod->rr_branch[rcat];
		lnL = Lk_At_Given_Edge(b,tree);

		if(lnL < worst_lnL) worst_lnL = lnL;
		if(lnL > best_lnL)  best_lnL  = lnL;
		post_prob[rcat] = lnL;

		tree->rates->br_r[b->num] = 1.0;
	}

	mid_lnL = worst_lnL + fabs(worst_lnL - best_lnL)/2.;

	/* exp(log(P(D|r))) is hard to compute. Try exp(log(K*P(D|r)))
     instead. The value of K is chosen such that the most accurate
     estimates of the posterior probabilities are obtained for the
     most likely rate classes.
	 */
	scale_int = 0;
	do scale_int++;
	while(best_lnL + scale_int * log2 < 20.0);

	norm_factor = .0;
	For(rcat,tree->mod->n_rr_branch)
	{
		/*       PhyML_Printf("\n. best_lnL=%f curr_lnL=%f %f %E", */
		/*  	     best_lnL, */
		/* 	     post_prob[rcat] , */
		/* 	     post_prob[rcat] + scale_int * log2, */
		/* 	     exp(post_prob[rcat] + scale_int * log2)); */

		post_prob[rcat] = exp(post_prob[rcat] + scale_int * log2);


		post_prob[rcat] *= tree->mod->p_rr_branch[rcat];
		norm_factor += post_prob[rcat];
	}

	sum = .0;
	For(rcat,tree->mod->n_rr_branch)
	{
		post_prob[rcat] /= norm_factor;
		/*       PhyML_Printf("%f ",post_prob[rcat]); */
		sum += post_prob[rcat];
	}

	if(sum < 0.999 || sum > 1.001)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	return post_prob;
}

/*********************************************************/

m3ldbl Lk_With_MAP_Branch_Rates(arbre *tree)
{
	int br,rcat,best_post_prob_cat;
	m3ldbl *post_prob;
	m3ldbl best_post_prob;
	edge *b;


	post_prob = (m3ldbl *)mCalloc(tree->mod->n_rr_branch,sizeof(m3ldbl));

	Lk(tree);
	Record_Br_Len(NULL,tree);

	For(br,2*tree->n_otu-3)
	{
		b = tree->t_edges[br];

		/* Compute the posterior probability of each rate class on edge b */
		post_prob = (m3ldbl *)Post_Prob_Rates_At_Given_Edge(b,post_prob,tree);

		/* Find the most probable class */
		best_post_prob = UNLIKELY;
		best_post_prob_cat = -1;
		For(rcat,tree->mod->n_rr_branch)
		{
			if(post_prob[rcat] > best_post_prob)
			{
				best_post_prob = post_prob[rcat];
				best_post_prob_cat = rcat;
			}
		}

		/* The relative rate on this branch corresponds to the most probable rate class */
		tree->rates->br_r[br] = tree->mod->rr_branch[best_post_prob_cat];
	}

	Lk(tree);

	For(br,2*tree->n_otu-3) tree->rates->br_r[br] = 1.0;

	Free(post_prob);

	return tree->c_lnL;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
