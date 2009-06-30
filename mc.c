/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */

#ifdef MC
#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "mc.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"
#include "mcmc.h"
#include "numeric.h"

/*********************************************************/

int MC_main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  option *io;
  arbre *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
  time_t t_beg,t_end;
  phydbl best_lnL,most_likely_size,tree_size;
  int r_seed;
  char *most_likely_tree;

  
#ifdef MPI
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    PhyML_Printf("\n. Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD,&Global_numTask);
  MPI_Comm_rank(MPI_COMM_WORLD,&Global_myRank);
#endif

#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif



  tree             = NULL;
  mod              = NULL;
  data             = NULL;
  most_likely_tree = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed); rand();
  Make_Model_Complete(io->mod);
  mod = io->mod;
  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  io->compress_seq = 1;

  mat = NULL;
  tree_line_number = 0;


  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  

  For(num_data_set,io->n_data_sets)
    {
      n_otu = 0;
      best_lnL = UNLIKELY;
      data = Get_Seq(io,0);

      if(data)
	{
	  if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
	  PhyML_Printf("\n. Compressing sequences...\n");
	  alldata = Compact_Seq(data,io);
	  Free_Seq(data,alldata->n_otu);
	  Check_Ambiguities(alldata,io->mod->datatype,io->mod->stepsize);

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  Init_Model(alldata,mod);

		  /* A BioNJ tree is built here */
		  if(!io->in_tree) tree = Dist_And_BioNJ(alldata,mod,io);
		  /* A user-given tree is used here instead of BioNJ */
		  else             tree = Read_User_Tree(alldata,mod,io);
		  

		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));


		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/

		  /* IMPORTANCE SAMPLING STUFF */

/* 		  phydbl *cov; */
		  int i,j;

		  phydbl *mean, *cov, *min, *max, *x;
		  int dim;

		  dim = 2;
		  
		  mean = (phydbl *)mCalloc(dim,    sizeof(phydbl));
		  cov  = (phydbl *)mCalloc(dim*dim,sizeof(phydbl));
		  min  = (phydbl *)mCalloc(dim,    sizeof(phydbl));
		  max  = (phydbl *)mCalloc(dim,    sizeof(phydbl));
		  
		  mean[0] = 0.; mean[1] = 10.;
		  
		  cov[0*dim+0] = 9. ; cov[0*dim+1] = 0.;
		  cov[1*dim+0] =  0.; cov[1*dim+1] = 9.;
		  
		  min[0] = -1.; max[0] = 1.;
		  min[1] = -10.; max[1] = 10.;

		  For(i,1000)
		    {
		      x = Rnorm_Multid_Trunc(mean,cov,min,max,dim);
		      For(j,dim) printf("%10lf ",x[j]);
		      printf("\n");
		    }
		  Exit("\n");

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = alldata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

		  Fill_Dir_Table(tree);
		  Update_Dirs(tree);
		  Make_Tree_4_Pars(tree,alldata,alldata->init_len);
		  Make_Tree_4_Lk(tree,alldata,alldata->init_len);
		  
		  Round_Optimize(tree,tree->data,ROUND_MAX);
		  cov = Hessian(tree);
		  
/* 		  RATES_Initialize_Mean_Rates(tree); */

		  Exit("\n");

		  /* END OF IMPORTANCE SAMPLING STUFF */

		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/




		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/

		  /* COMPOUND POISSON STUFF */

/* 		  int n_otu,i; */

/* 		  n_otu = 30; */

/* 		  tree = Generate_Random_Tree_From_Scratch(n_otu,1); */

/* 		  tree->mod         = mod; */
/* 		  tree->io          = io; */
/* 		  tree->data        = alldata; */
/* 		  tree->both_sides  = 1; */
/* 		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize; */

/* 		  For(i,tree->n_otu) strcpy(tree->noeud[i]->name,alldata->c_seq[i]->name); */

/* 		  Fill_Dir_Table(tree); */
/* 		  Update_Dirs(tree); */
/* 		  Make_Tree_4_Pars(tree,alldata,alldata->init_len); */
/* 		  Make_Tree_4_Lk(tree,alldata,alldata->init_len); */

/* 		  Evolve(tree->data,tree->mod,tree); */
/* 		  Init_Ui_Tips(tree); */
/* 		  Init_P_Pars_Tips(tree); */
/* 		  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree); */
/* 		  else Init_P_Lk_Tips_Int(tree); */

/* 		  For(i,2*tree->n_otu-1) tree->rates->true_t[i] = tree->rates->nd_t[i]; */

/* 		  printf("%s\n",Write_Tree(tree)); */


/* 		  /\***********************************\/ */
/*  		  tree->rates->approx = 1; */
/* 		  printf("\n. USING APPROX %d",tree->rates->approx); */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct(tree->mcmc); */
/* 		  tree->both_sides = 1; */
/* 		  tree->rates->model = COMPOUND_COR; */
/* 		  Lk(tree); */
/* 		  RATES_Lk_Rates(tree); */
/* 		  printf("\n. CORRELATED lnL_data = %f lnL_rate = %f\n",tree->c_lnL,tree->rates->c_lnL); */
/* 		  tree->rates->bl_from_rt = 1; */
/* 		  MCMC(tree); */
/* 		  /\***********************************\/ */

/* 		  /\***********************************\/ */
/*  		  tree->rates->approx = 2; */
/* 		  printf("\n. USING APPROX %d",tree->rates->approx); */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct(tree->mcmc); */
/* 		  tree->both_sides = 1; */
/* 		  tree->rates->model = COMPOUND_COR; */
/* 		  Lk(tree); */
/* 		  RATES_Lk_Rates(tree); */
/* 		  printf("\n. CORRELATED lnL_data = %f lnL_rate = %f\n",tree->c_lnL,tree->rates->c_lnL); */
/* 		  tree->rates->bl_from_rt = 1; */
/* 		  MCMC(tree); */
/* 		  /\***********************************\/ */

/* 		  /\***********************************\/ */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct(tree->mcmc); */
/* 		  tree->both_sides = 1; */
/* 		  tree->rates->model = COMPOUND_NOCOR; */
/* 		  Lk(tree); */
/* 		  RATES_Lk_Rates(tree); */
/* 		  printf("\n. NOT CORRELATED lnL_data = %f lnL_rate = %f\n",tree->c_lnL,tree->rates->c_lnL); */
/* 		  tree->rates->bl_from_rt = 1; */
/* 		  MCMC(tree); */
/* 		  /\***********************************\/ */

/* 		  /\***********************************\/ */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct(tree->mcmc); */
/* 		  tree->both_sides = 1; */
/* 		  tree->rates->model = EXPONENTIAL; */
/* 		  Lk(tree); */
/* 		  RATES_Lk_Rates(tree); */
/* 		  printf("\n. EXPONENTIAL lnL_data = %f lnL_rate = %f\n",tree->c_lnL,tree->rates->c_lnL); */
/* 		  tree->rates->bl_from_rt = 1; */
/* 		  MCMC(tree); */
/* 		  /\***********************************\/ */

/* 		  /\***********************************\/ */
/* 		  tree->mcmc = (tmcmc *)MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Init_MCMC_Struct(tree->mcmc); */
/* 		  tree->both_sides = 1; */
/* 		  tree->rates->model = GAMMA; */
/* 		  Round_Optimize(tree,tree->data,100); */
/* 		  Lk(tree); */
/* 		  RATES_Lk_Rates(tree); */
/* 		  printf("\n. GAMMA lnL_data = %f lnL_rate = %f\n",tree->c_lnL,tree->rates->c_lnL); */
/* 		  tree->rates->bl_from_rt = 1; */
/* 		  MCMC(tree); */
/* 		  /\***********************************\/ */


/* 		  Exit("\n"); */

		  /* END OF COMPOUND POISSON STUFF */

		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/
		  /************************************/


		  /* Print the tree estimated using the current random (or BioNJ) starting tree */
		  if(io->mod->s_opt->n_rand_starts > 1)
		    {
		      Br_Len_Involving_Invar(tree);
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		    }

		  /* Record most likely tree in a string of characters */
		  if(tree->c_lnL > best_lnL)
		    {
		      best_lnL = tree->c_lnL;
		      Br_Len_Involving_Invar(tree);
		      most_likely_tree = Write_Tree(tree);
		      most_likely_size = Get_Tree_Size(tree);
		    }

/* 		  JF(tree); */

		  time(&t_end);
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?
			       (num_rand_tree):(num_tree));
		  
		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree))
		    {
		      /* Do one more iteration in the loop, but don't randomize the tree */
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }
		  
		  Free_Spr_List(tree);
		  Free_One_Spr(tree->best_spr);
		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}

	      /* Launch bootstrap analysis */
	      if(tree->mod->bootstrap) 
		{
		  PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");

                  #ifdef MPI
		  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		  PhyML_Printf("\n. The bootstrap analysis will use %d CPUs.",Global_numTask);
		  #endif

		  Bootstrap_From_String(most_likely_tree,alldata,mod,io);
		}
	      else if(tree->io->ratio_test) 
		{
		  /* Launch aLRT */
		  PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
		  aLRT_From_String(most_likely_tree,alldata,mod,io);
		}

	      /* Print the most likely tree in the output file */
	      PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n",io->out_tree_file);
	      if(io->n_data_sets == 1) rewind(io->fp_out_tree);

	      PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);

	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(alldata);
	}
    }


  if(io->mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n\n. Best log likelihood : %f\n",best_lnL);

  Free_Model(mod);

  if(io->fp_in_seq)     fclose(io->fp_in_seq);
  if(io->fp_in_tree)    fclose(io->fp_in_tree);
  if(io->fp_out_lk)     fclose(io->fp_out_lk);
  if(io->fp_out_tree)   fclose(io->fp_out_tree);
  if(io->fp_out_trees)  fclose(io->fp_out_trees);
  if(io->fp_out_stats)  fclose(io->fp_out_stats);

  Free(most_likely_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

/*********************************************************/

void MC_Least_Square_Node_Times(edge *e_root, arbre *tree)
{

  /* Solve A.x = b, where x are the node time estimated
     under the least square criterion.

     A is a n x n matrix, with n being the number of
     nodes in a rooted tree (i.e. 2*n_otu-1).

   */

  phydbl *A, *b, *x;
  int n;
  int i,j;
  node *root;
  
  
  PhyML_Printf("\n. Making the tree molecular clock like.");
  
  n = 2*tree->n_otu-1;
  
  A = (phydbl *)mCalloc(n*n,sizeof(phydbl));
  b = (phydbl *)mCalloc(n,  sizeof(phydbl));
  x = (phydbl *)mCalloc(n,  sizeof(phydbl));
  
  (e_root)?(Add_Root(e_root,tree)):(Add_Root(tree->t_edges[0],tree));
  
  root = tree->n_root;

  MC_Least_Square_Node_Times_Pre(root,root->v[0],A,b,n,tree);
  MC_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);
  
  b[root->num] = tree->e_root->l/2.;
  
  A[root->num * n + root->num]       = 1.0;
  A[root->num * n + root->v[0]->num] = -.5;
  A[root->num * n + root->v[1]->num] = -.5;
    
  Matinv(A, n, n, NULL);

  For(i,n) x[i] = .0;
  For(i,n) For(j,n) x[i] += A[i*n+j] * b[j];

  For(i,n-1) { tree->rates->nd_t[tree->noeud[i]->num] = x[i]; }
  tree->rates->nd_t[root->num] = x[n-1];
  tree->n_root->l[0] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num];
  tree->n_root->l[1] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num];


  /* Rescale the node times such that the time at the root
     is -100. This constraint implies that the clock rate
     is fixed to the actual tree length divided by the tree
     length measured in term of differences of node times */

  phydbl scale_f,time_tree_length,tree_length;

  scale_f = -100./tree->rates->nd_t[root->num];
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= scale_f;

  time_tree_length = 0.0;
  For(i,2*tree->n_otu-3)
    if(tree->t_edges[i] != tree->e_root)
      time_tree_length +=
	fabs(tree->rates->nd_t[tree->t_edges[i]->left->num] -
	     tree->rates->nd_t[tree->t_edges[i]->rght->num]);
  time_tree_length += fabs(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num]);
  time_tree_length += fabs(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num]);
  
  tree_length = 0.0;
  For(i,2*tree->n_otu-3) tree_length += tree->t_edges[i]->l;

  tree->rates->clock_r = tree_length / time_tree_length;

  Free(A);
  Free(b);
  Free(x);

}

/*********************************************************/

void MC_Least_Square_Node_Times_Pre(node *a, node *d, phydbl *A, phydbl *b, int n, arbre *tree)
{
  if(d->tax)
    {
      A[d->num * n + d->num] = 1.;
      
      /* Set the time stamp at tip nodes to 0.0 */
/*       PhyML_Printf("\n. Tip node date set to 0"); */
      b[d->num] = 0.0;
      return;
    }
  else
    {
      int i;
      
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  MC_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      For(i,3)
	{
	  A[d->num * n + d->v[i]->num] = -1./3.;
	  if(d->v[i] != a) b[d->num] += d->b[i]->l;
	  else             b[d->num] -= d->b[i]->l;
	}
      b[d->num] /= 3.;
    }
}

/*********************************************************/

/* Adjust node times in order to have correct time stamp ranking with
 respect to the tree topology */

void MC_Adjust_Node_Times(arbre *tree)
{
  MC_Adjust_Node_Times_Pre(tree->n_root->v[0],tree->n_root->v[1],tree);
  MC_Adjust_Node_Times_Pre(tree->n_root->v[1],tree->n_root->v[0],tree);

  if(tree->rates->nd_t[tree->n_root->num] > MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
					     tree->rates->nd_t[tree->n_root->v[1]->num]))
    {
      tree->rates->nd_t[tree->n_root->num] = MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
					      tree->rates->nd_t[tree->n_root->v[1]->num]);
    }
}

/*********************************************************/

void MC_Adjust_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl min_height;

      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    MC_Adjust_Node_Times_Pre(d,d->v[i],tree);
	  }

      min_height = 0.0;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(tree->rates->nd_t[d->v[i]->num] < min_height)
		{
		  min_height = tree->rates->nd_t[d->v[i]->num];
		}
	    }
	}

      if(tree->rates->nd_t[d->num] > min_height) tree->rates->nd_t[d->num] = min_height;

      if(tree->rates->nd_t[d->num] < -100.) tree->rates->nd_t[d->num] = -100.;

    }
}

/*********************************************************/

  /* Multiply each time stamp at each internal 
     node by  'tree->time_stamp_mult'.
   */

void MC_Mult_Time_Stamps(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->noeud[i]->num] *= fabs(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] *= fabs(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

/* Divide each time stamp at each internal 
   node by  'tree->time_stamp_mult'.
*/
void MC_Div_Time_Stamps(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->noeud[i]->num] /= fabs(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] /= fabs(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

void MC_Bl_From_T(arbre *tree)
{
  phydbl mean_rate, branch_rate;

  /* Branch lengths are deduced from time stamps */
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Bl_From_T_Post(tree->n_root,tree->n_root->v[1],NULL,tree);

  mean_rate   = tree->rates->clock_r;
  branch_rate = tree->rates->br_r[tree->e_root->num];

  if(tree->rates->use_rates)
    tree->e_root->l = 
      mean_rate * branch_rate * (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num]) + 
      mean_rate * branch_rate * (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->rght->num]);
  else
    tree->e_root->l = 
      (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num]) + 
      (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->rght->num]);

  /* Actual formula =>  tree->e_root->l = 
     (tree->n_root->t - tree->e_root->left->t) + 
     (tree->n_root->t - tree->e_root->rght->t); */
  
  tree->n_root_pos = (tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->e_root->left->num])/tree->e_root->l;

}

/*********************************************************/

void MC_Bl_From_T_Post(node *a, node *d, edge *b, arbre *tree)
{

  if(b)
    {
      phydbl mean_rate, branch_rate;

      mean_rate   = tree->rates->clock_r;
      branch_rate = tree->rates->br_r[b->num];

      if(tree->rates->use_rates)
	{
	  b->l = (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]) * mean_rate * branch_rate;
	}
      else
	{
	  b->l = (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
	}

      if(b->l < 0.0)
	{
	  PhyML_Printf("\n. Correction failed.");
	  PhyML_Printf("\n. d->t = %f a->t = %f",tree->rates->nd_t[d->num],tree->rates->nd_t[a->num]);
	  PhyML_Printf("\n. a->num=%d d->num=%d",a->num,d->num);
	  Warn_And_Exit("\n");
	}
    }

  if(d->tax) return;
  else
    {
      int i;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) MC_Bl_From_T_Post(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void MC_Round_Optimize(arbre *tree)
{
  int n_round,each;
  phydbl lk_old, lk_new, tol;

  lk_new = UNLIKELY;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 5;
  tol = 1.e-2;

  tree->both_sides = 1;
  Lk(tree);

  while(n_round < ROUND_MAX)
    {
      if(tree->mod->s_opt->opt_bl)
	{
	  MC_Optimize_Node_Times_Serie(tree->n_root,tree->n_root->v[0],tree);
	  MC_Optimize_Node_Times_Serie(tree->n_root,tree->n_root->v[1],tree);
	  MC_Optimize_Root_Height(tree);
	}

      tree->both_sides = 1;
      Lk(tree);

      if(tree->mod->s_opt->print) Print_Lk(tree,"[Node times         ]");

      lk_new = tree->c_lnL;
      if((!each) ||
	 (fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global))
	{
	  each = 1;
	  MC_Optimize_Tree_Height(tree);
	  Optimiz_All_Free_Param(tree,tree->mod->s_opt->print);
	  tree->both_sides = 1;
	  Lk(tree);
	}

      lk_new = tree->c_lnL;

      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global*10.) 
	{
	  PhyML_Printf("\n. lk_new = %f lk_old = %f",lk_new,lk_old);
	  Warn_And_Exit("\n. Optimisation failed ! (Round_Optimize)\n");
	}
      if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global)  break;
/*       if(fabs(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_local)  break; */
      else lk_old  = lk_new;
      n_round++;
      each--;
    }
}



/*********************************************************/

void MC_Optimize_Node_Times_Serie(node *a, node *d, arbre *tree)
{
  int i;      

  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf;
      phydbl lk_init;
      
      lk_init = tree->c_lnL;
      
      v1 = v2 = NULL;
      For(i,3) if(d->v[i] != a) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MAX(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_sup = tree->rates->nd_t[a->num];

      if(t_sup < t_inf - MDBL_MAX)
	{
	  PhyML_Printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      else
	{
	  Node_Time_Brent(t_inf,tree->rates->nd_t[d->num],t_sup,
			  tree->mod->s_opt->min_diff_lk_local,
			  a,d,tree,
			  tree->mod->s_opt->brent_it_max);  
	}

      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local*10.)
/*       if(tree->c_lnL < lk_init - 1.E-03) */
	{
	  PhyML_Printf("\n. t-inf= %f t-sup=%f t-est=%f",t_inf,t_sup,tree->rates->nd_t[d->num]);
	  PhyML_Printf("\n. %f -- %f",lk_init,tree->c_lnL);
	  PhyML_Printf("\n. a->num = %d, d->num = %d",a->num,d->num);
	  Warn_And_Exit("\n. Err. in MC_Optimize_Node_Times_Serie.");
	}

/*       PhyML_Printf("\n. init_lnL = %f c_lnL = %f",lk_init,tree->c_lnL); */
      
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    MC_Optimize_Node_Times_Serie(d,d->v[i],tree);
	  }

      For(i,3) 
	if((d->v[i] == a) || (d->b[i] == tree->e_root)) 
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    break;
	  }
    }
}

/*********************************************************/

void MC_Print_Node_Times(node *a, node *d, arbre *tree)
{
  edge *b;
  int i;
  
  b = NULL;
  For(i,3) if((d->v[i]) && (d->v[i] == a)) {b = d->b[i]; break;}

  PhyML_Printf("\n. (%3d %3d) a->t = %f d->t = %f (#=%f) b->l = %f",
	 a->num,d->num,
	 tree->rates->nd_t[a->num],
	 tree->rates->nd_t[d->num],
	 tree->rates->nd_t[a->num]-tree->rates->nd_t[d->num],
	 (b)?(b->l):(-1.0));
  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  MC_Print_Node_Times(d,d->v[i],tree);
    }
}

/*********************************************************/

edge *MC_Find_Best_Root_Position(arbre *tree)
{
  int i;
  edge *best_edge;
  phydbl best_lnL,best_pos;

  Record_Br_Len(NULL,tree);
  best_pos = -1.;
  best_edge = NULL;
  best_lnL = UNLIKELY;
  For(i,2*tree->n_otu-3)
    {
      Restore_Br_Len(NULL,tree);
      PhyML_Printf("\n. Root positioned on edge %3d",i);
      MC_Least_Square_Node_Times(tree->t_edges[i],tree);      
      MC_Adjust_Node_Times(tree);
      MC_Round_Optimize(tree);
      if(tree->c_lnL > best_lnL)
	{
	  best_lnL = tree->c_lnL;
	  best_edge = tree->t_edges[i];
	  best_pos = tree->n_root_pos;
	}
    }
  tree->n_root_pos = best_pos; /* Set the root node to its best position */ 
  PhyML_Printf("\n. Best root position: edge %3d",best_edge->num);
  PhyML_Printf("\n. Best constrained lnL = %f",best_lnL);
  PhyML_Fprintf(tree->io->fp_out_stats,"\n. Best constrained lnL = %f",best_lnL);
  return best_edge;
}

/*********************************************************/

edge *MC_Find_Best_Root_Position_Approx(arbre *tree)
{
  int i;
  edge *best_edge;
  phydbl best_lnL,best_pos;

  Record_Br_Len(NULL,tree);
  best_pos  = -1.;
  best_edge = NULL;
  best_lnL  = UNLIKELY;
  For(i,2*tree->n_otu-3)
    {
      Restore_Br_Len(NULL,tree);
      PhyML_Printf("\n. Root positioned on edge %3d",i);
      MC_Least_Square_Node_Times(tree->t_edges[i],tree);      
      MC_Adjust_Node_Times(tree);
      Lk(tree);
      PhyML_Printf("\n. LnL = %f",tree->c_lnL);
      if(tree->c_lnL > best_lnL)
	{
	  best_lnL  = tree->c_lnL;
	  best_edge = tree->t_edges[i];
	  best_pos  = tree->n_root_pos;
	}
    }
  tree->n_root_pos = best_pos; /* Set the root node to its best position */ 
  PhyML_Printf("\n. Best root position: edge %3d",best_edge->num);
  PhyML_Printf("\n. Best constrained lnL = %f",best_lnL);
  PhyML_Fprintf(tree->io->fp_out_stats,"\n. Best constrained lnL = %f",best_lnL);
  return best_edge;
}

/*********************************************************/

void MC_Optimize_Tree_Height(arbre *tree)
{
  tree->mod->s_opt->tree_size_mult = 1.0;
  Time_Stamps_Mult_Brent(0.1,1.0,10.0,
			 tree->mod->s_opt->min_diff_lk_global,
			 tree,100);
}

/*********************************************************/

void MC_Optimize_Root_Height(arbre *tree)
{
  /*
           t_root
           --x---
l_2 / mu  |      | l_2 / mu
          |      |
          x      x t_r
l_1 / mu  |
          |
          x t_l

 l* = 2 l_2 + l_1 is the ML estimate of the root branch length.
 l_1 = (Max(t_l,t_r) - Min(t_l,t_r)) * mu 
 l_2 = (l* - l1)/2
 t_root = Max(t_l,t_r) + l_2 / mu 

  */
  
  phydbl l_1, l_2;
  phydbl mean_rate, branch_rate;

  mean_rate   = tree->rates->clock_r;
  branch_rate = tree->rates->br_r[tree->e_root->num];

  Br_Len_Brent_Default(tree->e_root,tree);

  l_1 = 
    (MAX(tree->rates->nd_t[tree->e_root->left->num],tree->rates->nd_t[tree->e_root->rght->num]) -
     MIN(tree->rates->nd_t[tree->e_root->left->num],tree->rates->nd_t[tree->e_root->rght->num])) *
    mean_rate * branch_rate;
  
  l_2 = (tree->e_root->l - l_1) / 2.;

  if(l_2 < 0.0) 
    {
      l_2 = 0.0;
      tree->e_root->l = l_1;
      Lk_At_Given_Edge(tree->e_root,tree);
    }

  tree->rates->nd_t[tree->n_root->num] = 
    MAX(tree->rates->nd_t[tree->e_root->left->num],tree->rates->nd_t[tree->e_root->rght->num]) +
    l_2 / (mean_rate * branch_rate);

/*  /\* Check that the optimal 'root branch' length is longer than the  */
/*      lower bound determined by time stamps on the left and */
/*      right handside of the root branch  */
/*   *\/ */
/*   if(tree->e_root->l <  */
/*      MAX(tree->e_root->left->t,tree->e_root->rght->t) -  */
/*      MIN(tree->e_root->left->t,tree->e_root->rght->t)) */
/*     { */
/*       tree->e_root->l =  */
/* 	MAX(tree->e_root->left->t,tree->e_root->rght->t) -  */
/* 	MIN(tree->e_root->left->t,tree->e_root->rght->t); */
/*       Lk_At_Given_Edge(tree->e_root,tree); */
/*     } */

/*   tree->n_root->t = .5 * (tree->e_root->l + tree->e_root->left->t + tree->e_root->rght->t); */
 
/* /\*   PhyML_Printf("\n. Root->t = %f left->t=%f rght->t=%f e_root->l=%f", *\/ */
/* /\* 		 tree->n_root->t, *\/ */
/* /\* 		 tree->e_root->left->t, *\/ */
/* /\* 		 tree->e_root->rght->t, *\/ */
/* /\* 		 tree->e_root->l); *\/ */

  if((tree->rates->nd_t[tree->n_root->num] < tree->rates->nd_t[tree->e_root->left->num]-1.E-4) ||
     (tree->rates->nd_t[tree->n_root->num] < tree->rates->nd_t[tree->e_root->rght->num]-1.E-4))
    {
      PhyML_Printf("\n. t_root = %f t_left = %f t_rght = %f",
	     tree->rates->nd_t[tree->n_root->num],
	     tree->rates->nd_t[tree->e_root->left->num],
	     tree->rates->nd_t[tree->e_root->rght->num]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
}

/*********************************************************/

void MC_Estimate_Branch_Rate_Parameter(arbre *tree)
{

  /* The tree should be clock-like already */
  Branch_Rate_Shape_Brent(0.3, tree->mod->rr_branch_alpha, 100.,
			  tree->mod->s_opt->min_diff_lk_global, 
			  &(tree->mod->rr_branch_alpha), 
			  tree,tree->mod->s_opt->brent_it_max);  
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments(arbre *tree)
{
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[0],NULL,tree);
  MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  if(d->tax) return;

  if(b)
    {
      phydbl t0, t1, t2, t3;
      phydbl mu1, mu2, mu3;
      phydbl K;

      t0 = tree->rates->nd_t[a->num];
      t1 = tree->rates->nd_t[d->num];

      mu1 = tree->rates->br_r[b->num];
      mu2 = -1.;
      mu3 = -1.;

      t2 = t3 = -1.;
      For(i,3)
	if(d->v[i] != a)
	  {
	    if(t2 < 0) 
	      {
		t2  = tree->rates->nd_t[d->v[i]->num];
		mu2 = tree->rates->br_r[d->b[i]->num];
	      }
	    if(t3 < 0) 
	      {
		t3  = tree->rates->nd_t[d->v[i]->num];
		mu3 = tree->rates->br_r[d->b[i]->num];
	      }
	  }

      PhyML_Printf("\n. t0=%f, t1=%f, t2=%f, t3=%f, mu1=%f, mu2=%f, mu3=%f",
	     t0,t1,t2,t3,
	     mu1,mu2,mu3);

      if(t2 < 1.E-10 && t3 < 1.E-10)
	{
	  K = t0 / t1;
	}
      else
	{
	  K = 
	    (pow(mu2,2)*t2)/(pow(t1-t2,2)) +
	    (pow(mu1,2)*t0)/(pow(t0-t1,2)) +
	    (pow(mu3,2)*t3)/(pow(t1-t3,2)) +
	    (pow(mu1,2)*t0)/(pow(t0-t1,2)) +
	    mu2 * mu1 * (t0 + t2) / ((t1-t2)*(t0-t1)) +
	    mu3 * mu1 * (t0 + t3) / ((t1-t3)*(t0-t1)) ;

	  K /= 
	    t1 *
	    (
	    (pow(mu2,2))/(pow(t1-t2,2))  +
	    (pow(mu1,2))/(pow(t0-t1,2))  +
	    (pow(mu3,2))/(pow(t1-t3,2))  +
	    (pow(mu1,2))/(pow(t0-t1,2))  +
	    2.*mu2*mu1/((t1-t2)*(t0-t1)) +
	    2.*mu3*mu1/((t1-t3)*(t0-t1)) 
	    );
	}
      PhyML_Printf("\n. K = %f",K); 
    }

  For(i,3) 
    if(d->v[i] != a)
      MC_Compute_Rates_And_Times_Least_Square_Adjustments_Post(d,d->v[i],d->b[i],tree);

}

/*********************************************************/

int MC_Check_MC(arbre *tree)
{
  int i;
  phydbl dist,eps;

  eps = 1.E-06;

  Dist_To_Root(tree->n_root,tree);
  dist = tree->noeud[0]->dist_to_root;
  
  For(i,tree->n_otu)
    {
      if(fabs(tree->noeud[i]->dist_to_root - dist) > eps)
	{
	  return 0;
	}
    }
  return 1;
}

/*********************************************************/


#endif
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
