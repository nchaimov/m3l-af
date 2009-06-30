/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for Markov-Modulated Markov Models (M4) */
#ifdef M4

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
#include "m4.h"
#include "mc.h"
#include "draw.h"
#ifdef MG
#include "mg.h"
#endif


int M4_main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  option *io;
  arbre *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
  m4 *m4mod;
  time_t t_beg,t_end;
  div_t hour,min;
  phydbl best_lnL;
  int bootstrap_this_tree;
  int r_seed;


#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif

  tree = NULL;
  mod  = NULL;
  data = NULL;
  bootstrap_this_tree = 1;
  best_lnL = UNLIKELY;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  Make_Model_Complete(io->mod);
  mod = io->mod;
  m4mod = mod->m4mod;
  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;


  if(io->mod->s_opt->random_input_tree) bootstrap_this_tree = 0;

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
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search == SPR_MOVE))
		    PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  Init_Model(alldata,mod);
		  if(io->m4_model) M4_Init_Model(m4mod,alldata,mod);		    

		  if(!io->in_tree)
		    {
		      PhyML_Printf("\n. Computing pairwise distances...\n");
		      mat = ML_Dist(alldata,mod);
		      Fill_Missing_Dist(mat);
		      PhyML_Printf("\n. Building BIONJ tree...\n");
		      mat->tree = Make_Tree_From_Scratch(alldata->n_otu,alldata);
		      Bionj(mat);
		      tree      = mat->tree;
		      tree->mat = mat;

		    }
		  else if(io->in_tree == 2)
		    {
		      if((io->n_trees == 1) || (!num_tree))
			{
			  rewind(io->fp_in_tree);
			  tree_line_number = 0;
			}

		      if(io->n_trees > 1) PhyML_Printf("\n. Reading tree [#%d]\n",tree_line_number+1);
		      else PhyML_Printf("\n. Reading tree...\n");
		      fflush(NULL);

		      tree = Read_Tree_File(io->fp_in_tree);
		      tree_line_number++;

		      if(!tree)
			{
			  PhyML_Printf("\n. Input tree not found...\n");
			  Exit("\n\n");
			}

		      if(!tree->has_branch_lengths)
			{
			  PhyML_Printf("\n. Computing branch length estimates...\n");
			  Order_Tree_CSeq(tree,alldata);
			  mat = ML_Dist(alldata,mod);
			  mat->tree = tree;
			  mat->method = 0;
			  Bionj_Br_Length(mat);
			  tree->mat = mat;
			}

		      tree->mod        = mod;
		      tree->io         = io;
		      tree->data       = alldata;
		      tree->both_sides = 1;
		      tree->n_pattern  = tree->data->crunch_len/tree->mod->stepsize;
		    }

		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));


		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = alldata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) 
		    {
//#ifndef BATCH
		      Check_Memory_Amount(tree);
//#endif
		    }

		  Order_Tree_CSeq(tree,alldata);

		  if((tree->mod->s_opt->random_input_tree) && (tree->mod->s_opt->topo_search == SPR_MOVE))
		    {
		      PhyML_Printf("\n. Randomising the tree...\n");
		      Random_Tree(tree);
		    }

		  Fill_Dir_Table(tree);
		  Update_Dirs(tree);
		  Make_Tree_4_Pars(tree,alldata,alldata->init_len);
		  Make_Tree_4_Lk(tree,alldata,alldata->init_len);
		  tree->triplet_struct = Make_Triplet_Struct(mod);
		  Br_Len_Not_Involving_Invar(tree);

 		  if((tree->mod->s_opt->topo_search == SPR_MOVE) ||
		     (tree->mod->s_opt->topo_search == NNI_MOVE))
		    {
		      Make_Spr_List(tree);
		      Make_Best_Spr(tree);
		    }

		  
		  if(tree->mod->s_opt->opt_topo)
		    {
		      if(tree->mod->s_opt->topo_search == NNI_MOVE)
			{
			  Simu(tree,1000);
			}
		      else
			{
			  if(tree->mod->s_opt->steph_spr)
			    {
			      Speed_Spr(tree);
			      Simu(tree,1000);
			    }
			  else
			    {
			      Init_SPR(tree);
			      Optim_SPR(tree,0,ALL);
			      Clean_SPR(tree);
			    }
			}
		    }
		  else
		    {
		      if(tree->mod->s_opt->opt_num_param || tree->mod->s_opt->opt_bl) 
			{
			  Round_Optimize(tree,tree->data);
			}
		      else
			{
			  Lk(tree);
			  Print_Lk(tree,"");
			}
		    }

		  if(tree->io->ratio_test) aLRT(tree);

		  Lk(tree);
		  PhyML_Printf("\n\n. Final log likelihood : %f",tree->c_lnL);
		  

		  /* */
		  M4_Compute_Proba_Hidden_States_On_Edges(tree);
		  /* */



		  if((tree->c_lnL > best_lnL) && (io->mod->s_opt->n_rand_starts > 1))
		    {
		      best_lnL = tree->c_lnL;
		      io->fp_out_trees = (FILE *)fopen(io->out_trees_file,"w");
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		      fclose(io->fp_out_trees);
		    }

		  if((tree->mod->bootstrap) && (bootstrap_this_tree))
		    {
		      if(num_rand_tree > 0) io->in_tree = 0;
		      Bootstrap(tree);
		      tree->mod->bootstrap = 0;
		    }

		  Br_Len_Involving_Invar(tree);
		  Print_Tree(io->fp_out_tree,tree);

		  Unconstraint_Lk(tree);
		  time(&t_end);
		  hour = div(t_end-t_beg,3600);
		  min  = div(t_end-t_beg,60  );
		  min.quot -= hour.quot*60;
		  PhyML_Printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
		  PhyML_Printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?(num_rand_tree):(num_tree));
		  
		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1)
		     && (io->mod->s_opt->n_rand_starts > 1)
		     && (tree->mod->s_opt->random_input_tree))
		    {
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }

		  if((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) &&
		     (!tree->mod->s_opt->random_input_tree) &&
		     (io->mod->s_opt->n_rand_starts > 1))
		    {
		      if(tree->mod->bootstrap)
			{
			  num_rand_tree--;
			  io->in_tree = 2;
			  io->fp_in_tree = io->fp_out_trees;
			  bootstrap_this_tree  = 1;
			  io->fp_in_tree = (FILE *)fopen(io->out_trees_file,"r");
			}
		      else
			{
			  io->fp_out_trees = (FILE *)fopen(io->out_trees_file,"w");
			  Print_Tree(io->fp_out_trees,tree);
			  fflush(NULL);
			  fclose(io->fp_out_trees);
			}
		    }

 		  if((tree->mod->s_opt->topo_search == SPR_MOVE) ||
		     (tree->mod->s_opt->topo_search == NNI_MOVE))
		    {
		      Free_Spr_List(tree);
		      Free_One_Spr(tree->best_spr);
		    }

		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}
	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(alldata);
	}
    }

  if(io->mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n\n. Best log likelihood : %f\n",best_lnL);

  Free_Model(mod);

  if(io->fp_in_seq)    fclose(io->fp_in_seq);
  if(io->fp_in_tree)   fclose(io->fp_in_tree);
  if(io->fp_out_lk)    fclose(io->fp_out_lk);
  if(io->fp_out_tree)  fclose(io->fp_out_tree);
  if(io->fp_out_stats) fclose(io->fp_out_stats);

  Free_Input(io);
  return 0;

}

/*********************************************************/

/* Allocate memory */
m4 *M4_Make_Light(int n_o)
{
  m4 *m4mod;

  m4mod = (m4 *)mCalloc(1,sizeof(m4));
  m4mod->n_o = n_o;
  m4mod->o_rr = (phydbl *)mCalloc(n_o*n_o,sizeof(phydbl));
  m4mod->o_fq = (phydbl *)mCalloc(n_o,sizeof(phydbl));
  M4_Set_M4mod_Default(m4mod);
  return m4mod;
}

/*********************************************************/

void M4_Set_M4mod_Default(m4 *m4mod)
{
  m4mod->use_cov_alpha = 1;
  m4mod->use_cov_alpha = 0;
}

/*********************************************************/
/* Allocate memory */
void M4_Make_Complete(int n_h, int n_o, m4 *m4mod)
{
  int i;

  m4mod->n_h = n_h;
  m4mod->n_o = n_o;
  m4mod->o_mats = (phydbl **)mCalloc(n_h,sizeof(phydbl *));
  For(i,n_h) m4mod->o_mats[i] = (phydbl *)mCalloc(n_o*n_o,sizeof(phydbl));
  m4mod->h_mat = (phydbl *)mCalloc(n_h*n_h,sizeof(phydbl));
  m4mod->h_rr = (phydbl *)mCalloc(n_h*n_h,sizeof(phydbl));
  m4mod->h_fq = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->multipl = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->multipl_unscaled = (phydbl *)mCalloc(n_h,sizeof(phydbl));
  m4mod->h_fq_unscaled = (phydbl *)mCalloc(n_h,sizeof(phydbl));
}

/*********************************************************/

/* Free memory */
void M4_Free_M4_Model(m4 *m4mod)
{
  int i;
  
  For(i,m4mod->n_h) Free(m4mod->o_mats[i]);
  Free(m4mod->o_mats);
  Free(m4mod->h_mat);
  Free(m4mod->o_rr);
  Free(m4mod->h_rr);
  Free(m4mod->o_fq);
  Free(m4mod->h_fq);
  Free(m4mod->multipl);
  Free(m4mod->multipl_unscaled);
  Free(m4mod->h_fq_unscaled);
  Free(m4mod);
}

/*********************************************************/

void M4_Init_Model(m4 *m4mod, allseq *data, model *mod)
{
  int i;
  phydbl fq;

  
  m4mod->n_o = (mod->datatype == NT)?(4):(20);
  mod->ns = m4mod->n_o * m4mod->n_h;
  For(i,m4mod->n_o) m4mod->o_fq[i] = data->b_frq[i];
  For(i,(int)(m4mod->n_h)) m4mod->multipl[i] = 1.;
  For(i,(int)(m4mod->n_o*(m4mod->n_o-1)/2)) m4mod->o_rr[i] = 1.;
  For(i,(int)(m4mod->n_h*(m4mod->n_h-1)/2)) m4mod->h_rr[i] = 1.;
  fq = (phydbl)(1./m4mod->n_h);

  if(mod->s_opt->opt_cov_delta) m4mod->delta = 1.0;
  if(mod->s_opt->opt_cov_alpha) m4mod->alpha = 1.0;
  For(i,m4mod->n_h) m4mod->h_fq[i] = fq;
  For(i,m4mod->n_h) m4mod->h_fq_unscaled[i] = 1.0;
  For(i,m4mod->n_h) m4mod->multipl[i] = i;
  For(i,m4mod->n_h) m4mod->multipl_unscaled[i] = i;

  mod->update_eigen = 1;
  M4_Update_Qmat(m4mod,mod);
}

/*********************************************************/

/* Fill the (big) rate matrix of the M4 model */ 
void M4_Update_Qmat(m4 *m4mod, model *mod)
{
  int i,j;
  int n_s, n_o, n_h;
  phydbl mr, sum;

  /* The number of states in M4 models is the product 
     of the number of hidden states (or classes) by the
     number of observable states 
   */
  n_s = mod->ns;
  n_o = m4mod->n_o;
  n_h = m4mod->n_h;
  
  /* Set the relative substitution rates */
  if(mod->m4mod->use_cov_alpha)
    {
      DiscreteGamma(m4mod->h_fq,m4mod->multipl,m4mod->alpha,m4mod->alpha,m4mod->n_h,m4mod->gamma_median);
    }
  else if(mod->m4mod->use_cov_free)
    {
      sum = .0;
      For(i,mod->m4mod->n_h) sum += fabs(mod->m4mod->h_fq_unscaled[i]);
      For(i,mod->m4mod->n_h) mod->m4mod->h_fq[i] = fabs(mod->m4mod->h_fq_unscaled[i])/sum;
      
      do
	{
	  sum = .0;
	  For(i,mod->m4mod->n_h)
	    {
	      if(mod->m4mod->h_fq[i] < 0.01) mod->m4mod->h_fq[i]=0.01;
	      if(mod->m4mod->h_fq[i] > 0.99) mod->m4mod->h_fq[i]=0.99;
	      sum += mod->m4mod->h_fq[i];
	    }
	  For(i,mod->m4mod->n_h) mod->m4mod->h_fq[i]/=sum;
	}
      while((sum > 1.01) || (sum < 0.99));


      /* Make sure the multipliers are centered around 1.0 */
      sum = .0;
      For(i,mod->m4mod->n_h) sum += fabs(mod->m4mod->multipl_unscaled[i]) * mod->m4mod->h_fq[i];
      For(i,mod->m4mod->n_h) mod->m4mod->multipl[i] = mod->m4mod->multipl_unscaled[i] / sum;
      

/*       mod->m4mod->h_fq[0] = 0.33; */
/*       mod->m4mod->h_fq[1] = 0.33; */
/*       mod->m4mod->h_fq[2] = 0.33; */

/*       mod->m4mod->multipl[0] =  0.1; */
/*       mod->m4mod->multipl[1] =  1.0; */
/*       mod->m4mod->multipl[2] = 10.0; */


      sum = 0;
      For(i,mod->m4mod->n_h) sum += mod->m4mod->multipl[i] * mod->m4mod->h_fq[i];
      if(sum < 0.99 || sum > 1.01)
	{
	  PhyML_Printf("\n. sum = %f",sum);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}

/*       PhyML_Printf("\n__ "); */
/*       For(i,mod->m4mod->n_h) PhyML_Printf("\n.%f %f %f", */
/* 				    mod->m4mod->h_fq[i], */
/* 				    mod->m4mod->h_fq_unscaled[i], */
/* 				    mod->m4mod->multipl[i]); */
    }


/*   PhyML_Printf("\n."); */
/*   PhyML_Printf("\n. M4 model parameters"); */
/*   PhyML_Printf("\n. Delta = %f",m4mod->delta); */
/*   For(i,mod->m4mod->n_h) PhyML_Printf("\n. multipl %d = %f",i,m4mod->multipl[i]); */
/*   For(i,mod->m4mod->n_h) PhyML_Printf("\n. fq %d = %f",i,m4mod->h_fq[i]); */


  /* Set up the stationary frequency vector */
  For(i,n_s) mod->pi[i] = m4mod->o_fq[i%n_o] * m4mod->h_fq[i/n_o];

  /* Fill the matrices of nucleotide substitution rates here */
  Update_Qmat_Generic(m4mod->o_rr, m4mod->o_fq, m4mod->n_o, m4mod->o_mats[0]);

  /* Multiply each of these matrix by a relative substitution rate */
  for(i=1;i<m4mod->n_h;i++) For(j,n_o*n_o) m4mod->o_mats[i][j] = m4mod->o_mats[0][j]*m4mod->multipl[i];
  For(j,n_o*n_o) m4mod->o_mats[0][j] *= m4mod->multipl[0];

  For(i,n_s*n_s) mod->qmat[i] = .0;

  /* Diagonal blocks (i.e, nucleotide substitutions), symmetric */
  For(i,n_s)
    {
      for(j=i+1;j<n_s;j++)
	{
	  if((int)(j/n_o) == (int)(i/n_o))
	    {
	      mod->qmat[i*n_s+j] = m4mod->o_mats[(int)(i/n_o)][(i%n_o)*n_o+j%n_o];
	      mod->qmat[j*n_s+i] = mod->qmat[i*n_s+j] * m4mod->o_fq[i%n_o] / m4mod->o_fq[j%n_o];
	    }
	}
    }

  /* Note: nucleotide equilibrium frequencies are already built in the o_mats matrices.
     No need to 'add' these frequencies later on. */

  /* Work out scaling factor */
  mr = .0;
  For(i,n_s)
    {
      sum = .0;
      For(j,n_s) sum += mod->qmat[i*n_s+j];
      mr += sum * m4mod->o_fq[i%n_o] * m4mod->h_fq[(int)(i/n_o)]; 
    }
  
  /* Scale the diagonal blocks */
  For(i,n_s*n_s) mod->qmat[i] /= mr;
  
  /* We are done with the diagonal blocks. Let's fill the non-diagonal ones now. */

  /* Fill the matrix of substitution rate across classes (switches) here */
  Update_Qmat_Generic(m4mod->h_rr, m4mod->h_fq, m4mod->n_h, m4mod->h_mat);

/*   Print_Square_Matrix_Generic(m4mod->n_h,m4mod->h_mat); */

  /* Multiply this matrix by the switching rate */
  For(i,n_h*n_h) m4mod->h_mat[i] *= m4mod->delta;

  /* Fill the non diagonal blocks */
  For(i,n_s)
    {
      for(j=i+1;j<n_s;j++)
	{
	  if((int)(j/n_o) != (int)(i/n_o))
	    {
	      if(i%n_o == j%n_o)
		{
		  mod->qmat[i*n_s+j] = m4mod->h_mat[(int)(i/n_o)*n_h+(int)(j/n_o)];
		  mod->qmat[j*n_s+i] = mod->qmat[i*n_s+j] * m4mod->h_fq[(int)(i/n_o)] / m4mod->h_fq[(int)(j/n_o)]; 
		}
	    }
	}
    }

  /* Note: class equilibrium frequencies are already built in the h_mat matrix.
     No need to 'add' these frequencies later on. */

  /* We are done with the non diagonal blocks */

  /* Diagonal cells */
  For(i,n_s)
    {
      sum = .0;
      For(j,n_s)
	{
	  if(j != i)
	    sum += mod->qmat[i*n_s+j];
	}
      mod->qmat[i*n_s+i] = -sum;
    }

  /* Print_Square_Matrix_Generic(n_s,mod->qmat); */
  For(i,n_s*n_s) mod->eigen->q[i] = mod->qmat[i];
}

/*********************************************************/

void M4_Init_P_Lk_Tips_Double(arbre *tree)
{
  int curr_site,i,j,k,l;
  
  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      For(i,tree->n_otu)
	{
	  for(j=1;j<tree->mod->m4mod->n_h;j++)
	    {
	      For(k,tree->mod->m4mod->n_o)
		tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][j*tree->mod->m4mod->n_o+k] = 
		tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][k];
	      
	      For(k,tree->mod->m4mod->n_o)
		for(l=1;l<tree->mod->n_catg;l++)
		  tree->noeud[i]->b[0]->p_lk_rght[curr_site][l][j*tree->mod->m4mod->n_o+k] = 
		  tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][j*tree->mod->m4mod->n_o+k];
	    }
	}
    }
}

/*********************************************************/

void M4_Init_P_Lk_Tips_Int(arbre *tree)
{
  int curr_site,i,j,k;

  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      For(i,tree->n_otu)
	{
	  for(j=1;j<tree->mod->m4mod->n_h;j++)
	    {
	      For(k,tree->mod->m4mod->n_o)
		tree->noeud[i]->b[0]->p_lk_tip_r[curr_site][j*tree->mod->m4mod->n_o+k] = 
		tree->noeud[i]->b[0]->p_lk_tip_r[curr_site][k];
	    }
	}
    }
}

/*********************************************************/

phydbl ****M4_Integral_Term_On_One_Edge(edge *b, arbre *tree)
{
  phydbl ****integral,***P1,***P2;  
  int ns;
  int g,i,j,k,l;
  int step;

  ns = tree->mod->ns;

  P1 = (phydbl ***)mCalloc(tree->mod->n_catg,sizeof(phydbl **));
  For(g,tree->mod->n_catg) 
    {
      P1[g] = (phydbl **)mCalloc(ns,sizeof(phydbl *));
      For(j,ns) P1[g][j] = (phydbl *)mCalloc(ns,sizeof(phydbl));
    }

  P2 = (phydbl ***)mCalloc(tree->mod->n_catg,sizeof(phydbl **));
  For(g,tree->mod->n_catg) 
    {
      P2[g] = (phydbl **)mCalloc(ns,sizeof(phydbl *));
      For(j,ns) P2[g][j] = (phydbl *)mCalloc(ns,sizeof(phydbl));
    }


  integral = (phydbl ****)mCalloc(tree->mod->n_catg,sizeof(phydbl ***));
  For(g,tree->mod->n_catg)
    {
      integral[g] = (phydbl ***)mCalloc(ns,sizeof(phydbl **));
      For(j,ns)
	{
	  integral[g][j] = (phydbl **)mCalloc(ns,sizeof(phydbl *));
	  For(k,ns) integral[g][j][k] = (phydbl *)mCalloc(ns,sizeof(phydbl));
	}
    }

  /* Integral calculation */
  step = 100;

  PhyML_Printf("\n. [");
  For(i,step)
    {
      For(g,tree->mod->n_catg)
	{
	  PMat(((phydbl)(i+0.5)/step)*b->l*tree->mod->gamma_rr[g],tree->mod,P1+g);
	  PMat(((phydbl)(step-i-0.5)/step)*b->l*tree->mod->gamma_rr[g],tree->mod,P2+g);

	  For(j,ns)
	    {
	      For(k,ns)
		{
		  For(l,ns)
		    {
		      integral[g][j][k][l] += P1[g][j][k] * P2[g][j][l]  / ((phydbl)(step));
		    }
		}
	    }      
	}
      PhyML_Printf("."); fflush(NULL);
    }
  PhyML_Printf("]\n");

  For(g,tree->mod->n_catg)
    {
      For(i,ns) Free(P1[g][i]); 
      Free(P1[g]);
    }
  Free(P1);

  For(g,tree->mod->n_catg)
    {
      For(i,ns) Free(P2[g][i]);
      Free(P2[g]);
    }
  Free(P2);

  return integral;
}

/*********************************************************/

void M4_Post_Prob_H_Class_Edge_Site(edge *b, phydbl ****integral, phydbl *postprob, arbre *tree)
{
  /* Calculation of the expected frequencies of each hidden
     class at a given site. */

  phydbl site_lk;
  int g,i,j,k,l;
  int ns,n_h;
  phydbl sum;

  ns = tree->mod->ns;
  n_h = tree->mod->m4mod->n_h; /* number of classes, i.e., number of hidden states */

  site_lk = (phydbl)exp(tree->site_lk[tree->curr_site]);

  if(b->rght->tax)
    {
      sum = .0;
      For(i,n_h)
	{
	  postprob[i] = .0;
	  For(j,tree->mod->m4mod->n_o)
	    {
	      For(g,tree->mod->n_catg)
		{
		  For(k,tree->mod->ns)
		    {
		      For(l,tree->mod->ns)
			{
			  postprob[i] +=
			    (1./site_lk) *
			    tree->mod->gamma_r_proba[g] *
			    tree->mod->m4mod->h_fq[i] *
			    tree->mod->m4mod->o_fq[j] *
			    b->p_lk_left[tree->curr_site][g][k] *
			    b->p_lk_tip_r[tree->curr_site][l] *
			    /* 			b->p_lk_rght[tree->curr_site][0][l] * */
			    integral[g][i*tree->mod->m4mod->n_o+j][k][l];
			}
		    }
		}
	    }
	  sum += postprob[i];
	}
      For(i,n_h) postprob[i] *= exp(b->sum_scale_f_left[tree->curr_site]); 
    }
  else
    {
      sum = .0;
      For(i,n_h)
	{
	  postprob[i] = .0;
	  For(j,tree->mod->m4mod->n_o)
	    {
	      For(g,tree->mod->n_catg)
		{
		  For(k,tree->mod->ns)
		    {
		      For(l,tree->mod->ns)
			{
			  postprob[i] +=
			    (1./site_lk) *
			    tree->mod->gamma_r_proba[g] *
			    tree->mod->m4mod->h_fq[i] *
			    tree->mod->m4mod->o_fq[j] *
			    b->p_lk_left[tree->curr_site][g][k] *
			    b->p_lk_rght[tree->curr_site][g][l] *
			    integral[g][i*tree->mod->m4mod->n_o+j][k][l];
			}
		    }
		}
	    }
	  sum += postprob[i];
	}
      For(i,n_h) postprob[i] *= exp(b->sum_scale_f_left[tree->curr_site] + b->sum_scale_f_rght[tree->curr_site]); 
    }

  For(i,n_h) 
    if((postprob[i] < -1.E-5) || (postprob[i] > 1.0+1.E-5))
      {
	PhyML_Printf("\n. Cat : %d Prob : %f\n",i,postprob[i]);
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("\n");
      }

  sum = 0.0;
  For(i,n_h) sum += postprob[i];

  if((sum > 1.0+1.E-2) || (sum < 1.0-1.E-2))
    {
      PhyML_Printf("\n. Sum = %f\n",sum);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return;
}

/*********************************************************/

phydbl ***M4_Compute_Proba_Hidden_States_On_Edges(arbre *tree)
{
  int i;
  phydbl ***post_probs, *dwell;
  phydbl ****integral;


  dwell = (phydbl *)mCalloc(tree->mod->m4mod->n_h,sizeof(phydbl));
  post_probs = (phydbl ***)mCalloc(2*tree->n_otu-3,sizeof(phydbl **));

  For(i,2*tree->n_otu-3)
    {
      post_probs[i] = (phydbl **)mCalloc(tree->n_pattern,sizeof(phydbl *));
      For(tree->curr_site,tree->n_pattern) 
	post_probs[i][tree->curr_site] = (phydbl *)mCalloc(tree->mod->m4mod->n_h,sizeof(phydbl));
    }


  /* Compute posterior probabilities of each hidden class (usually a rate class) 
     on each edge, at each site.
  */
  For(i,2*tree->n_otu-3) 
    {
      PhyML_Printf("\n. Edge %4d/%4d",i+1,2*tree->n_otu-3);

      integral = M4_Integral_Term_On_One_Edge(tree->t_edges[i],tree);

      For(tree->curr_site,tree->n_pattern)
	M4_Post_Prob_H_Class_Edge_Site(tree->t_edges[i],
				       integral,
				       post_probs[i][tree->curr_site],
				       tree);
      
      M4_Free_Integral_Term_On_One_Edge(integral,tree);
    }
  return post_probs;
}

/*********************************************************/

/* Estimate the (posterior) mean relative rate of substitution on each branch
   at each site. The posterior mean rates averaged over sites is also estimated
   for each edge. The corresponding trees are printed in a postscript file. Tree 0
   is the tree with posterior mean rates averaged over the sites. The following trees
   have posterior mean rates computed for each site.
*/
void M4_Compute_Posterior_Mean_Rates(phydbl ***post_probs, arbre *tree)
{
  char *s;
  int i;
  phydbl **mean_post_probs;
  phydbl *mrr;
  phydbl sum;
  int patt,br,rcat;
  phydbl *mean_br_len;
  int best_r,len_var;
  phydbl max_prob;

  mean_br_len = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  mean_post_probs = (phydbl **)mCalloc(2*tree->n_otu-3,sizeof(phydbl *));
  For(i,2*tree->n_otu-3) mean_post_probs[i] = (phydbl *)mCalloc(tree->mod->m4mod->n_h,sizeof(phydbl ));
  mrr = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  Record_Br_Len(NULL,tree);
  M4_Scale_Br_Len(tree);

  /* Compute the posterior mean relative rate on each branch averaged over the 
     whole set of patterns (sites) */
  len_var = 0;
  For(patt,tree->n_pattern) 
    {
      if(!Is_Invar(patt,1,NT,tree->data))
	{
	  For(br,2*tree->n_otu-3)
	    {
	      max_prob = -1.;
	      best_r = -1;
	      For(rcat,tree->mod->m4mod->n_h)
		{
		  if(post_probs[br][patt][rcat] > max_prob)
		    {
		      max_prob = post_probs[br][patt][rcat];
		      best_r = rcat;
		    }
		}

/* /\* 	      Add weight on each category, weight is proportional to the corresponding posterior probability *\/ */
/* 	      For(rcat,tree->mod->m4mod->n_h) */
/* 		{ */
/* 		  mean_post_probs[br][rcat] += post_probs[br][patt][rcat] * tree->data->wght[patt]; */
/* 		} */

	      /* Add weight on the most probable rate category only */
	      mean_post_probs[br][best_r] += tree->data->wght[patt];
	    }
	  len_var += tree->data->wght[patt];
	}
    }

  For(br,2*tree->n_otu-3) 
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  mean_post_probs[br][rcat] /= (phydbl)len_var;
	}
    }

  /* Compute the posterior mean relative rate and scale
     each branch length using this factor */
  For(br,2*tree->n_otu-3)
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  mrr[br] += mean_post_probs[br][rcat] * tree->mod->m4mod->multipl[rcat];
	}
      tree->t_edges[br]->l *= mrr[br];
    }

  PhyML_Fprintf(tree->io->fp_out_stats,"\n. Mean posterior probabilities & rates\n");
  For(rcat,tree->mod->m4mod->n_h) PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",tree->mod->m4mod->multipl[rcat]);
  PhyML_Fprintf(tree->io->fp_out_stats,"\n");
  For(br,2*tree->n_otu-3) 
    {
      For(rcat,tree->mod->m4mod->n_h)
	{
	  PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",mean_post_probs[br][rcat]);
	}
/*       PhyML_Fprintf(tree->io->fp_out_stats," -- %f -> %f x %f = %f",mrr[br],tree->t_edges[br]->l,mrr[br],tree->t_edges[br]->l*mrr[br]); */

      PhyML_Fprintf(tree->io->fp_out_stats," mrr=%f ",mrr[br]);

      if(mrr[br] > 1.) PhyML_Fprintf(tree->io->fp_out_stats,"FAST ");
      else             PhyML_Fprintf(tree->io->fp_out_stats,"SLOW ");
      
      PhyML_Fprintf(tree->io->fp_out_stats,"%s",tree->t_edges[br]->labels[0]);

      PhyML_Fprintf(tree->io->fp_out_stats,"\n");
    }

  /* Print the tree */
  PhyML_Fprintf(tree->io->fp_out_tree,"Constrained tree with corrected branch lengths = ");
  s = Write_Tree(tree);
  PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
  Free(s);
  tree->ps_tree = DR_Make_Tdraw_Struct(tree);
  DR_Print_Postscript_Header(tree->n_pattern,tree->io->fp_out_ps);
  tree->ps_page_number = 0;
  DR_Print_Tree_Postscript(tree->ps_page_number++,tree->io->fp_out_ps,tree);

  /* Go back to the initial scaled branch lengths */
  For(br,2*tree->n_otu-3) tree->t_edges[br]->l /= mrr[br];

  /* Compute the posterior mean relative rate at each site, for each branch
     and each rate category. Scale branch lengths using these factors and
     print each tree (i.e., on tree per site pattern) */
  For(patt,tree->n_pattern) 
    {
      For(br,2*tree->n_otu-3) 
	{
	  mrr[br] = .0;
	  max_prob = -1.;
	  best_r = -1;
	  For(rcat,tree->mod->m4mod->n_h) /* For each rate class */
	    {
	      mrr[br] += post_probs[br][patt][rcat] * tree->mod->m4mod->multipl[rcat];
	      if(post_probs[br][patt][rcat] > max_prob)
		{
		  max_prob = post_probs[br][patt][rcat];
		  best_r = rcat;
		}
	    }
/* 	  mrr[br] = tree->mod->m4mod->multipl[best_r]; /\* Use the most probable relative rate instead of mean *\/ */
	  tree->t_edges[br]->l *= mrr[br];
	}

      For(br,2*tree->n_otu-3) mean_br_len[br] += tree->t_edges[br]->l * tree->data->wght[patt];

      PhyML_Fprintf(tree->io->fp_out_stats,"\n. Posterior probabilities site %4d (weight=%d, is_inv=%d)\n",
	     patt,
	     tree->data->wght[patt],
	     Is_Invar(patt,1,NT,tree->data));

      For(rcat,tree->mod->m4mod->n_h) PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f ",tree->mod->m4mod->multipl[rcat]);
      PhyML_Fprintf(tree->io->fp_out_stats,"\n");
      For(br,2*tree->n_otu-3)
	{
	  PhyML_Fprintf(tree->io->fp_out_stats,"Edge %3d ",br);
	  max_prob = -1.0;
	  best_r = -1;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      if(post_probs[br][patt][rcat] > max_prob)
		{
		  max_prob = post_probs[br][patt][rcat];
		  best_r = rcat;
		}
	    }

	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      PhyML_Fprintf(tree->io->fp_out_stats,"%2.4f",post_probs[br][patt][rcat]);
	      if(rcat == best_r) PhyML_Fprintf(tree->io->fp_out_stats,"* ");
	      else               PhyML_Fprintf(tree->io->fp_out_stats,"  ");
	    }

/* 	  PhyML_Fprintf(tree->io->fp_out_stats," -- %f -> %f x %f = %f",mrr[br],tree->t_edges[br]->l,mrr[br],tree->t_edges[br]->l*mrr[br]); */
	  
	  if(mrr[br] > 1.01)      PhyML_Fprintf(tree->io->fp_out_stats," %s ","FAST");
	  else if(mrr[br] < 0.99) PhyML_Fprintf(tree->io->fp_out_stats," %s ","SLOW");
	  else 	                  PhyML_Fprintf(tree->io->fp_out_stats," %s ","MEDIUM");
	  PhyML_Fprintf(tree->io->fp_out_stats,"%s ",tree->t_edges[br]->labels[0]);
	  PhyML_Fprintf(tree->io->fp_out_stats,"\n");
	}

      PhyML_Fprintf(tree->io->fp_out_tree,"tree %d = ",patt+1);
      s = Write_Tree(tree);
      PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
      Free(s);
      DR_Print_Tree_Postscript(tree->ps_page_number++,tree->io->fp_out_ps,tree);

      /* Go back to the initial scaled branch lengths */
      For(br,2*tree->n_otu-3) tree->t_edges[br]->l /= mrr[br];

      For(br,2*tree->n_otu-3) 
	{
	  sum = .0;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
	      sum += post_probs[br][patt][rcat];
	    }
	  
	  if((sum < 0.99) || (sum > 1.01))
	    {
	      PhyML_Fprintf(tree->io->fp_out_stats,"\n. sum = %f\n",sum);
	      PhyML_Fprintf(tree->io->fp_out_stats,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }
  
  /* Mean branch lengths */
  For(br,2*tree->n_otu-3)
    {
      mean_br_len[br] /= (phydbl)tree->data->init_len;
      tree->t_edges[br]->l = mean_br_len[br];
    }
  PhyML_Fprintf(tree->io->fp_out_tree,"Mean branch lengths=");
  s = Write_Tree(tree);
  PhyML_Fprintf(tree->io->fp_out_tree,"%s\n",s);
  Free(s);
/*   DR_Print_Tree_Postscript(tree->ps_page_number++,tree->io->fp_out_ps,tree); */

  Restore_Br_Len(NULL,tree);

  DR_Print_Postscript_EOF(tree->io->fp_out_ps);

  For(br,2*tree->n_otu-3)
    {
      For(tree->curr_site,tree->n_pattern)
	Free(post_probs[br][tree->curr_site]);
      Free(post_probs[br]);
    }
  Free(post_probs);
  For(i,2*tree->n_otu-3) Free(mean_post_probs[i]);
  Free(mean_post_probs);
  Free(mrr);
  Free(mean_br_len);
}

/*********************************************************/

/* Classifiy each branch, at each site, among one of the rate classes */
phydbl **M4_Site_Branch_Classification(phydbl ***post_probs, arbre *tree)
{
  int patt, br, rcat, i;
  phydbl **best_probs;
  phydbl post_prob_fast, post_prob_slow;

  best_probs = (phydbl **)mCalloc(tree->n_pattern,sizeof(phydbl *));
  For(i,tree->n_pattern) best_probs[i] = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  tree->print_labels = 1;

  For(patt,tree->n_pattern) 
    {
      For(br,2*tree->n_otu-3) 
	{
	  post_prob_fast = .0;
	  post_prob_slow = .0;

	  For(rcat,tree->mod->m4mod->n_h) /* For each rate class */
	    {	      
	      if(tree->mod->m4mod->multipl[rcat] > 1.0) 
		post_prob_fast += post_probs[br][patt][rcat];
	      else
		post_prob_slow += post_probs[br][patt][rcat];
	    }

	  best_probs[patt][br] = (post_prob_fast > post_prob_slow)?(post_prob_fast):(post_prob_slow);

	  if(!(tree->t_edges[br]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(tree->t_edges[br]);

/* 	  if((post_prob_fast > post_prob_slow) && (best_probs[patt][br] > 0.95)) */
/* 	    strcpy(tree->t_edges[br]->labels[tree->t_edges[br]->n_labels],"FASTER"); */
/* 	  else if((post_prob_fast < post_prob_slow) && (best_probs[patt][br] > 0.95)) */
/* 	    strcpy(tree->t_edges[br]->labels[tree->t_edges[br]->n_labels],"SLOWER"); */
/* 	  else */
/* 	    strcpy(tree->t_edges[br]->labels[tree->t_edges[br]->n_labels],"UNKNOWN"); */

	  if(post_prob_fast > post_prob_slow)
	    strcpy(tree->t_edges[br]->labels[tree->t_edges[br]->n_labels],"FASTER");
	  else 
	    strcpy(tree->t_edges[br]->labels[tree->t_edges[br]->n_labels],"SLOWER");

	  tree->t_edges[br]->n_labels++;
	}
    }
  return best_probs;
}

/*********************************************************/

void M4_Site_Branch_Classification_Experiment(arbre *tree)
{
  allseq *ori_data,*cpy_data;
  short int **true_rclass, **est_rclass;
  phydbl **best_probs;
  int i,j;
  phydbl correct_class, mis_class, unknown;


  true_rclass = (short int **)mCalloc(tree->data->init_len, sizeof(short int *));
  est_rclass  = (short int **)mCalloc(tree->data->init_len, sizeof(short int *));
 
  For(i,tree->data->init_len)
    {
      true_rclass[i] = (short int *)mCalloc(2*tree->n_otu-3,sizeof(short int));
      est_rclass[i]  = (short int *)mCalloc(2*tree->n_otu-3,sizeof(short int));
    }

  ori_data = tree->data;

  cpy_data = Copy_Cseq(tree->data,
		       tree->data->init_len,
		       (tree->mod->datatype == NT)?(4):(20));

  /* Generate a simulated data set under H0, with the right sequence length. */
  PhyML_Printf("\n. Evolving sequences (delta=%f, alpha=%f) ...\n",tree->mod->m4mod->delta,tree->mod->m4mod->alpha);
  Evolve(cpy_data,tree->mod,tree);

  For(i,cpy_data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
	  if(!strcmp(tree->t_edges[j]->labels[i],"FASTER"))
	    {
	      true_rclass[i][j] = 1;
	    }
	  else if(!strcmp(tree->t_edges[j]->labels[i],"SLOWER"))
	    {
	      true_rclass[i][j] = 0;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }
  
  For(j,2*tree->n_otu-3) 
    {
      Free_Edge_Labels(tree->t_edges[j]);
      tree->t_edges[j]->n_labels = 0;
    }

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays
  */
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

  /* Estimate model parameters */
  PhyML_Printf("\n. Estimating model parameters...\n");
  tree->mod->s_opt->opt_cov_alpha = 1;
  tree->mod->s_opt->opt_cov_delta = 1;
  Round_Optimize(tree,tree->data);

  tree->both_sides = 1;
  Lk(tree);

  /* Classify branches */
  best_probs = M4_Site_Branch_Classification(M4_Compute_Proba_Hidden_States_On_Edges(tree),tree);

  For(i,tree->data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
	  if(!strcmp(tree->t_edges[j]->labels[i],"FASTER"))
	    {
	      est_rclass[i][j] = 1;
	    }
	  else if(!strcmp(tree->t_edges[j]->labels[i],"SLOWER"))
	    {
	      est_rclass[i][j] = 0;
	    }
	  else if(!strcmp(tree->t_edges[j]->labels[i],"UNKNOWN"))
	    {
	      est_rclass[i][j] = -1;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
    }

  unknown       = .0;
  correct_class = .0;
  mis_class     = .0;
  For(i,tree->data->init_len)
    {
      For(j,2*tree->n_otu-3)
	{
/* 	  PhyML_Printf("\n. Edge %3d %4d %4d - %f",j,true_rclass[i][j],est_rclass[i][j],best_probs[i][j]); */
	  if(est_rclass[i][j] == -1)
	    {
	      unknown += 1.;
	    }
	  else if(est_rclass[i][j] == true_rclass[i][j])
	    {
	      correct_class += 1.;
	    }
	  else if(est_rclass[i][j] != true_rclass[i][j])
	    {
	      mis_class += 1.;
	    }
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
/*       PhyML_Printf("\n"); */
    }

  correct_class /= ((tree->data->init_len * (2*tree->n_otu-3)) - unknown);
  mis_class     /= ((tree->data->init_len * (2*tree->n_otu-3)) - unknown);
  unknown       /= (tree->data->init_len  * (2*tree->n_otu-3));
  
  PhyML_Printf("\n. correct_class = %f mis_class = %f unknown = %f",
	 correct_class,
	 mis_class,
	 unknown);


  For(i,tree->data->init_len)
    {
      Free(true_rclass[i]);
      Free(est_rclass[i]);
      Free(best_probs[i]);
    }
  Free(true_rclass);
  Free(est_rclass);
  Free(best_probs);

}

/*********************************************************/

  /* Scale branch lengths such that they express expected number
     of nucleotide or amino-acid substitutions */

void M4_Scale_Br_Len(arbre *tree)
{
  phydbl scale_fact,mrs;
  int i,j;

  /* (1) Work out the relative mean rate of switches */
  mrs = .0;
  For(i,tree->mod->m4mod->n_h)
    {
      For(j,tree->mod->m4mod->n_h)
	{
	  if(j != i)
	    mrs += tree->mod->m4mod->h_fq[i] * tree->mod->m4mod->h_mat[i*tree->mod->m4mod->n_h+j];
	}
    }

  /* (2) scale_fact = (1 + delta x mrs) */
  scale_fact = 1.0 + tree->mod->m4mod->delta * mrs;

  /* (3) Scale branch lengths */
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= scale_fact;
}

/*********************************************************/

void M4_Free_Integral_Term_On_One_Edge(phydbl ****integral, arbre *tree)
{
  int g,i,j;

  For(g,tree->mod->n_catg)
    {
      For(i,tree->mod->m4mod->n_h)
	{
	  For(j,tree->mod->m4mod->n_h)
	    {
	      Free(integral[g][i][j]);
	    }
	  Free(integral[g][i]);
	}      
      Free(integral[g]);
    }
  Free(integral);
}

/*********************************************************/

void M4_Detect_Site_Switches_Experiment(arbre *tree)
{
  model *nocov_mod,*cov_mod,*ori_mod;
  allseq *ori_data,*cpy_data;
  int i,n_iter;
  phydbl *nocov_bl,*cov_bl;
  phydbl *site_lnl_nocov, *site_lnl_cov;

  nocov_bl       = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  cov_bl         = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
  site_lnl_nocov = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));
  site_lnl_cov   = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  ori_data = tree->data;
  ori_mod  = tree->mod;
  cpy_data = Copy_Cseq(tree->data,
		       tree->data->init_len,
		       (tree->mod->datatype == NT)?(4):(20));

  PhyML_Printf("\n. Estimate model parameters under non-switching substitution model.\n");
  Switch_From_M4mod_To_Mod(tree->mod);
  Simu_Loop(tree);
  nocov_mod = (model *)Copy_Model(tree->mod); /* Record model parameters */
  For(i,2*tree->n_otu-3) nocov_bl[i] = tree->t_edges[i]->l; /* Record branch lengths */
  For(i,tree->data->crunch_len) site_lnl_nocov[i] = tree->site_lk[i];
  Print_Lk(tree,"[LnL under non-switching substitution model]");
  
  PhyML_Printf("\n. Estimate model parameters under switching substitution model.\n");
  Switch_From_Mod_To_M4mod(tree->mod);
  Simu_Loop(tree);
  cov_mod = (model *)Copy_Model(tree->mod); /* Record model parameters */
  For(i,2*tree->n_otu-3) cov_bl[i] = tree->t_edges[i]->l; /* Record branch lengths */
  For(i,tree->data->crunch_len) site_lnl_cov[i] = tree->site_lk[i];
  Print_Lk(tree,"[LnL under switching substitution model]");
  

  PhyML_Printf("\n");
  For(i,tree->data->crunch_len) PhyML_Printf("TRUTH %f %f\n",site_lnl_nocov[i],site_lnl_cov[i]);

  /* Generate a simulated data set under H0, with the right sequence length. */
  tree->mod = nocov_mod;
  Evolve(cpy_data, nocov_mod, tree);

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays 
  */
  tree->mod = cov_mod;
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;
  tree->mod       = cov_mod;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

 
  n_iter = 0;
  do
    {
      /* Get the transition proba right to generate sequences */
      tree->mod = nocov_mod;
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = nocov_bl[i];
      For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
      
      /* Generate sequences */
      Evolve(cpy_data, nocov_mod, tree);
      tree->data = cpy_data;

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);
      
      tree->mod = nocov_mod;
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = nocov_bl[i];
      Lk(tree);
      For(i,tree->data->crunch_len) site_lnl_nocov[i] = tree->site_lk[i];
      Print_Lk(tree,"[CPY LnL under non-switching substitution model]");

      tree->mod = cov_mod;
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = cov_bl[i];
      Lk(tree);
      For(i,tree->data->crunch_len) site_lnl_cov[i] = tree->site_lk[i];
      Print_Lk(tree,"[CPY LnL under switching substitution model]");

      PhyML_Printf("\n");
      For(i,tree->data->crunch_len) PhyML_Printf("SYNTH %f %f\n",site_lnl_nocov[i],site_lnl_cov[i]);
    }
  while(++n_iter < 200);

  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  /* Back to the original data set to check that everything is ok */
  tree->data      = ori_data;
  tree->n_pattern = tree->data->crunch_len;

  Make_Tree_4_Pars(tree,ori_data,ori_data->init_len);
  Make_Tree_4_Lk(tree,ori_data,ori_data->init_len);

  tree->mod = nocov_mod;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l = nocov_bl[i];  
  Lk(tree);
  Print_Lk(tree,"[FINAL LnL under non-switching substitution model]");
  
  tree->mod = cov_mod;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l = cov_bl[i];  
  Lk(tree);
  Print_Lk(tree,"[FINAL LnL under switching substitution model]");

  tree->mod = ori_mod;

  Free_Model(cov_mod);
  Free_Model(nocov_mod);
  Free(site_lnl_cov);
  Free(site_lnl_nocov);

  Free_Cseq(cpy_data);
  Free(nocov_bl);
  Free(cov_bl);
}

/*********************************************************/

void M4_Posterior_Prediction_Experiment(arbre *tree)
{
  model *ori_mod;
  allseq *ori_data,*cpy_data;
  int i,n_iter,n_simul;
  FILE *fp_nocov,*fp_cov,*fp_obs;
  char *s;
  edge *best_edge;

  s = (char *)mCalloc(100,sizeof(char));

  strcpy(s,tree->io->in_seq_file);
  fp_nocov = Openfile(strcat(s,"_nocov"),1);
  strcpy(s,tree->io->in_seq_file);
  fp_cov   = Openfile(strcat(s,"_cov"),1);
  strcpy(s,tree->io->in_seq_file);
  fp_obs = Openfile(strcat(s,"_obs"),1);
  
  Free(s);

  Print_Diversity_Header(fp_nocov, tree);
  Print_Diversity_Header(fp_cov, tree);
  Print_Diversity_Header(fp_obs, tree);

  ori_data = tree->data;
  ori_mod  = tree->mod;

  cpy_data = Copy_Cseq(tree->data,
		       tree->data->init_len,
		       (tree->mod->datatype == NT)?(4):(20));

  /* Generate a simulated data set under H0, with the right sequence length. */
  Set_Model_Parameters(tree->mod);
  For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);
  Evolve(cpy_data,tree->mod,tree);

  /* Generate the memory needed for likelihood calculation because
     we will need bigger arrays
  */
  Free_Tree_Lk(tree);
  Free_Tree_Pars(tree);

  tree->data      = cpy_data;
  tree->n_pattern = tree->data->crunch_len;

  /* Allocate memory and initialize likelihood structure with
     simulated sequences (i.e., columns are not compressed)
  */
  Make_Tree_4_Pars(tree,cpy_data,cpy_data->init_len);
  Make_Tree_4_Lk(tree,cpy_data,cpy_data->init_len);

  /* Go back to the original data set */
  tree->data      = ori_data;
  tree->n_pattern = ori_data->crunch_len;
  
  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
  else Init_P_Lk_Tips_Int(tree);

  PhyML_Printf("\n. Estimate model parameters under non-switching substitution model.\n");
  Switch_From_M4mod_To_Mod(tree->mod);

  tree->bl_from_node_stamps = 1;
  best_edge = MC_Find_Best_Root_Position(tree);
  PhyML_Printf("\n. Put root on edge %3d",i);
  MC_Least_Square_Node_Times(best_edge,tree);
  MC_Adjust_Node_Times(tree);
  MC_Round_Optimize(tree);

/*   Round_Optimize(tree,tree->data); */
/*   Simu_Loop(tree); */
  Print_Lk(tree,"[LnL under non-switching substitution model]");
  Print_Tree(tree->io->fp_out_tree,tree);
  
  /* Print observed diversities */
  Init_Ui_Tips(tree);
  Site_Diversity(tree);
  Print_Diversity(fp_obs,tree);

  n_simul = 100;
  n_iter = 0;
  do
    {
      Evolve(cpy_data,tree->mod,tree);
      tree->data      = cpy_data;
      tree->n_pattern = cpy_data->init_len;

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);

      Lk(tree);

      Init_Ui_Tips(tree);
      Site_Diversity(tree);
      Print_Diversity(fp_nocov,tree);

      Print_Lk(tree,"[CPY under non-switching substitution model]");
    }while(++n_iter < n_simul);


  /* Go back to the original data set */
  tree->data      = ori_data;
  tree->n_pattern = ori_data->crunch_len;
  
  if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
  else Init_P_Lk_Tips_Int(tree);

  PhyML_Printf("\n. Estimate model parameters under switching substitution model.\n");
  Switch_From_Mod_To_M4mod(tree->mod);
  MC_Round_Optimize(tree);
/*   Round_Optimize(tree,tree->data); */
  /*   Simu_Loop(tree); */
  Print_Lk(tree,"[LnL under switching substitution model]");
  Print_Tree(tree->io->fp_out_tree,tree);

  n_iter = 0;
  do
    {
      Evolve(cpy_data,tree->mod,tree);
      tree->data      = cpy_data;
      tree->n_pattern = cpy_data->init_len;
      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(tree);
      else Init_P_Lk_Tips_Int(tree);

      Lk(tree);

      Init_Ui_Tips(tree);
      Site_Diversity(tree);
      Print_Diversity(fp_cov,tree);

      Print_Lk(tree,"[LnL under switching substitution model]");
    }while(++n_iter < n_simul);

  fclose(fp_obs);
  fclose(fp_nocov);
  fclose(fp_cov);
}

/*********************************************************/

m4 *M4_Copy_M4_Model(model *ori_mod, m4 *ori_m4mod)
{
  int i,j,n_h, n_o;
  m4 *cpy_m4mod;

  cpy_m4mod = (m4 *)M4_Make_Light((ori_mod->datatype == NT)?(4):(20));
  cpy_m4mod->n_h = ori_m4mod->n_h;

  M4_Make_Complete(cpy_m4mod->n_h,
		   cpy_m4mod->n_o,
		   cpy_m4mod);

  n_h = cpy_m4mod->n_h;
  n_o = cpy_m4mod->n_o;
  
  cpy_m4mod->n_h = ori_m4mod->n_h;
  cpy_m4mod->n_o = ori_m4mod->n_o;
  For(i,n_h) For(j,n_o*n_o) cpy_m4mod->o_mats[i][j] = ori_m4mod->o_mats[i][j];
  For(i,n_h) cpy_m4mod->multipl[i] = ori_m4mod->multipl[i];
  For(i,n_h) cpy_m4mod->multipl_unscaled[i] = ori_m4mod->multipl_unscaled[i];  
  For(i,n_o*n_o) cpy_m4mod->o_rr[i] = ori_m4mod->o_rr[i];
  For(i,n_h*n_h) cpy_m4mod->h_rr[i] = ori_m4mod->h_rr[i];
  For(i,n_h*n_h) cpy_m4mod->h_mat[i] = ori_m4mod->h_mat[i];
  For(i,n_o) cpy_m4mod->o_fq[i] = ori_m4mod->o_fq[i];
  For(i,n_h) cpy_m4mod->h_fq[i] = ori_m4mod->h_fq[i];
  For(i,n_h) cpy_m4mod->h_fq_unscaled[i] = ori_m4mod->h_fq_unscaled[i];
  cpy_m4mod->delta = ori_m4mod->delta;
  cpy_m4mod->alpha = ori_m4mod->alpha;

  return cpy_m4mod;
}

/*********************************************************/
#endif
