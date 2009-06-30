/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "simu.h"
#include "spr.h"
#include "pars.h"
#include "alrt.h"


/*********************************************************/

void Simu_Loop(arbre *tree)
{
  phydbl lk_old;

  tree->both_sides = 0;
  Lk(tree);

  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n. Maximizing likelihood (using NNI moves)...\n");
  
  do
    {
      lk_old = tree->c_lnL;
      Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
      if(!Simu(tree,10)) Check_NNI_Five_Branches(tree);
    }
  while(tree->c_lnL > lk_old + tree->mod->s_opt->min_diff_lk_global);


  do
    {
      if(!Check_NNI_Five_Branches(tree)) break;
    }while(1);
  /*****************************/
  
  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n");

}

/*********************************************************/

int Simu(arbre *tree, int n_step_max)
{
  phydbl old_loglk,n_iter,lambda;
  int i,n_neg,n_tested,n_without_swap,n_tot_swap,step,it_lim_without_swap;
  edge **sorted_b,**tested_b;
  int opt_free_param;
  int recurr;

  sorted_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  tested_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  
  old_loglk           = UNLIKELY;
  tree->c_lnL         = UNLIKELY;
  n_iter              = 1.0;
  it_lim_without_swap = (tree->mod->invar)?(8):(5);
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  lambda              = .75;
  n_tot_swap          = 0;
  opt_free_param      = 0;
  recurr              = 0;
  
  Update_Dirs(tree);
      
  if(tree->lock_topo)
    {
      PhyML_Printf("\n. The tree topology is locked.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  
  do
    {
      ++step;
                 
      old_loglk = tree->c_lnL;	    
      tree->both_sides = 1;
      Lk(tree);
      
      if(tree->c_lnL < old_loglk)
	{
	  if((tree->mod->s_opt->print) && (!tree->io->quiet)) printf("\n\n. Moving backward\n");
	  if(!Mov_Backward_Topo_Bl(tree,old_loglk,tested_b,n_tested))
	    Exit("\n. Err: mov_back failed\n");
	  if(!tree->n_swap) n_neg = 0;
	  
	  For(i,2*tree->n_otu-3) tree->t_edges[i]->l_old = tree->t_edges[i]->l;	    
	  tree->both_sides = 1;
	  Lk(tree);
	}

      if(step > n_step_max) break;

      if(tree->io->print_trace)
	{
	  PhyML_Fprintf(tree->io->fp_out_trace,"[%f]%s\n",tree->c_lnL,Write_Tree(tree)); fflush(tree->io->fp_out_trace);
	  if(tree->io->print_site_lnl) Print_Site_Lk(tree,tree->io->fp_out_lk); fflush(tree->io->fp_out_lk);
	}

      if((tree->mod->s_opt->print) && (!tree->io->quiet)) Print_Lk(tree,"[Topology           ]");
      
      if(((tree->c_lnL > old_loglk) &&
	  (fabs(old_loglk-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global)) ||
	 (n_without_swap > it_lim_without_swap)) break;
      

      Fill_Dir_Table(tree);
      Fix_All(tree);
      n_neg = 0;
      For(i,2*tree->n_otu-3)
	if((!tree->t_edges[i]->left->tax) && 
	   (!tree->t_edges[i]->rght->tax)) 
	  NNI(tree,tree->t_edges[i],0);
      
      
      Select_Edges_To_Swap(tree,sorted_b,&n_neg);	    	  
      Sort_Edges_NNI_Score(tree,sorted_b,n_neg);	    
      Optimiz_Ext_Br(tree);	  	    
      Update_Bl(tree,lambda);
	  	        
      n_tested = 0;
      For(i,(int)ceil((phydbl)n_neg*(lambda)))
	tested_b[n_tested++] = sorted_b[i];
      
      Make_N_Swap(tree,tested_b,0,n_tested);
            
      n_tot_swap += n_tested;
      
      if(n_tested > 0) n_without_swap = 0;
      else             n_without_swap++;
      
      n_iter+=1.0;
    }
  while(1);
    
/*   Round_Optimize(tree,tree->data); */

  Free(sorted_b);
  Free(tested_b);

  return n_tested;
}

/*********************************************************/

void Simu_Pars(arbre *tree, int n_step_max)
{
  phydbl old_pars,n_iter,lambda;
  int i,n_neg,n_tested,n_without_swap,n_tot_swap,step;
  edge **sorted_b,**tested_b;
  int each;
  
  sorted_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  tested_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  
  old_pars            = 0;
  tree->c_pars        = 0;
  n_iter              = 1.0;
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  each                = 4;
  lambda              = .75;
  n_tot_swap          = 0;
  
  Update_Dirs(tree);
  
  if((tree->mod->s_opt->print) && (!tree->io->quiet)) PhyML_Printf("\n. Starting simultaneous NNI moves (parsimony criterion)...\n");
  
  do
    {
      ++step;
      
      if(step > n_step_max) break;
      
      each--;
      
      tree->both_sides = 1;
      Pars(tree);
      
      if((tree->mod->s_opt->print) && (!tree->io->quiet)) 
	{
	  Print_Pars(tree);
	  if(step > 1) (n_tested > 1)?(printf("[%4d NNIs]",n_tested)):(printf("[%4d NNI ]",n_tested));
	}
      

      if(old_pars == tree->c_pars) break;
	      
      if((tree->c_pars > old_pars) && (step > 1))
	{
	  if((tree->mod->s_opt->print) && (!tree->io->quiet))
	    PhyML_Printf("\n\n. Moving backward (topology) \n");
	  if(!Mov_Backward_Topo_Pars(tree,old_pars,tested_b,n_tested))
	    Exit("\n. Err: mov_back failed\n");
	  if(!tree->n_swap) n_neg = 0;
	  
	  tree->both_sides = 1;
	  Pars(tree);
	}
      else 
	{
	  
	  old_pars = tree->c_pars;	    
	  Fill_Dir_Table(tree);

	  n_neg = 0;
	  For(i,2*tree->n_otu-3)
	    if((!tree->t_edges[i]->left->tax) && 
	       (!tree->t_edges[i]->rght->tax)) 
	      NNI_Pars(tree,tree->t_edges[i],0);
	  
	  Select_Edges_To_Swap(tree,sorted_b,&n_neg);	    	  
	  Sort_Edges_NNI_Score(tree,sorted_b,n_neg);	    
	  
	  n_tested = 0;
	  For(i,(int)ceil((phydbl)n_neg*(lambda)))
	    tested_b[n_tested++] = sorted_b[i];
	  
	  Make_N_Swap(tree,tested_b,0,n_tested);
	  
	  n_tot_swap += n_tested;
	  
	  if(n_tested > 0) n_without_swap = 0;
	  else             n_without_swap++;
	}
      n_iter+=1.0;
    }
  while(1);
  
  Free(sorted_b);
  Free(tested_b);
}

/*********************************************************/

void Select_Edges_To_Swap(arbre *tree, edge **sorted_b, int *n_neg)
{
  int i;
  edge *b;
  int min;
  phydbl best_score;

  *n_neg = 0;
  min = 0;
  
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      best_score = b->nni->score;

      if((!b->left->tax) && (!b->rght->tax) && (b->nni->score < -tree->mod->s_opt->min_diff_lk_move)) 
	{
	  Check_NNI_Scores_Around(b->left,b->rght,b,&best_score);
	  Check_NNI_Scores_Around(b->rght,b->left,b,&best_score);
	  if(best_score < b->nni->score) continue;
	  sorted_b[*n_neg] = b;
	  (*n_neg)++;
	}
    }
}

/*********************************************************/

void Update_Bl(arbre *tree, phydbl fact)
{
  int i;
  edge *b;

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      b->l = b->l_old + (b->nni->l0 - b->l_old)*fact;
    }
}

/*********************************************************/

void Make_N_Swap(arbre *tree,edge **b, int beg, int end)
{
  int i;

/*   PhyML_Printf("\n. Beg Actually performing swaps\n"); */
  tree->n_swap = 0;
  for(i=beg;i<end;i++)
    {
      /* we use t_dir here to take into account previous modifications of the topology */
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]],
	   b[i]->nni->swap_node_v2,
	   b[i]->nni->swap_node_v3,
	   b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]],
	   tree);

      b[i]->l = b[i]->nni->best_l;

      tree->n_swap++;
    }

/*   PhyML_Printf("\n. End Actually performing swaps\n"); */

}

/*********************************************************/

int Make_Best_Swap(arbre *tree)
{
  int i,j,return_value;
  edge *b,**sorted_b;
  

  sorted_b = (edge **)mCalloc(tree->n_otu-3,sizeof(edge *));
  
  j=0;
  For(i,2*tree->n_otu-3) if((!tree->t_edges[i]->left->tax) &&
			    (!tree->t_edges[i]->rght->tax))
                              sorted_b[j++] = tree->t_edges[i];

  Sort_Edges_NNI_Score(tree,sorted_b,tree->n_otu-3);

  if(sorted_b[0]->nni->score < -0.0)
    {
      b = sorted_b[0];
      return_value = 1;

      Swap(b->nni->swap_node_v2->v[tree->t_dir[b->nni->swap_node_v2->num][b->nni->swap_node_v1->num]],
	   b->nni->swap_node_v2,
	   b->nni->swap_node_v3,
	   b->nni->swap_node_v3->v[tree->t_dir[b->nni->swap_node_v3->num][b->nni->swap_node_v4->num]],
	   tree);


      b->l = b->nni->best_l;

/*       (b->nni->best_conf == 1)? */
/* 	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v1],tree)): */
/* 	(Swap(b->left->v[b->l_v2],b->left,b->rght,b->rght->v[b->r_v2],tree)); */
      
/*       b->l =  */
/* 	(b->nni->best_conf == 1)? */
/* 	(b->nni->l1): */
/* 	(b->nni->l2); */


    }
  else return_value = 0;

  Free(sorted_b);

  return return_value;
}

/*********************************************************/

int Mov_Backward_Topo_Bl(arbre *tree, phydbl lk_old, edge **tested_b, int n_tested)
{
  phydbl *l_init;
  int i,step,beg,end;
  edge *b;


  l_init = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  For(i,2*tree->n_otu-3) l_init[i] = tree->t_edges[i]->l;
  
  step = 2;
  tree->both_sides = 0;
  do
    {
      For(i,2*tree->n_otu-3) 
	{
	  b = tree->t_edges[i];
	  b->l = b->l_old + (1./step) * (l_init[i] - b->l_old);
	}

      beg = (int)floor((phydbl)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)floor((phydbl)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);
      
      if(!end) tree->n_swap = 0;
      
      tree->both_sides = 0;
      Lk(tree);
      
      step++;

    }while((tree->c_lnL < lk_old) && (step < 1000));

  
  if(step == 1000)
    {
      if(tree->n_swap)  Exit("\n. Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");

      For(i,2*tree->n_otu-3) 
	{
	  b = tree->t_edges[i];
	  b->l = b->l_old;
	}
      
      tree->both_sides = 0;
      Lk(tree);
    }

  Free(l_init);

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3) 
    {
      if(tree->t_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->t_edges[i]->nni->score = +1.0;
    }

  
  if(tree->c_lnL > lk_old)                                return  1;
  else if((tree->c_lnL > lk_old-tree->mod->s_opt->min_diff_lk_local) && 
	  (tree->c_lnL < lk_old+tree->mod->s_opt->min_diff_lk_local)) return -1;
  else                                                    return  0;
}

/*********************************************************/

int Mov_Backward_Topo_Pars(arbre *tree, int pars_old, edge **tested_b, int n_tested)
{
  int i,step,beg,end;
  
  step = 2;
  tree->both_sides = 0;
  do
    {
      beg = (int)floor((phydbl)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)floor((phydbl)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);
      
      if(!end) tree->n_swap = 0;
      
      tree->both_sides         = 0;
      Pars(tree);
      
      step++;

    }while((tree->c_pars > pars_old) && (step < 1000));

  
  if(step == 1000)
    {
      if(tree->n_swap)  Exit("\n. Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");

      tree->both_sides = 0;
      Pars(tree);
    }

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3) 
    {
      if(tree->t_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->t_edges[i]->nni->score = +1.0;
    }

  
  if(tree->c_pars < pars_old)       return  1;
  else if(tree->c_pars == pars_old) return -1;
  else                              return  0;
}

/*********************************************************/

void Unswap_N_Branch(arbre *tree, edge **b, int beg, int end)
{
  int i;

  if(end>beg)
    {
      for(i=beg;i<end;i++)
	{
	  
/* 	  PhyML_Printf("MOV BACK UNSWAP Edge %d Swap nodes %d(%d) %d %d %d(%d)\n", */
/* 		 b[i]->num, */
/* 		 b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]]->num, */
/* 		 b[i]->nni->swap_node_v4->num, */
/* 		 b[i]->nni->swap_node_v2->num, */
/* 		 b[i]->nni->swap_node_v3->num, */
/* 		 b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]]->num, */
/* 		 b[i]->nni->swap_node_v1->num */
/* 		 ); */

	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]],
	       tree);


/* 	  (b[i]->nni->best_conf == 1)? */
/* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v1],tree)): */
/* 	    (Swap(b[i]->left->v[b[i]->l_v2],b[i]->left,b[i]->rght,b[i]->rght->v[b[i]->r_v2],tree)); */

	  b[i]->l = b[i]->l_old;
	}
    }
  else
    {
      for(i=beg-1;i>=end;i--)
	{	
	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]],
	       tree);

	  b[i]->l = b[i]->l_old;
	}
    }
}

/*********************************************************/

void Swap_N_Branch(arbre *tree,edge **b, int beg, int end)
{
  int i;
  
  if(end>beg)
    {
      for(i=beg;i<end;i++)
	{
	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]],
	       tree);


	  b[i]->l = b[i]->nni->best_l;

	}
    }
  else
    {
      for(i=beg-1;i>=end;i--)
	{
	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num][b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num][b[i]->nni->swap_node_v4->num]],
	       tree);

	  b[i]->l = b[i]->nni->best_l;
	}
    }
}

/*********************************************************/

void Check_NNI_Scores_Around(node *a, node *d, edge *b, phydbl *best_score)
{

  int i;
  For(i,3)
    {
      if((d->v[i] != a) && (!d->v[i]->tax))
	{
	  if((d->b[i]->nni->score > *best_score-1.E-10) &&
	     (d->b[i]->nni->score < *best_score+1.E-10)) /* ties */
	     {
	       d->b[i]->nni->score = *best_score+1.;
	     }

	  if(d->b[i]->nni->score < *best_score) 
	    {
	      *best_score = d->b[i]->nni->score;
	    }
	}
    }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
