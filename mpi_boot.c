/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "mpi_boot.h"
#include "utilities.h"
#include "bionj.h"
#include "lk.h"
#include "pars.h"
#include "free.h"
#include "models.h"
#include "simu.h"
#include "spr.h"

#ifdef MPI

/*********************************************************/

void Bootstrap_MPI(arbre *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position, init_len, nbRep;

  allseq *boot_data;
  arbre *boot_tree;
  model *boot_mod;
  matrix *boot_mat;
  char *s;

  MPI_Status Stat;
  int randomRecv, bootRecv, nbElem, i;
  int *score_par, *score_tot;
  char *bootStr, *t;
  
  randomRecv = nbElem = bootRecv = 0;

  tree->print_boot_val       = 1;
  tree->print_alrt_val       = 0;
  boot_tree                  = NULL;

  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));

  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

  n_site = 0;
  For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
    {
      site_num[n_site] = j;
      n_site++;
    }

  boot_data = Copy_Cseq(tree->data, tree->data->crunch_len, tree->mod->ns);

  if (Global_myRank == 0)
    printf("\n. Non parametric bootstrap analysis \n");
  
  //number of bootstraps for each process
  if (tree->mod->bootstrap%Global_numTask != 0) 
    {
      nbRep = (tree->mod->bootstrap / Global_numTask) + 1;
      tree->mod->bootstrap = nbRep * Global_numTask;
      if (Global_myRank == 0) {
        PhyML_Printf("\n. The number of replicates is not a multiple of %d CPUs.\n", Global_numTask);
        PhyML_Printf("\n. Will run %d replicates analysis.\n", tree->mod->bootstrap);
      }
    }
  else
    nbRep = tree->mod->bootstrap/Global_numTask;
  
  //Bip score
  if (Global_myRank == 0) 
    {
      score_tot = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
      For(i,2*tree->n_otu-3)
	score_tot[i] = 0;
    }
  else
    score_tot = NULL;

  score_par = (int *)mCalloc((2*tree->n_otu - 3),sizeof(int));
  For(i,2*tree->n_otu-3)
    score_par[i] = 0;

  if (Global_myRank == 0)
    PhyML_Printf("\n  [");

  For(replicate, nbRep)
    {
      For(j,boot_data->crunch_len) boot_data->wght[j] = 0;

      // Send random data to other process
      if (Global_myRank == 0) {
        // Compute number of data to send
        if (tree->mod->bootstrap - randomRecv > Global_numTask)
          nbElem = Global_numTask;
        else
          nbElem = tree->mod->bootstrap - randomRecv;

        For(i,nbElem) {
          For(j,boot_data->crunch_len) boot_data->wght[j] = 0;
          init_len = 0;
          // Create random data
          For(j,boot_data->init_len)
	    {
	      position = Rand_Int(0,(int)(tree->data->init_len-1.0));
	      boot_data->wght[site_num[position]] += 1;
	      init_len++;
	    }
	    
          if (init_len != tree->data->init_len) {
            MPI_Finalize();
            Warn_And_Exit("\n. Pb when copying sequences\n");
          }
          // Send random data to other process, not to current process
	  if (i < nbElem-1) {
	    MPI_Ssend (boot_data->wght, boot_data->crunch_len, MPI_INT, i+1, Global_myRank, MPI_COMM_WORLD);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, sending random to %d done\n", Global_myRank, i+1);
fflush(stderr);
#endif
	  }
	  randomRecv++;
        }
      }
      else {
        MPI_Recv (boot_data->wght, boot_data->crunch_len, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, receiving random from task %d done\n", Global_myRank, Stat.MPI_SOURCE);
fflush(stderr);
#endif
      }

      init_len = 0;
      For(j,boot_data->crunch_len) init_len += boot_data->wght[j];

      if(init_len != tree->data->init_len) {
        MPI_Finalize();
        Warn_And_Exit("\n. Pb when copying sequences\n");
      }

      (tree->mod->datatype == NT)?
	(Get_Base_Freqs(boot_data)):
	(Get_AA_Freqs(boot_data));

      if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

      boot_mod = Copy_Model(tree->mod);
      Init_Model(boot_data,boot_mod);

      if(tree->io->in_tree == 2)
	{
	  rewind(tree->io->fp_in_tree);
	  boot_tree = Read_Tree_File(tree->io->fp_in_tree);
	}
      else
	{
	  boot_mat = ML_Dist(boot_data,boot_mod);
	  boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
	  Fill_Missing_Dist(boot_mat);
	  Bionj(boot_mat);
	  boot_tree = boot_mat->tree;
	  boot_tree->mat = boot_mat;
	}

      boot_tree->mod                = boot_mod;
      boot_tree->io                 = tree->io;
      boot_tree->data               = boot_data;
      boot_tree->both_sides         = 1;
      boot_tree->mod->s_opt->print  = 0;
      boot_tree->n_pattern          = boot_tree->data->crunch_len/
	                              boot_tree->mod->stepsize;
      boot_tree->io->print_site_lnl = 0;
      boot_tree->io->print_trace    = 0;

      if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
      Order_Tree_CSeq(boot_tree,boot_data);
      Share_Lk_Struct(tree,boot_tree);
      Share_Spr_Struct(tree,boot_tree);
      Share_Pars_Struct(tree,boot_tree);
      Fill_Dir_Table(boot_tree);
      Update_Dirs(boot_tree);

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
      else                         Init_P_Lk_Tips_Int(boot_tree);
      Init_Ui_Tips(boot_tree);
      Init_P_Pars_Tips(boot_tree);
      Br_Len_Not_Involving_Invar(boot_tree);
      
      if(boot_tree->mod->s_opt->opt_topo)
	{
	  if(boot_tree->mod->s_opt->topo_search == NNI_MOVE) 
	    {
	      Simu_Loop(boot_tree);
	    }
	  else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
		  (boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
	    {
	      Speed_Spr_Loop(boot_tree);
	    }
	}
      else
	{
	  if(boot_tree->mod->s_opt->opt_num_param || boot_tree->mod->s_opt->opt_bl)
	    Round_Optimize(boot_tree,boot_tree->data,ROUND_MAX);
	  else
	    Lk(boot_tree);
	}

      Alloc_Bip(boot_tree);

      Get_Bip(boot_tree->noeud[0],
	      boot_tree->noeud[0]->v[0],
	      boot_tree);

      if(!tree->io->collapse_boot) 
	Compare_Bip(tree,boot_tree);
      else
	Compare_Bip_On_Existing_Edges(1,tree,boot_tree);

      Br_Len_Involving_Invar(boot_tree);

      if(tree->io->print_boot_trees)
	{
	  s = Write_Tree(boot_tree);
          t=(char *)mCalloc(T_MAX_LINE,sizeof(char));
          Print_Fp_Out_Lines_MPI(boot_tree, tree->io, replicate+1, t);
          
          // Get bootstrap trees from other process and write to boot file
          if (Global_myRank == 0) {
            fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
            fprintf(tree->io->fp_out_boot_stats,"%s\n",t);
            bootRecv++;
            PhyML_Printf(".");
            // Compute number of bootstraps to receive
            if (tree->mod->bootstrap - bootRecv > Global_numTask)
              nbElem = Global_numTask;
            else
              nbElem = tree->mod->bootstrap - bootRecv + 1;
              
            bootStr=(char *)mCalloc(T_MAX_LINE,sizeof(char));
            for (i=1; i<nbElem; i++) {
              MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, receiving bootstrap from task %d tag %d done\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
fflush(stderr);
#endif
              if (Stat.MPI_TAG == BootTreeTag)
                fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
              if (Stat.MPI_TAG == BootStatTag)
                fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);

              MPI_Recv (bootStr, T_MAX_LINE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Stat);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, receiving bootstrap from task %d tag %d done\n", Global_myRank, Stat.MPI_SOURCE, Stat.MPI_TAG);
fflush(stderr);
#endif
              if (Stat.MPI_TAG == BootTreeTag)
                fprintf(tree->io->fp_out_boot_tree,"%s\n", bootStr);
              if (Stat.MPI_TAG == BootStatTag)
                fprintf(tree->io->fp_out_boot_stats,"%s\n", bootStr);

              bootRecv++;
              PhyML_Printf(".");
              if(!((bootRecv)%20)) {
	        printf("] %4d/%4d\n  ",bootRecv,tree->mod->bootstrap);
	        if(bootRecv != tree->mod->bootstrap) printf("[");
	      }
            }
            Free(bootStr);
          }
          else {
             MPI_Ssend (s, T_MAX_LINE, MPI_CHAR, 0, BootTreeTag, MPI_COMM_WORLD);
             MPI_Ssend (t, T_MAX_LINE, MPI_CHAR, 0, BootStatTag, MPI_COMM_WORLD);
#ifdef MPI_DEBUG
fprintf (stderr, "\ntask %d, sending bootstraps done\n", Global_myRank);
fflush(stderr);
#endif
          }
          Free (t);
	  Free(s);
	}

#ifndef QUIET
fflush(stdout);
#endif

      if(boot_tree->mat) Free_Mat(boot_tree->mat);
      Free_Tree(boot_tree);      
      Free_Model(boot_mod);

      //Each process computes the Bip score sum for all its bootstrap trees
      For(i,2*tree->n_otu-3) score_par[i] = tree->t_edges[i]->bip_score;

      //Each process sends its Bip score sum. The sums are summed.
      MPI_Reduce(score_par, score_tot, 2*tree->n_otu - 3, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }

  
  if (Global_myRank == 0) 
    {
      For(i,2*tree->n_otu-3)
	tree->t_edges[i]->bip_score = score_tot[i];
      Free (score_tot);
    }
  Free (score_par);

  if (Global_myRank == 0)
    if(((bootRecv)%20)) printf("] %4d/%4d\n ",bootRecv,tree->mod->bootstrap);

  tree->lock_topo = 1; /* Topology should not be modified afterwards */

  if(tree->io->print_boot_trees)
    {
      fclose(tree->io->fp_out_boot_tree);
      fclose(tree->io->fp_out_boot_stats);
    }

  Free_Cseq(boot_data);
  Free(site_num);
}

/*********************************************************/

void Print_Fp_Out_Lines_MPI(arbre *tree, option *io, int n_data_set, char *bootStr)
{
  char *s, *tmp;

  // Build a string to be sent to the writing process
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  tmp=(char *)mCalloc(T_MAX_LINE,sizeof(char));
  
  if (Global_myRank == 0 && n_data_set == 1) {
    snprintf(tmp, T_MAX_LINE, ". Sequence file : [%s]\n\n", io->in_seq_file); strncat (s, tmp, T_MAX_LINE);
    
    (tree->mod->datatype == NT)?
	(snprintf(tmp, T_MAX_LINE, ". Model of nucleotides substitution : %s\n\n", io->mod->modelname)):
	(snprintf(tmp, T_MAX_LINE, ". Model of amino acids substitution : %s\n\n", io->mod->modelname));
    strncat (s, tmp, T_MAX_LINE);
    
    switch(io->in_tree)
      {
      case 0: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [BioNJ]\n\n");               break; }
      case 1: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [parsimony]\n\n");           break; }
      case 2: { snprintf(tmp, T_MAX_LINE, ". Initial tree : [%s]\n\n",io->in_tree_file); break; }
      }
    strncat (s, tmp, T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);
    
    /*headline 1*/
    strncat (s, ". Data\t", T_MAX_LINE);

    strncat (s, "Nb of \t", T_MAX_LINE);

    strncat (s, "Likelihood\t", T_MAX_LINE);

    strncat (s, "Discrete   \t", T_MAX_LINE);

    if(tree->mod->n_catg > 1)
      strncat (s, "Number of \tGamma shape\t", T_MAX_LINE);

    strncat (s, "Proportion of\t", T_MAX_LINE);

    if(tree->mod->whichmodel <= 6)
      strncat (s, "Transition/ \t", T_MAX_LINE);
    
    strncat (s, "Nucleotides frequencies               \t", T_MAX_LINE);

    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      strncat (s, "Instantaneous rate matrix              \t", T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);

    /*headline 2*/
    strncat (s, "  set\t", T_MAX_LINE);

    strncat (s, "taxa\t", T_MAX_LINE);

    strncat (s, "loglk     \t", T_MAX_LINE);

    strncat (s, "gamma model\t", T_MAX_LINE);

    if(tree->mod->n_catg > 1)
      strncat (s, "categories\tparameter  \t", T_MAX_LINE);

    strncat (s, "invariant    \t", T_MAX_LINE);

    if(tree->mod->whichmodel <= 6)
      strncat (s, "transversion\t", T_MAX_LINE);

    strncat (s, "f(A)      f(C)      f(G)      f(T)    \t", T_MAX_LINE);

    if((tree->mod->whichmodel == GTR) ||
       (tree->mod->whichmodel == CUSTOM))
      strncat (s, "[A---------C---------G---------T------]\t", T_MAX_LINE);

    strncat (s, "\n", T_MAX_LINE);

    /*headline 3*/
    if(tree->mod->whichmodel == TN93) {
      strncat (s, "    \t      \t          \t           \t", T_MAX_LINE);
      if(tree->mod->n_catg > 1)
        strncat (s, "         \t         \t", T_MAX_LINE);
      strncat (s, "             \t", T_MAX_LINE);
      strncat (s, "purines pyrimid.\t", T_MAX_LINE);
      strncat (s, "\n", T_MAX_LINE);
    }

    strncat (s, "\n", T_MAX_LINE);
  }

  /*line items*/

  snprintf(tmp, T_MAX_LINE, "  #%d\t", (((n_data_set-1)*Global_numTask)+Global_myRank+1));
  strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%d   \t",tree->n_otu); strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%.5f\t",tree->c_lnL); strncat (s, tmp, T_MAX_LINE);
  
  snprintf(tmp, T_MAX_LINE, "%s        \t",
	  (tree->mod->n_catg>1)?("Yes"):("No ")); strncat (s, tmp, T_MAX_LINE);
  
  if(tree->mod->n_catg > 1)
    {
      snprintf(tmp, T_MAX_LINE, "%d        \t",tree->mod->n_catg); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%.3f    \t",tree->mod->alpha); strncat (s, tmp, T_MAX_LINE);
    }
  
  snprintf(tmp, T_MAX_LINE, "%.3f    \t",tree->mod->pinvar); strncat (s, tmp, T_MAX_LINE);
  
  if(tree->mod->whichmodel <= 5)
    {
      snprintf(tmp, T_MAX_LINE, "%.3f     \t",tree->mod->kappa); strncat (s, tmp, T_MAX_LINE);
    }
  else if(tree->mod->whichmodel == TN93)
    {
      snprintf(tmp, T_MAX_LINE, "%.3f   ",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%.3f\t",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
      strncat (s, tmp, T_MAX_LINE);
    }

  if(tree->mod->datatype == NT)
    {
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->pi[0]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->pi[1]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->pi[2]); strncat (s, tmp, T_MAX_LINE);
      snprintf(tmp, T_MAX_LINE, "%8.5f\t",tree->mod->pi[3]); strncat (s, tmp, T_MAX_LINE);
    }
    
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      For(i,4)
	{
	  if (i!=0) {
	    /*format*/
	    snprintf(tmp, T_MAX_LINE, "      \t     \t          \t           \t");
	    strncat (s, tmp, T_MAX_LINE);
	    if(tree->mod->n_catg > 1) {
	      snprintf(tmp, T_MAX_LINE, "          \t           \t");
	      strncat (s, tmp, T_MAX_LINE);
	    }
	    snprintf(tmp, T_MAX_LINE, "             \t                                      \t");
	    strncat (s, tmp, T_MAX_LINE);
	  }
	  For(j,4) {
	    snprintf(tmp, T_MAX_LINE, "%8.5f  ",tree->mod->qmat[i*4+j]);
	    strncat (s, tmp, T_MAX_LINE);
	  }
	  if (i<3) {
	    snprintf(tmp, T_MAX_LINE, "\n");
	    strncat (s, tmp, T_MAX_LINE);
	  }
	}
    }
    
    Free (tmp);
    strncpy (bootStr, s, T_MAX_LINE);
    Free (s);

  return;
}
#endif

