/*
 *
 * Hello World!!!-John St John

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

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
#include "annealing.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef PHYML

int main(int argc, char **argv)
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
  m3ldbl best_lnL,most_likely_size,tree_size;
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
	  if(!io->quiet) PhyML_Printf("\n. Compressing sequences...\n");
	  alldata = Compact_Seq(data,io);

	  Free_Seq(data,alldata->n_otu);
	  Check_Ambiguities(alldata,io->mod->datatype,io->mod->stepsize);

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    if(!io->quiet) PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);

		  /**
		   * JSJ: we need to initialize the model with some branch lengths
		   */
//		  int n_l = 2;
		  int i;
//		  m3ldbl *props = (m3ldbl *)mCalloc(n_l,sizeof(m3ldbl));
//		  For(i,n_l){
//			  props[i] = 1.0/n_l;
//		  }
		  Init_Model(alldata,mod);

		  switch(io->in_tree)
		    {
		    case 0 : case 1 : { tree = Dist_And_BioNJ(alldata,mod,io);    break; }
		    case 2 :
		    {
		    	tree = Read_User_Tree(alldata,mod,io);
		    	//JSJ: Make sure that the user hasn't defined starting
		    	//parameters that we need to copy into the tree.

		    	if(io->user_props != 0){//this value has been set! lets swap the trees value's with this one
		    		//PhyML_Printf("JSJ: Size of io props array is: %i\n[",((int)(sizeof(io->props)/sizeof(m3ldbl))));
		    		int tmp = io->n_l;
		    		if(tmp == 1) io->n_l = tree->n_l; //if io->n_l hasn't been initialized, set it...
		    		else if(tmp != tree->n_l){
		    			//it has been initialized, but not to the correct value!
		    			PhyML_Printf("Warning, the starting tree you supplied doesn't have the same number of branch length sets that you requested.\n");
		    			PhyML_Printf("Requested %i branch length sets. Found %i branch length sets on the starting tree.\n",tmp,tree->n_l);
		    			PhyML_Printf("\n. Err. in file %s on line %d\n",__FILE__,__LINE__);
		    			Warn_And_Exit("\n");
		    		}
		    		Normalize_Props_IO(io);
		    		For(i,tree->n_l){
//		    			if(i != (tree->n_l - 1)) PhyML_Printf(" %lf,",io->props[i]);
//		    			else PhyML_Printf(" %lf ]",io->props[i]);
		    			tree->props[i] = io->props[i];
		    		}

		    	}

		    	break;
		    }
		    }

		  /**
		   * JSJ: now lets fill the model with the values stored in the tree
		   * for n_l and props
		   */

		  if(!tree) continue;
		  /**
		   * JSJ: print out the current settings...
		   * It no longer makes sense to print this out right away, since some of these settings
		   * are dependent on the user's input tree. The output is thus a little sloppier because
		   * the user is identified that the input tree is read, the sequence is read, all before
		   * the analysis output is printed to the screen.
		   */
		  Print_Settings(io);

		  time(&t_beg);
		  time(&(tree->t_beg));

		  tree->mod         = mod;
//		  tree->mod->n_l    = tree->n_l;
//		  For(i,tree->n_l) tree->mod->props[i] = tree->props[i];
		  tree->io          = io;
		  tree->data        = alldata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

		  Prepare_Tree_For_Lk(tree);

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

		  if(io->in_tree == 1) Spr_Pars(tree);

		  if(tree->mod->s_opt->opt_topo)
		    {
			  switch(tree->mod->s_opt->topo_search){
			  case NNI_MOVE:
				  Simu_Loop(tree);
				  break;
			  case SPR_MOVE:
				  Speed_Spr_Loop(tree);
				  break;
			  case BEST_OF_NNI_AND_SPR:
				  Best_Of_NNI_And_SPR(tree);
				  break;
			  case SIMULATED_THERMAL_ANNEALING:
				  Thermal_anneal_all_free_params(tree, 1);//JSJ: 1 for verbose...
				  break;
			  case SIMULATED_QUANTUM_ANNEALING:
				  Quantum_anneal_all_free_params(tree, 1); //JSJ: 1 for verbose...
				  break;
			  default:
				  PhyML_Printf("\n The topology search option was not recognized...");
				  PhyML_Printf("\n. Err. in file %s on line %d\n",__FILE__,__LINE__);
				  Warn_And_Exit("\n");
			  }

//		      if(tree->mod->s_opt->topo_search      == NNI_MOVE) Simu_Loop(tree);
//		      else if(tree->mod->s_opt->topo_search == SPR_MOVE) Speed_Spr_Loop(tree);
//		      else if(tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) Best_Of_NNI_And_SPR(tree);
//		      else if
		    }
		  else
		    {
		      if(tree->mod->s_opt->opt_num_param ||
			 tree->mod->s_opt->opt_bl)                       Round_Optimize(tree,tree->data,ROUND_MAX);
		      else                                               Lk(tree);
		    }

		  tree->both_sides = 1;
		  Lk(tree);
		  Pars(tree);
		  Get_Tree_Size(tree);
		  PhyML_Printf("\n. Log likelihood of the current tree: %f.\n",tree->c_lnL);

		  /* Print the tree estimated using the current random (or BioNJ) starting tree */
		  if(io->mod->s_opt->n_rand_starts > 1)
		    {
		      Br_Len_Involving_Invar(tree);
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		    }

		  /* Record the most likely tree in a string of characters */
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
	      if(mod->bootstrap)
		{
		  if(!io->quiet) PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");

                  #ifdef MPI
		  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		  if(!io->quiet)  PhyML_Printf("\n. The bootstrap analysis will use %d CPUs.",Global_numTask);
		  #endif

		  most_likely_tree = Bootstrap_From_String(most_likely_tree,alldata,mod,io);
		}
	      else if(io->ratio_test)
		{
		  /* Launch aLRT */
		  if(!io->quiet) PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
		  most_likely_tree = aLRT_From_String(most_likely_tree,alldata,mod,io);
		}

	      /* Print the most likely tree in the output file */
	      if(!io->quiet) PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n", Basename(io->out_tree_file));
	      if(io->n_data_sets == 1) rewind(io->fp_out_tree);
	      PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);

	      Free(most_likely_tree);

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

  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

#elif(MG)
#include "mg.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#elif(M4)
#include "m4.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#elif(MC)
#include "mc.h"
int main(int argc, char **argv)
{
  MC_main(argc, argv);
  return 1;
}

#elif(RF)
int main(int argc, char **argv)
{
  arbre *tree1, *tree2;
  FILE *fp_tree1, *fp_tree2;
  int i,j;

  fp_tree1 = (FILE *)fopen(argv[1],"r");
  fp_tree2 = (FILE *)fopen(argv[2],"r");

  tree1 = Read_Tree_File(fp_tree1);
  tree2 = Read_Tree_File(fp_tree2);


  For(i,tree1->n_otu)
    {
      For(j,tree2->n_otu)
	{
	  if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
	    {
	      break;
	    }
	}

      if(j == tree2->n_otu)
	{
	  Prune_Subtree(tree1->noeud[i]->v[0],tree1->noeud[i],NULL,NULL,tree1);
	}
    }

  For(i,tree2->n_otu)
    {
      For(j,tree1->n_otu)
	{
	  if(!strcmp(tree2->noeud[i]->name,tree1->noeud[j]->name))
	    {
	      break;
	    }
	}

      if(j == tree1->n_otu)
	{
	  Prune_Subtree(tree2->noeud[i]->v[0],tree2->noeud[i],NULL,NULL,tree2);
	}
    }

  printf("%s\n",Write_Tree(tree1));
  printf("%s\n",Write_Tree(tree2));



/*   arbre *tree1, *tree2; */
/*   FILE *fp_tree1, *fp_tree2; */
/*   int i,j,rf,n_edges,n_common,bip_size; */
/*   m3ldbl thresh; */
/*   edge *b; */


/*   fp_tree1 = (FILE *)fopen(argv[1],"r"); */
/*   fp_tree2 = (FILE *)fopen(argv[2],"r"); */
/*   thresh = (m3ldbl)atof(argv[3]); */

/*   tree1 = Read_Tree_File(fp_tree1); */
/*   tree2 = Read_Tree_File(fp_tree2); */

/*   Get_Rid_Of_Prefix('_',tree1); */

/* /\*   Find_Common_Tips(tree1,tree2); *\/ */

/*   Alloc_Bip(tree1); */
/*   Alloc_Bip(tree2); */

/*   Get_Bip(tree1->noeud[0],tree1->noeud[0]->v[0],tree1); */
/*   Get_Bip(tree2->noeud[0],tree2->noeud[0]->v[0],tree2); */

/* /\*   PhyML_Printf("\n. rf=%f\n",Compare_Bip_On_Existing_Edges(thresh,tree1,tree2)); *\/ */
/*   For(i,2*tree1->n_otu-3) tree1->t_edges[i]->bip_score = 0; */
/*   For(i,2*tree2->n_otu-3) tree2->t_edges[i]->bip_score = 0; */

/*   rf = 0; */
/*   n_edges = 0; */

/*   /\* First tree *\/ */
/*   For(i,2*tree1->n_otu-3)  */
/*     { */
/*       /\* Consider the branch only if the corresponding bipartition has size > 1 *\/ */
/*       b = tree1->t_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

/*       if(bip_size > 1) */
/* 	{ */
/* 	  /\* with non-zero length *\/ */
/* 	  if(tree1->t_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This edge is not found in tree2 *\/ */
/* 	      if(!tree1->t_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */


/*   /\* Second tree *\/ */
/*   For(i,2*tree2->n_otu-3)  */
/*     { */
/*       b = tree2->t_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

/*       if(bip_size > 1) */
/* 	{ */
/* 	  if(tree2->t_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This edge is not found in tree1 *\/ */
/* 	      if(!tree2->t_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */

/*   if(!n_edges) */
/*     { */
/*       Exit("\n. No comparable internal edges were found.\n"); */
/*     } */
/*   else */
/*     { */
/*       PhyML_Printf("\n. Robinson and Foulds distance: %f.",(double)rf/(n_edges)); */
/* /\*       PhyML_Printf("\n. %d internal edges were processed (%d in the first tree, %d in the second).\n",n_edges,n_edges_t1,n_edges-n_edges_t1); *\/ */
/*       PhyML_Printf("\n"); */
/*     } */

  return 1;
}

#endif
