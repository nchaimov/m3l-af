/*
 *
 * Hello World

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
#include "unittests.h"
#include "eb.h"
#include "modeltest.h"

#ifdef MPI
#include "mpi_boot.h"
#endif


int main(int argc, char **argv)
{
#ifdef RUN_TESTS
	Run_tests(argc, argv);
	exit(1);
#endif

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
	  data = Get_Seq(io);
	  Make_Model_Complete(io->mod);
	  mod = io->mod;
	  Print_Settings(io);

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

				  Init_Model(alldata,mod);

				  switch(io->in_tree)
				  {
					case 0 : case 1 : {
						tree = Dist_And_BioNJ(alldata,mod,io);
						break; }
					case 2 : { 	tree = Read_User_Tree(alldata,mod,io); break; }
				  }

				  if(!tree) continue;

				  time(&t_beg);
				  time(&(tree->t_beg));

				  tree->mod         = mod;
				  tree->io          = io;
				  tree->data        = alldata;
				  tree->both_sides  = 1;
				  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

          if(tree->io->mod->datatype == AA){
            AIC(tree);
          } else {
            HLRT(tree);
          }

				  Prepare_Tree_For_Lk(tree);

				  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

				  if(io->in_tree == 2) Spr_Pars(tree); /* minimize parsimony, if it was specified */

				  //
				  // VHS: If subalignment compression is enabled, then the following precondition is true
				  // as we enter this switch block: (1) the tree contains allocated space for
				  // red. arrays, (2) the red arrays are initialized, (3) the patterns have been collapsed according
				  // to the phylogeny.
				  //
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
						  Thermal_Anneal_All_Free_Params(tree, (io->quiet)?(1):(0));
						  break;
					  case SIMULATED_QUANTUM_ANNEALING:
						  Quantum_Anneal_All_Free_Params(tree, (io->quiet)?(1):(0));
						  break;
					  case EMPIRICAL_BAYES:
						  Empirical_Bayes(tree);
						  break;
					  default:
						  PhyML_Printf("\n The topology search option was not recognized...");
						  PhyML_Printf("\n. Err. in file %s on line %d\n",__FILE__,__LINE__);
						  Warn_And_Exit("\n");
					  }
				  }
				  else
				  {
					  if(tree->mod->s_opt->opt_num_param)
					  {
						  Round_Optimize(tree,tree->data,ROUND_MAX);
					  }
					  //else
					  //{
					//	  Lk(tree);
					  //}
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

				  if(tree->c_lnL > best_lnL)
				  {
					  best_lnL = tree->c_lnL;
					  Br_Len_Involving_Invar(tree);
					  most_likely_tree = Write_Tree(tree);
					  most_likely_size = Get_Tree_Size(tree);
				  }

				  time(&t_end);
				  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
						   io,num_data_set+1,
						   (tree->mod->s_opt->n_rand_starts > 1)?
						   (num_rand_tree):(num_tree));

				  if(tree->io->print_site_lnl)
				  {
					  Print_Site_Lk(tree,io->fp_out_lk);
				  }

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
			  else if(io->post_probs > 0) /* Launch posterior probability measurement, using MCMC sample */
			  {
				  if(!io->quiet) PhyML_Printf("\n. Compute posterior probabilities of clades on the most likely tree...\n");
				  most_likely_tree = PostProb_From_String(most_likely_tree,alldata,mod,io);
			  }


			  /* Print the most likely tree in the output file */
			  if(!io->quiet) PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n", Basename(io->out_tree_file));
			  if(io->n_data_sets == 1) rewind(io->fp_out_tree);
			  PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);

			  Free(most_likely_tree);

			  if(io->n_trees > 1 && io->n_data_sets > 1) break;
	  } //end if(data)
	  Free_Cseq(alldata);
	 } //end For
  }


  if(io->mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n\n. Best log likelihood : %f\n",best_lnL);

  /* OK, we're done.  Now cleanup. . . */

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

