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

#ifdef COMPRESS_SUBALIGNMENTS
#include "compress.h"
#endif


// s=(char *)mCalloc(T_MAX_LINE,sizeof(char));
// append = 0 = no append
// append = 1 = yes append
void Write_File(char *filename, char *s, int append)
{
	FILE *file;
	if (append == 0)
	{
		file = fopen(filename,"w+");
	}
	else{
		file = fopen(filename, "a+");
	}
	fprintf(file, "%s", s); /*writes*/
	fclose(file); /*done!*/
}

void Test_alignment_read(option *io, seq **data, int testid)
{
	PhyML_Printf("Test #%d: alignment read\n", testid);
	if (testid == 1)
	{
		if (data[0]->len != 10)
		{
			PhyML_Printf("Test %d failed, the alignment should be 10 sites long.\n", testid);
			exit(1);
		}
		else
		{
			PhyML_Printf("OK. Alignment length = %d sites.\n", data[0]->len);
		}
		int j;
		For(j,io->mod->n_otu)
		{
			if (!data[j])
			{
				PhyML_Printf("Test %d failed, I cannot find taxon #%d of %d.\n", testid, j, io->mod->n_otu);
				exit(1);
			}
		}
		PhyML_Printf("OK. Alignment contains %d taxa.\n", io->mod->n_otu);
	}
}

void Run_tests(int argc, char **argv)
{
	/*
	 * Test case #1:
	 */
	//char align[] = " 4 10\nta\nAAAAAAAAAA\ntb\nDDDDDDDDDD\nta\nEEEEEEEEEE\nta\nLLLLLLLLLL\n";
	//Write_File("align.phy", align, 0);
	//int x = Test_main(argc, argv, 1);
}

// This is a faux main for unit tests.
int Test_main(int argc, char **argv, int testid)
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
	  data = Get_Seq(io);

	  //
	  Test_alignment_read(io, data, testid);

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
			case 0 : case 1 : { tree = Dist_And_BioNJ(alldata,mod,io);    break; }
			case 2 : { 	tree = Read_User_Tree(alldata,mod,io); break; }
		  }

// VHS: this is code from JSJ, which I think we can scrap:

			//JSJ: Make sure that the user hasn't defined starting
				//parameters that we need to copy into the tree.
//				if(tree->mod->n_l > 1) io->mod->s_opt->opt_five_branch = 0;
//				if(io->user_topo == 0 && tree->mod->n_l > 1){ //if more than one bl set default to simulated thermal annealing.
//					mod->s_opt->topo_search = SIMULATED_THERMAL_ANNEALING;
//				}
//				if(io->user_props != 0){
//					int tmp = io->mod->n_l;
//					//if(tmp == 1) io->mod->n_l = tree->mod->n_l; //if io->n_l hasn't been initialized, set it...
//					if(tmp != tree->mod->n_l){
//						PhyML_Printf("Warning, the starting tree you supplied doesn't have the same number of branch length sets that you requested.\n");
//						PhyML_Printf("Requested %i branch length sets. Found %i branch length sets on the starting tree.\n",tmp,tree->mod->n_l);
//						PhyML_Printf("\n. Err. in file %s on line %d\n",__FILE__,__LINE__);
//						Warn_And_Exit("\n");
//					}

//					For(i,tree->mod->n_l){
//						tree->mod->bl_props[i] = io->mod->bl_props[i];
//					}
//				}

		  if(!tree) continue;

		  Print_Settings(io);

		  time(&t_beg);
		  time(&(tree->t_beg));

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = alldata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

		  Prepare_Tree_For_Lk(tree);

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);

		  if(io->in_tree == 1) Spr_Pars(tree);

		  //
		  // VHS: precondition entering all these cases: the tree contains allocated space for
		  // red. arrays, the red arrays are initialized, the patterns have been collapsed according
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
			  default:
				  PhyML_Printf("\n The topology search option was not recognized...");
				  PhyML_Printf("\n. Err. in file %s on line %d\n",__FILE__,__LINE__);
				  Warn_And_Exit("\n");
			  }

		  }
		  else
		  {
			  if(tree->mod->s_opt->opt_num_param || tree->mod->s_opt->opt_bl)
			  {
				  Round_Optimize(tree,tree->data,ROUND_MAX);
			  }
			  else
			  {
				  Lk(tree);
			  }
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
