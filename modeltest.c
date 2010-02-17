#include "modeltest.h"

//
// You'll create a for-loop or a switch structure, to iterate
// through different models.
//

//
// Here is how you find the likelihood, given one model
//

//
// 1. Set model-specific values in the 'mod' struct.
// See the specifications in the spreadsheet I emailed to you.
//
// i.e., make changes to tree->mod
//
void assignModel(arbre* tree,int model){
  switch (model){
    case JC69:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 0;
      tree->mod->s_opt->opt_lambda = 0;
      tree->mod->s_opt->opt_alpha = 0;
      tree->mod->invar = 0;
      tree->mod->n_l = 1;
      tree->mod->s_opt->opt_state_freq = 0;
      //tree->mod->modelname = "JC69";
      break;
    case F81:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 0;
      tree->mod->s_opt->opt_lambda = 0;
      tree->mod->s_opt->opt_alpha = 0;
      tree->mod->invar = 0;
      tree->mod->n_l = 1;
      tree->mod->s_opt->opt_state_freq = 0;
      //tree->mod->modelname = "F81";
      break;
    default:
      return;
  }
  tree->mod->whichmodel = model;
  Init_Model(tree->data, tree->mod);
}
void AIC(arbre* tree)
{
  printf("\n\n\nChecking tree. WOOOOOOOOOOOOOOO\n\n\n");
				float Pjc69 = likelihood(tree, JC69);
				float Pf81 = likelihood(tree, F81);
				
				if (Pjc69>Pf81) {
          assignModel(tree,JC69);
				}
				char outfilename[1000];
				sprintf(outfilename, "%s.modeltest",	tree->io->out_tree_file);
				FILE* outfile = fopen(outfilename,"w");
				//fprintf(outfile, "Probability of JC69:  %f. /n Probability of F81:  %f /n Therefore tree->mod->modelname = %s", Pjc69, Pf81, tree->mod->modelname);
				printf("Probability of JC69:  %f. /n Probability of F81:  %f /n Therefore tree->mod->modelname = %s", Pjc69, Pf81, tree->mod->modelname);
				fclose(outfile);
}

//likelihood method
float likelihood(arbre* tree, int mod)
{
				// 2. Call this method to initialize several important values in the 'mod' struct.

        assignModel(tree,mod);
				
				//
				// 3. Prepare for the likelihood calculation 
				//
				
				Prepare_Tree_For_Lk(tree);
				//if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);
				
				//
				// 4. Optimize all the free parameters, using maximum likelihood
				//
				if(tree->mod->s_opt->opt_topo)
				{
								switch(tree->mod->s_opt->topo_search)
								{
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
																Thermal_Anneal_All_Free_Params(tree, 1); //(io->quiet)?(1):(0)
																break;
												case SIMULATED_QUANTUM_ANNEALING:
																Quantum_Anneal_All_Free_Params(tree, 1); //(io->quiet)?(1):(0)
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
				}
				
				//
				// 5. Calculate the likelihood of the tree, given the model
				//
				tree->both_sides = 1;
				Lk(tree);
				Pars(tree);
				Get_Tree_Size(tree);
				//PhyML_Printf("\n. Log likelihood of the current tree: %f.\n",tree->c_lnL);
				
				//
				// 6. Record the likelihood of tree, given this model
				// (this step is for you to create)
				
				return tree->c_lnL;
}
