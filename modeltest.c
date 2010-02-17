include "modeltest.h";
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
AIC(arbre* tree)
{
				model JayCeeSixNine;
				model EffEightyOne;
				//__Model DAYHOFF;
				//__Model JTT;
				//__Model HIVW;
				//__Model HIVB;
				
				JayCeeSixNine->datatype = NT;
				JayCeeSixNine->n_catg = 1;
				JayCeeSixNine->s_opt->opt_kappa = 0;
				JayCeeSixNine->s_opt->opt_lambda = 0;
				//JayCeeSixNine->alpha;
				JayCeeSixNine->s_opt->opt_alpha = 0;
				JayCeeSixNine->invar = 0;
				//JayCeeSixNine->pinvar;
				//JayCeeSixNine->s_opt->opt_pinvar;
				JayCeeSixNine->n_l = 1;
				JayCeeSixNine->s_opt->opt_state_freq = 0;
				JayCeeSixNine->modelname = "JC69";
				//
				// --MODEL F81 --
				//
				EffEightyOne->datatype = NT;
				EffEightyOne->n_catg = 1;
				EffEightyOne->s_opt->opt_kappa = 0;
				EffEightyOne->s_opt->opt_lambda = 0;
				//EffEightyOne->alpha;
				EffEightyOne->s_opt->opt_alpha = 0;
				EffEightyOne->invar = 0;
				//EffEightyOne->pinvar;
				//EffEightyOne->s_opt->opt_pinvar;
				EffEightyOne->n_l = 1;
				EffEightyOne->s_opt->opt_state_freq = 0;
				EffEightyOne->modelname = "F81";

				float Pjc69 = likelihood(tree, JayCeeSixNine*);
				float Pf81 = likelihood(tree, EffEightyOne*);
				
				if (Pjc69>Pf81) {
								Init_Model(tree->data, JayCeeSixNine);
								tree->mod = JayCeeSixNine;
				}
				char outfilename[1000];
				sprintf(outfilename, "%s.modeltest",	tree->io->out_tree_file);
				FILE* outfile = fopen(outfilename,"w");
				fprintf(outfile, "Probability of JC69:  %f. /n Probability of F81:  %f /n Therefore tree->mod->modelname = %s", Pjc69, Pf81, tree->mod->modelname);
				fclose(outfile);
}

//likelihood method
likelihood(arbre* tree, model* mod)
{
				// 2. Call this method to initialize several important values in the 'mod' struct.

				tree->mod->datatype = mod->datatype;
				tree->mod->ncat_g = mod->ncat_g;
				tree->mod->s_opt->opt_kappa = mod->s_opt->opt_kappa;
				tree->mod->s_opt->opt_lambda = mod->s_opt->opt_lambda;
				tree->mod->s_opt->opt_alpha = mod->s_opt->opt_alpha;
				tree->mod->invar = mod->invar;
				tree->mod->pinvar = mod->pinvar;
			 tree->mod->s_opt->opt_pinvar = mod->s_opt->opt_pinvar;
				tree->mod->n_l = mod->n_l;
				tree->mod->s_opt->opt_state_freq = mod->s_opt->opt_state_freq;
				tree->mod->modelname = mod->modelname;
				
				
				Init_Model(tree->data, mod);
				
				//
				// 3. Prepare for the likelihood calculation 
				//
				
				Prepare_Tree_For_Lk(tree);
				if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);
				
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




