#include "modeltest.h"

struct Node{
  int i;
  int j;
  struct Node * left;
  struct Node * right;
};

void AIC(arbre* tree)
{
  struct Node * root = constructTree();
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
  destructTree(root);
}

void destructTree(struct Node* n){
  if(n->left != NULL)
    destructTree(n->left);
  
  if(n->right != NULL)
    destructTree(n->right);

  free(n);
}

struct Node * constructTree(){
  struct Node * root = (struct Node *)malloc(sizeof(struct Node));
  root->i = JC69;
  root->j = F81;
  root->left = (struct Node *)malloc(sizeof(struct Node));
  root->right = (struct Node *)malloc(sizeof(struct Node));

  root->left->i = JC69;
  root->left->j = K80;

  root->right->i = F81;
  root->right->j = HKY85;
  
  return root;
}

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
      break;
   case K80:
   	tree->mod->datatype = NT;
	tree->mod->n_catg = 1;
	tree->mod->s_opt->opt_kappa = 0;
	tree->mod->s_opt->opt_lambda = 0;
	tree->mod->s_opt->opt_alpha = 0;
	tree->mod->invar = 0;
	tree->mod->n_1 = 1;
	tree->mod->s_opt->opt_state_freq = 0;
	break;
  case HKY85:
	tree->mod->datatype = NT;
	tree->mod->n_catg = 1;
	tree->mod->s_opt->opt_kappa = 1;
	tree->mod->s_opt->opt_lambda = 1;
	tree->mod->s_opt->opt_alpha = 0;
	tree->mod->invar = 0;
	tree->mod->n_1 = 1;
	tree->mod->s_opt->opt_state_freq = 0;
	break;
    default:
      return;
  }
  tree->mod->whichmodel = model;
  Init_Model(tree->data, tree->mod);
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


ubreakk
											Case SIMULATED_THERMAL_ANNEALING:
qe 
g;lkqerl'kjg 'lsfd jkg:
`
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
