#include "modeltest.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_log.h>
#include <float.h>

void AIC(arbre* tree){
	int models[] = {LG, JTT, WAG, DAYHOFF, BLOSUM62, MTREV, RTREV, CPREV, DCMUT, VT, MTMAM, MTART, HIVW, HIVB};
	double bestScore = DBL_MAX;
	int bestModel = -1;
	int i;
	for(i = 0; i < 14; ++i) {
		double logLikelihood = likelihood(tree, models[i]);
		int params = getParams(models[i], tree);
		double aic = 2.0*params - 2.0*gsl_sf_log(logLikelihood);
		if(aic < bestScore) {
			bestScore = aic;
			bestModel = models[i];
		}
	}
	assignModel(tree, bestModel);
}

void HLRT(arbre* tree)
{
  Node * root = constructTree();
  //wikipedia.org/wiki/Likelihood-ratio_test
  int mod = runTests(tree,root,likelihood(tree,F81),F81);
  assignModel(tree,mod);
  destructTree(root);

  /*
  char outfilename[1000];
  sprintf(outfilename, "%s.modeltest",	tree->io->out_tree_file);
  FILE* outfile = fopen(outfilename,"w");
  printf("Probability :%f\nTherefore tree->whichmodel = %i\n", likelihood(tree,tree->mod->whichmodel), tree->mod->whichmodel);
  fclose(outfile);*/
}

void destructTree(Node* n){
  if(n->left)
    destructTree(n->left);
  
  if(n->right)
    destructTree(n->right);

  free(n);
}

int runTests(arbre* tree,Node* n, double previousLikelihood,int previousMod){
  float thisLikelihood = likelihood(tree,n->mod); 
  //float j = likelihood(tree,n->j); 

  double D = (-2)*(thisLikelihood/previousLikelihood);
  double x = 0.05; //significant P-value
  double df = getParams(previousMod,tree) - getParams(n->mod,tree); //degress of freedom. TODO: make a method to derive this
  if(df<0){
    df *= -1;
  } // absolute value
  double c = gsl_cdf_chisq_Qinv(x,df);

  if(c>D){
    return (!n->left) ? n->mod : runTests(tree,n->left,thisLikelihood,n->mod);
  } else {
    return (!n->right) ? previousMod : runTests(tree,n->right,previousLikelihood,previousMod);
  }
}

Node * constructTree(){
  Node * root = (Node *)malloc(sizeof(Node));
  root->mod = JC69;
  root->left = (Node *)malloc(sizeof(Node));
  root->right = (Node *)malloc(sizeof(Node));

  root->left->mod = K80;
  root->left->left = NULL;
  root->left->right = NULL;

  root->right->mod = HKY85;
  root->right->left = NULL;
  root->right->right = (Node *)malloc(sizeof(Node));

  root->right->right->mod = GTR;
  root->right->right->left = NULL;
  root->right->right->right = NULL;
  
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
      tree->mod->n_l = 1;
      tree->mod->s_opt->opt_state_freq = 0;
      break;
    case HKY85:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 1;
      tree->mod->s_opt->opt_lambda = 1;
      tree->mod->s_opt->opt_alpha = 0;
      tree->mod->invar = 0;
      tree->mod->n_l = 1;
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

int getParams(int mod, arbre* tree){
  switch(mod){
	case JC69:
	  return 1 + 2*tree->n_otu-3;
	  break;
	case K80:
	  return 2 + 2*tree->n_otu-3;
	  break;
	case F81:
	  return 5 + 2*tree->n_otu-3;
	  break;
	case F84:
	  return 6 + 2*tree->n_otu-3;
	  break;
	case HKY:
	  return 6 + 2*tree->n_otu-3;
	  break;
	case TN93:
	  return 6 + 2*tree->n_otu-3;
	  break;
	case GTR:
	  return 10 + 2*tree->n_otu-3;
	  break;
    default:
      return 2*tree->n_otu-3;
	  break;
  }
}
