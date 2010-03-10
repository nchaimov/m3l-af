#include "modeltest.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_log.h>
#include <float.h>

FILE* outfile;

void modeltest(arbre* tree){
  char outfilename[1000];
  sprintf(outfilename, "%s.modeltest", tree->io->out_tree_file);
  outfile = fopen(outfilename,"w");

  if(tree->io->mod->datatype == AA){
    fprintf(outfile,"Entering AIC:\n");
    AIC(tree);
  } else {
    fprintf(outfile,"Entering HLRT:\n");
    HLRT(tree);
  }

  double finalLikelihood = likelihood(tree,tree->mod->whichmodel);

  char modelname[10];
  getName(modelname,tree->mod->whichmodel);
  fprintf(outfile,"\n  Best model: %s with likelihood %f\n",modelname,finalLikelihood);

  fclose(outfile);
}

void AIC(arbre* tree){
	int models[] = {LG, JTT, WAG, DAYHOFF, BLOSUM62, MTREV, RTREV, CPREV, DCMUT, VT, MTMAM, MTART, HIVW, HIVB};
	double bestScore = DBL_MAX;
	int bestModel = -1;
	int i;
	for(i = 0; i < 14; ++i) {
		double logLikelihood = likelihood(tree, models[i]);
    char modelname[10];
    getName(modelname,models[i]);
    fprintf(outfile,"  Likelihood of %s: %f\n",modelname,logLikelihood);
		int params = getParams(models[i], tree);
    //fprintf(outfile," --------- params: %i -----------",params);
		//double aic = 2.0*params - 2.0*gsl_sf_log(logLikelihood);
		double aic = 2.0*params - 2.0*logLikelihood;
		if(aic < bestScore) {
			bestScore = aic;
			bestModel = models[i];
		}
	}
	assignModel(tree, bestModel);
  testOpts(tree,bestScore,bestModel);
}


void HLRT(arbre* tree)
{
  Node * root = constructTree();
  int mod = runTests(tree,root,likelihood(tree,JC69),JC69);
  destructTree(root);

  assignModel(tree,mod);
  testOpts(tree,likelihood(tree,tree->mod->whichmodel),tree->mod->whichmodel);

  
  //char outfilename[10];
  //sfprintf(outfile,outfilename, "%s.modeltest",	tree->io->out_tree_file);
  //FILE* outfile = fopen(outfilename,"w");
  //ffprintf(outfile,outfile,"Probability :%f\nTherefore tree->whichmodel = %i\n", likelihood(tree,tree->mod->whichmodel), tree->mod->whichmodel);
  //fclose(outfile);
}

void testOpts(arbre* tree, double bestscore, int bestModel)
{
  fprintf(outfile,"\n  Model testing completed, now testing the Gamma things\n");

	double fscore;
	tree->mod->s_opt->opt_state_freq = 1;
	fscore = likelihood(tree,bestModel);
	double D = 2.0 * (fscore - bestscore);
	double df = (tree->io->mod->datatype == AA)? 19 : 3;
	double c =tree->mod->s_opt->opt_pinvar = 0;
			tree->mod->s_opt->opt_num_param = 0;
			 gsl_cdf_chisq_Qinv(0.5,df);
	if(D >= c)
	{
		bestscore = fscore;
		double gscore;
		tree->mod->s_opt->opt_alpha = 1;
		tree->mod->n_catg = 2;
		gscore = likelihood(tree,bestModel);
		D = 2.0 * (gscore - bestscore);
		df = 1;
		c = gsl_cdf_chisq_Qinv(0.5,df);
		if(D >= c)
		{
			bestscore = gscore;
			tree->mod->n_catg = 4;
			gscore = likelihood(tree,bestModel);
			D = 2.0 * (gscore - bestscore);
			df = 2;
			c = gsl_cdf_chisq_Qinv(0.5,df);
			if(D >= c)
			{
				bestscore = gscore;
				tree->mod->n_catg = 8;
				gscore = likelihood(tree,bestModel);
				D = 2.0 * (gscore - bestscore);
				df = 4;
				c = gsl_cdf_chisq_Qinv(0.5,df);
				if(D >= c)
				{
					bestscore = gscore;
				}
			}
			double iscore;
			tree->mod->invar = 1;
			tree->mod->s_opt->opt_pinvar = 1;
			tree->mod->s_opt->opt_num_param = 1;
			iscore = likelihood(tree,bestModel);
			D = 2.0 * (gscore - bestscore);
			df = 1;
			c = gsl_cdf_chisq_Qinv(0.5,df);
			if(D >= c)
			{
							bestscore = iscore;
							return; //+F+G+I is the best
			}
			tree->mod->invar = 0;
			tree->mod->s_opt->opt_pinvar = 0;
			tree->mod->s_opt->opt_num_param = 0;
			return; //+F+G is the best
		}
	}
	tree->mod->s_opt->opt_alpha = 0;
	tree->mod->n_catg = 1;
	tree->mod->s_opt->opt_state_freq = 1;
	double gscore;
	tree->mod->s_opt->opt_alpha = 1;
	tree->mod->n_catg = 2;
	gscore = likelihood(tree,bestModel);
	D = 2.0 * (gscore - bestscore);
	df = 1;
	c = gsl_cdf_chisq_Qinv(0.5,df);
	if(D >= c)
	{
		bestscore = gscore;
		tree->mod->n_catg = 4;
		gscore = likelihood(tree,bestModel);
		D = 2.0 * (gscore - bestscore);
		df = 2;
		c = gsl_cdf_chisq_Qinv(0.5,df);
		if(D >= c)
		{
			bestscore = gscore;
			tree->mod->n_catg = 8;
			gscore = likelihood(tree,bestModel);
			D = 2.0 * (gscore - bestscore);
			df = 4;
			c = gsl_cdf_chisq_Qinv(0.5,df);
			if(D >= c)
			{
				bestscore = gscore;
			}
		}
		double iscore;
		tree->mod->invar = 1;
		tree->mod->s_opt->opt_pinvar = 1;
		tree->mod->s_opt->opt_num_param = 1;
		iscore = likelihood(tree,bestModel);
		D = 2.0 * (gscore - bestscore);
		df = 1;
		c = gsl_cdf_chisq_Qinv(0.5,df);
		if(D >= c)
		{
			bestscore = iscore;
			return; //+G+I is the best
		}
		tree->mod->invar = 0;
		tree->mod->s_opt->opt_pinvar = 0;
		tree->mod->s_opt->opt_num_param = 0;
		return; //+G is the best
	}
	tree->mod->s_opt->opt_alpha = 0;
	tree->mod->n_catg = 1;
	//The original was the best

  fprintf(outfile,"  Gamma tests complete. Results:\n");
  if(tree->mod->s_opt->opt_state_freq == 1)
    fprintf(outfile,"    +F is enabled\n");
  if(tree->mod->s_opt->opt_alpha == 1){
    fprintf(outfile,"    +G is enabled; n_catg=%i\n",tree->mod->n_catg);
    if(tree->mod->invar == 1)
      fprintf(outfile,"    +I is enabled\n");
  }

  if(tree->mod->s_opt->opt_state_freq!=1 && tree->mod->s_opt->opt_alpha!=1)
    fprintf(outfile,"    nothing enabled\n");

}

void getName(char* ret, int mod){
  if(ret==NULL){
    fprintf(outfile,"WTF!!!!");
  }
  switch(mod){
    case JC69:
      sprintf(ret,"jc69");
      break;
    case K80:
      sprintf(ret,"k80");
      break;
    case F81:
      sprintf(ret,"f81");
      break;
    case HKY85:
      sprintf(ret,"hky85");
      break;
    case F84:
      sprintf(ret,"f84");
      break;
    case TN93:
      sprintf(ret,"tn93");
      break;
    case GTR:
      sprintf(ret,"gtr");
      break;
    case LG:
      sprintf(ret,"lg");
      break;
    case JTT:
      sprintf(ret,"jtt");
      break;
    case WAG:
      sprintf(ret,"wag");
      break;
    case DAYHOFF:
      sprintf(ret,"dayhoff");
      break;
    case BLOSUM62:
      sprintf(ret,"blosum62");
      break;
    case MTREV:
      sprintf(ret,"mtrev");
      break;
    case RTREV:
      sprintf(ret,"rtrev");
      break;
    case CPREV:
      sprintf(ret,"cprev");
      break;
    case DCMUT:
      sprintf(ret,"dcmut");
      break;
    case VT:
      sprintf(ret,"vt");
      break;
    case MTMAM:
      sprintf(ret,"mtmam");
      break;
    case MTART:
      sprintf(ret,"mtart");
      break;
    case HIVW:
      sprintf(ret,"hivw");
      break;
    case HIVB:
      sprintf(ret,"hivb");
      break;
    default:
      sprintf(ret,"UNKNOWN(%i)",mod);
      break;
  }
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

  char x1[10];
  getName(x1,n->mod);
  char x2[10];
  getName(x2,previousMod);

  fprintf(outfile,"  checking %s against %s ...\n",x1,x2);

  double D = (2.0)*(thisLikelihood-previousLikelihood);
  double x = 0.05; //significant P-value
  double df = getParams(previousMod,tree) - getParams(n->mod,tree); //degress of freedom. TODO: make a method to derive this
  if(df<0){
    df = -1.0*df;
  } // absolute value
  double c = gsl_cdf_chisq_Qinv(x,df);

  fprintf(outfile,"    likelihood of %s: %f\n",x1,thisLikelihood);
  fprintf(outfile,"    likelihood of %s: %f\n",x2,previousLikelihood);
  fprintf(outfile,"    D=%f, x=%f, df=%f, so therefore c=%f\n",D,x,df,c);

  if(c>D){
    fprintf(outfile,"    taking a left; %f>%f\n",c,D);
    return (!n->left) ? previousMod : runTests(tree,n->left,previousLikelihood,previousMod);
  } else {
    fprintf(outfile,"    taking a right; %f>%f\n",c,D);
    return (!n->right) ? n->mod : runTests(tree,n->right,thisLikelihood,n->mod);
  }
}

Node * constructTree(){
  Node * root = (Node *)malloc(sizeof(Node));
  root->mod = F81;
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
    case F84:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 1;
      tree->mod->s_opt->opt_lambda = 1;
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
    case TN93:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 1;
      tree->mod->s_opt->opt_lambda = 1;
      tree->mod->s_opt->opt_alpha = 1;
      tree->mod->invar = 0;
      tree->mod->n_l = 1;
      tree->mod->s_opt->opt_state_freq = 0;
      break;
    case GTR:
      tree->mod->datatype = NT;
      tree->mod->n_catg = 1;
      tree->mod->s_opt->opt_kappa = 1;
      tree->mod->s_opt->opt_lambda = 1;
      tree->mod->s_opt->opt_alpha = 0;
      tree->mod->invar = 0;
      tree->mod->n_l = 1;
      tree->mod->s_opt->opt_state_freq = 0;
      break;
    case DAYHOFF:
    case JTT:
    case MTREV:
    case LG:
    case WAG:
    case DCMUT:
    case CPREV:
    case VT:
    case BLOSUM62:
    case MTMAM:
    case MTART:
    case HIVW:
    case HIVB:
      tree->mod->datatype = AA;
      tree->mod->n_catg = 1;
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
  /*
    +g = +(n_catg-1)
    +f = +3 for DNA, +19 for AA
    +i = +1
  */
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
	case HKY85:
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
