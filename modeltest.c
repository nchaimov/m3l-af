AIC(arbre* tree)
{

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


	// 2. Call this method to initialize several important values in the 'mod' struct.
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
	PhyML_Printf("\n. Log likelihood of the current tree: %f.\n",tree->c_lnL);

	//
	// 6. Record the likelihood of tree, given this model
	// (this step is for you to create)

}