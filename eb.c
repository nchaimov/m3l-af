#include "eb.h"
#include "utilities.h"
#include "optimiz.h"
#include "spr.h"
#include "lk.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

void Empirical_Bayes(arbre* tree)
{
	int i;
	
	char outfilestr[1000];
	sprintf(outfilestr, "%s.eb", tree->io->out_tree_file);
	
	FILE* outfile = fopen(outfilestr, "w");
	
	if(!outfile)
	{
		fprintf(stderr, "ERROR: cant open output file for EB: %s\n",
			    outfilestr);
		exit(1);	
	}
	
	/* create random number generator */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, (int)time(NULL));
	
	PhyML_Printf("\nstarting empirical Bayes MCMC of %d generations.\n", 
				 tree->io->eb_n_gens);
	PhyML_Printf("Sampled trees will be written to %s\n\n", outfilestr);
	
	/* optimize starting tree and sample it */
	Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	fprintf(outfile, "[%f] ", tree->c_lnL);
	Print_Tree(outfile, tree);
	
	/* create storage for trees and models */
	arbre *best_tree = Make_Tree(tree->n_otu);
	Init_Tree(best_tree,tree->n_otu);
	Make_All_Tree_Nodes(best_tree, tree->mod->n_l);
	Make_All_Tree_Edges(best_tree, tree->mod->n_l);
	best_tree->mod = Copy_Model(tree->mod);
	Copy_Tree(tree, best_tree);
	Record_Model(tree->mod, best_tree->mod);

	arbre *last_tree = Make_Tree(tree->n_otu);
	Init_Tree(last_tree,tree->n_otu);
	Make_All_Tree_Nodes(last_tree, tree->mod->n_l);
	Make_All_Tree_Edges(last_tree, tree->mod->n_l);
	last_tree->mod = Copy_Model(best_tree->mod);
	Copy_Tree(tree, last_tree);
	Record_Model(tree->mod, last_tree->mod);

	m3ldbl lnL_best = tree->c_lnL;
	m3ldbl lnL_current = tree->c_lnL;
	
	/* MCMC routine */
	for(i=0; i<tree->io->eb_n_gens; i++)
	{		
		/* modify current tree */
		Random_Spr(1, tree);
		
		/* optimize current tree */
		Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
		
		/* check acceptance */
		//printf("LR %f/%f [%f]\n", tree->c_lnL, lnL_current, exp(tree->c_lnL - lnL_current));
		
		if(tree->c_lnL >= lnL_current || gsl_rng_uniform(rng) < exp(tree->c_lnL - lnL_current))
		{			
			Copy_Tree(tree, last_tree);
			Record_Model(tree->mod, last_tree->mod);
			lnL_current = tree->c_lnL;
			
			/* check if better than best so far */
			if(tree->c_lnL > lnL_best)
			{
				Copy_Tree(tree, best_tree);
				Record_Model(tree->mod, best_tree->mod);
				lnL_best = tree->c_lnL;
			}
		}
		
		else
		{
			/* copy old tree back to current tree */
			Copy_Tree(last_tree, tree);
			Record_Model(last_tree->mod, tree->mod);
			tree->c_lnL = lnL_current;
		}
		
		/* print sampled tree */
		fprintf(outfile, "[%f] ", tree->c_lnL);
		Print_Tree(outfile, tree);
	}
		
	/* copy best tree found back into tree */
	Copy_Tree(best_tree, tree);
	Record_Model(best_tree->mod, tree->mod);
	tree->c_lnL = lnL_best;
	
	printf("BEST: %f\n", tree->c_lnL);
	
	/* clean up and done */
	fclose(outfile);
	gsl_rng_free(rng);
	PhyML_Printf("\ndone empirical Bayes MCMC.\n");
}

