#include "eb.h"
#include "utilities.h"
#include "optimiz.h"
#include "spr.h"

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

	/* create storage in tree */
	tree->best_tree = Make_Tree(tree->n_otu);
	Init_Tree(tree->best_tree,tree->n_otu);
	Make_All_Tree_Nodes(tree->best_tree, tree->mod->n_l);
	Make_All_Tree_Edges(tree->best_tree, tree->mod->n_l);
	tree->best_tree->mod = Copy_Model(tree->mod);

	tree->old_tree = Make_Tree(tree->n_otu);
	Init_Tree(tree->old_tree,tree->n_otu);
	Make_All_Tree_Nodes(tree->old_tree, tree->mod->n_l);
	Make_All_Tree_Edges(tree->old_tree, tree->mod->n_l);
	tree->old_tree->mod = Copy_Model(tree->mod);
	
	PhyML_Printf("\nstarting empirical Bayes MCMC of %d generations.\n", 
				 tree->io->eb_n_gens);
	PhyML_Printf("Sampled trees will be written to %s\n\n", outfilestr);
	
	/* optimize starting tree and sample it */
	Round_Optimize(tree,tree->data,ROUND_MAX);
	fprintf(outfile, "[%f] ", tree->c_lnL);
	Print_Tree(outfile, tree);

	/* current tree is best so far */
	Copy_Tree(tree, tree->best_tree);
	Record_Model(tree->mod, tree->best_tree->mod);
	tree->best_lnL = tree->c_lnL;

	m3ldbl last_lnL;
	
	/* MCMC routine */
	for(i=0; i<tree->io->eb_n_gens; i++)
	{
		/* copy current tree into old tree */
		Copy_Tree(tree, tree->old_tree);
		Record_Model(tree->mod, tree->old_tree->mod);
		last_lnL = tree->c_lnL;
		
		/* modify current tree */
		Random_Spr(1, tree);
		
		/* optimize current tree */
		Round_Optimize(tree,tree->data,ROUND_MAX);
		
		/* check acceptance */
		printf("LR %f/%f [%f]\n", tree->c_lnL, last_lnL, exp(tree->c_lnL - last_lnL));
		
		if(tree->c_lnL >= last_lnL || gsl_rng_uniform(rng) < exp(tree->c_lnL - last_lnL))
		{			
			printf("OKAY!\n");
			
			/* check if better than best so far */
			if(tree->c_lnL > tree->best_lnL)
			{
				Copy_Tree(tree, tree->best_tree);
				Record_Model(tree->mod, tree->best_tree->mod);
				tree->best_lnL = tree->c_lnL;
			}
		}
		
		else
		{
			/* copy old tree back to current tree */
			Copy_Tree(tree->old_tree, tree);
			Record_Model(tree->old_tree->mod, tree->mod);
			tree->c_lnL = last_lnL;
		}
		
		/* print sampled tree */
		fprintf(outfile, "[%f] ", tree->c_lnL);
		Print_Tree(outfile, tree);
	}
		
	/* copy best tree found back into tree */
	Copy_Tree(tree->best_tree, tree);
	Record_Model(tree->best_tree->mod, tree->mod);
	tree->c_lnL = tree->best_lnL;
	
	printf("BEST: %f\n", tree->c_lnL);
	
	/* clean up and done */
	fclose(outfile);
	gsl_rng_free(rng);
	PhyML_Printf("\ndone empirical Bayes MCMC.\n");
}

