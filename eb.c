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

	PhyML_Printf("\n . Starting empirical Bayes MCMC of %d generations.\n",
				 tree->io->eb_n_gens);
	PhyML_Printf(" . Sampled trees will be written to %s\n\n", outfilestr);

	/* Record clock time, to estimate time-to-completion */
	time_t start = time(NULL);
	time_t now;

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
	for(i=0; i<tree->io->eb_n_gens-1; i++)
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

		/* print estimated remaining time*/
		if (i%DISPLAY_TIME_REMAINING_PERIOD == 0)
		{
			now = time(NULL);
			Print_time_remaining(now, start, i, tree->io->eb_n_gens);
		}
	}

	/* copy best tree found back into tree */
	Copy_Tree(best_tree, tree);
	Record_Model(best_tree->mod, tree->mod);
	tree->c_lnL = lnL_best;

	printf("\n . Empirical Bayes MCMC, lnL of best sample = %f\n", tree->c_lnL);

	/* clean up and done */
	fclose(outfile);
	gsl_rng_free(rng);
	PhyML_Printf("\n . Done with empirical Bayes MCMC.\n");

	//Calculate_PP(tree);
}

// Calculate posterior probability values for clades,
// given the *.eb file written by the Empricial_Bayes method.
void Calculate_PP(arbre* tree)
{
	PhyML_Printf("\n . Calculating posterior probabilities of clades.\n");

	char outfilestr[1000];
	sprintf(outfilestr, "%s.eb", tree->io->out_tree_file);
	FILE* ebfile = fopen(outfilestr, "r");

	int i,j,k;

	int count_lines = 0; // how many lines (i.e. samples) have we seen?
	struct __Node ***parts; // parts[i] = a list_of_reachable_tips for partition i
	int MAX_PARTS = 10*(2*tree->n_otu-2);
	parts = (node ***)mCalloc(MAX_PARTS,sizeof(node **));
	For(i,MAX_PARTS) parts[i] = (node **)mCalloc(tree->n_otu,sizeof(node *));
	int *size_of_parts; // how big is each partition?
	size_of_parts = (int *)mCalloc(MAX_PARTS,sizeof(int));
	For(i,MAX_PARTS) size_of_parts[i] = 0;
	int *freq_of_parts; // how often have we seen each partition?
	freq_of_parts = (int *)mCalloc(MAX_PARTS,sizeof(int));
	For(i,MAX_PARTS) freq_of_parts[i] = 0;
	int count_parts = 0; // how many partitions exist?

	char line[T_MAX_LINE];
	while( fgets( line, T_MAX_LINE, ebfile) != NULL ){
		char *tokens = strtok( line, " ;");
		tokens = strtok( NULL, " ;");
		//PhyML_Printf("tree=%s\n", tokens);
		count_lines++;
		PhyML_Printf(" . parsing MCMC sample %d\n", count_lines);
		arbre* sampled_tree = Read_Tree(tokens);
		For(j, 2*tree->n_otu-2) // for every node in the tree
		{	if (!tree->noeud[j]->tax) // if the node is not terminal
			{	int found_parti = 0;
				For(i,count_parts) // for each partition
				{	PhyML_Printf("debug: partition %d\n", i);
					// First look to see if the partition i exists in the tree:
					For(k,3)
					{	PhyML_Printf("debug: direction %d\n", k);
						int score = Compare_List_Of_Reachable_Tips(parts[i],
								size_of_parts[i],
								tree->noeud[j]->list_of_reachable_tips[k],
								tree->noeud[j]->n_of_reachable_tips[k]);

						if(score == size_of_parts[i])
						{	PhyML_Printf("debug: found partition %d\n", i);
							found_parti = 1;
							// partition i is present within the sampled_tree
							// increment our frequency count for partition i
							freq_of_parts[i]++;
							break;

						}
					}
				}
				if (found_parti == 0) // this node has not been seen in previous samples.
				{	PhyML_Printf("debug: adding partition\n");
					// ... then add the partition to our list of partitions
					parts[count_parts] = tree->noeud[j]->list_of_reachable_tips[0];
					size_of_parts[count_parts] = tree->noeud[j]->n_of_reachable_tips[0];
					freq_of_parts[count_parts]++;
					count_parts++;
				}
			}
		}
	}

	// debugging code:
	For(i,count_parts)
	{
		PhyML_Printf(" partition %d: size = %d, freq = %d\n", i, size_of_parts[i], freq_of_parts[i]);
	}
	// end debugging code

	fclose(ebfile);
}

