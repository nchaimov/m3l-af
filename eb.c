#include "eb.h"
#include "utilities.h"
#include "optimiz.h"
#include "spr.h"
#include "lk.h"
#include "free.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>

void Empirical_Bayes(arbre* tree)
{
	int CHAIN_THINNING = 100;
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

	/* optimize starting tree */
	Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	//fprintf(outfile, "[%f] ", tree->c_lnL);
	//Print_Tree(outfile, tree);

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

	//Prepare_Tree_For_Lk(last_tree);
	//Prepare_Tree_For_Lk(best_tree);


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
			//Prepare_Tree_For_Lk(last_tree);


			/* check if better than best so far */
			if(tree->c_lnL > lnL_best)
			{
				Copy_Tree(tree, best_tree);
				Record_Model(tree->mod, best_tree->mod);
				lnL_best = tree->c_lnL;
				//Prepare_Tree_For_Lk(best_tree);
			}
		}
		else
		{
			/* copy old tree back to current tree */
			Copy_Tree(last_tree, tree);
			Record_Model(last_tree->mod, tree->mod);
			tree->c_lnL = lnL_current;
			//Prepare_Tree_For_Lk(tree);
		}

		/* print sampled tree */
		if(i % CHAIN_THINNING == 0)
		{
			fprintf(outfile, "[%f] ", tree->c_lnL);
			Print_Tree(outfile, tree);
		}

		/* print estimated remaining time*/
		if (i%DISPLAY_TIME_REMAINING_PERIOD == 0)
		{	now = time(NULL);
			Print_time_remaining(now, start, i, tree->io->eb_n_gens);
		}
	}

	/* copy best tree found back into tree */
	Copy_Tree(best_tree, tree);
	Record_Model(best_tree->mod, tree->mod);
	//Prepare_Tree_For_Lk(tree);
	tree->c_lnL = lnL_best;

	/* clean up and done */
	fclose(outfile);
	gsl_rng_free(rng);
	PhyML_Printf("\n. Done with empirical Bayes MCMC.\n");

	Calculate_PP(tree);

	Free_Tree(best_tree);
	Free_Tree(last_tree);
}

// Calculate posterior probabilities of clades,
// given the *.eb file written by the method named "Empricial_Bayes"
void Calculate_PP(arbre* tree)
{
	PhyML_Printf("\n. Calculating posterior probabilities of clades.");

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

	/*
	 * 1. Parse the MCMC samples and count the frequency of topological partitions among samples
	 */
	char line[T_MAX_LINE];
	while( fgets( line, T_MAX_LINE, ebfile) != NULL ){
		char *tokens = strtok( line, " ;");
		tokens = strtok( NULL, " ;");
		count_lines++;
		PhyML_Printf(".");

		arbre* sampled_tree = Read_Tree(tokens); // Read the Newick-formatted tree from the *.eb file and create an arbre structure. . .
		sampled_tree->mod         = tree->mod;
		sampled_tree->io          = tree->io;
		sampled_tree->data        = tree->data;
		sampled_tree->both_sides  = tree->both_sides;
		sampled_tree->n_pattern   = tree->n_pattern;
		Fill_Dir_Table(sampled_tree);
		Update_Dirs(sampled_tree);

		For(j,2*sampled_tree->n_otu-3) // for every edge in the tree...
		{
			/*
			 * Does edge j connect to a terminal taxa?
			 * If so, then this partition will *always* have posterior probability = 1.0.
			 * Therefore, let's not bother computing the PP for this edge.
			 */
			if (sampled_tree->t_edges[j]->rght->tax || sampled_tree->t_edges[j]->left->tax)
			{	continue;
			}

			/*
			 * Here we find that largest partition attached to edge j.
			 * The tips in the largest partition will be saved in tiplist.
			 */
			node **tiplist;
			int proposed_part_size = 0;
			For(k,3)
			{
				if (sampled_tree->t_edges[j]->left->n_of_reachable_tips[k] > proposed_part_size)
				{	tiplist = sampled_tree->t_edges[j]->left->list_of_reachable_tips[k];
					proposed_part_size = sampled_tree->t_edges[j]->left->n_of_reachable_tips[k];
				}
				if (sampled_tree->t_edges[j]->rght->n_of_reachable_tips[k] > proposed_part_size)
				{	tiplist = sampled_tree->t_edges[j]->rght->list_of_reachable_tips[k];
					proposed_part_size = sampled_tree->t_edges[j]->rght->n_of_reachable_tips[k];
				}
			}

			int found_parti = 0;
			For(i,count_parts) // for each partition...
			{	// do partition i and subtree k have the same number of descendants?
				if (size_of_parts[i] == proposed_part_size)
				{	int score = Compare_List_Of_Reachable_Tips_version2(parts[i],
							size_of_parts[i],
							tiplist,
							proposed_part_size);
					// do partition i and subtree k contain the exact same set of descendant taxa?
					if(score == size_of_parts[i])
					{	// ... then partition i is present within the sampled_tree.
						// Increment the observed frequency for partition i:
						found_parti = 1;
						freq_of_parts[i]++;
						break;
					}
					// a case we must deal with: if the partition size = 1/2 the number of terminal taxa, and if
					// the match score = 0, then we found the taxa in tiplist are the inverse of the taxa
					// in partition i.  In other words, we found the partition in sampled_tree.
					if (score == 0 && proposed_part_size == sampled_tree->n_otu * 0.5)
					{
						found_parti = 1;
						freq_of_parts[i]++;
						break;
					}
				}
				if(found_parti == 1) break;
			}
			if (found_parti == 0) // if the partition of node j has not been observed in previous samples...
			{	// ... then add the partition to our collection of observed partitions.
				parts[count_parts] = tiplist;
				size_of_parts[count_parts] = proposed_part_size;
				freq_of_parts[count_parts]++; // ... and increment the observed frequency of this partition.
				count_parts++;
			}
		}
	}
	PhyML_Printf("\n");
	fclose(ebfile);

	/*
	 * 2. Scale the raw frequency counts to posterior probabilities.
	 */
	double *pp_of_parts; // how often have we seen each partition?
	pp_of_parts = (double *)mCalloc(count_parts,sizeof(double));
	For(i,count_parts) pp_of_parts[i] = 0.0;
	For(i,count_parts)
	{	pp_of_parts[i] = (double)freq_of_parts[i] / (double)count_lines;
	}

	/*
	 * 3. Write posterior probabilities onto the branches of the given tree
	 */
	For(j,2*tree->n_otu-3) // for every edge in the tree...
	{
		/*
		 * Does edge j connect to a terminal taxa?
		 * If so, then this partition will *always* have posterior probability = 1.0.
		 * Therefore, let's not bother computing the PP for this edge.
		 */
		if (tree->t_edges[j]->rght->tax || tree->t_edges[j]->left->tax)
		{	tree->t_edges[j]->post_prob = 1.0;
		}

		/*
		 * Determine which partition matches edge j
		 */
		node **tiplist;
		int proposed_part_size = 0;
		For(k,3)
		{
			if (tree->t_edges[j]->left->n_of_reachable_tips[k] > proposed_part_size)
			{	tiplist = tree->t_edges[j]->left->list_of_reachable_tips[k];
				proposed_part_size = tree->t_edges[j]->left->n_of_reachable_tips[k];
			}
			if (tree->t_edges[j]->rght->n_of_reachable_tips[k] > proposed_part_size)
			{	tiplist = tree->t_edges[j]->rght->list_of_reachable_tips[k];
				proposed_part_size = tree->t_edges[j]->rght->n_of_reachable_tips[k];
			}
		}
		For(i,count_parts)
		{	if (size_of_parts[i] == proposed_part_size)
			{	int score = Compare_List_Of_Reachable_Tips_version2(parts[i],
						size_of_parts[i],
						tiplist,
						proposed_part_size);
				if( (score == size_of_parts[i]) ||
						(score == 0 && proposed_part_size == tree->n_otu * 0.5) )
				{	tree->t_edges[j]->post_prob = pp_of_parts[i];
					break;
				}
			}
		}
	}
	tree->print_pp_val = 1;

	/*
	 * debugging code: spew PP results
	 */
	/*
	Print_Tree_Screen(tree);
	For(i,count_parts)
	{
		PhyML_Printf(" partition %d: size = %d, freq = %d, pp = %f\n", i, size_of_parts[i], freq_of_parts[i], pp_of_parts[i]);
		For(j,size_of_parts[i])
		{	PhyML_Printf("%s ", parts[i][j]->name );
		}
		PhyML_Printf("\n");
	}
	*/
	// end debugging code
}

