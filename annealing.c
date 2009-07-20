#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"
#include "rates.h"
#include "annealing.h"


// This method returns the scaled and estimated acceptance ratio for thermal annealing.
// The variable "tree" will be modified by this method, such that tree will be in a good
// starting position for thermal annealing after this method concludes.
m3ldbl Scale_acceptance_ratio(arbre *tree){

	/* psuedo-code:

	lnL_best = tree's lnL
	lnL_proposed = max_lnL

	fsum = 0.0 // function mean
	fsqsum  0.0 // function standard deviation
	n = 0 // number of points tried

	tree_proposed = copy of tree

	for iter in range(0, iters_per_temp):
		Get_TA_neighbor_proposition(tree_proposed)
		lnL_proposed = tree_proposed's likelihood
		if (lnL_proposed > -DBL_MAX):
			n++
			fsum += working_lnL
			fsqsum += working_lnL * working_lnL
		if (lnL_proposed > lnL_best):
			lnL_best = lnL_proposed
			tree = copy of tree_proposed

	fsum = fsum / n
	fsqsum = fsqsum / n
	fsqsum = sqrt( fsqsum - fsum * fsum)
	aratio = accept_ratio * fsqsum / start_temp
	return aratio
	*/
}

// This method modifies some or all parameters of "tree", including:
// * topology
// * branch lengths
// * branch length proportions
// * alpha (for gamma distributed ASRV)
// * gamma proportions
// * mu (evolutionary rate)
void Get_TA_neighbor_proposition(arbre *tree){
	// For now, we'll always perturb every parameter.
	// In the future, we'll do something more sophisticated, where the probability of perturbing
	// any particular parameter will be drawn from a probability distribution.
	Step_brlen_proportion(tree);
	Step_branch_lengths(tree);
	Step_gamma(tree);
	Step_topology(tree);
}

// helper for "Get_TA_neighbor_proposition"
void Step_brlen_proportion(tree){
	// For now, do something simple and stupid, like stepping branch lengths by +/- 0.001.
	// We'll eventually do something more sophisticated, like selected the new branch
	// lengths from a Dirichlet distribution.
}

// helper for "Get_TA_neighbor_proposition"
void Step_gamma(tree){
	// For now, do something stupid like incrementing gamma +/- 0.001.
	// In reality, we'll want to perturb the gamma parameter based on a normal distribution
	// with mean equal to the current value of gamma and standard deviation equal to some
	// user-specified parameter sigma.
}

// helper for "Get_TA_neighbor_proposition"
void Step_branch_lengths(tree){
	// For now, do something stupid like incrementing lengths by +/- 0.001.
}

// helper for "Get_TA_neighbor_proposition"
void Step_topology(tree){
	// invoke NNI or SPR to create a new topology.
}


void Get_QA_neighbor_proposition(arbre *tree){
}

// INPUT: Given the likelihood (lnl_curr) of our current state, the likelihood (lnl_new) of our new proposed state,
// and the temperature (temperature) of the system, calculate the probability of the new state,
// where the probability is drawn from a Boltzmann distribution.
//
// OUTPUT: a floating point decimal, between 0.0 and 1.0
m3ldbl Boltzmann_P(m3ldbl lnl_curr, m3ldbl lnl_new, m3ldbl temperature){
	m3ldbl result = 1.0;
	/* psuedocode:
	if (lnl_new < lnl_curr)
	{
		result = exp( (lnl_new - lnl_curr)/(accept_ratio * temperature) );
	}
	*/
	return result;
}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
//
// At the end of this method, the object "tree" will contain the best-found topology, branch lengths,
// and model parameters.
m3ldbl Thermal_anneal_all_free_params(arbre *tree, int verbose){
	/* Here is psuedocode for the algorithm:
	 *
	 * m3ldbl temp = start_temp // start_temp is a global
	 * m3ldbl tempmult = exp(log(end_temp / start_temp) / double(temp_count - 1) )
	 *
	 * aratio = Scale_acceptance_ratio(tree) // accept_ratio is a global value
	 *
	 * lnL_best = tree's lnL
	 * tree_current = copy of tree
	 * lnL_current = tree_current's lnL
	 * tree_proposed = copy of tree
	 * lnL_proposed = tree_proposed's lnL
	 *
	 * for (int itemp = 0; itemp < temp_count; itemp++): // temp_count is a global
	 *
	 *		// Here we recenter the search at each temperature, using the best values found thus far
	 * 		tree_current = copy of tree
	 * 		lnL_current = lnL_best
	 *
	 *		steps_tries = 0
	 * 		steps_accepted = 0
	 *
	 *		for (int iter = 0; iter < iters_per_temp; iter++):
	 *			steps_tried++;
	 *			tree_proposed = copy of tree_current
	 *			Get_TA_neighbor_proposition(tree_proposed)
	 *			lnL_proposed = tree_proposed's likelihood
	 *			if lnL_proposed > lnL_best:
	 *				lnL_best = lnL_proposed
	 *				tree = copy of tree_proposed
	 *				iter -= set_back // keep going at this temperature if we are improving the likelihood
	 *				if (iter < 0):
	 *					iter = 0
	 *			acc_prob = Boltzmann_P(lnL_proposed, lnL_current, itemp)
	 *			if acc_prob > random(0,1):
	 *				steps_accepted++
	 *				lnL_current = lnL_proposed
	 *				tree_current = copy of tree_proposed
	 *
	 *		temp *= tempmult // reduce the temperature with geometric descent
	 */

}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
m3ldbl Quantum_anneal_all_free_params(arbre *tree, int verbose){
	m3ldbl result = 1.0;

	return result;

}

