#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"
#include "rates.h"
#include "annealing.h"


//do we want to make a new tree here??
//otherwise we should return a void and just do work on the input tree
void Get_TA_neighbor_proposition(arbre *tree){


}

void Get_QA_neighbor_proposition(arbre *tree){


}

// INPUT: Given that we're on iteration k, what is the scheduled temperature for this iteration?
//
// OUTPUT: a temperature value, as a floating-point decimal
m3ldbl Temperature(m3ldbl k){
	m3ldbl result = 1.0;

	return result;
	// Here is where we program the annealing schedule.
}

// INPUT: Given the energy (e) of our current state, the energy (enew) of our new proposed state,
// and the temperature (temperature) of the system, calculate the probability of the new state having
// energy equal to enew, where the probability is drawn from a Boltzmann distribution.
//
// OUTPUT: a floating point decimal, between 0.0 and 1.0
m3ldbl Boltzmann_P(m3ldbl e, m3ldbl enew, m3ldbl temperature){
	m3ldbl result = 1.0;
	return result;
}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
m3ldbl Thermal_anneal_all_free_params(arbre *tree, int verbose){
	m3ldbl result = 1.0;
	/* Here is psudocode for the algorithm:
	 * s = input tree with input parameters
	 * e = the likelihood of s
	 * sbest = s
	 * ebest = e
	 * iter = 0
	 * while (iter < MAX_QA_ITERS):
	 * 		snew = Get_TA_neighbor_proposition(s)
	 * 		enew = likelihood of snew
	 * 		if enew < ebest:
	 * 			sbest = snew
	 * 			ebest = enew
	 * 		else if Boltzmann_P(e, enew, Temperature(iter / MAX_QA_ITERS) ) > random(0,1):
	 * 			s = snew
	 * 			e = enew
	 * 		iter = iter + 1
	 * return sbest
	 */

	return result;
}

// INPUT: a tree structure, with parameters to optimize specified in tree->mod->s_opt
// OUTPUT: the likelihood of the best found tree
m3ldbl Quantum_anneal_all_free_params(arbre *tree, int verbose){
	m3ldbl result = 1.0;

	return result;

}

