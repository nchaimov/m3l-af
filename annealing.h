void Get_TA_Neighbor_Proposition(arbre *tree);
void Get_QA_Neighbor_Proposition(arbre *tree);
m3ldbl Boltzmann_P(m3ldbl e, m3ldbl enew, m3ldbl temperature);
m3ldbl Thermal_Anneal_All_Free_Params(arbre *tree, int verbose);
m3ldbl Quantum_Anneal_All_Free_Params(arbre *tree, int verbose);
m3ldbl Scale_Acceptance_Ratio(arbre *tree);
void Step_Brlen_Proportion(arbre *tree);
void Swap_Tree_And_Mod(arbre *tree1, arbre *tree2);
void Step_Gamma(arbre *tree);
void Step_Branch_Lengths(arbre *tree);
void Step_Topology(arbre *tree);
void Step_Pi(arbre * tree);
void Step_Lambda(arbre * tree);
void Step_Kappa(arbre * tree);
void Step_RR(arbre * tree);
void Step_Pinvar(arbre *tree);

//
// Global parameters.
// We should provide the option to set these parameters at runtime,
// but we should also provide default values.
//
typedef struct __ANNEALING{
	int temp_count; 	// the number of temperate values to try
	m3ldbl start_temp;	// initial temperature
	m3ldbl end_temp;	// order of magnitude of desired accuracy
	int iters_per_temp;	// number of tries at each temperature
	int set_back;	// number of steps to step back if improving
	m3ldbl accept_ratio;	// acceptance ratio
}annealing;
