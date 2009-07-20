void Get_TA_Neighbor_Proposition(arbre *tree);
void Get_QA_Neighbor_Proposition(arbre *tree);
m3ldbl Boltzmann_P(m3ldbl e, m3ldbl enew, m3ldbl temperature);
m3ldbl Thermal_Anneal_All_Free_Params(arbre *tree, int verbose);
m3ldbl Quantum_Anneal_All_Free_Params(arbre *tree, int verbose);
m3ldbl Scale_Acceptance_Ratio(arbre *tree);
void Step_Brlen_Proportion(arbre *tree);
void Step_Gamma(arbre *tree);
void Step_Branch_Lengths(arbre *tree);
void Step_Topology(arbre *tree);

//
// Global parameters.
// We should provide the option to set these parameters at runtime,
// but we should also provide default values.
//
m3ldbl temp_count; 	// the number of temperate values to try
m3ldbl start_temp;	// initial temperature
m3ldbl end_temp;	// order of magnitude of desired accuracy
m3ldbl iters_per_temp;	// number of tries at each temperature
m3ldbl set_back;	// amount to step back if improving
m3ldbl accept_ratio;	// acceptance ratio

