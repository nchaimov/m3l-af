void Get_TA_neighbor_proposition(arbre *tree);
void Get_QA_neighbor_proposition(arbre *tree);
m3ldbl Boltzmann_P(m3ldbl e, m3ldbl enew, m3ldbl temperature);
m3ldbl Thermal_anneal_all_free_params(arbre *tree, int verbose);
m3ldbl Quantum_anneal_all_free_params(arbre *tree, int verbose);
m3ldbl Scale_acceptance_ratio(arbre *tree);

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

