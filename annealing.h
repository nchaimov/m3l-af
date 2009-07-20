void Get_TA_neighbor_proposition(arbre *tree);
void Get_QA_neighbor_proposition(arbre *tree);
m3ldbl Temperature(m3ldbl k);
m3ldbl Boltzmann_P(m3ldbl e, m3ldbl enew, m3ldbl temperature);
m3ldbl Thermal_anneal_all_free_params(arbre *tree, int verbose);
m3ldbl Quantum_anneal_all_free_params(arbre *tree, int verbose);

