//
// This file contains code relating to Edge Migration (EMIG)
//


m3ldbl emig_sigma; // for selecting a perturbation distance (a float from a Gaussian distribution)

void Emig_Swap(arbre *tree, edge *a, node *v_a, edge *b, node *v_b);
void Migrate_One_Edge(arbre *tree);
