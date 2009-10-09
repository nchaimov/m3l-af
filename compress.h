#ifdef COMPRESS_SUBALIGNMENTS
void Compute_Red_Arrays(node *a, node *d, arbre *tree);
void Post_Order_Lk_Red(node *a, node *d, arbre *tree, int site);
void Pre_Order_Lk_Red(node *a, node *d, arbre *tree, int site);
void Update_P_Lk_Red(arbre *tree, edge *b, node *n, int site);
void Make_All_Nodes_Red(arbre *tree);
void Init_All_Nodes_Red(arbre *tree);
#endif
