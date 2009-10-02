/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef ML_H
#define ML_H


void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, plkflt *p_lk);
void Init_Tips_At_One_Site_AA_Float(char aa, int pos, plkflt *p_lk);
void Get_All_Partial_Lk_Scale(arbre *tree,edge *b_fcus,node *d);
void Post_Order_Lk(node *pere, node *fils, arbre *tree);
void Pre_Order_Lk(node *pere, node *fils, arbre *tree);
void Lk(arbre *tree);
void debug_Lk_nocompress(arbre *tree);
void Site_Lk(arbre *tree, int site);
m3ldbl Lk_At_Given_Edge(edge *b_fcus,arbre *tree);
m3ldbl Return_Lk(arbre *tree);
m3ldbl Return_Abs_Lk(arbre *tree);
matrix *ML_Dist(allseq *data,model *mod);
m3ldbl Lk_Given_Two_Seq(allseq *data,int numseq1,int numseq2,m3ldbl dist,model *mod,m3ldbl *loglk);
void Unconstraint_Lk(arbre *tree);
void Update_P_Lk(arbre *tree,edge *b_fcus,node *n);
void Make_Tree_4_Lk(arbre *tree,allseq *alldata,int n_site);
void Init_P_Lk_Tips_Double(arbre *tree);
void Init_P_Lk_Tips_Int(arbre *tree);
void Init_P_Lk_At_One_Node(node *a, arbre *tree);
void Update_PMat_At_Given_Edge(edge *b_fcus, arbre *tree);
void Sort_Sites_Based_On_Lk(arbre *tree);
void Get_Partial_Lk_Scale(arbre *tree, edge *b_fcus, node *a, node *d);
void Get_Partial_Lk(arbre *tree, edge *b_fcus, node *a, node *d);
void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars);
void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars);
void Update_P_Lk_Along_A_Path(node **path, int path_length, arbre *tree);
m3ldbl Lk_Dist(m3ldbl *F, m3ldbl *dist, model *mod, arbre *tree);
m3ldbl Lk_Dist_No_Bl(m3ldbl *F, m3ldbl dist, model *mod);
m3ldbl Update_Lk_At_Given_Edge(edge *b_fcus, arbre *tree);
void Update_P_Lk_Greedy(arbre *tree, edge *b_fcus, node *n);
void Get_All_Partial_Lk_Scale_Greedy(arbre *tree, edge *b_fcus, node *a, node *d);
m3ldbl Lk_Core(edge *b, arbre *tree, int site);
m3ldbl Lk_Triplet(node *a, node *d, arbre *tree);
void Print_Lk_Given_Edge_Recurr(node *a, node *d, edge *b, arbre *tree);
m3ldbl *Post_Prob_Rates_At_Given_Edge(edge *b, m3ldbl *post_prob, arbre *tree);
m3ldbl Lk_With_MAP_Branch_Rates(arbre *tree);



#endif






