/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef ALRT_H
#define ALRT_H

void aLRT(arbre *tree);
int Check_NNI_Five_Branches(arbre *tree);
int Compute_Likelihood_Ratio_Test(edge *tested_edge, arbre *tree);
int NNI_Neigh_BL(edge *b_fcus, arbre *tree);
void Make_Target_Swap(arbre *tree, edge *b_fcus, int swaptodo);
m3ldbl Statistics_To_Probabilities(m3ldbl in);
m3ldbl Statistics_To_RELL(arbre *tree);
m3ldbl Statistics_To_SH(arbre *tree);
m3ldbl Update_Lk_At_Given_Edge_Excluding(edge *b_fcus, arbre *tree, node *exclude);

#endif
