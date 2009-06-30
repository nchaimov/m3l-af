/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef RATES_H
#define RATES_H

void RATES_Monte_Carlo_Mean_Rates(arbre *tree);
void RATES_Monte_Carlo_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl curr_rate, arbre *tree);
void RATES_Print_Rates(arbre *tree);
void RATES_Print_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
trate *RATES_Make_Rate_Struct(int n_otu);
void RATES_Init_Rate_Struct(trate *rates, int n_otu);
void RATES_Classify_Branches(arbre *tree);
void RATES_Adjust_Rates(arbre *tree);
void RATES_Adjust_Rates_Local_Pre(node *a, node *d, edge *b, arbre *tree);
void RATES_Adjust_Rates_Local(node *a, node *d, edge *b1, arbre *tree);
void RATES_Record_T(arbre *tree);
void RATES_Restore_T(arbre *tree);
void RATES_Monte_Carlo_Mean_Rates_Core(phydbl t_lim_sup, phydbl t_lim_inf, phydbl *curr_rate, phydbl *mean_rate, phydbl lexp, phydbl alpha);
phydbl RATES_Lk_Rates(arbre *tree);
void RATES_Lk_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
void RATES_Fill_Node_Rates_Pre(node *a, node *d, edge *b, phydbl *node_r, arbre *tree);
void RATES_Fill_Node_Rates(phydbl *node_r, arbre *tree);
void RATES_Optimize_Node_Times_Serie_Fixed_Br_Len(node *a, node *d, arbre *tree);
void RATES_Optimize_Lexp(arbre *tree);
void RATES_Round_Optimize(arbre *tree);
void RATES_Optimize_Lexp(arbre *tree);
void RATES_Optimize_Alpha(arbre *tree);
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens);
phydbl RATES_Dr_X_Dx(phydbl r, phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y_Trpzd(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp, 
			       int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu_Given_Y_Std(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y_Romb(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_Given_Y(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dy_Given_Mu(phydbl mu, phydbl y, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(phydbl mu1, phydbl mu2, phydbl y, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Trpz(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp,
				 int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu2_And_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Romb(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl a, phydbl b, phydbl lexp);
void RATES_Random_Branch_Lengths(arbre *tree);
void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param);
void RATES_Set_Node_Times(arbre *tree);
void RATES_Set_Node_Times_Pre(node *a, node *d, arbre *tree);
void RATES_Randomize_Node_Times(arbre *tree);
void RATES_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree);
void RATES_Optimize_Node_Times(arbre *tree);
phydbl RATES_Exp_Y(phydbl mu1, phydbl mu2, phydbl dt1, phydbl lexp);
phydbl RATES_Dmu2_Given_Mu1_Bis(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp);
void RATES_Replace_Br_Lengths_By_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
void RATES_Replace_Br_Lengths_By_Rates(arbre *tree);

void RATES_Get_Mean_Rates(arbre *tree);
void RATES_Get_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl *r_a, phydbl alpha, phydbl lexp, arbre *tree);
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl *r_beg, phydbl alpha, phydbl lexp, int *n_jumps, phydbl *mean_r);
void RATES_Optimize_Clock_Rate(arbre *tree);
phydbl RATES_Dmu1_Given_Lbda_And_Mu2(phydbl lbda, phydbl mu1, phydbl mu2, phydbl alpha, phydbl beta);
phydbl RATES_Dmu1_And_Mu2_One_Jump_Trpz(phydbl mu1, phydbl mu2, phydbl a, phydbl b,
					  int nsteps, phydbl beg, phydbl end, phydbl prevs);
phydbl RATES_Dmu1_And_Mu2_One_Jump_One_Interval(phydbl mu1, phydbl mu2, phydbl a, phydbl b);
phydbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl dt1, phydbl dt2, phydbl mu1, phydbl mu2, phydbl a, phydbl b);

phydbl RATES_Dmu1_And_Mu2_One_Jump_Old(phydbl mu1, phydbl mu2, phydbl a, phydbl b);

phydbl RATES_Dmu2_And_Min_N_Given_Mu1(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_And_Mu1_Given_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Lk_Rates_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, arbre *tree);
void RATES_Init_Triplets(arbre *tree);
phydbl RATES_Lk_Change_One_Time(node *n, phydbl new_t, arbre *tree);
void RATES_Update_Triplet(node *n, arbre *tree);
void RATES_Print_Triplets(arbre *tree);
phydbl RATES_Lk_Change_One_Rate(node *d, phydbl new_rate, arbre *tree);
phydbl RATES_Dmu2_And_Mu1_Given_Min_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu2_And_Mu1_Given_N_Normal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Coeff_Corr(phydbl alpha, phydbl beta, int n1, int n2);
phydbl RATES_Dmu2_And_Mu1_Given_N_Full(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu1_Given_V_And_N(phydbl mu1, phydbl v, int n, phydbl dt1, phydbl a, phydbl b);
phydbl RATES_Yule(arbre *tree);
phydbl RATES_Check_Mean_Rates(arbre *tree);
void RATES_Check_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl *sum, arbre *tree);
void RATES_Adjust_Clock_Rate(arbre *tree);
void RATES_Discretize_Rates(arbre *tree);
void RATES_Discretize_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
phydbl RATES_Dmu_Given_V_And_MinN(phydbl mu, phydbl dt, phydbl v, int minn, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Dmu_One(phydbl mu, phydbl dt, phydbl a, phydbl b, phydbl lexp);
phydbl RATES_Compound_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx);
void RATES_Record_Rates(arbre *tree);
void RATES_Reset_Rates(arbre *tree);
void RATES_Record_Times(arbre *tree);
void RATES_Reset_Times(arbre *tree);
void RATES_Update_T_Rates_Pre(node *a, node *d, arbre *tree);
void RATES_Update_T_Rates(arbre *tree);
void RATES_Get_Br_Len(arbre *tree);
void RATES_Get_Rates_From_Bl(arbre *tree);
phydbl RATES_Compound_Core_Joint(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, 
				 phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx);
phydbl RATES_Dmu_Joint(phydbl mu, int n, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n);
phydbl RATES_Compound_Core_Marginal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, 
				    phydbl beta, phydbl lexp, phydbl eps, int approx);
phydbl RATES_Lk_Jumps(arbre *tree);
void RATES_Set_Rates_Prior_Mean(arbre *tree);
void RATES_Set_Rates_Prior_Mean_Pre(node *a, node *d, arbre *tree);
void RATES_Set_Rates_Post_Mean_And_Cov(arbre *tree);

#endif
