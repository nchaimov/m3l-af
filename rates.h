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
void RATES_Monte_Carlo_Mean_Rates_Pre(node *a, node *d, edge *b, m3ldbl curr_rate, arbre *tree);
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
void RATES_Monte_Carlo_Mean_Rates_Core(m3ldbl t_lim_sup, m3ldbl t_lim_inf, m3ldbl *curr_rate, m3ldbl *mean_rate, m3ldbl lexp, m3ldbl alpha);
m3ldbl RATES_Lk_Rates(arbre *tree);
void RATES_Lk_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
void RATES_Fill_Node_Rates_Pre(node *a, node *d, edge *b, m3ldbl *node_r, arbre *tree);
void RATES_Fill_Node_Rates(m3ldbl *node_r, arbre *tree);
void RATES_Optimize_Node_Times_Serie_Fixed_Br_Len(node *a, node *d, arbre *tree);
void RATES_Optimize_Lexp(arbre *tree);
void RATES_Round_Optimize(arbre *tree);
void RATES_Optimize_Lexp(arbre *tree);
void RATES_Optimize_Alpha(arbre *tree);
m3ldbl RATES_Dmu(m3ldbl mu, int n_jumps, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp, int min_n, int jps_dens);
m3ldbl RATES_Dr_X_Dx(m3ldbl r, m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu_Given_Y_Trpzd(m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp, 
			       int nsteps, m3ldbl beg, m3ldbl end, m3ldbl prevs);
m3ldbl RATES_Dmu_Given_Y_Std(m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu_Given_Y_Romb(m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu_Given_Y(m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dy_Given_Mu(m3ldbl mu, m3ldbl y, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_Given_Y_X_Dy_Given_Mu1(m3ldbl mu1, m3ldbl mu2, m3ldbl y, m3ldbl dt1, m3ldbl dt2, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_Given_Mu1_Trpz(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl a, m3ldbl b, m3ldbl lexp,
				 int nsteps, m3ldbl beg, m3ldbl end, m3ldbl prevs);
m3ldbl RATES_Dmu2_And_Mu1(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_Given_Mu1(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_Given_Mu1_Romb(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl a, m3ldbl b, m3ldbl lexp);
void RATES_Random_Branch_Lengths(arbre *tree);
void RATES_Bracket_N_Jumps(int *up, int *down, m3ldbl param);
void RATES_Set_Node_Times(arbre *tree);
void RATES_Set_Node_Times_Pre(node *a, node *d, arbre *tree);
void RATES_Randomize_Node_Times(arbre *tree);
void RATES_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree);
void RATES_Optimize_Node_Times(arbre *tree);
m3ldbl RATES_Exp_Y(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl lexp);
m3ldbl RATES_Dmu2_Given_Mu1_Bis(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl alpha, m3ldbl beta, m3ldbl lexp);
void RATES_Replace_Br_Lengths_By_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
void RATES_Replace_Br_Lengths_By_Rates(arbre *tree);

void RATES_Get_Mean_Rates(arbre *tree);
void RATES_Get_Mean_Rates_Pre(node *a, node *d, edge *b, m3ldbl *r_a, m3ldbl alpha, m3ldbl lexp, arbre *tree);
void RATES_Expect_Number_Subst(m3ldbl t_beg, m3ldbl t_end, m3ldbl *r_beg, m3ldbl alpha, m3ldbl lexp, int *n_jumps, m3ldbl *mean_r);
void RATES_Optimize_Clock_Rate(arbre *tree);
m3ldbl RATES_Dmu1_Given_Lbda_And_Mu2(m3ldbl lbda, m3ldbl mu1, m3ldbl mu2, m3ldbl alpha, m3ldbl beta);
m3ldbl RATES_Dmu1_And_Mu2_One_Jump_Trpz(m3ldbl mu1, m3ldbl mu2, m3ldbl a, m3ldbl b,
					  int nsteps, m3ldbl beg, m3ldbl end, m3ldbl prevs);
m3ldbl RATES_Dmu1_And_Mu2_One_Jump_One_Interval(m3ldbl mu1, m3ldbl mu2, m3ldbl a, m3ldbl b);
m3ldbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(m3ldbl dt1, m3ldbl dt2, m3ldbl mu1, m3ldbl mu2, m3ldbl a, m3ldbl b);

m3ldbl RATES_Dmu1_And_Mu2_One_Jump_Old(m3ldbl mu1, m3ldbl mu2, m3ldbl a, m3ldbl b);

m3ldbl RATES_Dmu2_And_Min_N_Given_Mu1(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, int n_min, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_And_Mu1_Given_N(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, int n, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Lk_Rates_Core(m3ldbl mu1, m3ldbl mu2, int n1, int n2, m3ldbl dt1, m3ldbl dt2, arbre *tree);
void RATES_Init_Triplets(arbre *tree);
m3ldbl RATES_Lk_Change_One_Time(node *n, m3ldbl new_t, arbre *tree);
void RATES_Update_Triplet(node *n, arbre *tree);
void RATES_Print_Triplets(arbre *tree);
m3ldbl RATES_Lk_Change_One_Rate(node *d, m3ldbl new_rate, arbre *tree);
m3ldbl RATES_Dmu2_And_Mu1_Given_Min_N(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, int n_min, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu2_And_Mu1_Given_N_Normal(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, int n, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Coeff_Corr(m3ldbl alpha, m3ldbl beta, int n1, int n2);
m3ldbl RATES_Dmu2_And_Mu1_Given_N_Full(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, int n, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu1_Given_V_And_N(m3ldbl mu1, m3ldbl v, int n, m3ldbl dt1, m3ldbl a, m3ldbl b);
m3ldbl RATES_Yule(arbre *tree);
m3ldbl RATES_Check_Mean_Rates(arbre *tree);
void RATES_Check_Mean_Rates_Pre(node *a, node *d, edge *b, m3ldbl *sum, arbre *tree);
void RATES_Adjust_Clock_Rate(arbre *tree);
void RATES_Discretize_Rates(arbre *tree);
void RATES_Discretize_Rates_Pre(node *a, node *d, edge *b, arbre *tree);
m3ldbl RATES_Dmu_Given_V_And_MinN(m3ldbl mu, m3ldbl dt, m3ldbl v, int minn, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Dmu_One(m3ldbl mu, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp);
m3ldbl RATES_Compound_Core(m3ldbl mu1, m3ldbl mu2, int n1, int n2, m3ldbl dt1, m3ldbl dt2, m3ldbl alpha, m3ldbl beta, m3ldbl lexp, m3ldbl eps, int approx);
void RATES_Record_Rates(arbre *tree);
void RATES_Reset_Rates(arbre *tree);
void RATES_Record_Times(arbre *tree);
void RATES_Reset_Times(arbre *tree);
void RATES_Update_T_Rates_Pre(node *a, node *d, arbre *tree);
void RATES_Update_T_Rates(arbre *tree);
void RATES_Get_Br_Len(arbre *tree);
void RATES_Get_Rates_From_Bl(arbre *tree);
m3ldbl RATES_Compound_Core_Joint(m3ldbl mu1, m3ldbl mu2, int n1, int n2, m3ldbl dt1, m3ldbl dt2, 
				 m3ldbl alpha, m3ldbl beta, m3ldbl lexp, m3ldbl eps, int approx);
m3ldbl RATES_Dmu_Joint(m3ldbl mu, int n, m3ldbl dt, m3ldbl a, m3ldbl b, m3ldbl lexp, int min_n);
m3ldbl RATES_Compound_Core_Marginal(m3ldbl mu1, m3ldbl mu2, m3ldbl dt1, m3ldbl dt2, m3ldbl alpha, 
				    m3ldbl beta, m3ldbl lexp, m3ldbl eps, int approx);
m3ldbl RATES_Lk_Jumps(arbre *tree);
void RATES_Set_Rates_Prior_Mean(arbre *tree);
void RATES_Set_Rates_Prior_Mean_Pre(node *a, node *d, arbre *tree);
void RATES_Set_Rates_Post_Mean_And_Cov(arbre *tree);

#endif
