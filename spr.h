/*
** spr.h: Header file for the SPR routines.
**
** Wim Hordijk   Last modified: 28 August 2006
*/

#ifndef _SPR_H_
#define _SPR_H_

#include "utilities.h"

#define ALL   1
#define BEST  2
#define ONE   3

/*
** _move_: Structure for holding the relevant information for candidate SPR moves.
*/

//JSJ: will probably have to modify the move data structure if we want to use SPR in our analyses

typedef struct
{
  node   *v_prune, *u_prune, *v_n, *v_nx1, *u_n, **path;
  edge   *e_prune, *e_regraft;
  phydbl  l_connect, l_est[3], delta_lk, d_L, d_up_v, d_un_v;
  int     dist, rgrft_rank, optim_rank, globl_rank;
} _move_;



void Init_SPR          (arbre *tree);
void Clean_SPR         (arbre *tree);
void Optim_SPR         (arbre *tree, int max_size, int method);
int  Perform_SPR_Moves (arbre *tree, int max_size);
int  Perform_Best_SPR  (arbre *tree, int max_size);
int  Perform_One_SPR   (arbre *tree, int max_size);

void Calc_Tree_Length (edge *e_prune, node *v_prune, arbre *tree);
void Tree_Length      (node *v_prune, node *u_prune, node *v_n, node *v_n_1,
		       node *v_nx1, node *v_0, node *u_n, phydbl d_up_v_1,
		       phydbl d_uu, phydbl d_L_1, int n, arbre *tree);
int  Est_Lk_Change    (edge *e_prune, node *v_prune, arbre *tree);
int  Best_Lk_Change   (edge *e_prune, node *v_prune, arbre *tree);
void Make_Move        (_move_ *move, int type, arbre *tree);
int  Find_Optim_Local (arbre *tree);
int  Find_Optim_Globl (arbre *tree);
void Prune            (edge *e, node *v, edge **e_connect, edge **e_avail,
		       arbre *tree);
void Regraft          (edge *e, node *v, edge *avail, arbre *tree);
void PostOrder_v      (arbre *tree, node *v, edge *e);
void PostOrder_w      (arbre *tree, node *v, edge *v_e, node *w, edge *e);





void Speed_Spr(arbre *tree, int max_cycles);
void Speed_Spr_Loop(arbre *tree);
void Make_Spr_List(arbre *tree);
void Init_One_Spr(spr *a_spr);
spr *Make_One_Spr(arbre *tree);
int Spr(phydbl init_lnL, arbre *tree);
int Spr_Recur(node *a, node *d, arbre *tree);
int Test_All_Spr_Targets(edge *pulled, node *link, arbre *tree);
void Randomize_Spr_List(arbre *tree);
void Test_One_Spr_Target_Recur(node *a, node *d, edge *pulled, node *link, edge *residual, int *best_found, arbre *tree);
phydbl Test_One_Spr_Target(edge *target, edge *arrow, node *link, edge *residual, arbre *tree);
void Apply_Spr_Moves_One_By_One(arbre *tree);
int Try_One_Spr_Move_Triple(spr *move, arbre *tree);
int Try_One_Spr_Move_Full(spr *move, arbre *tree);
void Make_Best_Spr(arbre *tree);
void Random_Spr(int n_moves, arbre *tree);
void Include_One_Spr_To_List_Of_Spr(spr *move, arbre *tree);
void Reset_Spr_List(arbre *tree);
int Evaluate_List_Of_Regraft_Pos_Triple(spr **spr_list, int list_size, arbre *tree);
void Best_Spr(arbre *tree);
int Check_Spr_Move_Validity(spr *this_spr_move, arbre *tree);
void Spr_Subtree(edge *b, node *link, arbre *tree);
void Spr_Pars(arbre *tree);



#endif  /* _SPR_H_ */


/*
** EOF: spr.h
*/
