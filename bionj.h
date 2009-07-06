/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef NJ_H
#define NJ_H

#include "utilities.h"
#include "optimiz.h"
/*#include "tools.h"*/

void   Bionj(matrix *mat);
void   Finish(matrix *mat);
void   Compute_Sx(matrix *mat);
m3ldbl Sum_S(matrix *mat, int i);
m3ldbl Dist(matrix *mat, int x, int y);
m3ldbl Q_Agglo(matrix *mat, int x, int y);
m3ldbl Variance(matrix *mat, int x, int y);
m3ldbl Br_Length(matrix *mat, int x, int y);
void   Update_Dist(matrix *mat, int x, int y);
m3ldbl Lamda(matrix *mat, int x, int y, m3ldbl vxy);
void   Best_Pair(matrix *mat, int *x, int *y, m3ldbl *score);
m3ldbl Var_Red(matrix *mat, int x, int y, int i, m3ldbl lamda, m3ldbl vxy);
void   Update_Tree(matrix *mat, int x, int y, m3ldbl lx, m3ldbl ly, m3ldbl score);
void   Update_Mat(matrix *mat, int x, int y, 
		  m3ldbl lx, m3ldbl ly, m3ldbl vxy, m3ldbl lamda);
m3ldbl Dist_Red(matrix *mat, int x, m3ldbl lx, int y, 
		m3ldbl ly, int i, m3ldbl lamda);
int    Bionj_Br_Length_Post(node *a, node *d, matrix *mat);
void   Bionj_Br_Length(matrix *mat);

#endif
