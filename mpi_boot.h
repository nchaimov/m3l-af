/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "utilities.h"

#ifdef MPI

#include "mpi.h"
#define BootTreeTag 0
#define BootStatTag 1
int Global_numTask, Global_myRank;


void Bootstrap_MPI(arbre *tree);
void Print_Fp_Out_Lines_MPI(arbre *tree, option *io, int n_data_set, char *bootStr);

#endif  // MPI
