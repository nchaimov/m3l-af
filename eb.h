#ifndef __EB__
#define __EB__

#include "utilities.h"

void Empirical_Bayes(arbre* tree);
char *PostProb_From_String(char *s_tree, allseq *alldata, model *mod, option *io);

#endif

