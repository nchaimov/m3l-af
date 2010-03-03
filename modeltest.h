/*
 *  modeltest.h
 *  
 *
 *  Created by David Elliott on 2/15/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef modeltest_h
#define modeltest_h

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "annealing.h"
#include "unittests.h"
#include "eb.h"

typedef struct _Node{
  int mod;
  struct _Node * left;
  struct _Node * right;
} Node;

void AIC(arbre* tree);
void testOpts(arbre* tree, double bestscore, int bestModel);
void HLRT(arbre* tree);
void printName(int mod);
int runTests(arbre* tree,Node* n, double previousLikelihood, int previousMod);
void destructTree(Node* n);
Node * constructTree();
void assignModel(arbre* tree, int model);
float likelihood(arbre* tree, int mod);
int getParams(int mod, arbre* tree);

#endif //modeltest_h
