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

void assignModel(arbre* tree, int model);
void AIC(arbre* tree);
float likelihood(arbre* tree, int mod);

#endif //modeltest_h
