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

#include "utilities.h"

void AIC(arbre* tree);
float likelihood(arbre* tree, model* mod);

#endif //modeltest_h
