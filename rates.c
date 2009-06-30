/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */
#if defined (MC) ||  defined (RWRAPPER)

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
#include "mc.h"
#include "m4.h"
#include "draw.h"
#include "rates.h"
#include "numeric.h"

#ifdef RWRAPPER
#include <R.h>
#endif



/*********************************************************/

phydbl RATES_Lk_Rates(arbre *tree)
{
  
  RATES_Set_Node_Times(tree);
  RATES_Init_Triplets(tree);
/*   if(tree->rates->model == EXPONENTIAL) */
/*   RATES_Adjust_Clock_Rate(tree); */

  tree->rates->c_lnL = .0;
  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Lk_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

  return tree->rates->c_lnL;
}

/*********************************************************/

void RATES_Lk_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{
  int i,n1,n2;
  phydbl dens,mu1,mu2,dt1,dt2;
      

  if(d->anc != a)
    {
      printf("\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(a != tree->n_root)
    {
      dt1 = fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[a->anc->num]);
      mu1 = tree->rates->nd_r[a->num];
      n1  = tree->rates->n_jps[a->num];

      dt2 = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
      mu2 = tree->rates->nd_r[d->num];
      n2  = tree->rates->n_jps[d->num];

      dens = RATES_Lk_Rates_Core(mu1,mu2,n1,n2,dt1,dt2,tree);
    }
  else
    {
      dt2 = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);
      mu2 = tree->rates->nd_r[d->num];
      n2  = tree->rates->n_jps[d->num];

      dens = RATES_Dmu(mu2,n2,dt2,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
    }

  tree->rates->c_lnL += log(dens);

  tree->rates->triplet[a->num] += log(dens);

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Lk_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Rate(node *d, phydbl new_rate, arbre *tree)
{
  tree->rates->nd_r[d->num] = new_rate;
  RATES_Update_Triplet(d,tree);
  RATES_Update_Triplet(d->anc,tree);
  return(tree->rates->c_lnL);
}

/*********************************************************/

phydbl RATES_Lk_Change_One_Time(node *n, phydbl new_t, arbre *tree)
{  
  if(n == tree->n_root)
    {
      printf("\n. Moving the time of the root node is not permitted.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i;
      
      tree->rates->nd_t[n->num] = new_t;

      RATES_Update_Triplet(n,tree);
      
      For(i,3)
	{
	  if(n->b[i] != tree->e_root) RATES_Update_Triplet(n->v[i],tree);
	  else RATES_Update_Triplet(tree->n_root,tree);
	}
    }
  return(tree->rates->c_lnL);
}

/*********************************************************/

void RATES_Update_Triplet(node *n, arbre *tree)
{
  phydbl curr_triplet,new_triplet;
  phydbl dt0,dt1,dt2;
  phydbl mu1_mu0,mu2_mu0;
  phydbl mu0,mu1,mu2;
  int n0,n1,n2;
  int i;
  node *v1,*v2;

  if(n->tax) return;

  curr_triplet = tree->rates->triplet[n->num];

  dt0 = dt1 = dt2 = -100.0;

  if(n == tree->n_root)
    {
      phydbl dens;
      
      dens = 0.0;

      dt0 = tree->rates->nd_t[tree->n_root->v[0]->num] - tree->rates->nd_t[tree->n_root->num];
      dt1 = tree->rates->nd_t[tree->n_root->v[1]->num] - tree->rates->nd_t[tree->n_root->num];
      
      mu0 = tree->rates->nd_r[tree->n_root->v[0]->num];
      mu1 = tree->rates->nd_r[tree->n_root->v[1]->num];
      
      n0  = tree->rates->n_jps[tree->n_root->v[0]->num];
      n1  = tree->rates->n_jps[tree->n_root->v[1]->num];

      if((tree->rates->model == COMPOUND_COR) ||
	 (tree->rates->model == COMPOUND_NOCOR))
	{
	  dens  = RATES_Dmu(mu0,n0,dt0,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
	  dens *= RATES_Dmu(mu1,n1,dt1,tree->rates->alpha,1./tree->rates->alpha,tree->rates->lexp,0,1);
	}
      else
	{
/* 	  dens = Dexp(mu0,tree->rates->nu) * Dexp(mu1,tree->rates->nu); */
	  printf("\n TO DOOOO");
	}

      new_triplet = +1.;
      new_triplet = log(dens);
    }
  else
    {
      mu0 = mu1 = mu2 = -1.;
      n0 = n1 = n2 = -1;

      mu0 = tree->rates->nd_r[n->num];
      dt0 = fabs(tree->rates->nd_t[n->num] - tree->rates->nd_t[n->anc->num]);
      n0  = tree->rates->n_jps[n->num];

      v1 = v2 = NULL;
      For(i,3)
	{
	  if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	    {
	      if(!v1)
		{
		  v1  = n->v[i]; 
		  mu1 = tree->rates->nd_r[v1->num];
		  dt1 = fabs(tree->rates->nd_t[v1->num] - tree->rates->nd_t[n->num]);
		  n1  = tree->rates->n_jps[v1->num];
		}
	      else
		{
		  v2  = n->v[i]; 
		  mu2 = tree->rates->nd_r[v2->num];
		  dt2 = fabs(tree->rates->nd_t[v2->num] - tree->rates->nd_t[n->num]);
		  n2  = tree->rates->n_jps[v2->num];
		}
	    }
	}
 
      mu1_mu0 = RATES_Lk_Rates_Core(mu0,mu1,n0,n1,dt0,dt1,tree);
      mu2_mu0 = RATES_Lk_Rates_Core(mu0,mu2,n0,n2,dt0,dt2,tree);
      
      new_triplet = log(mu1_mu0) + log(mu2_mu0);
    }

  tree->rates->c_lnL = tree->rates->c_lnL + new_triplet - curr_triplet;
  tree->rates->triplet[n->num] = new_triplet;
}

/*********************************************************/
/* Returns f(mu2;mu1) */
phydbl RATES_Lk_Rates_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, arbre *tree)
{
  phydbl dens;
  phydbl alpha, beta, lexp;

  lexp = tree->rates->lexp;
  alpha = tree->rates->alpha;
  beta = 1./alpha;
  dens = UNLIKELY;

  if(mu1 < tree->rates->min_rate) mu1 = tree->rates->min_rate;
  if(mu1 > tree->rates->max_rate) mu1 = tree->rates->max_rate;
  
  if(mu2 < tree->rates->min_rate) mu2 = tree->rates->min_rate;
  if(mu2 > tree->rates->max_rate) mu2 = tree->rates->max_rate;
 
  if(dt1 < MIN_DT) dt1 = MIN_DT;
  if(dt2 < MIN_DT) dt2 = MIN_DT;
  
  switch(tree->rates->model)
    {
    case COMPOUND_COR:
      {       
	dens = RATES_Compound_Core(mu1,mu2,n1,n2,dt1,dt2,alpha,beta,lexp,tree->rates->step_rate,tree->rates->approx);
	break;
      }
      
    case COMPOUND_NOCOR :
      {
	dens = RATES_Dmu(mu2,n2,dt2,alpha,beta,lexp,0,1);
	break;
      }
      
    case EXPONENTIAL :
      {
	dens = Dexp(mu2,tree->rates->nu);
	break;
      }
      
    case GAMMA :
      {
	dens = Dgamma(mu2,tree->rates->alpha,1./tree->rates->alpha);
	break;
      }
      
    default : 
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
      }
    }
  return dens;
}

/*********************************************************/

phydbl RATES_Compound_Core(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  if((n1 > -1) && (n2 > -1))
    {
      return RATES_Compound_Core_Joint(mu1,mu2,n1,n2,dt1,dt2,alpha,beta,lexp,eps,approx);
    }
  else
    {
      return RATES_Compound_Core_Marginal(mu1,mu2,dt1,dt2,alpha,beta,lexp,eps,approx);
    }
}

/*********************************************************/

phydbl RATES_Compound_Core_Marginal(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  phydbl p0,p1,v0,v1,v2;
  phydbl dmu1;
  int    equ;
  phydbl dens;
  
  v0 = v1 = v2 = 0.0;
  
  /* Probability of 0 and 1 jumps */
  p0   = Dpois(0,lexp*(dt2+dt1));       
  p1   = Dpois(1,lexp*(dt2+dt1));
  
  dmu1 = RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,0,0);
  
  /* Are the two rates equal ? */
  equ = 0;
  if(fabs(mu1-mu2) < eps) equ = 1;
  
  /* No jump */
  if(equ)
    {
      v0 = 1.0*Dgamma(mu1,alpha,beta)/dmu1;
      /*       Rprintf("\n. mu1=%f mu2=%f",mu1,mu2); */
    }
  else
    {
      v0 = 1.E-100;
    }

  /* One jump */
  v1 = RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(mu1,mu2,dt1,dt2,alpha,beta);
  v1 /= dmu1;
    
  /* Two jumps and more (approximation) */
  if(approx == 1)
    {
      v2 =
	RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,0,0)*RATES_Dmu(mu2,-1,dt2,alpha,beta,lexp,0,0) -
	Dpois(0,lexp*dt1) * Dpois(0,lexp*dt2) *
	Dgamma(mu1,alpha,beta) * Dgamma(mu2,alpha,beta);
    }
  else
    {
      v2 = 
	RATES_Dmu_One(mu1,dt1,alpha,beta,lexp) * 
	RATES_Dmu_One(mu2,dt2,alpha,beta,lexp);

      v2 += Dpois(0,lexp*dt1)*Dgamma(mu1,alpha,beta)*RATES_Dmu(mu2,-1,dt2,alpha,beta,lexp,1,0);
      v2 += Dpois(0,lexp*dt2)*Dgamma(mu2,alpha,beta)*RATES_Dmu(mu1,-1,dt1,alpha,beta,lexp,1,0);

    }
/*   printf("\n. %f %f %f %f %f ",mu1,mu2,dt1,dt2,v2); */
  v2 /= dmu1;

  dens = p0*v0 + p1*v1 + v2;
/*   dens = p1*v1 + v2; */
  /*       dens = p1*v1 + v2; */
  /*   dens = v0; */
/*   dens *= dmu1; */
  
  if(dens < MDBL_MIN)
    {
      printf("\n. dens=%12G mu1=%12G mu2=%12G dt1=%12G dt2=%12G lexp=%12G alpha=%f v0=%f v1=%f v2=%f p0=%f p1=%f p2=%f",
	     dens,
	     mu1,mu2,dt1,dt2,
	     lexp,
	     alpha,
	     v0,v1,v2,
	     p0,p1,1.-p0-p1);
    }
  
  return dens;

}

/**********************************************************/

phydbl RATES_Compound_Core_Joint(phydbl mu1, phydbl mu2, int n1, int n2, phydbl dt1, phydbl dt2, 
				 phydbl alpha, phydbl beta, phydbl lexp, phydbl eps, int approx)
{
  phydbl density;
  phydbl dmu1;
  

  if(n1 < 0 || n2 < 0)
    {
      printf("\n. n1=%d n2=%d",n1,n2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dmu1 = RATES_Dmu(mu1,n1,dt1,alpha,beta,lexp,0,0);
  
  if((n1 < 0) || (n2 < 0))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if((n1 == 0) && (n2 == 0))
    {
      if(fabs(mu1-mu2) < eps) { density = Dgamma(mu1,alpha,beta); }
      else                    { density = 1.E-70; }
    }
  else if((n1 == 0) && (n2 == 1))
    {
      density = 
	Dgamma(mu1,alpha,beta) *
	RATES_Dmu1_Given_V_And_N(mu2,mu1,1,dt2,alpha,beta);
    }
  else if((n1 == 1) && (n2 == 0))
    {
      density = 
	Dgamma(mu2,alpha,beta) *
	RATES_Dmu1_Given_V_And_N(mu1,mu2,1,dt1,alpha,beta);
    }
  else /* independent */
    {
      density = 
	RATES_Dmu(mu1,n1,dt1,alpha,beta,lexp,0,0) * 
	RATES_Dmu(mu2,n2,dt2,alpha,beta,lexp,0,0);
    }

  density /= dmu1;

  density *= Dpois(n2,dt2*lexp);
  
  if(density < 1.E-70) density = 1.E-70;

/*   printf("\n. density = %15G mu1=%3.4f mu2=%3.4f dt1=%3.4f dt2=%3.4f n1=%2d n2=%2d",density,mu1,mu2,dt1,dt2,n1,n2); */
  return density;
}

/**********************************************************/

void RATES_Print_Triplets(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) printf("\n. Node %3d t=%f",i,tree->rates->triplet[i]);
}


/**********************************************************/

void RATES_Print_Rates(arbre *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

/*********************************************************/

void RATES_Print_Rates_Pre(node *a, node *d, edge *b, arbre *tree)
{

  printf("\n. a=%d d=%d rate=%f n_jps=%d t_left=%f t_rght=%f",
	 a->num,d->num,
	 tree->rates->nd_r[d->num],
	 tree->rates->n_jps[d->num],
	 tree->rates->nd_t[a->num],tree->rates->nd_t[d->num]);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Print_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Check_Mean_Rates(arbre *tree)
{
  phydbl sum;
  int i;
  
  sum = 0.0;
  For(i,2*tree->n_otu-2) sum += tree->rates->nd_r[i];
  return(sum/(phydbl)(2*tree->n_otu-2));
}

/*********************************************************/


trate *RATES_Make_Rate_Struct(int n_otu)
{
  trate *rates;
  
  rates               = (trate *)mCalloc(1,sizeof(trate));
  rates->br_r         = (phydbl *)mCalloc(2*n_otu-3,sizeof(phydbl));
  rates->old_r        = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->nd_r         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->old_t        = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->nd_t         = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->true_t       = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->dens         = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->triplet      = (phydbl *)mCalloc(2*n_otu-1,sizeof(phydbl));
  rates->n_jps        = (int *)mCalloc(2*n_otu-2,sizeof(int));
  rates->t_jps        = (int *)mCalloc(2*n_otu-2,sizeof(int));
  rates->prior_r_mean = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->prior_r_cov  = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));
  rates->post_r_mean  = (phydbl *)mCalloc(2*n_otu-2,sizeof(phydbl));
  rates->post_r_cov   = (phydbl *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(phydbl));

  return rates;
}

/*********************************************************/

void RATES_Init_Rate_Struct(trate *rates, int n_otu)
{
  int i;

  rates->model         = COMPOUND_COR;
  rates->clock_r       = 1.E-3;
  rates->c_lnL         = -INFINITY;
  rates->c_lnL_jps     = -INFINITY;
  rates->adjust_rates  = 0;
  rates->use_rates     = 1;
  rates->lexp          = 0.05;
  rates->alpha         = 2.0;
  rates->birth_rate    = 0.001;
  rates->max_rate      = 100.;
  rates->min_rate      = 0.01;
  rates->step_rate     = 1.E-4;
  rates->nu            = 1.0;
  rates->approx        = 1;
  rates->bl_from_rt    = 0;

  For(i,2*n_otu-2) 
    {
      rates->old_r[i]  = 1.0;
      rates->nd_r[i]   = 1.0;
      rates->br_r[i]   = 1.0;
      rates->n_jps[i]  = -1;
      rates->t_jps[i]  = -1;
    }

  For(i,2*n_otu-1) 
    {
      rates->nd_t[i]   = 0.0;
      rates->true_t[i] = 0.0;
    }
}

/*********************************************************/

void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param)
{
  phydbl cdf,eps,a,b,c;
  int step;

  step = 10;
  eps = 1.E-10;
  cdf = 0.0;
  c = 1;
  
  while(cdf < 1.-eps)
    {
      c = (int)floor(c * step);
      cdf = CDF_Pois(c,param);      
    }
  
  a = 0.0;
  b = (c-a)/2.;
  step = 0;
  do
    {
      step++;
      cdf = CDF_Pois(b,param);
      if(cdf < eps) a = b;
      else 
	{
	  break;
	}
      b = (c-a)/2.;
    }
  while(step < 1000);
  
  if(step == 1000)
    {
      PhyML_Printf("\n. a=%f b=%f c=%f param=%f",a,b,c,param);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  *up = c;
  *down = a;
}

/*********************************************************/
/* 
   mu   : average rate of the time period dt
   dt   : time period to be considered
   a    : rate at a given time point is gamma distributed. a is the shape parameter
   b    : rate at a given time point is gamma distributed. b is the scale parameter
   lexp : the number of rate switches is Poisson distributed with parameter lexp * dt
*/ 
/* compute f(mu;dt,a,b,lexp), the probability density of mu. We need to integrate over the
   possible number of jumps (n) during the time interval dt */
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens)
{
  if(n_jumps < 0) /* Marginal, i.e., the number of jumps is not fixed */
    {
      phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt,*suminv,b2;
      int n,up,down;
      
      var          = 0.0;
      cumpoissprob = 0.0;
      dens         = 0.0;
      n            = 0;
      mean         = a*b;
      ab2          = a*b*b;
      lexpdt       = lexp*dt;  
      b2           = b*b;
      suminv       = NULL;
      
      RATES_Bracket_N_Jumps(&up,&down,lexpdt);
      For(n,MAX(down,min_n)-1) cumpoissprob += Dpois(n,lexpdt);
      
      for(n=MAX(down,min_n);n<up+1;n++)
	{
	  poissprob    = Dpois(n,lexpdt); /* probability of having n jumps */      
	  var          = (2./(n+2.))*ab2; /* var(mu|n) = var(mu|n=0) * 2 / (n+2) */
	  gammadens    = Dgamma_Moments(mu,mean,var);
	  dens         += poissprob * gammadens;
	  cumpoissprob += poissprob;
	  if(cumpoissprob > 1.-1.E-04) break;
	}
      
      if(dens < 1.E-70) dens = 1.E-70;

      return(dens);      
    }
  else /* Joint, i.e., return P(mu | dt, n_jumps) */
    {
      phydbl mean, var, density;


      mean = 1.0;
      var = (2./(n_jumps+2.))*a*b*b;

      if(jps_dens)
	density = Dgamma_Moments(mu,mean,var) * Dpois(n_jumps,dt*lexp);
      else
	density = Dgamma_Moments(mu,mean,var);
      
      if(density < 1.E-70) density = 1.E-70;

      return density;
    }
}
  
/*********************************************************/
  
phydbl RATES_Dmu_One(phydbl mu, phydbl dt, phydbl a, phydbl b, phydbl lexp)
{
  phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt,*suminv,b2;
  int n,up,down;
  
  var          = 0.0;
  cumpoissprob = 0.0;
  dens         = 0.0;
  n            = 0;
  mean         = a*b;
  ab2          = a*b*b;
  lexpdt       = lexp*dt;  
  b2           = b*b;
  suminv       = NULL;
  
  if(dt < 0.0)
    {
      PhyML_Printf("\n. dt=%f",dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  if(lexpdt < MDBL_MIN)
    {
      PhyML_Printf("\n. lexpdt=%G lexp=%G dt=%G",lexpdt,lexp,dt);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(mu < 1.E-10)
    {
      PhyML_Printf("\n. mu=%G",mu);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");      
    }
  
  RATES_Bracket_N_Jumps(&up,&down,lexpdt);

  For(n,MAX(1,down)-1) cumpoissprob += Dpois(n,lexpdt);

  for(n=MAX(1,down);n<up+1;n++) /* WARNING: we are considering that at least one jump occurs in the interval */
    {
      poissprob    = Dpois(n,lexpdt); /* probability of having n jumps */      
      var          = (n/((n+1)*(n+1)*(n+2)))*(pow(1-a*b,2) + 2/(n+1)*ab2) + 2*n*n*ab2/pow(n+1,3);      
      gammadens    = Dgamma_Moments(mu,mean,var);
      dens         += poissprob * gammadens;
      cumpoissprob += poissprob;
      if(cumpoissprob > 1.-1.E-06) break;
    }

  return(dens);
}

/*********************************************************/

/* Given the times of nodes a (ta) and d (td), the shape of the gamma distribution of instantaneous
   rates, the parameter of the exponential distribution of waiting times between rate jumps and the 
   instantaneous rate at node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl *r_beg, phydbl alpha, phydbl lexp, int *n_jumps, phydbl *mean_r)
{
  phydbl curr_r, curr_t, next_t;

  curr_r = *r_beg;
  *mean_r = *r_beg;
      
  curr_t = t_beg + Rexp(lexp); /* Exponentially distributed waiting times */
  next_t = curr_t;
  
  *n_jumps = 0;
  while(curr_t < t_end)
    {
      curr_r = Rgamma(alpha,1./alpha); /* Gamma distributed random instantaneous rate */
      
      (*n_jumps)++;
      
      next_t = curr_t + Rexp(lexp);
	  
      if(next_t < t_end)
	{
	  *mean_r = (1./(next_t - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (next_t - curr_t));
	}
      else
	{
	  *mean_r = (1./(t_end - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (t_end - curr_t));
	}
      curr_t = next_t;
    }

/*   printf("\n. [%3d %f %f]",*n_jumps,*mean_r,*r_beg); */

  *r_beg = curr_r;  
}

/*********************************************************/


void RATES_Get_Mean_Rates_Pre(node *a, node *d, edge *b, phydbl *r_a, phydbl alpha, phydbl lexp, arbre *tree)
{
  phydbl a_t,d_t;
  phydbl r_d;
  phydbl mean_r;
  int n_jumps;

  a_t = tree->rates->nd_t[a->num];
  d_t = tree->rates->nd_t[d->num];
      
  RATES_Expect_Number_Subst(a_t,d_t,r_a,alpha,lexp,&n_jumps,&mean_r);
  
  tree->rates->nd_r[d->num] = mean_r;
/*   tree->rates->n_jps[d->num] = n_jumps; */
  tree->rates->t_jps[d->num] = n_jumps;

/*   printf("\n. a=%3d d=%3d mu=%15f nj=%3d r_a=%f", */
/* 	 a->num,d->num,mean_r,n_jumps,*r_a); */
  

/*   /\*******\/ */
/*   printf("\n. No correlation"); */
/*   *r_a = Rgamma(alpha,1./alpha); /\* Gamma distributed random instantaneous rate *\/ */
/*   /\*******\/ */

  r_d = *r_a;

  /* Move to the next branches */
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      *r_a = r_d;
	      RATES_Get_Mean_Rates_Pre(d,d->v[i],d->b[i],r_a,alpha,lexp,tree);
	    }
	}
    }
}

/*********************************************************/

void RATES_Random_Branch_Lengths(arbre *tree)
{
  phydbl r0,r_root,alpha,lexp;
   
  alpha = tree->rates->alpha;
  lexp  = tree->rates->lexp;
  
  printf("\n. alpha = %f lexp = %f",alpha,lexp);

  r0 = Rgamma(alpha,1./alpha);
  r_root = r0;
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[0],NULL,&r_root,alpha,lexp,tree);
  r_root = r0;
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,&r_root,alpha,lexp,tree);

  RATES_Adjust_Clock_Rate(tree);
  RATES_Get_Br_Len(tree);

/*   RATES_Print_Rates(tree); */

}

/*********************************************************/

void RATES_Randomize_Node_Times(arbre *tree)
{
  RATES_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Randomize_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Set_Node_Times(arbre *tree)
{
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[0],tree);
  RATES_Set_Node_Times_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void RATES_Init_Triplets(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->triplet[i] = 0.0;
}
/*********************************************************/

void RATES_Set_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf;
      int i;

      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_sup = tree->rates->nd_t[a->num];

      if(t_sup > t_inf)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
            
      if(tree->rates->nd_t[d->num] > t_inf)      tree->rates->nd_t[d->num] = t_inf;
      else if(tree->rates->nd_t[d->num] < t_sup) tree->rates->nd_t[d->num] = t_sup;

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Set_Node_Times_Pre(d,d->v[i],tree);	      
	    }
	}
    }
}

/*********************************************************/

void RATES_Randomize_Node_Times_Pre(node *a, node *d, arbre *tree)
{
  if(d->tax) return;
  else
    {
      node *v1, *v2; /* the two sons of d */
      phydbl t_sup, t_inf, u;
      int i;

      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
	  
      t_inf = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_sup = tree->rates->nd_t[a->num];

      if(t_sup > t_inf)
	{
	  printf("\n. t_sup = %f t_inf = %f",t_sup,t_inf);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      u = rand();
      u /= RAND_MAX;

      tree->rates->nd_t[d->num] = t_sup + u * fabs(t_sup-t_inf);

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Randomize_Node_Times_Pre(d,d->v[i],tree);	      
	    }
	}
    }
}

/*********************************************************/

phydbl RATES_Dmu1_And_Mu2_One_Jump_Two_Intervals(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, phydbl alpha, phydbl beta)
{
  phydbl dens;

  if(mu2 < 1.E-10)
    {
      printf("\n. mu2=%G",mu2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(mu1 < 1.E-10)
    {
      printf("\n. mu2=%G",mu1);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dens =
    ((dt1/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu1,mu2,1,dt1,alpha,beta) * Dgamma(mu2,alpha,beta)) +
    ((dt2/(dt1+dt2)) * RATES_Dmu1_Given_V_And_N(mu2,mu1,1,dt2,alpha,beta) * Dgamma(mu1,alpha,beta));

  return dens;
}

/*********************************************************/

phydbl RATES_Dmu1_Given_V_And_N(phydbl mu1, phydbl v, int n, phydbl dt1, phydbl a, phydbl b)
{
  phydbl lbda,dens,h,u;
  phydbl mean,var;
  int n_points,i;
  phydbl ndb;
  phydbl end, beg;

  n_points = 20;

  end = MIN(mu1/v-0.01,0.99);
  beg = 0.01;
  
  dens = 0.0;
  
  if(end > beg)
    {
      mean = a*b;
      var = a*b*b*2./(n+1.);
      ndb = (phydbl)n/dt1;
      
      h = (end - beg) / (phydbl)n_points;
      
      lbda = beg;
      For(i,n_points-1) 
	{
	  lbda += h;
	  u = (mu1 - lbda*v)/(1.-lbda);
	  
	  if(u < 1.E-10)
	    {
	      printf("\n. u = %G",u);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  
	  dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
	}
      dens *= 2.;
      
      lbda = beg;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u < 1.E-10)
	{
	  printf("\n. mu1 = %f lambda = %f v = %f u = %G beg = %f end = %f",mu1,lbda,v,u,beg,end);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
      
      lbda = end;
      u = (mu1 - lbda*v)/(1.-lbda);
      if(u < 1.E-10)
	{
	  printf("\n. u = %G",u);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      dens += Dgamma_Moments(u,mean,var) / (1.-lbda) * ndb * pow(1.-lbda,n-1);
      
      dens *= (h/2.);
      dens *= dt1;
    }

  return(dens);
}

/*********************************************************/
/* Joint density of mu1 and a minimum number of jumps occuring in the interval dt1+dt2 given mu1. 
   1 jump occurs at the junction of the two intervals, which makes mu1 and mu2 independant */
phydbl RATES_Dmu2_And_Mu1_Given_Min_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n_min, phydbl a, phydbl b, phydbl lexp)
{
  phydbl density, lexpdt,cumpoiss,poiss;
  int i;
  int up,down;

  density = 0.0;
  lexpdt = lexp * (dt1+dt2);
  cumpoiss = 0.0;

  RATES_Bracket_N_Jumps(&up,&down,lexpdt);

  For(i,MAX(up,n_min)-1)
    {
      poiss = Dpois(i,lexpdt);
      cumpoiss = cumpoiss + poiss;
    }

  for(i=MAX(up,n_min);i<up;i++)
    {
/*       poiss = Dpois(i-1,lexpdt); /\* Complies with the no correlation model *\/ */
      poiss = Dpois(i,lexpdt);
      cumpoiss = cumpoiss + poiss;

      density = density + poiss * RATES_Dmu2_And_Mu1_Given_N(mu1,mu2,dt1,dt2,i-1,a,b,lexp);
/*       density = density + poiss * RATES_Dmu2_And_Mu1_Given_N_Full(mu1,mu2,dt1,dt2,i,a,b,lexp); */
      if(cumpoiss > 1.-1.E-6) break;
    }

  if(density < 0.0)
    {
      printf("\n. density=%f cmpoiss = %f i=%d n_min=%d mu1=%f mu2=%f dt1=%f dt2=%f",
	     density,cumpoiss,i,n_min,mu1,mu2,dt1,dt2);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return(density);
}

/*********************************************************/
/* Joint density of mu1 and mu2 given the number of jumps (n) in the interval dt1+dt2, which are considered
   as independant. Hence, for n jumps occuring in dt1+dt2, the number of jumps occuring in dt1 and dt2 are
   n and 0, or n-1 and 1, or n-2 and 2 ... or 0 and n. This function sums over all these scenarios, with
   weights corresponding to the probability of each partitition.
 */

phydbl RATES_Dmu2_And_Mu1_Given_N(phydbl mu1, phydbl mu2, phydbl dt1, phydbl dt2, int n, phydbl a, phydbl b, phydbl lexp)
  {
    phydbl density,lexpdt1,lexpdt2,texpn,cumpoiss,poiss,abb,ab,lognf,logdt1,logdt2,nlogdt;
    int i;

    density  = 0.0;
    lexpdt1  = lexp*dt1;
    lexpdt2  = lexp*dt2;
    texpn    = pow(dt1+dt2,n);
    abb      = a*b*b;
    ab       = a*b;
    cumpoiss = 0.0;
    poiss    = 0.0;
    lognf    = LnFact(n);
    logdt1   = log(dt1);
    logdt2   = log(dt2);
    nlogdt   = n*log(dt1+dt2);

    For(i,n+1)
      {
        poiss = lognf - LnFact(i) - LnFact(n-i) + i*logdt1 + (n-i)*logdt2 - nlogdt;
	poiss = exp(poiss);
	cumpoiss = cumpoiss + poiss;

	if(mu2 < 1.E-10)
	  {
	    printf("\n. mu2=%f",mu2);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }

	if(mu1 < 1.E-10)
	  {
	    printf("\n. mu1=%f",mu1);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }


	density = density + poiss * Dgamma_Moments(mu2,ab,2./((n-i)+2.)*abb) * Dgamma_Moments(mu1,ab,2./(i+2.)*abb);
        if(cumpoiss > 1.-1.E-6) break;
      }

    return(density);
  }

/*********************************************************/

/* Logarithm of Formula (9) in Rannala and Yang (1996) */
phydbl RATES_Yule(arbre *tree)
{
  phydbl sumti,density,lambda;
  int n,i;

  sumti = 0.0;
  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) sumti += tree->rates->nd_t[i];
  sumti -= tree->rates->nd_t[i];

  lambda = tree->rates->birth_rate;
  n = tree->n_otu;
  
  density = 
    (n-1.)*log(2.) + 
    (n-2.)*log(lambda) - 
    lambda*sumti - 
    Factln(n) - 
    log(n-1.) - 
    (n-2.)*log(1.-exp(-lambda));
  
  /* Ad-hoc */
  if(tree->rates->nd_t[tree->n_root->num] > -100.) density = -INFINITY;

  return density;
}

/*********************************************************/

/* Set the clock rate such that the relative rates are centered on 1.0 */
void RATES_Adjust_Clock_Rate(arbre *tree)
{
  int i;
  phydbl mean;

  mean = 0.0;
  For(i,2*tree->n_otu-2) mean += tree->rates->nd_r[i];
  mean /= (phydbl)(2*tree->n_otu-2);
  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] /= mean;
  tree->rates->clock_r *= mean;
}

/*********************************************************/

void RATES_Record_Rates(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->old_r[i] = tree->rates->nd_r[i];
}

/*********************************************************/

void RATES_Reset_Rates(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_r[i] = tree->rates->old_r[i];
}

/*********************************************************/

void RATES_Record_Times(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->old_t[i] = tree->rates->nd_t[i];
}

/*********************************************************/

void RATES_Reset_Times(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] = tree->rates->old_t[i];
}

/*********************************************************/

void RATES_Get_Br_Len(arbre *tree)
{
  phydbl dt,rr,cr;
  node *left, *rght;
  int i;

  dt = rr = -1.0;
  cr = tree->rates->clock_r;


  if(tree->n_root)
    {
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num]);
      rr = tree->rates->nd_r[tree->n_root->v[0]->num];
      tree->e_root->l = dt*rr*cr;
      tree->n_root->l[0] = dt*rr*cr;
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]);
      rr = tree->rates->nd_r[tree->n_root->v[1]->num];
      tree->e_root->l += dt*rr*cr;
      tree->n_root->l[1] = dt*rr*cr;
    }
  
  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  left = tree->t_edges[i]->left;
	  rght = tree->t_edges[i]->rght;
	  dt = fabs(tree->rates->nd_t[left->num] - tree->rates->nd_t[rght->num]);	  
	  rr = (left->anc == rght)?(tree->rates->nd_r[left->num]):(tree->rates->nd_r[rght->num]);
	  tree->t_edges[i]->l = dt*rr*cr;
	}
    }
}

/*********************************************************/

void RATES_Get_Rates_From_Bl(arbre *tree)
{
  phydbl dt,cr;
  node *left, *rght;
  int i;

  dt = -1.0;
  cr = tree->rates->clock_r;

  if(tree->n_root)
    {
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[0]->num]);
      tree->rates->nd_r[tree->n_root->v[0]->num] = 0.5 * tree->e_root->l / (dt*cr);
      dt = fabs(tree->rates->nd_t[tree->n_root->num] - tree->rates->nd_t[tree->n_root->v[1]->num]);
      tree->rates->nd_r[tree->n_root->v[1]->num] = 0.5 * tree->e_root->l / (dt*cr);
    }
  

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i] != tree->e_root)
	{
	  left = tree->t_edges[i]->left;
	  rght = tree->t_edges[i]->rght;
	  dt = fabs(tree->rates->nd_t[left->num] - tree->rates->nd_t[rght->num]);	  
	  
	  if(left->anc == rght) tree->rates->nd_r[left->num] = tree->t_edges[i]->l / (dt*cr);
	  else                  tree->rates->nd_r[rght->num] = tree->t_edges[i]->l / (dt*cr);
	}
    }

  RATES_Adjust_Clock_Rate(tree);
}

/*********************************************************/

phydbl RATES_Lk_Jumps(arbre *tree)
{
  int i,n_jps;
  phydbl dens,dt,lexp;
  node *n;

  n = NULL;
  lexp = tree->rates->lexp;
  n_jps = 0;
  dt = 0.0;
  dens = 0.0;

  For(i,2*tree->n_otu-2)
    {
      n = tree->noeud[i];
      dt = fabs(tree->rates->nd_t[n->num]-tree->rates->nd_t[n->anc->num]);
      n_jps = tree->rates->n_jps[n->num];
      dens += log(Dpois(n_jps,lexp*dt));
    }

  tree->rates->c_lnL_jps = dens;
  
  return dens;
}

/*********************************************************/

/* void RATES_Set_Prior_Mean_Rates(arbre *tree) */
/* { */
/*   tree->rates->nd_r[tree->n_root->num] = 1.0; */
/*   RATES_Set_Rates_Prior_Mean_Pre(tree->n_root,tree->n_root->v[0],tree); */
/*   RATES_Set_Rates_Prior_Mean_Pre(tree->n_root,tree->n_root->v[1],tree); */
/* } */

/* /\*********************************************************\/ */

/* void RATES_Set_Prior_Cov_Rates(arbre *tree) */
/* { */
/*   RATES_Set_Rates_Prior_Cov_Pre(tree->n_root,tree->n_root->v[0],tree); */
/*   RATES_Set_Rates_Prior_Cov_Pre(tree->n_root,tree->n_root->v[1],tree); */
/* } */

/* /\*********************************************************\/ */

/* void RATES_Set_Prior_Mean_Rates_Pre(node *a, node *d, arbre *tree) */
/* { */
/*   tree->rates->prior_mean_rates[d->num] = tree->rates->nd_r[a->num] * tree->rates->clock_r; */
  
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       int i; */
/*       For(i,3) if(d->v[i] != a) RATES_Set_Prior_Mean_Rates_Pre(d,d->v[i],tree);      */
/*     } */
/* } */

/* /\*********************************************************\/ */


/* void RATES_Set_Prior_Cov_Rates_Pre(node *a, node *d, arbre *tree) */
/* { */
/*   tree->rates->prior_cov_rates[d->num*(2*tree->n_otu-2)+d->num] =  */
/*     fabs(tree->rates->nd_t[a->num] - tree->rates->nd_t[d->num]) * tree->rates->autocor; */
  
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       int i; */
/*       For(i,3) if(d->v[i] != a) RATES_Set_Prior_Cov_Rates_Pre(d,d->v[i],tree); */
/*     } */
/* } */

/* /\*********************************************************\/ */

/* void RATES_Set_Post_Mean_And_Cov_Bl(arbre *tree) */
/* { */
/*   int i,j; */
/*   int dim; */
/*   phydbl *inv_prior_cov, *inv_like_cov; */
/*   phydbl *post_mean, *post_cov; */
/*   phydbl *A, *B, *buff; */

/*   dim = 2*tree->n_otu-3; */

/*   inv_prior_cov = (phydbl *)mCalloc(dim*dim,sizeof(phydbl)); */
/*   buff          = (phydbl *)mCalloc(dim*dim,sizeof(phydbl)); */

/*   inv_like_cov  = tree->hessian; */
/*   post_cov      = tree->rates->post_cov_bl; */
/*   post_mean     = tree->rates->post_mean_bl; */

/*   For(i,dim*dim) inv_prior_cov[i] = tree->rates->prior_r_cov[i]; */

/*   Matinv(inv_prior_r_cov,dim,dim,buff); */
  
/*   /\* sum of the inverse of the covariance matrix for the likelihood (i.e., the hessian) and */
/*      the covariance matrix for the prior *\/ */
/*   For(i,dim) For(j,dim) post_cov[i*dim+j] = inv_like_cov[i*dim+j] + inv_prior_cov[i*dim+j]; */


/*   A = Matrix_Mult( */






  
  
/*   Matinvc(tree->rates->post_r_cov,dim,dim,buff); */
  
  
  
  
  
  



  
/*   Free(inv_prior_r_cov); */
/*   Free(A); */
/*   Free(B); */
/* } */

/* /\*********************************************************\/ */



#endif
