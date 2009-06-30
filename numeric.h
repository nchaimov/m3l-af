/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef NUMERIC_H
#define NUMERIC_H

phydbl *Covariance_Matrix(arbre *tree);
phydbl *Hessian(arbre *tree);
void Recurr_Hessian(node *a, node *b, int plus_minus, phydbl eps, phydbl *res, arbre *tree);
double stdnormal_inv(double p);
double Uni();
int    Rand_Int(int min, int max);
double Ahrensdietergamma(double alpha);
double Rgamma(double shape, double scale);
double Rexp(double lambda);
phydbl Bico(int n, int k);
phydbl Factln(int n);
phydbl Gammln(phydbl xx);
phydbl Pbinom(int N, int ni, phydbl p);
phydbl LnGamma (phydbl alpha);
phydbl IncompleteGamma(phydbl x, phydbl alpha, phydbl ln_gamma_alpha);
phydbl PointChi2 (phydbl prob, phydbl v);
phydbl Bivariate_Normal_Density(phydbl x, phydbl y, phydbl mux, phydbl muy, phydbl sdx, phydbl sdy, phydbl rho);
phydbl PointNormal (phydbl prob);
int    DiscreteGamma (phydbl freqK[], phydbl rK[],phydbl alfa, phydbl beta, int K, int median);
phydbl CDF_Normal(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dnorm(phydbl x, phydbl mean, phydbl sd);
phydbl CDF_Gamma(phydbl x, phydbl shape, phydbl scale);
phydbl Dgamma_Moments(phydbl x, phydbl mean, phydbl var);
phydbl Dgamma(phydbl x, phydbl shape, phydbl scale);
phydbl LnFact(int n);
int    Choose(int n, int k);
phydbl CDF_Pois(phydbl x, phydbl param);
phydbl Dexp(phydbl x, phydbl param);
phydbl Dpois(phydbl x, phydbl param);
phydbl Rand_Normal_Deviate(phydbl mean, phydbl sd);
phydbl Rnorm(phydbl mean, phydbl sd);
phydbl *Rnorm_Multid(phydbl *mu, phydbl *cov, int dim);
phydbl Rnorm_Trunc(phydbl mean, phydbl sd, phydbl min, phydbl max);
phydbl *Rnorm_Multid_Trunc(phydbl *mean, phydbl *cov, phydbl *min, phydbl *max, int dim);


#endif NUMERIC_H
