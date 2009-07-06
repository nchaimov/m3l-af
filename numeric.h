/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef NUMERIC_H
#define NUMERIC_H

m3ldbl *Covariance_Matrix(arbre *tree);
m3ldbl *Hessian(arbre *tree);
void Recurr_Hessian(node *a, node *b, int plus_minus, m3ldbl eps, m3ldbl *res, arbre *tree);
double stdnormal_inv(double p);
double Uni();
int    Rand_Int(int min, int max);
double Ahrensdietergamma(double alpha);
double Rgamma(double shape, double scale);
double Rexp(double lambda);
m3ldbl Bico(int n, int k);
m3ldbl Factln(int n);
m3ldbl Gammln(m3ldbl xx);
m3ldbl Pbinom(int N, int ni, m3ldbl p);
m3ldbl LnGamma (m3ldbl alpha);
m3ldbl IncompleteGamma(m3ldbl x, m3ldbl alpha, m3ldbl ln_gamma_alpha);
m3ldbl PointChi2 (m3ldbl prob, m3ldbl v);
m3ldbl Bivariate_Normal_Density(m3ldbl x, m3ldbl y, m3ldbl mux, m3ldbl muy, m3ldbl sdx, m3ldbl sdy, m3ldbl rho);
m3ldbl PointNormal (m3ldbl prob);
int    DiscreteGamma (m3ldbl freqK[], m3ldbl rK[],m3ldbl alfa, m3ldbl beta, int K, int median);
m3ldbl CDF_Normal(m3ldbl x, m3ldbl mean, m3ldbl var);
m3ldbl Dnorm_Moments(m3ldbl x, m3ldbl mean, m3ldbl var);
m3ldbl Dnorm(m3ldbl x, m3ldbl mean, m3ldbl sd);
m3ldbl CDF_Gamma(m3ldbl x, m3ldbl shape, m3ldbl scale);
m3ldbl Dgamma_Moments(m3ldbl x, m3ldbl mean, m3ldbl var);
m3ldbl Dgamma(m3ldbl x, m3ldbl shape, m3ldbl scale);
m3ldbl LnFact(int n);
int    Choose(int n, int k);
m3ldbl CDF_Pois(m3ldbl x, m3ldbl param);
m3ldbl Dexp(m3ldbl x, m3ldbl param);
m3ldbl Dpois(m3ldbl x, m3ldbl param);
m3ldbl Rand_Normal_Deviate(m3ldbl mean, m3ldbl sd);
m3ldbl Rnorm(m3ldbl mean, m3ldbl sd);
m3ldbl *Rnorm_Multid(m3ldbl *mu, m3ldbl *cov, int dim);
m3ldbl Rnorm_Trunc(m3ldbl mean, m3ldbl sd, m3ldbl min, m3ldbl max);
m3ldbl *Rnorm_Multid_Trunc(m3ldbl *mean, m3ldbl *cov, m3ldbl *min, m3ldbl *max, int dim);


#endif NUMERIC_H
