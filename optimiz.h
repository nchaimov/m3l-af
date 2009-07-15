/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef OPTIMIZ_H
#define OPTIMIZ_H

void      Optimiz_Ext_Br(arbre *tree);
void      Optimize_Alpha(arbre *tree);
void      Optimize_Kappa(arbre *tree);
void      Optimize_Lambda(arbre *tree);
void      Optimize_Param_Parall(arbre *tree);
m3ldbl    Optimize_Branch_Quad(arbre *tree, allseq *alldata, edge *b_fcus);
void      Optimize_After_Hide(arbre *tree, allseq *alldata, node *h);
void      Round_Optimize(arbre *tree, allseq *data, int n_round_max);
int       Dist_Seq_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
			m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
			allseq *data, int num1, int num2, model *mod);
m3ldbl    Dist_Seq_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			 m3ldbl *xmin, allseq *data, 
			 int num1, int num2, model *mod);
m3ldbl    Kappa_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
		       m3ldbl *xmin, arbre *tree, allseq *alldata);
m3ldbl    Lambda_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			m3ldbl *xmin, arbre *tree, allseq *alldata);
m3ldbl    Alpha_Golden_Br_Opt(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			      m3ldbl *xmin, arbre *tree, allseq *alldata, 
			      int n_opt, m3ldbl *init_l);
m3ldbl    Alpha_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,m3ldbl *xmin, 
		       arbre *tree, allseq *alldata);
m3ldbl    Br_Len_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			m3ldbl *xmin, edge *b_fcus, arbre *tree);
m3ldbl    Br_Len_Brent(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, m3ldbl tol,
		       edge *b_fcus, arbre *tree, int n_iter_max, int quickdirty);
int       Br_Len_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		      m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
		      edge *b_fcus, arbre *tree);
m3ldbl    Optimize_Path_Length(model *mod, allseq *alldata, edge *a, 
			       int lra, edge *b, int lrb, m3ldbl i_len);
void      Optimize_Param_Serie(node *a, node *d, edge *b_fcus, arbre *tree, 
			       allseq *alldata, int n_passes);
m3ldbl    Optimize_Dist(model *mod, m3ldbl init, allseq *twoseqs);
m3ldbl    Pinvar_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			m3ldbl *xmin, arbre *tree, allseq *alldata, int n_iter_max);
void      Optimize_Pinvar(arbre *tree);
int       Lambda_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		      m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
		      arbre *tree);
int       Kappa_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		      m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
		      arbre *tree);
int       Alpha_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		      m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
		      arbre *tree);
int       Pinvar_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		      m3ldbl *fa, m3ldbl *fb, m3ldbl *fc, 
		      arbre *tree);
void Optimiz_All_Free_Param(arbre *tree, int verbose);
void      Optimiz_RRparam_GTR(arbre *tree, int num_param);
m3ldbl    RRparam_GTR_Golden(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
		   	     m3ldbl *xmin, arbre *tree, allseq *alldata, m3ldbl *param, int n_iter_max);

int Powell_GTR_Param(arbre *tree, m3ldbl *p, int n, m3ldbl ftol);
m3ldbl Linmin_GTR_Param(arbre *tree,m3ldbl *p, m3ldbl *xi, int n);
m3ldbl F1dim(arbre *tree, m3ldbl x, m3ldbl *p, m3ldbl *xi, m3ldbl n);
int Mnbrak_1dim(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		m3ldbl *fa, m3ldbl *fb, m3ldbl *fc,
		arbre *tree,
		m3ldbl *p,  m3ldbl *xi, m3ldbl n);
m3ldbl Brent_1dim(m3ldbl ax, m3ldbl bx, m3ldbl cx, 
		  m3ldbl tol, m3ldbl *xmin,
		  arbre *tree,
		  m3ldbl *p, m3ldbl *xi, m3ldbl n);

int Min_With_Derivatives(arbre *tree, m3ldbl *p, int n, m3ldbl ftol, m3ldbl step_size, 
			 m3ldbl (*func) (), void (*dfunc)(), m3ldbl (*linmin)());
void BFGS(arbre *tree, m3ldbl *p, int n, m3ldbl gtol, m3ldbl step_size,
	  m3ldbl(*func)(), void (*dfunc)(), void (*lnsrch)(),int *failed);
void Lnsrch_RR_Param(arbre *tree, int n, m3ldbl *xold, m3ldbl fold, m3ldbl *g, m3ldbl *p, m3ldbl *x,
		     m3ldbl *f, m3ldbl stpmax, int *check);
void Optimize_Single_Param_Generic(arbre *tree, m3ldbl *param, m3ldbl lim_inf, m3ldbl lim_sup, m3ldbl tol, int n_max_iter, int quickdirty);
int Generic_Brak(m3ldbl *param,
		 m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, 
		 m3ldbl *fa, m3ldbl *fb, m3ldbl *fc,
		 m3ldbl lim_inf, m3ldbl lim_sup,
		 arbre *tree);
m3ldbl Generic_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,
		     m3ldbl *xmin, arbre *tree, int n_iter_max,int quickdirty);
void Optimize_Br_Len_Serie(node *a, node *d, edge *b_fcus, arbre *tree,allseq *alldata);
void Lnsrch_Nucleotide_Frequencies(arbre *tree, int n, m3ldbl *xold, 
				   m3ldbl fold, m3ldbl *g, m3ldbl *p, m3ldbl *x,
				   m3ldbl *f, m3ldbl stpmax, int *check);

void Optimize_Global_Rate(arbre *tree);
m3ldbl Br_Len_Brent_Default(edge *b_fcus, arbre *tree);
m3ldbl Br_Len_Brent_Iter(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,
		    edge *b_fcus, arbre *tree, int n_iter_max, int quickdirty, int lnum);

void EM_Dist(model *mod, allseq *data);
m3ldbl Dist_F_Brent(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, m3ldbl tol, int n_iter_max,
		    m3ldbl *param, m3ldbl *F, model *mod, arbre *tree);
m3ldbl Dist_F_Brent_No_Bl(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, int n_iter_max,
		    m3ldbl *param, m3ldbl *F, model *mod);
int Dist_F_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, m3ldbl *F, m3ldbl *param, model *mod);
void Opt_Dist_F(m3ldbl *dist, m3ldbl *F, model *mod, arbre *tree);
void Opt_Dist_F_No_Bl(m3ldbl *dist, m3ldbl *F, model *mod);
m3ldbl Missing_Dist_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, int n_iter_max, 
			  int x, int y, matrix *mat);
int Missing_Dist_Brak(m3ldbl *ax, m3ldbl *bx, m3ldbl *cx, int x, int y, matrix *mat);
void Opt_Missing_Dist(int x, int y, matrix *mat);
int Optimiz_Alpha_And_Pinv(arbre *tree);
void Lnsrch_RR_Cov_Param(arbre *tree, int n, m3ldbl *xold, m3ldbl fold, 
			 m3ldbl *g, m3ldbl *p, m3ldbl *x,
			 m3ldbl *f, m3ldbl stpmax, int *check);
m3ldbl Node_Time_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,
		       node *anc, node *des, arbre *tree, int n_iter_max);
m3ldbl Time_Stamps_Mult_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,
			      arbre *tree, int n_iter_max);
m3ldbl Branch_Rate_Shape_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, 
			       m3ldbl *xmin, arbre *tree, int n_iter_max);
m3ldbl Node_Time_Brent_Fixed_Br_Len(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol,
				    node *n, arbre *tree, int n_iter_max);

m3ldbl Rates_Generic_Brent(m3ldbl ax, m3ldbl bx, m3ldbl cx, m3ldbl tol, m3ldbl *param, arbre *tree, int n_iter_max);

#endif

