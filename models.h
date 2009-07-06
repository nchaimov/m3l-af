/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef MODELS_H
#define MODELS_H

void  PMat(m3ldbl l, model *mod, int pos, double *Pij);
void  PMat_K80(m3ldbl l,m3ldbl kappa, int pos, double *Pij);
void  PMat_TN93(m3ldbl l, model *mod, int pos, double *Pij);
void  PMat_Empirical(m3ldbl l, model *mod, int pos, double *Pij);
void PMat_Zero_Br_Len(model  *mod, int pos, double *Pij);
void PMat_Gamma(m3ldbl l, model *mod, int pos, double *Pij);

int GetDaa (m3ldbl *daa, m3ldbl *pi, char *file_name);
int Matinv (double *x, int n, int m, double *space);
void Init_Model(allseq *data, model *mod);
void Update_Qmat_GTR(double *rr, m3ldbl *rr_val, int *rr_num, double *pi, double *qmat);
void Update_Qmat_HKY(double kappa, double *pi, double *qmat);
void Update_Qmat_Generic(double *rr, double *pi, int ns, double *qmat);
void Translate_Custom_Mod_String(model *mod);
void Set_Model_Parameters(model *mod);
m3ldbl GTR_Dist(m3ldbl *F, m3ldbl alpha, eigen *eigen_struct);
m3ldbl General_Dist(m3ldbl *F, model *mod, eigen *eigen_struct);

int Init_Qmat_WAG(double *daa, m3ldbl *pi);
int Init_Qmat_Dayhoff(double *daa, m3ldbl *pi);
int Init_Qmat_JTT(double *daa, m3ldbl *pi);
int Init_Qmat_RtREV(double *daa, m3ldbl *pi);
int Init_Qmat_CpREV(double *daa, m3ldbl *pi);
int Init_Qmat_VT(double *daa, m3ldbl *pi);
int Init_Qmat_Blosum62(double *daa, m3ldbl *pi);
int Init_Qmat_MtMam(double *daa, m3ldbl *pi);
int Init_Qmat_MtArt(double *daa, double *pi); // Added by Federico Abascal
int Init_Qmat_HIVb(double *daa, double *pi); //added by Federic Abascal
int Init_Qmat_HIVw(double *daa, double *pi); //added by Federico Abascal
void Switch_From_Mod_To_M4mod(model *mod);
void Switch_From_M4mod_To_Mod(model *mod);

#endif
