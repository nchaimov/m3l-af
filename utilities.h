/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.
*/

#ifndef UTILITIES_H
#define UTILITIES_H


//#define USE_OPENMP 1 // use openmp? if no, then remove this line.

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#define FALSE 0
#define TRUE 1

//#define USE_OPENMP 1 // comment/uncomment this line to use OpenMP parallelization
#ifdef USE_OPENMP
#include <omp.h>
#endif

//#define MEASURE 1 // for timing and speedup measurements

//#define COMPRESS_SUBALIGNMENTS 1 // this pragma enables compression of sub-alignments according to the phylogeny.
								 // although this optimization can increase runtime, it costs more in memory.

#define VERSION "v3.0 (179M)"

#define For(i,n)                     for(i=0; i<n; i++)
#define Fors(i,n,s)                  for(i=0; i<n; i+=s)
#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define SHFT2(a,b,c)                 (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d)               (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b)                     ((a)>(b)?(a):(b))
#define MIN(a,b)                     ((a)<(b)?(a):(b))
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);

#ifndef isnan
# define isnan(x)						 \
  (sizeof (x) == sizeof (long double) ? isnan_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isnan_d (x)		 \
   : isnan_f (x))
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x)						 \
  (sizeof (x) == sizeof (long double) ? isinf_ld (x)		 \
   : sizeof (x) == sizeof (double) ? isinf_d (x)		 \
   : isinf_f (x))
static inline int isinf_f  (float       x) { return isnan (x - x); }
static inline int isinf_d  (double      x) { return isnan (x - x); }
static inline int isinf_ld (long double x) { return isnan (x - x); }
#endif



#define  NNI_MOVE                    0
#define  SPR_MOVE                    1
#define  BEST_OF_NNI_AND_SPR         2
#define  SIMULATED_THERMAL_ANNEALING 3
#define  EMPIRICAL_BAYES             4
#define  SIMULATED_QUANTUM_ANNEALING 5


#define  PI      3.141593
#define  SQRT2PI 2.506628

#define  YES 1
#define  NO  0

#define  NT 0 /* nucleotides */
#define  AA 1 /* amino acids */

#define  ACGT 0 /* A,G,G,T encoding */
#define  RY   1 /* R,Y     encoding */

#define INTERFACE_DATA_TYPE      0
#define INTERFACE_MULTIGENE      1
#define INTERFACE_MODEL          2
#define INTERFACE_TOPO_SEARCH    3
#define INTERFACE_BRANCH_SUPPORT 4
#define INTERFACE_MBL_MODEL		 5

#ifndef INFINITY
#define INFINITY HUGE
#endif

#define  N_MAX_OPTIONS        100
#define  MIN_DT              0.01
#define  H_MCMC_RATES         0.5
#define  H_MCMC_LEXP          0.5
#define  H_MCMC_NU            1.0
#define  MAX_BL_SET            12

#define  T_MAX_FILE           500
#define  T_MAX_LINE       2000000
#define  T_MAX_NAME           500
#define  T_MAX_SEQ        2000000
#define  T_MAX_OPTION         100
#define  T_MAX_LABEL           10
#define  N_MAX_LABEL           10
#define  BLOCK_LABELS         100

#define  NODE_DEG_MAX          50
#define  BRENT_ITMAX       100000
#define  BRENT_CGOLD    0.3819660
#define  BRENT_ZEPS        1.e-10
#define  MNBRAK_GOLD     1.618034
#define  MNBRAK_GLIMIT      100.0
#define  MNBRAK_TINY       1.e-20
#define  ALPHA_MIN           0.04
#define  ALPHA_MAX            100
#define  BL_MIN            1.e-10
#define  BL_START          1.e-04
#define  BL_MAX             100.0
#define  GOLDEN_R      0.61803399
#define  GOLDEN_C  (1.0-GOLDEN_R)
#define  N_MAX_INSERT          20
#define  N_MAX_OTU           4000
#define  UNLIKELY          -1.e10
#define  NJ_SEUIL             0.1
#define  ROUND_MAX            100
#define  DIST_MAX            2.00
#define  AROUND_LK           50.0
#define  PROP_STEP            1.0
#define  T_MAX_ALPHABET        30
#define  MDBL_MIN   2.225074E-308
#define  MDBL_MAX   1.797693E+308
#define  POWELL_ITMAX         200
#define  LINMIN_TOL       2.0E-04
#define  LIM_SCALE              3
#define  LIM_SCALE_VAL     1.E-50
/* #define  LIM_SCALE           3000 */
/* #define  LIM_SCALE_VAL    1.E-500 */
#define  DEFAULT_SIZE_SPR_LIST  20
#define  OUTPUT_TREE_FORMAT  0 /* 0-->Newick; 1-->Nexus */
#define  MAX_PARS        1000000000
#define  MIN_DIFF_LK_LOCAL	 1.E-04 // VHS: changed from 1.E-03
#define  MIN_DIFF_LK_GLOBAL	 1.E-03
#define  MIN_DIFF_LK_MOVE    1.E-02

#define  MIN_CLOCK_RATE   1.E-10

#define JC69       1
#define K80        2
#define F81        3
#define HKY85      4
#define F84        5
#define TN93       6
#define GTR        7
#define CUSTOM     8

#define WAG       11
#define DAYHOFF   12
#define JTT       13
#define BLOSUM62  14
#define MTREV     15
#define RTREV     16
#define CPREV     17
#define DCMUT     18
#define VT        19
#define MTMAM     20
#define MTART     21
#define HIVW      22
#define HIVB      23
#define CUSTOMAA  24
#define LG        25

#define COMPOUND_COR   0
#define COMPOUND_NOCOR 1
#define EXPONENTIAL    2
#define GAMMA          3

#define DISPLAY_TIME_REMAINING_PERIOD 10 // Display the "estimated time remaining" counter every DISPLAY_TIME_REMAINING_PERIOD generations during EB MCMC

typedef	double m3ldbl;
typedef double plkflt;

/*********************************************************/

typedef struct __Node {
  struct __Node                       **v; /* table of pointers to neighbor nodes. Dimension = 3 */
  struct __Node               ***bip_node; /* three lists of pointer to tip nodes. One list for each direction */
  struct __Edge                       **b; /* table of pointers to neighbor branches */
  struct __Node ***list_of_reachable_tips; /* list of tip nodes that can be reached in each direction from that node */
  struct __Node                      *anc; /* direct ancestor node (for rooted tree only) */

  int                *n_of_reachable_tips; /* sizes of the list_of_reachable_tips (in each direction) */
  int                           *bip_size; /* Size of each of the three lists from bip_node */
  int                                 num; /* node number */
  int                                 tax; /* tax = 1 -> external node, else -> internal node */
  int                        check_branch; /* check_branch=1 is the corresponding branch is labelled with '*' */
  char                        ***bip_name; /* three lists of tip node names. One list for each direction */
  char                              *name; /* taxon name (if exists) */

  m3ldbl                           *score; /* score used in BioNJ to determine the best pair of nodes to agglomerate */
  m3ldbl                               *l; /* lengths of the (three or one) branch length sets connected this node */
  m3ldbl                     dist_to_root; /* distance to the root node */

  short int                        common;
}node;


/*********************************************************/

typedef struct __Edge {
  /*
    syntax :  (node) [edge]
(left_1) .                   .(right_1)
          \ (left)  (right) /
           \._____________./
           /    [b_fcus]   \
          /                 \
(left_2) .                   .(right_2)

  */

  struct __Node               *left,*rght; /* node on the left/right side of the edge */
  int         l_r,r_l,l_v1,l_v2,r_v1,r_v2;
  /* these are directions (i.e., 0, 1 or 2): */
  /* l_r (left to right) -> left[b_fcus->l_r] = right */
  /* r_l (right to left) -> right[b_fcus->r_l] = left */
  /* l_v1 (left node to first node != from right) -> left[b_fcus->l_v1] = left_1 */
  /* l_v2 (left node to secnd node != from right) -> left[b_fcus->l_v2] = left_2 */
  /* r_v1 (right node to first node != from left) -> right[b_fcus->r_v1] = right_1 */
  /* r_v2 (right node to secnd node != from left) -> right[b_fcus->r_v2] = right_2 */

  struct __NNI                       *nni;


  int                                 num; /* branch number */
  m3ldbl                               *l;//[MAX_BL_SET]; /* JSJ: branch lengths */
  m3ldbl                           *l_old;//[MAX_BL_SET]; /* JSJ: old branch lengths */

  int                           bip_score; /* score of the bipartition generated by the corresponding edge
					      bip_score = 1 iif the branch is found in both trees to be compared,
					      bip_score = 0 otherwise. */

  plkflt            *p_lk_left,*p_lk_rght; /* likelihoods of the subtree on the left and
					      right side (for each site and each relative rate category) */
  short int      *p_lk_tip_r, *p_lk_tip_l;
  short int           *div_post_pred_left; /* posterior prediction of nucleotide/aa diversity (left-hand subtree) */
  short int           *div_post_pred_rght; /* posterior prediction of nucleotide/aa diversity (rght-hand subtree) */

  double                          *Pij_rr; /* set of matrices of change probabilities and its first and second derivates
										    * VHS: usage = Pij_rr[gamma_cat*dimd + BL_set*dimaa + state_i*dima + state_j]
										    * = transition probability from state i to state j on this edge,
										    * where dimd = number of branch length sets times the number of states squared
										    * where dimaa = number of states squared
										    * where dima = number of states
										    */
  int                     *pars_l,*pars_r; /* parsimony of the subtree on the left and right sides (for each site) */
  unsigned int               *ui_l, *ui_r; /* union - intersection vectors used in Fitch's parsimony algorithm */
  int                *p_pars_l, *p_pars_r; /* conditional parsimony vectors */

  int                         num_st_left; /* number of the subtree on the left side */
  int                         num_st_rght; /* number of the subtree on the right side */


  /* Below are the likelihood scaling factors (used in functions
     `Get_All_Partial_Lk_Scale' in lk.c */
  int                          scale_left;
  int                          scale_rght;
  plkflt                *sum_scale_f_left;
  plkflt                *sum_scale_f_rght;

  m3ldbl                          bootval; /* bootstrap value (if exists) */

  short int                      is_alive; /* is_alive = 1 if this edge is used in a tree */

  m3ldbl                   dist_btw_edges;
  int                 topo_dist_btw_edges;

  int							      n_l; /* the number of branch length sets on this edge */
  int                    has_zero_br_len[MAX_BL_SET]; /* JSJ: branch length in given set is zero? */

#ifdef COMPRESS_SUBALIGNMENTS
  /*
   * Notes for red_left and red_right:
   * red_left and red_right will have one entry for each sequence site.
   * if red_left[i] == j, then the state pattern for site j and site i are the same down the left side
   * of this branch.
   * Similarly, if red_right[i] == j, then the state patterns for site j and i are the same down the
   * right side of this branch.
   * if red_right[i] == -1, then this site is unique, or is the primary pattern from which other sites
   * are identical.
   */
  int								*red_left; /* VHS: a list of redundant sites.  For examples if red[i] = j, then site j contains the same sub-alignment pattern as site i for this node.*/
  int							    *red_right;
 #endif

  m3ldbl                       ratio_test; /* approximate likelihood ratio test */
  m3ldbl                   alrt_statistic; /* aLRT statistic */
  int          num_tax_left, num_tax_rght; /* number of taxa in subtrees          */
  m3ldbl     avg_dist_left, avg_dist_rght; /* average taxon distance in subtrees  */

  m3ldbl						post_prob; /* posterior probability, see methods in eb.c and eb.h */

  int                       is_p_lk_l_u2d; /* is the conditional likelihood vector on the left up
					      to data ? */
  int                       is_p_lk_r_u2d; /* is the conditional likelihood vector on the right up
					      to data ? */

  char                           **labels; /* string of characters that labels the corresponding edge */
  int                            n_labels; /* number of labels */
  int                             n_jumps; /* number of jumps of substitution rates */
}edge;

/*********************************************************/

typedef struct __Arbre {
  struct __Node                       *n_root; /* root node */
  struct __Edge                       *e_root; /* edge on which lies the root */
  struct __Node                       **noeud; /* array of nodes that defines the tree topology */
  struct __Edge                     **t_edges; /* array of edges */
  struct __Arbre                    *old_tree; /* old copy of the tree */
  struct __Arbre                   *best_tree; /* best tree found so far */
  struct __Model                         *mod; /* substitution model */
  struct __AllSeq                       *data; /* sequences */
  struct __Option                         *io; /* input/output */
  struct __Matrix                        *mat; /* pairwise distance matrix */
  struct __Node                   **curr_path; /* list of nodes that form a path in the tree */
  struct __SPR                     **spr_list;
  struct __SPR                      *best_spr;
  struct __Tdraw                     *ps_tree; /* structure for drawing trees in postscript format */
  struct __Trate                       *rates; /* structure for handling rates of evolution */
  struct __Tmcmc                        *mcmc;
  //m3ldbl							   props[MAX_BL_SET]; /* JSJ: the proportion of sites in each branch length set */
  int                          ps_page_number; /* when multiple trees are printed, this variable give the current page number */
  int                         depth_curr_path; /* depth of the node path defined by curr_path */
  int                                 has_bip; /*if has_bip=1, then the structure to compare
						 tree topologies is allocated, has_bip=0 otherwise */
  int                                   n_otu; /* number of taxa */
  int                               curr_site; /* current site of the alignment to be processed */
  int                               curr_catg; /* current class of the discrete gamma rate distribution */
  int                                  n_swap; /* number of NNIs performed */
  int                               n_pattern; /* number of distinct site patterns */
  int                      has_branch_lengths; /* =1 iff input tree displays branch lengths */
  int                          print_boot_val; /* if print_boot_val=1, the bootstrap values are printed */
  int                          print_alrt_val; /* if print_alrt_val=1, the aLRT values are printed */
  int                            print_pp_val; /* if print_pp_val=1, the posterior probabilities are printed */
  int                              both_sides; /* both_sides=1 -> a pre-order and a post-order tree
												  traversal are required to compute the likelihood
												  of every subtree in the phylogeny*/
  int               num_curr_branch_available; /*gives the number of the next cell in t_edges that is free to receive a pointer to a branch */
  int                                 **t_dir;
  int                          n_improvements;
  int                                 n_moves;

  int                                      dp; /* Data partition */
  int                               s_mod_num; /* Substitution model number */
  int                      number_of_lk_calls;
  int               number_of_branch_lk_calls;
  int                               lock_topo; /* = 1 any subsequent topological modification will be banished */
  int                            print_labels;

  m3ldbl                              init_lnL;
  m3ldbl                              best_lnL; /* highest value of the loglikelihood found so far */
  int                                best_pars; /* highest value of the parsimony found so far */
  m3ldbl                                 c_lnL; /* loglikelihood */
  m3ldbl                         *c_lnL_sorted; /* used to compute c_lnL by adding sorted terms to minimize CPU errors */
  m3ldbl                              *site_lk; /* vector of likelihoods at individual sites */
  m3ldbl                     **log_site_lk_cat; /* loglikelihood at individual sites and for each class of rate*/
  m3ldbl                       unconstraint_lk; /* unconstrained (or multinomial) likelihood  */
  m3ldbl             prop_of_sites_to_consider;
  m3ldbl                        **log_lks_aLRT; /* used to compute several branch supports */
  m3ldbl                            n_root_pos; /* position of the root on its edge */
  m3ldbl                                  size; /* tree size */
  int                              *site_pars;
  int                                  c_pars;
  int                               *step_mat;

  int                           size_spr_list;
  int                  perform_spr_right_away;

  time_t                                t_beg;
  time_t                            t_current;

  struct __Triplet            *triplet_struct;

  int                     bl_from_node_stamps; /* == 1 -> Branch lengths are determined by node times */

  int					   red_arrays_invalid; /* == 1 -> the red arrays are dirty, == 0 -> the red arrays are up-to-date. */

}arbre;

/*********************************************************/
//JSJ: not dealing with this yet...
typedef struct __Super_Arbre {
  struct __Arbre                           *tree;
  struct __List_Arbre                  *treelist; /* list of trees. One tree for each data set to be processed */
  struct __AllSeq              *data_of_interest;
  struct __Option                   **optionlist; /* list of pointers to input structures (used in supertrees) */

  struct __Node           ***match_st_node_in_gt;
  /*  match_st_in_gt_node[subdataset number][supertree node number]
   *  gives the node in tree estimated from 'subdataset number' that corresponds
   *  to 'supertree node number' in the supertree
   */

  struct __Node           *****map_st_node_in_gt;
  /*  mat_st_gt_node[gt_num][st_node_num][direction] gives the
   *  node in gt gt_num that maps node st_node_num in st.
   */

  struct __Edge             ***map_st_edge_in_gt;
  /*  map_st_gt_br[gt_num][st_branch_num] gives the
   *  branch in gt gt_num that maps branch st_branch_num
   *  in st.
   */

  struct __Edge            ****map_gt_edge_in_st;
  /*  mat_gt_st_br[gt_num][gt_branch_num][] is the list of
   *  branches in st that map branch gt_branch_num
   *  in gt gt_num.
   */

  int                   **size_map_gt_edge_in_st;
  /*  size_map_gt_st_br[gt_num][gt_branch_num] gives the
   *  size of the list map_gt_st_br[gt_num][gt_branch_num][]
   */


  struct __Edge             ***match_st_edge_in_gt;
  /* match_st_edge_in_gt[gt_num][st_branch_num] gives the
   * branch in gt gt_num that matches branch st_branch_num
   */

  struct __Edge             ***match_gt_edge_in_st;
  /* match_gt_edge_in_st[gt_num][gt_branch_num] gives the
   * branch in st that matches branch gt_branch_num
   */

  struct __Node                  ****closest_match;
  /* closest_match[gt_num][st_node_num][dir] gives the
   * closest node in st that matches a node in gt gt_num
   */

  int                              ***closest_dist;
  /* closest_dist[gt_num][st_node_num][dir] gives the
   * number of edges to traverse to get to node
   * closest_match[gt_num][st_node_num][dir]
   */

  int                                         n_gt;
  /* number of trees */

  m3ldbl                                      **bl;
  /* bl[gt_num][gt_branch_num] gives the length of
   * branch gt_branch_num
   */

  m3ldbl                                  **bl_cpy;
  /* copy of bl */

  m3ldbl                                     **bl0;
  /* bl estimated during NNI (original topo)
   * See Mg_NNI.
   */

  m3ldbl                                     **bl1;
  /* bl estimated during NNI (topo conf 1)
   * See Mg_NNI.
   */

  m3ldbl                                     **bl2;
  /* bl estimated during NNI (topo conf 2)
   * See Mg_NNI.
   */

  int                                *bl_partition;
  /* partition[gt_num] gives the edge partition number
   * gt_num belongs to.
   */
  int                               n_bl_partition;

  struct __Model                          **s_mod; /* substitution model */

  int                                     n_s_mod;
  int                                 lock_br_len;

}superarbre;

/*********************************************************/

typedef struct __List_Arbre { /* a list of trees */
  struct __Arbre   **tree;
  int           list_size;                /* number of trees in the list */
}arbrelist;

/*********************************************************/

typedef struct __Seq {
  char           *name; /* sequence name */
  int              len; /* sequence length */
  char          *state; /* sequence itself */
  short int *is_ambigu; /* is_ambigu[site] = 1 if state[site] is an ambiguous character. 0 otherwise */
}seq;

/*********************************************************/


typedef struct __AllSeq {
  struct __Seq **c_seq;             /* compressed sequences      */
  m3ldbl        *b_frq;             /* observed state frequencies */
  short int     *invar;             /* 1 -> states are identical, 0 states vary */
  int            *wght;             /* # of each site in c_seq */
  short int    *ambigu;             /* ambigu[i]=1 is one or more of the sequences at site
				       i display an ambiguous character */
  m3ldbl    obs_pinvar;
  int            n_otu;             /* number of taxa */
  int        clean_len;             /* uncrunched sequences lengths without gaps */
  int       crunch_len;             /* crunched sequences lengths */
  int         init_len;             /* length of the uncompressed sequences */
  int        *sitepatt;             /* this array maps the position of the patterns in the
				       compressed alignment to the positions in the uncompressed
				       one */
}allseq;

/*********************************************************/

typedef struct __Matrix { /* mostly used in BIONJ */
  m3ldbl    **P,**Q,**dist; /* observed proportions of transition, transverion and  distances
			       between pairs of  sequences */

  arbre              *tree; /* tree... */
  int              *on_off; /* on_off[i]=1 if column/line i corresponds to a node that has not
			       been agglomerated yet */
  int                n_otu; /* number of taxa */
  char              **name; /* sequence names */
  int                    r; /* number of nodes that have not been agglomerated yet */
  struct __Node **tip_node; /* array of pointer to the leaves of the tree */
  int             curr_int; /* used in the NJ/BIONJ algorithms */
  int               method; /* if method=1->NJ method is used, BIONJ otherwise */
}matrix;

/*********************************************************/

typedef struct __Model {
  struct __Optimiz  *s_opt; /* pointer to parameters to optimize */
  struct __Eigen    *eigen;
  struct __M4       *m4mod;

  char          *modelname;
  char  *custom_mod_string; /* string of characters used to define custom models of substitution */
  int              *rr_num;
  int        *n_rr_per_cat; /* number of rate parameters in each category */
  int            n_diff_rr; /* number of different relative substitution rates in the custom model */
  int         update_eigen; /* update_eigen=1-> eigen values/vectors need to be updated */
  int           whichmodel;
  int             datatype; /* 0->DNA, 1->AA */
  int               n_catg; /* number of categories in the discrete gamma distribution */
  int 				   n_l; /* VHS: the number of branch length sets*/
  m3ldbl 		 *bl_props; /* branch length set proportions */
  int                invar; /* =1 iff the substitution model takes into account invariable sites */
  int                   ns; /* number of states (4 for ADN, 20 for AA) */
  int              seq_len; /* sequence length */
  int             stepsize; /* stepsize=1 for nucleotide models, 3 for codon models */
  int                n_otu; /* number of taxa */
  int            bootstrap; /* Number of bootstrap replicates (0 : no bootstrap analysis is launched) */
  int            use_m4mod; /* Use a Makrkov modulated Markov model ? */
  int         gamma_median; /* 1: use the median of each bin in the discrete gamma distribution. 0: the mean is used */


  m3ldbl               *pi; /* states frequencies */
  m3ldbl      *pi_unscaled; /* states frequencies (unscaled) */

  m3ldbl    *gamma_r_proba; /* probabilities of the substitution rates defined by the discrete gamma distribution */
  m3ldbl         *gamma_rr; /* substitution rates defined by the discrete gamma distribution */
  m3ldbl             kappa; /* transition/transversion rate */
  m3ldbl            lambda; /* parameter used to define the ts/tv ratios in the F84 and TN93 models */
  m3ldbl             alpha; /* gamma shapa parameter */
  m3ldbl            pinvar; /* proportion of invariable sites */
  m3ldbl         alpha_old;
  m3ldbl         kappa_old;
  m3ldbl        lambda_old;
  m3ldbl        pinvar_old;

  m3ldbl               *rr; /* relative rate parameters of the GTR or custom model (given by rr_val[rr_num[i]]) */
  m3ldbl           *rr_val; /* relative rate parameters of the GTR or custom model */
  double           *Pij_rr; /* matrix of change probabilities */
  m3ldbl                mr; /* mean rate = branch length/time interval  mr = -sum(i)(vct_pi[i].mat_Q[ii]) */
  m3ldbl      *user_b_freq; /* user-defined nucleotide frequencies */
  m3ldbl             *qmat;
  m3ldbl        *qmat_buff;

  m3ldbl        *rr_branch; /* relative substitution rates on each branch, for the whole set of sites */
  m3ldbl      *p_rr_branch; /* corresponding frequencies */
  int          n_rr_branch; /* number of classes */
  m3ldbl   rr_branch_alpha; /* Shape of the gamma distribution that defines the rr_branch and p_rr_branch values */
}model;


/*********************************************************/

typedef struct __Eigen{
  int              size;
  double             *q; /* matrix which eigen values and vectors are computed */
  double         *space;
  int        *space_int;
  double         *e_val; /* eigen values (vector), real part. */
  double      *e_val_im; /* eigen values (vector), imaginary part */
  double      *r_e_vect; /* right eigen vector (matrix), real part */
  double   *r_e_vect_im; /* right eigen vector (matrix), imaginary part */
  double      *l_e_vect; /* left eigen vector (matrix), real part */
}eigen;

/*********************************************************/

typedef struct __Option { /* mostly used in 'options.c' */
  struct __Model                *mod; /* pointer to a substitution model */
  struct __Arbre               *tree; /* pointer to the current tree */
  struct __Seq                **data; /* pointer to the uncompressed sequences */
  struct __AllSeq           *alldata; /* pointer to the compressed sequences */
  struct __Super_Arbre           *st; /* pointer to supertree */
  int                    interleaved; /* interleaved or sequential sequence file format ? */
  int                        in_tree; /* =1 if a user input tree is used as input */
  int				      user_props; /* JSJ: 0-> proportions are not user defined. */
  int 					   user_topo; /* JSJ: 0-> user has not defined a topo search method */
  int					 fixed_props;
  int 					  num_anneal_stages; /* JSJ: the temp count for use in Simulated Annealing */
  int 					   iters_per_stage; /* JSJ: the number of iterations per temp in Simulated Annealing*/
  double					temp_end; /* JSJ: the following, until the break, are options for simulated annealing*/
  double				  temp_start;
  double				   tau_start; // VHS: the starting value of the kinetic energy term for quantum annealing
  double					 tau_end; // VHS: the final value of the kinetic energy term for quantum annelaing
  double 				   acc_ratio; /* JSJ: the acceptance ratio, defaults to -1.0, used if set positive */
  int						set_back;
  double                   max_alpha;
  double                 brlen_sigma;
  double                pinvar_sigma;
  double                 gamma_sigma;
  double                  emig_sigma;
  double                    prob_NNI;
  double                    prob_SPR;
  double                  prob_brlen;
  double                  prob_gamma;
  double                  prob_kappa;
  double                 prob_lambda;
  double                     prob_rr;
  double                     prob_pi;
  double        prob_rate_proportion;
  double               prob_topology;
  double                 prob_pinvar;
  double 				   prob_emig;
  int                      eb_n_gens;


  char                  *in_seq_file; /* sequence file name */
  FILE                    *fp_in_seq; /* pointer to the sequence file */

  char                 *in_tree_file; /* input tree file name */
  FILE                   *fp_in_tree; /* pointer to the input tree file */

  char                *out_tree_file; /* name of the tree file */
  FILE                  *fp_out_tree;

  char               *out_trees_file; /* name of the tree file */
  FILE                 *fp_out_trees; /* pointer to the tree file containing all the trees estimated using random starting trees */

  char           *out_boot_tree_file; /* name of the tree file */
  FILE             *fp_out_boot_tree; /* pointer to the bootstrap tree file */

  char          *out_boot_stats_file; /* name of the tree file */
  FILE            *fp_out_boot_stats; /* pointer to the statistics file */

  char               *out_stats_file; /* name of the statistics file */
  FILE                 *fp_out_stats;

  char               *out_trace_file; /* name of the file in which the likelihood of the model is written */
  FILE                 *fp_out_trace;

  char                  *out_lk_file; /* name of the file in which the likelihood of the model is written */
  FILE                    *fp_out_lk;

  char                  *out_ps_file; /* name of the file in which tree(s) is(are) written */
  FILE                    *fp_out_ps;

  int               print_boot_trees; /* =1 if the bootstrapped trees are printed in output */
  int       out_stats_file_open_mode; /* opening file mode for statistics file */
  int        out_tree_file_open_mode; /* opening file mode for tree file */
  int                    n_data_sets; /* number of data sets to be analysed */
  int                        n_trees; /* number of trees */
  int                        seq_len; /* sequence length */
  int               n_data_set_asked; /* number of bootstrap replicates */
  char                     *nt_or_cd; /* nucleotide or codon data ? (not used) */
  int                      multigene; /* if=1 -> analyse several partitions. */
  int               config_multigene;
  int                           n_gt; /* number of gene trees */
  int                        curr_gt;
  int                     ratio_test; /* from 1 to 4 for specific branch supports, 0 of not */
  int                    ready_to_go;

  int					  post_probs; /* 1 = calculate PPs using the default *.eb file produce by EB MCMC, 2 = calculate PPs using a specified existing *.eb file, 0 = no PPs */
  char                   *in_eb_file; /* pathname of *.eb file, to be used to calculate PPs */
  FILE                     *fp_in_eb; /* pointer to the *.eb file */

  int                 curr_interface;
  int                         r_seed; /* random seed */
  int                  collapse_boot; /* 0 -> branch length on bootstrap trees are not collapsed if too small */
  int          random_boot_seq_order; /* !0 -> sequence order in bootstrapped data set is random */
  int                    print_trace;
  int                 print_site_lnl;
  int                       m4_model;
  int                      rm_ambigu; /* 0 is the default. 1: columns with ambiguous characters are discarded prior further analysis */
  int                   compress_seq;
  int                  append_run_ID;
  char                *run_id_string;
  int                          quiet; /* 0 is the default. 1: no interactive question (for batch mode) */

}option;

/*********************************************************/

typedef struct __Optimiz { /* parameters to be optimised (mostly used in 'optimiz.c') */
  int                 print; /* =1 -> verbose mode  */

  int             opt_alpha; /* =1 -> the gamma shape parameter is optimised */
  int             opt_kappa; /* =1 -> the ts/tv ratio parameter is optimised */
  int            opt_lambda; /* =1 -> the F84|TN93 model specific parameter is optimised */
  int            opt_pinvar; /* =1 -> the proportion of invariants is optimised */
  int        opt_state_freq; /* =1 -> the nucleotide frequencies are optimised */
  int                opt_rr; /* =1 -> the relative rate parameters of the GTR or the customn model are optimised */
  int         opt_num_param; /* if opt_topo=0 and opt_num_param=1 -> the numerical parameters of the
				model are optimised. if opt_topo=0 and opt_free_param=0 -> no parameter is
				optimised */
  int         opt_cov_delta;
  int         opt_cov_alpha;
  int    opt_cov_free_rates;


  int                opt_bl; /* =1 -> the branch lengths are optimized */
  int             opt_props; /* =1 -> the proportion of sites described by each member of
							          the branch length sets are optimized */
  int              opt_topo; /* =1 -> the tree topology is optimized */
  int           topo_search;
  m3ldbl            init_lk; /* initial loglikelihood value */
  int              n_it_max; /* maximum number of iterations during an optimization step */
  int              last_opt; /* =1 -> the numerical parameters are optimized further while the
			       tree topology remains fixed */
  int     random_input_tree; /* boolean */
  int         n_rand_starts; /* number of random starting points */
  int          brent_it_max;
  int             steph_spr;
  int       user_state_freq;
  int       opt_five_branch;
  int           pars_thresh;
  int         hybrid_thresh;

  m3ldbl     tree_size_mult; /* tree size multiplier */
  m3ldbl  min_diff_lk_local;
  m3ldbl min_diff_lk_global;
  m3ldbl   min_diff_lk_move;
  m3ldbl p_moves_to_examine;
  int              fast_nni;
  int                greedy;
  int          general_pars;
  int            quickdirty;
  int              spr_pars;
  int               spr_lnL;
  int        max_depth_path;
  int        min_depth_path;
  int          deepest_path;
  m3ldbl  max_delta_lnL_spr;



  int           wim_n_rgrft;
  int           wim_n_globl;
  int          wim_max_dist;
  int           wim_n_optim;
  int            wim_n_best;
  int        wim_inside_opt;



}optimiz;

/*********************************************************/
//JSJ: we will probably have to modify the __NNI and __SPR structs so that they take
// arrays of l
typedef struct __NNI{

  struct __Node         *left;
  struct __Node         *rght;
  struct __Edge            *b;

  int					  n_l; //JSJ: stores the number of branch length sets
  m3ldbl                score;
  m3ldbl              *init_l; //JSJ: array of initial branch lengths
  m3ldbl              init_lk;
  m3ldbl              *best_l; //JSJ: array of best branch lengths
  m3ldbl          lk0,lk1,lk2;
  m3ldbl          *l0,*l1,*l2; //JSJ: array of bls extending in three directions.

  struct __Node *swap_node_v1;
  struct __Node *swap_node_v2;
  struct __Node *swap_node_v3;
  struct __Node *swap_node_v4;

  int       best_conf;   /* best topological configuration :
			    ((left_1,left_2),right_1,right_2) or
			    ((left_1,right_2),right_1,left_2) or
			    ((left_1,right_1),right_1,left_2)  */
}nni;

/*********************************************************/
//JSJ: should probably modify this as well..
typedef struct __SPR{
  struct __Node         *n_link;
  struct __Node  *n_opp_to_link;
  struct __Edge  *b_opp_to_link;
  struct __Edge       *b_target;
  struct __Edge  *b_init_target;
  struct __Node          **path;
  m3ldbl         init_target_l[MAX_BL_SET]; //JSJ: an array of init_target_lengths
  m3ldbl            l0[MAX_BL_SET],l1[MAX_BL_SET],l2[MAX_BL_SET]; //JSJ: an array of edge direction lengths
  int						n_l; //JSJ: store the number of branch length sets
  m3ldbl                    lnL;
  int                depth_path;
  int                      pars;
  int                      dist;
}spr;

/*********************************************************/
//JSJ: this may require some examination
typedef struct __Triplet{
  int    size;
  m3ldbl *F_bc;
  m3ldbl *F_cd;
  m3ldbl *F_bd;
  m3ldbl ****core;
  m3ldbl ***p_one_site;
  m3ldbl ***sum_p_one_site;
  m3ldbl *pi_bc;
  m3ldbl *pi_cd;
  m3ldbl *pi_bd;
  struct __Eigen *eigen_struct;
  struct __Model *mod;
}triplet;

/*********************************************************/

typedef struct __Pnode{
  struct __Pnode **next;
  int weight;
  int num;
}pnode;

/*********************************************************/

typedef struct __M4 {
  int                  n_h; /* number of hidden states */
  int                  n_o; /* number of observable states  */
  int        use_cov_alpha;
  int         use_cov_free;

  m3ldbl          **o_mats; /* set of matrices of substitution rates across observable states */
  m3ldbl          *multipl; /* vector of values that multiply each o_mats matrix */
  m3ldbl             *o_rr; /* relative rates (symmetric) of substitution between observable states */
  m3ldbl             *h_rr; /* relative rates (symmetric) of substitution between hidden states */
  m3ldbl            *h_mat; /* matrix that describes the substitutions between hidden states (aka switches) */
  m3ldbl             *o_fq; /* equilibrium frequencies for the observable states */
  m3ldbl             *h_fq; /* equilibrium frequencies for the hidden states */
  m3ldbl    *h_fq_unscaled; /* unscaled equilibrium frequencies for the hidden states */
  m3ldbl *multipl_unscaled; /* unscaled  vector of values that multiply each o_mats matrix */

  m3ldbl             delta; /* switching rate */
  m3ldbl             alpha; /* gamma shape parameter */
}m4;

/*********************************************************/

typedef struct __Tdraw {
  int             *xcoord; /* node coordinates on the x axis */
  int             *ycoord; /* node coordinates on the y axis */
  int          page_width;
  int         page_height;
  int      tree_box_width;

  double max_dist_to_root;
}tdraw;

/*********************************************************/

typedef struct __Trate {
  m3ldbl clock_r; /* Mean substitution rate, i.e., 'molecular clock' rate */
  m3ldbl  *br_r; /* Relative substitution rate, i.e., multiplier of mean_r on each branch */
  m3ldbl  lexp; /* Parameter of the exponential distribution that governs the rate at which substitution between rate classes ocur */
  m3ldbl alpha;
  m3ldbl *true_t;
  int *n_jps;
  int *t_jps;
  m3ldbl *dens; /* Probability densities of mean substitution rates at the nodes */
  m3ldbl c_lnL; /* Prob(Br len | time stamps, model of rate evolution) */
  m3ldbl c_lnL_jps; /* Prob(# Jumps | time stamps, rates, model of rate evolution) */
  int adjust_rates; /* if = 1, branch rates are adjusted such that a modification of a given node time
		       does not modify any branch lengths */
  int use_rates; /* if = 0, branch lengths are expressed as differences between node times */
  m3ldbl *triplet;
  m3ldbl less_likely;
  m3ldbl birth_rate;
  m3ldbl min_rate;
  m3ldbl max_rate;
  m3ldbl step_rate;
  m3ldbl  *nd_r;  /* Current rates at nodes and the corresponding incoming edges */
  m3ldbl  *old_r; /* Old node rates */
  m3ldbl  *nd_t; /* Current node times */
  m3ldbl  *old_t; /* Old node times */


  int bl_from_rt; /* if =1, branch lengths are obtained as the product of cur_r and t */
  int approx;
  int model; /* Model number */
  m3ldbl nu; /* Parameter of the Exponential distribution for the corresponding model */

  m3ldbl *prior_r_mean;
  m3ldbl *prior_r_cov;
  m3ldbl *post_r_mean;
  m3ldbl *post_r_cov;

  m3ldbl autocor;
}trate;

/*********************************************************/

typedef struct __Tmcmc {
  int run;
  int sample_interval;
  int acc_lexp;
  int acc_rates;
  int acc_times;
  int acc_nu;
  int n_rate_jumps;

  m3ldbl *dt_prop;
  m3ldbl *p_no_jump;
  m3ldbl *t_rate_jumps;
  int    *t_rank;
  m3ldbl *r_path;
}tmcmc;

/*********************************************************/
/*********************************************************/

m3ldbl Rand_Normal_Deviate(m3ldbl mean, m3ldbl sd);
m3ldbl bico(int n,int k);
m3ldbl factln(int n);
m3ldbl gammln(m3ldbl xx);
m3ldbl Pbinom(int N,int ni,m3ldbl p);
void Plim_Binom(m3ldbl pH0,int N,m3ldbl *pinf,m3ldbl *psup);
m3ldbl LnGamma(m3ldbl alpha);
m3ldbl IncompleteGamma(m3ldbl x,m3ldbl alpha,m3ldbl ln_gamma_alpha);
m3ldbl PointChi2(m3ldbl prob,m3ldbl v);
m3ldbl PointNormal(m3ldbl prob);
int DiscreteGamma(m3ldbl freqK[],m3ldbl rK[],m3ldbl alfa,m3ldbl beta,int K,int median);
arbre *Read_Tree(char *s_tree);
void Restore_Tree_From_String(char *s_tree, arbre* tree);
void Make_All_Edges_Light(node *a,node *d);
void Make_All_Edges_Lk(node *a,node *d,arbre *tree);
void R_rtree(char *s_tree_a, char *s_tree_d, node *a, arbre *tree, int *n_int, int *n_ext);
void Clean_Multifurcation(char **subtrees,int current_deg,int end_deg);
char **Sub_Trees(char *tree,int *degree);
int Next_Par(char *s,int pos);
int Next_Brac(char *s, int pos);
char *Write_Tree(arbre *tree);
void R_wtree(node *pere,node *fils,char *s_tree,arbre *tree);
void Init_Tree(arbre *tree, int n_otu);
edge *Make_Edge_Light(node *a, node *d, int num, int n_l);
void Init_Edge_Light(edge *b, int num);
void Make_Edge_Dirs(edge *b,node *a,node *d);
void Make_Edge_Lk(edge *b, arbre *tree);

node *Make_Node_Light(int num, int n_l);
void Make_Node_Lk(node *n);
seq **Get_Seq(option *input);
seq **Read_Seq_Sequential(FILE *in,int *n_otu);
seq **Read_Seq_Interleaved(FILE *in,int *n_otu);
int Read_One_Line_Seq(seq ***data,int num_otu,FILE *in);
void Uppercase(char *ch);
allseq *Compact_Seq(seq **data,option *input);
allseq *Compact_CSeq(allseq *data,model *mod);
void Get_Base_Freqs(allseq *data);
void Get_AA_Freqs(allseq *data);
arbre *Read_Tree_File(FILE *fp_input_tree);
void Connect_Edges_To_Nodes(node *a,node *d,arbre *tree,int *cur);
void Exit(char *message);
void *mCalloc(int nb,size_t size);
void *mRealloc(void *p,int nb,size_t size);
/* arbre *Make_Light_Tree_Struct(int n_otu); */
int Sort_Phydbl_Decrease(const void *a, const void *b);
void Qksort(m3ldbl *A, m3ldbl *B, int ilo,int ihi);
void Print_P_Lk(plkflt *p_lk, int site, arbre *tree);
void Print_Site(allseq *alldata,int num,int n_otu,char *sep,int stepsize);
void Print_Seq(seq **data,int n_otu);
void Print_CSeq(FILE *fp,allseq *alldata);
void Print_Mat(matrix *mat);
void Print_Pij(double *Pij, model *mod);
void Print_Plk(double *plk, arbre *tree, int is_tax);
void Order_Tree_Seq(arbre *tree,seq **data);
void Order_Tree_CSeq(arbre *tree,allseq *data);
matrix *Make_Mat(int n_otu);
void Init_Mat(matrix *mat,allseq *data);
void Print_Dist(matrix *mat);
void Print_Node(node *a,node *d,arbre *tree);
void Share_Lk_Struct(arbre *t_full,arbre *t_empt);
void Init_Constant();

int Sort_Edges_NNI_Score(arbre *tree, edge **sorted_edges, int n_elem);
void NNI(arbre *tree, edge *b_fcus, int do_swap);
void Swap(node *a,node *b,node *c,node *d,arbre *tree);
void Update_All_Partial_Lk(edge *b_fcus,arbre *tree);
void Update_SubTree_Partial_Lk(edge *b_fcus,node *a,node *d,arbre *tree);
#ifdef MEASURE
void Count_Mean_Compressability(arbre *tree);
int Post_Order_Count_Compressability(node *a, node *d, arbre *tree, int site);
#endif MEASURE
allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names);
allseq *Copy_Cseq(allseq *ori, int len, int ns);
optimiz *Alloc_Optimiz();
int Filexists(char *filename);
FILE *Openfile(char *filename,int mode);
void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *input, int n_data_set, int num_rand_tree);
void Print_Fp_Out_Lines(FILE *fp_out,time_t t_beg,time_t t_end,arbre *tree,option *input,int n_data_set);
matrix *K80_dist(allseq *data,m3ldbl g_shape);
matrix *JC69_Dist(allseq *data,model *mod);
matrix *Hamming_Dist(allseq *data,model *mod);
int Is_Ambigu(char *state,int datatype,int stepsize);
void Check_Ambiguities(allseq *data,int datatype,int stepsize);
int Assign_State(char *c,int datatype,int stepsize);
void Bootstrap(arbre *tree);
void Br_Len_Involving_Invar(arbre *tree);
void Br_Len_Not_Involving_Invar(arbre *tree);
void Getstring_Stdin(char *file_name);
void Print_Freq(arbre *tree);
m3ldbl Num_Derivatives_One_Param(m3ldbl(*func)(arbre *tree),arbre *tree,m3ldbl f0,m3ldbl *param,m3ldbl stepsize,m3ldbl *err,int precise);
void Num_Derivative_Several_Param(arbre *tree,m3ldbl *param,int n_param,m3ldbl stepsize,m3ldbl(*func)(arbre *tree),m3ldbl *derivatives);
int Compare_Two_States(char *state1,char *state2,int state_size);
void Copy_One_State(char *from,char *to,int state_size);
model *Make_Model_Basic();
void Make_Model_Complete(model *mod);
model *Copy_Model(model *ori);
void Set_Defaults_Input(option *input);
void Set_Defaults_Model(model *mod);
void Set_Defaults_Optimiz(optimiz *s_opt);
void Copy_Optimiz(optimiz *ori,optimiz *cpy);
void Get_Bip(node *a,node *d,arbre *tree);
void Alloc_Bip(arbre *tree);
int Sort_Phydbl_Increase(const void *a,const void *b);
int Sort_String(const void *a,const void *b);
void Compare_Bip(arbre *tree1,arbre *tree2);
void Test_Multiple_Data_Set_Format(option *input);
int Are_Compatible(char *statea,char *stateb,int stepsize,int datatype);
void Hide_Ambiguities(allseq *data);
void Print_Site_Lk(arbre *tree, FILE *fp);
arbrelist *Make_Tree_List(int n_trees);
option *Make_Input();
arbre *Make_Tree();
void Make_All_Tree_Nodes(arbre *tree, int n_l);
void Make_All_Tree_Edges(arbre *tree, int n_l);
void Copy_Tax_Names_To_Tip_Labels(arbre *tree, allseq *data);
arbre *Make_And_Init_Tree(allseq *data);
void Connect_Edges_To_Nodes_Recur(node *a, node *d, arbre *tree);
void Connect_One_Edge_To_Two_Nodes(node *a, node *d, edge *b, arbre *tree);
arbre *Make_Tree_From_Scratch(int n_otu, allseq *data, int n_l);
arbrelist *Make_Treelist(int list_size);
void Put_Subtree_In_Dead_Objects(node *a, node *d, arbre *tree);
void Prune_Subtree(node *a, node *d, edge **target, edge **residual, arbre *tree);
void Reassign_Node_Nums(node *a, node *d, int *curr_ext_node, int *curr_int_node, arbre *tree);
void Copy_Tree(arbre *ori, arbre *cpy);
void Reassign_Edge_Nums(node *a, node *d, int *curr_br, arbre *tree);
void Init_Node_Light(node *n, int num);
void Make_All_Edge_Dirs(node *a, node *d, arbre *tree);
void Get_List_Of_Reachable_Tips(arbre *tree);
void Get_List_Of_Reachable_Tips_Post(node *a, node *d, arbre *tree);
void Get_List_Of_Reachable_Tips_Pre(node *a, node *d, arbre *tree);
void Make_List_Of_Reachable_Tips(arbre *tree);
void Graft_Subtree(edge *target, node *link, edge *add_br, arbre *tree);
int Get_Subtree_Size(node *a, node *d);
void Pull_Subtree_From_Dead_Objects(node *a, node *d, arbre *tree);
void Make_Edge_NNI(edge *b);
nni *Make_NNI(int n_l); //JSJ: so we can make it with appropriate arrays.
void Init_NNI(nni *a_nni, int n_l); //JSJ: to help initialize values
void Insert(arbre *tree);
void Connect_Two_Nodes(node *a, node *d);
void Get_List_Of_Target_Edges(node *a, node *d, edge **list, int *list_size, arbre *tree);
void Fix_All(arbre *tree);
void Record_Br_Len(m3ldbl **where, arbre *tree);
void Restore_Br_Len(m3ldbl **from, arbre *tree);
void Get_Dist_Btw_Edges(node *a, node *d, arbre *tree);
void Detect_Polytomies(edge *b, m3ldbl l_thresh, arbre *tree);
int Compare_List_Of_Reachable_Tips(node **list1, int size_list1, node **list2, int size_list2);
int Compare_List_Of_Reachable_Tips_version2(node **list1, int size_list1, node **list2, int size_list2);
void Find_Mutual_Direction(node *n1, node *n2, int *dir_n1_to_n2, int *dir_n2_to_n1);
void Fill_Dir_Table(arbre *tree);
void Get_List_Of_Nodes_In_Polytomy(node *a, node *d, node ***list, int *size_list);
void NNI_Polytomies(arbre *tree, edge *b_fcus, int do_swap);
void Check_Path(node *a, node *d, node *target, arbre *tree);
void Get_List_Of_Adjacent_Targets(node *a, node *d, node ***node_list, edge ***edge_list, int *list_size);
void Sort_List_Of_Adjacent_Targets(edge ***list, int list_size);
node *Common_Nodes_Btw_Two_Edges(edge *a, edge *b);
void Make_Site_Lk_Backup(arbre *tree);
int KH_Test(m3ldbl *site_lk_m1, m3ldbl *site_lk_M2, arbre *tree);
void Store_P_Lk(m3ldbl ****ori, m3ldbl ****cpy, arbre *tree);
m3ldbl Triple_Dist(node *a, arbre *tree, int approx);
void Make_Symmetric(m3ldbl **F, int n);
void Round_Down_Freq_Patt(m3ldbl **F, arbre *tree);
m3ldbl Get_Sum_Of_Cells(m3ldbl *F, arbre *tree);
void Divide_Cells(m3ldbl **F, m3ldbl div, arbre *tree);
void Divide_Mat_By_Vect(m3ldbl **F, m3ldbl *vect, int size);
void Multiply_Mat_By_Vect(m3ldbl **F, m3ldbl *vect, int size);
void Triple_Dist_Recur(node *a, node *d, arbre *tree);
triplet *Make_Triplet_Struct(model *mod);
void Fast_Br_Len(edge *b, arbre *tree, int approx);
void Fast_Br_Len_Recur(node *a, node *d, edge *b, arbre *tree);
void Print_Tree(FILE *fp, arbre *tree);
void Found_In_Subtree(node *a, node *d, node *target, int *match, arbre *tree);
void Random_Tree(arbre *tree);
void Copy_Dist(m3ldbl **cpy, m3ldbl **orig, int n);
eigen *Make_Eigen_Struct(model *mod);
void Random_NNI(int n_moves, arbre *tree);
void Make_Edge_Pars(edge *b, arbre *tree);
void Make_Tree_Path(arbre *tree);
void Share_Pars_Struct(arbre *t_full, arbre *t_empt);
void Share_Spr_Struct(arbre *t_full, arbre *t_empt);
void Share_List_Of_Reachable_Tips_Struct(arbre *t_full, arbre *t_empt);
void Clean_Tree_Connections(arbre *tree);
void Print_Settings(option *input);
void Fill_Missing_Dist(matrix *mat);
void Fill_Missing_Dist_XY(int x, int y, matrix *mat);
m3ldbl Least_Square_Missing_Dist_XY(int x, int y, m3ldbl dxy, matrix *mat);
void Update_Dirs(arbre *tree);
void Print_Banner(FILE *fp);
void Qksort_matrix(m3ldbl **A, int col, int ilo, int ihi);
void Check_Memory_Amount(arbre *tree);
int Get_State_From_P_Lk(m3ldbl *p_lk, int pos, arbre *tree);
int Get_State_From_P_Pars(short int *p_pars, int pos, arbre *tree);
void Unroot_Tree(char **subtrees);
void Print_Lk(arbre *tree, char *string);
void Print_Pars(arbre *tree);
void Print_Lk_And_Pars(arbre *tree);
void Check_Dirs(arbre *tree);
void Warn_And_Exit(char *s);
void Print_Data_Set_Number(option *input, FILE *fp);
m3ldbl Compare_Bip_On_Existing_Edges(m3ldbl thresh_len, arbre *tree1, arbre *tree2);
void NNI_Pars(arbre *tree, edge *b_fcus, int do_swap);
void Evaluate_One_Regraft_Pos_Triple(spr *move, arbre *tree);
int Get_State_From_Ui(int ui, int datatype);
void Read_Qmat(double *daa, m3ldbl *pi, FILE *fp);
void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, seq **data, option *input, pnode *n);
pnode *Create_Pnode(int size);
int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize);
void Randomize_Sequence_Order(allseq *data);
void Dist_To_Node_Pre(node *a, node *d, edge *b, arbre *tree);
void Add_Root(edge *target, arbre *tree);
int Is_Invar(int patt_num, int stepsize, int datatype, allseq *data);
void Update_Root_Pos(arbre *tree);
void Read_Branch_Label(char *sub_part, char *full_part, edge *b);
void Read_Branch_Lengths(char *s_d, char *s_a, arbre *tree);
void Read_Node_Name(node *d, char *s_tree_d, arbre *tree);
arbre *Generate_Random_Tree_From_Scratch(int n_otu, int rooted);
void Random_Lineage_Rates(node *a, node *d, edge *b, m3ldbl stick_prob, m3ldbl *rates, int curr_rate, int n_rates, arbre *tree);
edge *Find_Edge_With_Label(char *label, arbre *tree);
void Print_Square_Matrix_Generic(int n, m3ldbl *mat);
int Pick_State(int n, m3ldbl *prob);
char Reciproc_Assign_State(int i_state, int datatype);
void Evolve(allseq *data, model *mod, arbre *tree);
int Pick_State(int n, m3ldbl *prob);
void Evolve_Recur(node *a, node *d, edge *b, int a_state, int r_class, int site_num, allseq *gen_data, model *mod, arbre *tree);
void Site_Diversity(arbre *tree);
void Site_Diversity_Post(node *a, node *d, edge *b, arbre *tree);
void Site_Diversity_Pre(node *a, node *d, edge *b, arbre *tree);
void Subtree_Union(node *n, edge *b_fcus, arbre *tree);
void Binary_Decomposition(int value, int *bit_vect, int size);
void Print_Diversity(FILE *fp, arbre *tree);
void Print_Diversity_Header(FILE *fp, arbre *tree);
void Print_Diversity_Pre(node *a, node *d, edge *b, FILE *fp, arbre *tree);
void Make_New_Edge_Label(edge *b);
void Print_Qmat_AA(double *daa, m3ldbl *pi);
trate *Make_Rate_Struct(arbre *tree);
void Init_Rate_Struct(trate *rates, arbre *tree);
m3ldbl CDF_Normal(m3ldbl x, m3ldbl mean, m3ldbl var);
m3ldbl CDF_Gamma(m3ldbl x, m3ldbl mean, m3ldbl var);
double Uni();
double Ahrensdietergamma(double alpha);
double Rgamma(double shape, double scale);
double Rexp(double lambda);
m3ldbl Univariate_Kernel_Density_Estimate(m3ldbl where, m3ldbl *x, int n);
m3ldbl Var(m3ldbl *x, int n);
m3ldbl Mean(m3ldbl *x, int n);
m3ldbl Multivariate_Kernel_Density_Estimate(m3ldbl *where, m3ldbl **x, int sample_size, int vect_size);
m3ldbl Dgamma(m3ldbl x, m3ldbl shape, m3ldbl scale);
m3ldbl Dgamma_Moments(m3ldbl x, m3ldbl mean, m3ldbl var);
m3ldbl Dpois(m3ldbl x, m3ldbl param);
void Record_Model(model *ori, model *cpy);
void Best_Of_NNI_And_SPR(arbre *tree);
int Polint(m3ldbl *xa, m3ldbl *ya, int n, m3ldbl x, m3ldbl *y, m3ldbl *dy);
void Reset_Model_Parameters(model *mod);
void Print_Banner_Small(FILE *fp);
void JF(arbre *tree);
arbre *Dist_And_BioNJ(allseq *alldata, model *mod, option *io);
void Add_BioNJ_Branch_Lengths(arbre *tree, allseq *alldata, model *mod, option *io);
arbre *Read_User_Tree(allseq *alldata, model *mod, option *io);
void Print_Time_Info(time_t t_beg, time_t t_end);
void Prepare_Tree_For_Lk(arbre *tree);
char *Bootstrap_From_String(char *s_tree, allseq *alldata, model *mod, option *io);
char *aLRT_From_String(char *s_tree, allseq *alldata, model *mod, option *io);
int Rand_Int(int min, int max);
void PhyML_Printf(char *format, ...);
void PhyML_Fprintf(FILE *fp, char *format, ...);
m3ldbl CDF_Pois(m3ldbl x, m3ldbl param);
m3ldbl Dnorm_Moments(m3ldbl x, m3ldbl mean, m3ldbl var);
int Choose(int n, int k);
m3ldbl Dnorm(m3ldbl x, m3ldbl mean, m3ldbl sd);
m3ldbl LnFact(int n);
void Update_Ancestors(node *a, node *d, arbre *tree);
void Find_Common_Tips(arbre *tree1, arbre *tree2);
m3ldbl Get_Tree_Size(arbre *tree);
int Find_Bipartition(char **target_bip, int bip_size, arbre *tree);
int Find_Duplication_Node(char **tax_set, int n_tax, arbre *tree);
void Get_Rid_Of_Prefix(char delim, arbre *tree);
m3ldbl Bivariate_Normal_Density(m3ldbl x, m3ldbl y, m3ldbl mux, m3ldbl muy, m3ldbl sdx, m3ldbl sdy, m3ldbl rho);
void Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree);
void Dist_To_Root(node *n_root, arbre *tree);
int Is_Duplication_Node(node *n, char **tax_set, int n_tax, arbre *tree);
m3ldbl Dexp(m3ldbl x, m3ldbl param);
int Sort_Edges_Depth(arbre *tree, edge **sorted_edges, int n_elem);
char *Basename(char *path);
m3ldbl Rnorm(m3ldbl mean, m3ldbl sd);
m3ldbl *Rnorm_Multid(m3ldbl *mu, m3ldbl *cov, int dim);
m3ldbl *Matrix_Mult(m3ldbl *A, m3ldbl *B, int nra, int nca, int nrb, int ncb);
m3ldbl *Matrix_Transpose(m3ldbl *A, int dim);
void Print_time_remaining(time_t now, time_t start, int i, int total);
void Normalize_Props(model *mod);
void Update_Default_Props(option *io);
void Print_Tree_Screen(arbre *tree);

#endif


