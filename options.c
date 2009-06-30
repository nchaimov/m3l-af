/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include <stdlib.h>
#include "cl.h"
#include "utilities.h"
#include "options.h"
#include "models.h"
#include "free.h"
#include "interface.h"
#ifdef MG
#include "mg.h"
#endif

/* int  T_MAX_FILE; */
/* phydbl MDBL_MIN; */
/* phydbl UNLIKELY; */

/*********************************************************/

void Usage()
{

  char *BOLD=(char *)mCalloc(10,sizeof(char));
  char *FLAT=(char *)mCalloc(10,sizeof(char));
  char *LINE=(char *)mCalloc(10,sizeof(char));
  char *cha;


  cha =getenv("OS");

  if(cha!=NULL)
    {
      strcpy(BOLD, "");
      strcpy(FLAT, "");
      strcpy(LINE, "");
    }
  else
    {
      strcpy(BOLD, "\033[00;01m");
      strcpy(FLAT, "\033[00;00m");
      strcpy(LINE, "\033[00;04m");
    }

  PhyML_Printf("%sNAME\n"
	 "%s\t- PhyML %s - \n\n"
	 "%s\t\''A simple, fast, and accurate algorithm to estimate\n"
	 "%s\tlarge phylogenies by maximum likelihood\''\n\n"
	 "%s\tStephane Guindon and Olivier Gascuel,\n"
	 "%s\tSystematic Biology 52(5):696-704, 2003.\n\n"
	 "%s\tPlease cite this paper if you use this software in your publications.\n",BOLD,FLAT,VERSION,FLAT,FLAT,FLAT,FLAT,FLAT);

  PhyML_Printf("%s\nSYNOPSIS:\n\n"
	 "%s\tphyml %s[command args]\n",BOLD,BOLD,BOLD);
  PhyML_Printf("%s\n\tAll the options below are optional (except '%s-i%s' if you want to use the command-line interface).\n\n",FLAT,BOLD,FLAT);

  PhyML_Printf("%s\nCommand options:\n%s",BOLD,FLAT);

  PhyML_Printf("\n\t%s-i (or --input) %sseq_file_name%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\t%sseq_file_name%s is the name of the nucleotide or amino-acid sequence file in PHYLIP format.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-d (or --datatype) ""%sdata_type%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("%s\t\tdata_type%s is 'nt' for nucleotide (default) and 'aa' for amino-acid sequences.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-q (or --sequential)\n",BOLD);
  PhyML_Printf("%s\t\tChanges interleaved format (default) to sequential format.\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-n (or --multiple) ""%snb_data_sets%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("%s\t\tnb_data_sets%s is an integer corresponding to the number of data sets to analyse.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-p (or --pars)%s\n",BOLD,FLAT);
  PhyML_Printf("%s\t\tUse a minimum parsimony starting tree. This option is taken into account when the '-u' option\n",FLAT);
  PhyML_Printf("%s\t\tis absent and when tree topology modifications are to be done.\n",FLAT);

  PhyML_Printf("%s\n\t-b (or --bootstrap) %sint%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\t%sint%s >  0 : %sint%s is the number of bootstrap replicates.\n",LINE,FLAT,LINE,FLAT);
  PhyML_Printf("\t\t%sint%s =  0 : neither approximate likelihood ratio test nor bootstrap values are computed.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sint%s = -1 : approximate likelihood ratio test returning aLRT statistics.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sint%s = -2 : approximate likelihood ratio test returning Chi2-based parametric branch supports.\n",LINE,FLAT);
/*   PhyML_Printf("\t\t%sint%s = -3 : minimum of Chi2-based parametric and SH-like branch supports.\n",LINE,FLAT); */
  PhyML_Printf("\t\t%sint%s = -4 : SH-like branch supports alone.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-m (or --model) %smodel%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tmodel%s : substitution model name.\n",FLAT);
  PhyML_Printf("\t\t%s- %sNucleotide%s-based models : %sHKY85%s (default) | %sJC69%s | %sK80%s | %sF81%s | %sF84%s | %sTN93%s | %sGTR%s | %scustom (*)%s\n",
	 FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT,LINE,FLAT);
  PhyML_Printf("\t\t(*) : for the custom option, a string of six digits identifies the model. For instance, 000000\n");
  PhyML_Printf("\t\t corresponds to F81 (or JC69 provided the distribution of nucleotide frequencies is uniform).\n");
  PhyML_Printf("\t\t 012345 corresponds to GTR. This option can be used for encoding any model that is a nested within GTR.\n");
  PhyML_Printf("\n");
  PhyML_Printf("\t\t%s- %sAmino-acid%s based models : %sLG%s (default) | %sWAG%s | %sJTT%s | %sMtREV%s | %sDayhoff%s | %sDCMut%s | %sRtREV%s | %sCpREV%s | %sVT%s\n",
	 FLAT,LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT);
  PhyML_Printf("\t\t %sBlosum62%s | %sMtMam%s | %sMtArt%s | %sHIVw%s |  %sHIVb%s | %scustom%s\n",
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT,
	 LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-f %se%s, %sm%s, or %s\"fA fC fG fT\"%s\n",BOLD,LINE,BOLD,LINE,BOLD,LINE,FLAT);
  PhyML_Printf("\t\t%se%s : the character frequencies are determined as follows : \n",LINE,FLAT);
  PhyML_Printf("%s\t\t- %sNucleotide%s sequences: (Empirical) the equilibrium base frequencies are estimated by counting\n"
	 "\t\t the occurence of the different bases in the alignment.\n",FLAT,LINE,FLAT);
  PhyML_Printf("%s\t\t- %sAmino-acid%s sequences: (Empirical) the equilibrium amino-acid frequencies are estimated by counting\n"
"\t\t the occurence of the different amino-acids in the alignment.\n",FLAT,LINE,FLAT);
  PhyML_Printf("\n");
  PhyML_Printf("\t\t%sm%s : the character frequencies are determined as follows : \n",LINE,FLAT);
  PhyML_Printf("%s\t\t- %sNucleotide%s sequences: (ML) the equilibrium base frequencies are estimated using maximum likelihood \n",FLAT,LINE,FLAT);
  PhyML_Printf("%s\t\t- %sAmino-acid%s sequences: (Model) the equilibrium amino-acid frequencies are estimated using\n"
"\t\t the frequencies defined by the substitution model.\n",FLAT,LINE,FLAT);
  PhyML_Printf("\n");
  PhyML_Printf("\t\t%s\"fA fC fG fT\"%s : only valid for nucleotide-based models. fA, fC, fG and fT are floating numbers that \n",LINE,FLAT);
  PhyML_Printf("\t\t correspond to the frequencies of A, C, G and T respectively.\n");

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-t (or --ts/tv) %sts/tv_ratio%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tts/tv_ratio%s : transition/transversion ratio. DNA sequences only.\n",FLAT);
  PhyML_Printf("\t\tCan be a fixed positive value (ex:4.0) or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-v (or --pinv) %sprop_invar%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tprop_invar%s : proportion of invariable sites.\n",FLAT);
  PhyML_Printf("\t\tCan be a fixed value in the [0,1] range or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-c (or --nclasses) %snb_subst_cat%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tnb_subst_cat%s : number of relative substitution rate categories. Default : %snb_subst_cat%s=1.\n",
	 FLAT,LINE,FLAT);
  PhyML_Printf("\t\tMust be a positive integer.\n");

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-a (or --alpha) %sgamma%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tgamma%s : distribution of the gamma distribution shape parameter.\n",FLAT);
  PhyML_Printf("\t\tCan be a fixed positive value or %se%s to get the maximum likelihood estimate.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-s (or --search) %smove%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tTree topology search operation option.\n");
  PhyML_Printf("\t\tCan be either %sNNI%s (default, fast) or %sSPR%s (a bit slower than NNI) or %sBEST%s (best of NNI and SPR search).\n",LINE,FLAT,LINE,FLAT,LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-u (or --inputtree) %suser_tree_file%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tuser_tree_file%s : starting tree filename. The tree must be in Newick format.\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t-o %sparams%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tThis option focuses on specific parameter optimisation.\n");
  PhyML_Printf("\t\t%sparams%s=tlr : tree topology (t), branch length (l) and rate parameters (r) are optimised.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sparams%s=tl  : tree topology and branch length are optimised.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sparams%s=lr  : branch length and rate parameters are optimised.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sparams%s=l   : branch length are optimised.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sparams%s=r   : rate parameters are optimised.\n",LINE,FLAT);
  PhyML_Printf("\t\t%sparams%s=n   : no parameter is optimised.\n",LINE,FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--rand_start%s\n",BOLD,FLAT);
  PhyML_Printf("\t\tThis option sets the initial tree to random.\n");
  PhyML_Printf("\t\tIt is only valid if SPR searches are to be performed.\n");

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--n_rand_starts %snum%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tnum%s is the number of initial random trees to be used.\n",FLAT);
  PhyML_Printf("\t\tIt is only valid if SPR searches are to be performed.\n");

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--r_seed %snum%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\tnum%s is the seed used to initiate the random number generator.\n",FLAT);
  PhyML_Printf("\t\tMust be an integer.\n");

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--print_site_lnl%s\n",BOLD,FLAT);
  PhyML_Printf("\t\t%sPrint the likelihood for each site in file *_phyml_lk.txt.\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--print_trace%s\n",BOLD,FLAT);
  PhyML_Printf("\t\t%sPrint each phylogeny explored during the tree search process\n",FLAT);
  PhyML_Printf("\t\t%sin file *_phyml_trace.txt.\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--run_id %sID_string%s\n",BOLD,LINE,FLAT);
  PhyML_Printf("\t\t%sAppend the string %sID_string%s at the end of each PhyML output file.\n",FLAT,LINE,FLAT);
  PhyML_Printf("\t\t%sThis option may be useful when running simulations involving PhyML.\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%s\n\t--quiet%s\n",BOLD,FLAT);
  PhyML_Printf("\t\t%sNo interactive question (for running in batch mode).\n",FLAT);

  PhyML_Printf("\n");

  PhyML_Printf("%sPHYLIP-LIKE INTERFACE\n""%s\n\tYou can also use PhyML with no argument, in this case change the value of\n",BOLD,FLAT);
  PhyML_Printf("%s\ta parameter by typing its corresponding character as shown on screen.\n\n",FLAT);

  PhyML_Printf("%sEXAMPLES\n\n"
	 "%s\tDNA interleaved sequence file, default parameters : ""%s  ./phyml -i seqs1"
	 "%s\n\tAA interleaved sequence file, default parameters :  ""%s  ./phyml -i seqs2 -d aa"
	 "%s\n\tAA sequential sequence file, with customization :   ""%s  ./phyml -i seqs3 -q -d aa -m JTT -c 4 -a e%s\n",BOLD,FLAT,BOLD,FLAT,BOLD,FLAT,BOLD,FLAT);
  Exit("");
}

/*********************************************************/

#define N_SEQUENCEFILE 1
#define N_DATATYPE 2
#define N_FORMAT 3
#define N_DATASETS 4
#define N_BOOTSTRAPSETS 5
#define N_MODELNAME 6
#define N_KAPPA 7
#define N_PROPORTIONINVAR 7 /*same as kappa*/
#define N_NBCATG 8
#define N_ALPHA 9
#define N_STARTINGTREE 10
#define N_OPT_TOPO 11
#define N_OPT_LENGTHSRATES 12

#define N_NB_PARAMS_DNA 13
#define N_NB_PARAMS_AA 12

option *Get_Input(int argc, char **argv)
{

  option *io;

  io = (option *)Make_Input();
  Set_Defaults_Input(io);
  Set_Defaults_Model(io->mod);
  Set_Defaults_Optimiz(io->mod->s_opt);

  putchar('\n');

  switch (argc)
    {
    case 1:
/* #if defined(MG) || defined(USE_OLD_INTERFACE) */
/*       Get_Input_Interactive(io); */
/* #else */
      Launch_Interface(io);
/* #endif */

      break;
    case 2:
      Usage();
      break;
    default:
      Read_Command_Line(io,argc,argv);
    }

  Print_Settings(io);
  return io;
}

/*********************************************************/


