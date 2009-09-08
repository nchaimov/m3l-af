/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

 */

#include <unistd.h>
#include <getopt.h>
#include "utilities.h"
#include "options.h"
#include "cl.h"
#include "models.h"
#include "free.h"
#include "interface.h"
#include "mg.h"
#include "m4.h"


/*********************************************************/
/**
* Fill the Option fields, with the argc array
*/
void Read_Command_Line(option *io, int argc, char **argv)
{
	int c;
	int open_ps_file;
	int use_gamma;
	int writemode;
	int fixed_props, proportions, numcatg;
	fixed_props = 0;
	proportions = 0;
	numcatg = 0;

	struct option longopts[] =
	{
			{"n_rgrft",           required_argument,NULL,0},
			{"n_globl",           required_argument,NULL,1},
			{"max_dist",          required_argument,NULL,2},
			{"n_optim",           required_argument,NULL,3},
			{"n_best",            required_argument,NULL,4},
			{"model",             required_argument,NULL,5},
			{"search",            required_argument,NULL,6},
			{"datatype",          required_argument,NULL,7},
			{"multiple",          required_argument,NULL,8},
			{"input",             required_argument,NULL,9},
			{"bootstrap",         required_argument,NULL,10},
			{"ts/tv",             required_argument,NULL,11},
			{"nclasses",          required_argument,NULL,12},
			{"pinv",              required_argument,NULL,13},
			{"alpha",             required_argument,NULL,14},
			{"inputtree",         required_argument,NULL,15},
			{"min_diff_lk_local", required_argument,NULL,16},
			{"min_diff_lk_global",required_argument,NULL,17},
			{"steph_spr",         no_argument,NULL,18},
			{"brent_it_max",      required_argument,NULL,19},
			{"rand_start",        no_argument,NULL,20},
			{"n_rand_starts",     required_argument,NULL,21},
			{"sequential",        no_argument,NULL,22},
			{"inside_opt",        no_argument,NULL,23},
			{"p_moves",           required_argument,NULL,24},
			{"fast_nni",          no_argument,NULL,25},
			{"g_pars",            no_argument,NULL,26},
			{"r_seed",            required_argument,NULL,27},
			{"collapse_boot",     required_argument,NULL,28},
			{"random_boot",       required_argument,NULL,29},
			{"print_trace",       no_argument,NULL,30},
			{"print_site_lnl",    no_argument,NULL,31},
			{"cov",               no_argument,NULL,32},
			{"cov_delta",         required_argument,NULL,33},
			{"cov_alpha",         required_argument,NULL,34},
			{"cov_ncats",         required_argument,NULL,35},
			{"ps",                no_argument,NULL,36},
			{"cov_free",          no_argument,NULL,37},
			{"no_gap",            no_argument,NULL,38},
			{"n_rr_branch",       required_argument,NULL,39},
			{"append",            no_argument,NULL,40},
			{"no_five_branch",    no_argument,NULL,41},
			{"pars_thresh",       required_argument,NULL,42},
			{"min_diff_lk_move",  required_argument,NULL,43},
			{"hybrid",            no_argument,NULL,44},
			{"use_median",        no_argument,NULL,45},
			{"run_id",            required_argument,NULL,46},
			{"pars",              no_argument,NULL,47},
			{"quiet",             no_argument,NULL,48},
			{"blprops",           required_argument,NULL,49},
			{"fixblprops",        no_argument,NULL,50},
			{"blclasses",         required_argument,NULL,51},
			{"iters_per_stage",         required_argument,NULL,52},
			{"num_anneal_stages",        required_argument,NULL,53},
			{"acc_ratio",         required_argument,NULL,54},
			{"temp_end",          required_argument,NULL,55},
			{"set_back",          required_argument,NULL,56},
			{"temp_start",        required_argument,NULL,57},
			{"max_alpha",         required_argument,NULL,58},
			{"brlen_sigma",       required_argument,NULL,59},
			{"pinvar_sigma",      required_argument,NULL,60},
			{"gamma_sigma",       required_argument,NULL,61},
			{"emig_sigma",        required_argument,NULL,62},
			{"prob_nni",          required_argument,NULL,63},
			{"prob_spr",          required_argument,NULL,64},
			{"prob_brlen",        required_argument,NULL,65},
			{"prob_gamma",        required_argument,NULL,66},
			{"prob_kappa",        required_argument,NULL,67},
			{"prob_lambda",       required_argument,NULL,68},
			{"prob_rr",           required_argument,NULL,69},
			{"prob_rate_proportion",required_argument,NULL,70},
			{"prob_topology",     required_argument,NULL,71},
			{"prob_pinvar",       required_argument,NULL,72},
			{"prob_pi",           required_argument,NULL,73},
			{"prob_emig",         required_argument,NULL,74},
			{"tau_start",		  required_argument,NULL,75},
			{"tau_end",			  required_argument,NULL,76},

			{0,0,0,0}
	};

	writemode = 1;
	open_ps_file = 0;
	use_gamma = 0;
	while((c = getopt_long(argc,argv,"qi:d:m:b:n:t:f:zk:v:c:a:u:l:w:ho:s:x:g:ep:z",longopts,NULL)) != -1)
	{
		switch(c)
		{
		case 48 :
		{
			io->quiet = 1;
			break;
		}
		case 'p' : case 47 :
		{
			io->in_tree = 1;
			break;
		}
		case 46 :
		{
			io->append_run_ID = 1;
			strcpy(io->run_id_string,optarg);
			break;
		}
		case 45 :
		{
			io->mod->gamma_median = 1;
			break;
		}
		case 44 :
		{
			io->mod->s_opt->hybrid_thresh = 0;
			break;
		}
		case 43 :
		{
			io->mod->s_opt->min_diff_lk_move = atof(optarg);
			if(io->mod->s_opt->min_diff_lk_move < 0)
			{
				char choix;
				PhyML_Printf("\n. Min_diff_lk_move must be a double greater than 0.\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			break;
		}
		case 42 :
		{
			io->mod->s_opt->pars_thresh = (int)atoi(optarg);
			/* 	      if(io->mod->s_opt->pars_thresh < 0) */
			/* 		{ */
			/* 		  PhyML_Printf("\n. The parsimony threshold must be an integer greater than 0.\n"); */
			/* 		  PhyML_Printf("\n. Type any key to exit.\n"); */
			/* 		  Exit("\n"); */
			/* 		} */
			break;
		}
		case 41 :
		{
			io->mod->s_opt->opt_five_branch = 0;
			break;
		}
		case 40 :
		{
			writemode = 2;
			break;
		}
		case 39 :
		{
			char choix;
			io->mod->n_rr_branch = (int)atoi(optarg);
			if(io->mod->n_rr_branch < 1)
			{
				PhyML_Printf("\n. The number of classes must be an integer greater than 0.\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			break;
		}
		case 38 :
		{
			io->rm_ambigu = 1;
			break;
		}
		case 37 :
		{
			io->mod->s_opt->opt_cov_free_rates = 1;
			break;
		}
		case 36 :
		{
			open_ps_file = 1;
			break;
		}
		case 35 :
		{
#ifdef M4
			io->m4_model = YES;
			if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));
			io->mod->m4mod->n_h = (int)atoi(optarg);

			if(io->mod->m4mod->n_h < 1)
			{
				char choix;
				PhyML_Printf("\n. The number of classes must be greater than 0.\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
#endif
			break;
		}
		case 34 :
		{
#ifdef M4
			io->m4_model = YES;
			if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));
			io->mod->m4mod->use_cov_alpha = 1;
			io->mod->m4mod->use_cov_free  = 0;

			if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
					!strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
			{
				io->mod->s_opt->opt_cov_alpha = 1;
				io->mod->m4mod->alpha         = 1.0;
			}
			else
			{
				io->mod->m4mod->alpha = (m3ldbl)atof(optarg);

				if(io->mod->m4mod->alpha < 1.E-5)
				{
					char choix;
					PhyML_Printf("\n. The value of alpha must be greater than 1.E-5.\n");
					PhyML_Printf("\n. Type any key to exit.\n");
					if(!scanf("%c",&choix)) Exit("\n");
					Exit("\n");
				}
			}
#endif
			break;
		}

		case 33 :
		{
#ifdef M4
			io->m4_model = YES;
			if(!io->mod->m4mod) io->mod->m4mod = (m4 *)M4_Make_Light((io->mod->datatype == NT)?(4):(20));

			if(!strcmp(optarg,"e") || !strcmp(optarg,"E") ||
					!strcmp(optarg,"estimated") || !strcmp(optarg,"ESTIMATED"))
			{
				io->mod->s_opt->opt_cov_delta = 1;
				io->mod->m4mod->delta         = 1.0;
			}
			else
			{
				io->mod->m4mod->delta = (m3ldbl)atof(optarg);

				if(atof(optarg) < 1.E-10)
				{
					char choix;
					PhyML_Printf("\n. The value of delta must be larger than 1.E-10.\n");
					PhyML_Printf("\n. Type any key to exit.\n");
					if(!scanf("%c",&choix)) Exit("\n");
					Exit("\n");
				}
			}
#endif
			break;
		}
		case 32 :
		{
#ifdef M4
			io->m4_model = YES;
#endif
			break;
		}
		case 31 :
		{
			io->print_site_lnl = 1;
			break;
		}
		case 30 :
		{
			io->print_trace = 1;
			break;
		}
		case 29 :
		{
			io->random_boot_seq_order = (int)atoi(optarg);
			break;
		}
		case 28 :
		{
			io->collapse_boot = (int)atoi(optarg);
			break;
		}
		case 27 :
		{
			io->r_seed = (int)atoi(optarg);
			break;
		}
		case 26 :
		{
			io->mod->s_opt->general_pars = 1;
			break;
		}
		case 25 :
		{
			io->mod->s_opt->fast_nni = 1;
			break;
		}
		case 24 :
		{
			io->mod->s_opt->p_moves_to_examine = (double)atof(optarg);
			break;
		}
		case 23 :
		{
			io->mod->s_opt->wim_inside_opt = 1;
			break;
		}
		case 0 :
		{
			io->mod->s_opt->wim_n_rgrft = atoi(optarg);
			break;
		}
		case 1 :
		{
			io->mod->s_opt->wim_n_globl = atoi(optarg);
			break;
		}
		case 2 :
		{
			io->mod->s_opt->wim_max_dist = atoi(optarg);
			break;
		}
		case 3 :
		{
			io->mod->s_opt->wim_n_optim = atoi(optarg);
			break;
		}
		case 4 :
		{
			io->mod->s_opt->wim_n_best = atoi(optarg);
			break;
		}
		case 16 :
		{
			io->mod->s_opt->min_diff_lk_local = atof(optarg);
			break;
		}

		case 17 :
		{
			io->mod->s_opt->min_diff_lk_global = atof(optarg);
			break;
		}
		case 18 :
		{
			io->mod->s_opt->steph_spr = 0;
			io->mod->s_opt->greedy    = 1;
			break;
		}
		case 19 :
		{
			io->mod->s_opt->brent_it_max = atoi(optarg);
			break;
		}
		case 20 :
		{
			io->mod->s_opt->random_input_tree = 1;
			break;
		}
		case 21 :
		{
			io->mod->s_opt->random_input_tree = 1;
			io->mod->s_opt->n_rand_starts = atoi(optarg);
			if(io->mod->s_opt->n_rand_starts < 1) Exit("\n. Nunmber of random starting trees must be > 0.\n\n");
		}
		case 's':case 6:
		{
			if((!strcmp(optarg,"spr")) || (!strcmp(optarg,"SPR")))
			{
				io->mod->s_opt->topo_search = SPR_MOVE;
				io->user_topo 						= 1;
				io->mod->s_opt->greedy      = (io->mod->s_opt->steph_spr)?(0):(1);
			}
			else if((!strcmp(optarg,"nni")) || (!strcmp(optarg,"NNI")))
			{
				io->mod->s_opt->topo_search         = NNI_MOVE;
				io->user_topo 						= 1;
				io->mod->s_opt->random_input_tree   = 0;
			}
			else if((!strcmp(optarg,"best")) || (!strcmp(optarg,"BEST")))
			{
				io->mod->s_opt->topo_search = BEST_OF_NNI_AND_SPR;
				io->user_topo 						= 1;
				io->mod->s_opt->greedy      = (io->mod->s_opt->steph_spr)?(0):(1);
			}
			else if((!strcmp(optarg,"sta")) || (!strcmp(optarg,"STA")))
			{
				io->mod->s_opt->topo_search = SIMULATED_THERMAL_ANNEALING;
				io->user_topo 						= 1;
			}
			else if((!strcmp(optarg,"sqa")) || (!strcmp(optarg,"SQA")))
			{
				io->mod->s_opt->topo_search = SIMULATED_QUANTUM_ANNEALING;
				io->user_topo 						= 1;
			}
			break;
		}

		case 'd':case 7:
			if(!strcmp(optarg,"nt"))
				/*         if(atoi(optarg) == NT) */
			{
				io->mod->datatype = NT;
				io->mod->stepsize = 1;
				io->mod->ns = 4;

				if(
						(io->mod->whichmodel == LG)       ||
						(io->mod->whichmodel == WAG)       ||
						(io->mod->whichmodel == DAYHOFF)   ||
						(io->mod->whichmodel == JTT)       ||
						(io->mod->whichmodel == BLOSUM62)  ||
						(io->mod->whichmodel == MTREV)     ||
						(io->mod->whichmodel == RTREV)     ||
						(io->mod->whichmodel == CPREV)     ||
						(io->mod->whichmodel == DCMUT)     ||
						(io->mod->whichmodel == VT)        ||
						(io->mod->whichmodel == MTMAM)     ||
						(io->mod->whichmodel == MTART)     ||
						(io->mod->whichmodel == HIVW)      ||
						(io->mod->whichmodel == HIVB)      ||
						(io->mod->whichmodel == CUSTOMAA)
				)
				{
					io->mod->whichmodel = HKY85;
					strcpy(io->mod->modelname, "HKY85\0");
				}
			}
			else if (!strcmp(optarg,"aa"))
				/*         else if (atoi(optarg) == AA) */
			{
				io->mod->datatype         = AA;
				io->mod->stepsize         = 1;
				io->mod->s_opt->opt_kappa = 0;
				io->mod->ns               = 20;
				if(
						(io->mod->whichmodel == JC69)   ||
						(io->mod->whichmodel == K80)    ||
						(io->mod->whichmodel == F81)    ||
						(io->mod->whichmodel == HKY85)  ||
						(io->mod->whichmodel == F84)    ||
						(io->mod->whichmodel == TN93)   ||
						(io->mod->whichmodel == GTR)    ||
						(io->mod->whichmodel == CUSTOM)
				)
				{
					io->mod->whichmodel = LG;
					strcpy(io->mod->modelname, "LG\0");
				}
			}
			else
			{
				char choix;
				PhyML_Printf("\n. Unknown argument to -d option: please use `nt' for DNA or `aa' for Amino-Acids\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}

			break;

		case 'm': case 5 :
		{
			if(!isalpha(optarg[0]))
			{
				strcpy(io->mod->custom_mod_string,optarg);

				if(strlen(io->mod->custom_mod_string) != 6)
				{
					Warn_And_Exit("\n. The string should be of length 6.\n");
				}
				else
				{
					Translate_Custom_Mod_String(io->mod);
				}

				io->mod->datatype = NT;
				io->mod->whichmodel = CUSTOM;
				strcpy(io->mod->modelname, "custom");
				io->mod->s_opt->opt_kappa     = 0;
				io->mod->s_opt->opt_rr        = 1;
				io->mod->s_opt->opt_num_param = 1;
			}

			if (strcmp(optarg, "JC69") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = JC69;
				strcpy(io->mod->modelname, "JC69");
				io->mod->s_opt->opt_kappa  = 0;
			}
			else if(strcmp(optarg, "K80") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = K80;
				strcpy(io->mod->modelname, "K80");
			}
			else if(strcmp(optarg, "F81") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = F81;
				strcpy(io->mod->modelname, "F81");
				io->mod->s_opt->opt_kappa  = 0;
			}
			else if (strcmp(optarg, "HKY85") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = HKY85;
				strcpy(io->mod->modelname, "HKY85");
			}
			else if(strcmp(optarg, "F84") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = F84;
				strcpy(io->mod->modelname, "F84");
			}
			else if (strcmp (optarg,"TN93") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = TN93;
				strcpy(io->mod->modelname, "TN93");
				if(io->mod->s_opt->opt_kappa)
				{
					io->mod->s_opt->opt_lambda = 1;
				}
			}
			else if(strcmp (optarg, "GTR") == 0)
			{
				io->mod->datatype = NT;
				io->mod->whichmodel = GTR;
				strcpy(io->mod->modelname, "GTR");
				io->mod->s_opt->opt_kappa = 0;
				io->mod->s_opt->opt_rr    = 1;
			}
			else if(strcmp(optarg, "Dayhoff") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = DAYHOFF;
				strcpy(io->mod->modelname, "Dayhoff");
			}
			else if(strcmp (optarg, "JTT") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = JTT;
				strcpy(io->mod->modelname, "JTT");
			}
			else if(strcmp(optarg, "MtREV") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = MTREV;
				strcpy(io->mod->modelname,"MtREV");
			}
			else if(strcmp (optarg, "LG") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = LG;
				strcpy(io->mod->modelname, "LG");
			}
			else if(strcmp (optarg, "WAG") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = WAG;
				strcpy(io->mod->modelname, "WAG");
			}
			else if(strcmp(optarg, "DCMut") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = DCMUT;
				strcpy(io->mod->modelname, "DCMut");
			}
			else if(strcmp (optarg, "RtREV") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = RTREV;
				strcpy(io->mod->modelname, "RtREV");
			}
			else if(strcmp(optarg, "CpREV") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = CPREV;
				strcpy(io->mod->modelname, "CpREV");
			}
			else if(strcmp(optarg, "VT") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = VT;
				strcpy(io->mod->modelname, "VT");
			}
			else if(strcmp(optarg, "Blosum62") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = BLOSUM62;
				strcpy(io->mod->modelname,"Blosum62");
			}
			else if(strcmp(optarg, "MtMam") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = MTMAM;
				strcpy(io->mod->modelname, "MtMam");
			}
			else if (strcmp(optarg,"MtArt") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = MTART;
				strcpy(io->mod->modelname, "MtArt");
			}
			else if (strcmp(optarg,"HIVw") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = HIVW;
				strcpy(io->mod->modelname, "HIVw");
			}
			else if(strcmp(optarg, "HIVb") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = HIVB;
				strcpy(io->mod->modelname, "HIVb");
			}
			else if (strcmp(optarg, "custom") == 0)
			{
				io->mod->datatype = AA;
				io->mod->whichmodel = CUSTOMAA;
				strcpy(io->mod->modelname, "Read from file");
			}
			break;
		}

		case 'a':case 14 :
		{
			use_gamma = 1;
			if ((strcmp (optarg, "e") == 0) ||
					(strcmp (optarg, "E") == 0) ||
					(strcmp (optarg, "estimated") == 0) ||
					(strcmp (optarg, "ESTIMATED") == 0))
			{
				io->mod->s_opt->opt_alpha     = 1;
				io->mod->s_opt->opt_num_param = 1;
			}
			else if ((!atof(optarg)) || (atof(optarg) < 1.E-10))
			{
				char choix;
				PhyML_Printf("\n. Alpha must be > 1.E-10.\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else
			{
				io->mod->alpha = (m3ldbl)atof(optarg);
				io->mod->s_opt->opt_alpha  = 0;
			}
			break;
		}
		case 'b':case 10:
		{
			if (atoi(optarg) < -4)
			{
				char choix;
				PhyML_Printf("\n. Branch test value must be a positive integer for bootstrap, or between -1 and -4 for aLRT branch test\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else
			{
				if((int)atoi(optarg) > 0)
				{
					io->ratio_test       = 0;
					io->mod->bootstrap   = (int)atoi(optarg);
					io->print_boot_trees = 1;

					if(io->n_data_sets > 1)
					{
						char choix;
						PhyML_Printf("\n. Bootstrap option is not allowed with multiple data sets\n");
						PhyML_Printf("\n. Type any key to exit.\n");
						if(!scanf("%c",&choix)) Exit("\n");
						Exit("\n");
					}
				}
				else if (atoi(optarg)==0)
				{
					io->mod->bootstrap = 0;
					io->ratio_test     = 0;
				}
				else
				{
					io->mod->bootstrap = 0;
					io->ratio_test     = -(int)atoi(optarg);
				}
			}
			break;
		}
		case 'c':case 12:
		{
			if ((!atoi(optarg)) || (atoi(optarg) < 0))
			{
				char choix;
				PhyML_Printf("\n. Unknown argument to -c option: the number of categories must be a positive integer\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else io->mod->n_catg = atoi(optarg);
			break;
		}
		case 'f':
		{
			if(!strcmp(optarg,"e"))
			{
				if (io->mod->datatype == NT)
					io->mod->s_opt->opt_state_freq = 0;
				else
					io->mod->s_opt->opt_state_freq = 1;

				if((io->mod->whichmodel == JC69) ||
						(io->mod->whichmodel == K80))
				{
					char choix;
					PhyML_Printf("\n. Invalid model settings (option '-f').\n");
					PhyML_Printf("\n. Type any key to exit.\n");
					if(!scanf("%c",&choix)) Exit("\n");
					Exit("\n");
				}
			}
			else if(!strcmp(optarg,"m"))
			{
				if (io->mod->datatype == NT)
					io->mod->s_opt->opt_state_freq = 1;
				else
					io->mod->s_opt->opt_state_freq = 0;
			}
			else if(!isalpha(optarg[0]))
			{
				m3ldbl sum;

				io->mod->s_opt->opt_state_freq  = 0;
				io->mod->s_opt->user_state_freq = 1;

				sum =
						(io->mod->user_b_freq[0] +
								io->mod->user_b_freq[1] +
								io->mod->user_b_freq[2] +
								io->mod->user_b_freq[3]);

				io->mod->user_b_freq[0] /= sum;
				io->mod->user_b_freq[1] /= sum;
				io->mod->user_b_freq[2] /= sum;
				io->mod->user_b_freq[3] /= sum;

				if(io->mod->user_b_freq[0] < .0 ||
						io->mod->user_b_freq[1] < .0 ||
						io->mod->user_b_freq[2] < .0 ||
						io->mod->user_b_freq[3] < .0 ||
						io->mod->user_b_freq[0] > 1. ||
						io->mod->user_b_freq[1] > 1. ||
						io->mod->user_b_freq[2] > 1. ||
						io->mod->user_b_freq[3] > 1.)
				{
					Warn_And_Exit("\n. Invalid base frequencies.\n");
				}
			}
			break;
		}

		case 'h':
		{
			Usage();
			break;
		}

		case 'i':case 9:
		{
			char *tmp;
			tmp = (char *) mCalloc (T_MAX_FILE, sizeof(char));
			if (strlen (optarg) > T_MAX_FILE -16)
			{
				char choix;
				strcpy (tmp, "\n. The file name'");
				strcat (tmp, optarg);
				strcat (tmp, "' is too long.\n");
				PhyML_Printf("%s",tmp);
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}

			else if (!Filexists (optarg))
			{
				char choix;
				strcpy (tmp, "\n. The file '");
				strcat (tmp, optarg);
				strcat (tmp, "' does not exist.\n");
				PhyML_Printf("%s",tmp);
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else
			{
				strcpy(io->in_seq_file, optarg);
				io->fp_in_seq = Openfile(io->in_seq_file,0);
				strcpy(io->out_tree_file,optarg);
#ifdef PHYML
				strcat(io->out_tree_file,"_phyml_tree");
#else
				strcat(io->out_tree_file,"_mc_tree.txt");
#endif
				strcpy(io->out_stats_file,optarg);
#ifdef PHYML
				strcat(io->out_stats_file,"_phyml_stats");
#else
				strcat(io->out_stats_file,"_mc_stats.txt");
#endif
			}
			Free (tmp);
			break;
		}

		case 't':case 11:
		{
			if ((io->mod->whichmodel != JC69) &&
					(io->mod->whichmodel != F81)  &&
					(io->mod->whichmodel != GTR))
			{
				if ((strcmp(optarg, "e") == 0) ||
						(strcmp(optarg, "E") == 0) ||
						(strcmp(optarg, "estimated") == 0) ||
						(strcmp(optarg, "ESTIMATED") == 0))
				{
					io->mod->kappa                 = 4.0;
					io->mod->s_opt->opt_num_param  = 1;
					io->mod->s_opt->opt_kappa      = 1;
					if (io->mod->whichmodel == TN93)
						io->mod->s_opt->opt_lambda   = 1;
				}
				else
				{
					if ((!atof(optarg)) || (atof(optarg) < .0))
					{
						char choix;
						PhyML_Printf("\n. The ts/tv ratio must be a positive number\n");
						PhyML_Printf("\n. Type any key to exit.\n");
						if(!scanf("%c",&choix)) Exit("\n");
						Exit("\n");
					}
					else
					{
						io->mod->kappa = (m3ldbl)atof(optarg);
						io->mod->s_opt->opt_kappa  = 0;
						io->mod->s_opt->opt_lambda = 0;
					}
				}
			}
			break;
		}
		case 'n':case 8:
		{
			if ((!atoi(optarg)) || (atoi(optarg) < 0))
			{
				char choix;
				PhyML_Printf("\n. The number of alignments must be a positive integer\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else io->n_data_sets = atoi (optarg);
			break;
		}
		case 'q':case 22:
		{
			io->interleaved = 0;
			break;
		}










		/**
		* The next three cases are for obtaining information about the mixed branch
		* length model.
		*/
		case 'w': case 49: //JSJ: proportions "blprops"
		{
			proportions = 0; // use this as the counter for the numbr of props read in.
			// doubles to let one know if any props have been read in...
			io->user_props = 1;
			if(optarg[0] != '['){
				char choix;
				PhyML_Printf ("\n. The proportions must be in the form [NUM1,NUM2,...,NUMN] with no spaces and surrounded by '[' and ']'.");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			int i,j = 0;
			char aprop[100];
			for(i=1; i < strlen(optarg); i++){
				if(optarg[i] != ',' && optarg[i] != ']'){
					aprop[j] = optarg[i];
					j++;
				}else{
					aprop[j+1] = '\0'; // end the string
					j = 0; //start j over
					io->props[proportions] = atof(aprop);
					proportions++; //increment the props by one
				}
			}
			break;
		}
		case 'z': case 50: //JSJ: fixed proportions "fixblprops"
		{
			fixed_props = 1; //set the flag to let other functions know that this has been set
			io->fixed_props = 1;
			io->mod->s_opt->opt_props = io->fixed_props;
			break;
		}
		case 'l': case 51: //JSJ: number branch length categories "blclasses"
		{
			io->n_l = atoi(optarg);
			numcatg = io->n_l;
			if(io->n_l > 1 && io->user_topo == 0){
				io->mod->s_opt->topo_search = SIMULATED_THERMAL_ANNEALING;
			}
			if(fixed_props == 0 && io->n_l > 1){ //if the fprops option hasn't been flagged, and more than one catg, default to optimize
				io->fixed_props = 0;
				io->mod->s_opt->opt_props = io->fixed_props;
			}
			if(io->n_l < 1 || io->n_l > MAX_BL_SET){
				PhyML_Printf("\n. The number of branch length categories (--nclasses or -l) must be at least 1 and no more than %i",MAX_BL_SET);
				PhyML_Printf("\n. Type any key to exit.\n");
				Exit("\n");
			}
			if(io->n_l > 1){
				io->mod->s_opt->opt_five_branch = 0;
				if(proportions == 0){
					Update_Default_Props(io);
				}
			}
			break;
		}
		case 52: //JSJ: iters_per_stage
		{
			int tmp = io->iters_per_stage;
			io->iters_per_stage = atoi(optarg);
			if(io->iters_per_stage < 1) io->iters_per_stage = tmp;
			break;
		}
		case 53: //JSJ: num_anneal_stages
		{
			int tmp = io->num_anneal_stages;
			io->num_anneal_stages = atoi(optarg);
			//printf("(cl.c) FOUND io->num_anneal_stages = %s, %d", optarg,  io->num_anneal_stages);
			if(io->num_anneal_stages < 1) io->num_anneal_stages = tmp;
			break;
		}
		case 54: //JSJ: acc_ratio
		{
			double tmp = io->acc_ratio;
			io->acc_ratio = atof(optarg);
			if(io->acc_ratio < 0.0) io->acc_ratio = tmp;
			break;
		}
		case 55:
		{
			double tmp = io->temp_end;
			io->temp_end = atof(optarg);
			if(io->temp_end < 0.0) io->temp_end = tmp;
			break;
		}
		//		{"temp_end",          required_argument,NULL,55},
		case 56:
		{
			int tmp = io->set_back;
			io->set_back = atoi(optarg);
			if(io->set_back < 0) io->set_back = tmp;
			break;
		}
		//		{"set_back",          required_argument,NULL,56},
		case 57:
		{
			double tmp = io->temp_start;
			io->temp_start = atof(optarg);
			if(io->temp_start < 0.0) io->temp_start = tmp;
			break;
		}
		//		{"temp_start",        required_argument,NULL,57},
		case 58:
		{
			double tmp = io->max_alpha;
			io->max_alpha = atof(optarg);
			if(io->max_alpha < 0.0) io->max_alpha = tmp;
			break;
		}
		//		{"max_alpha",         required_argument,NULL,58},
		case 59:
		{
			double tmp = io->brlen_sigma;
			io->brlen_sigma = atof(optarg);
			if(io->brlen_sigma < 0.0) io->brlen_sigma = tmp;
			break;
		}
		//		{"brlen_sigma",       required_argument,NULL,59},
		case 60:
		{
			double tmp = io->pinvar_sigma;
			io->pinvar_sigma = atof(optarg);
			if(io->pinvar_sigma < 0.0) io->pinvar_sigma = tmp;
			break;
		}
		//		{"pinvar_sigma",      required_argument,NULL,60},
		case 61:
		{
			double tmp = io->gamma_sigma;
			io->gamma_sigma = atof(optarg);
			if(io->gamma_sigma < 0.0) io->gamma_sigma = tmp;
			break;

		}
		//		{"gamma_sigma",       required_argument,NULL,61},
		case 62:
		{
			double tmp = io->emig_sigma;
			io->emig_sigma = atof(optarg);
			if(io->emig_sigma < 0.0) io->emig_sigma = tmp;
			break;
		}
		//		{"emig_sigma",        required_argument,NULL,62},
		case 63:
		{
			double tmp = io->prob_NNI;
			io->prob_NNI = atof(optarg);
			if(io->prob_NNI < 0.0) io->prob_NNI = tmp;
			break;
		}
		//		{"prob_nni",          required_argument,NULL,63},
		case 64:
		{
			double tmp = io->prob_SPR;
			io->prob_SPR = atof(optarg);
			if(io->prob_SPR < 0.0) io->prob_SPR = tmp;
			break;
		}
		//		{"prob_spr",          required_argument,NULL,64},
		case 65:
		{
			double tmp = io->prob_brlen;
			io->prob_brlen = atof(optarg);
			if(io->prob_brlen < 0.0) io->prob_brlen = tmp;
			break;
		}
		//		{"prob_brlen",        required_argument,NULL,65},
		case 66:
		{
			double tmp = io->prob_gamma;
			io->prob_gamma = atof(optarg);
			if(io->prob_gamma < 0.0) io->prob_gamma = tmp;
			break;
		}
		//		{"prob_gamma",        required_argument,NULL,66},
		case 67:
		{
			double tmp = io->prob_kappa;
			io->prob_kappa = atof(optarg);
			if(io->prob_kappa < 0.0) io->prob_kappa = tmp;
			break;

		}
		//		{"prob_kappa",        required_argument,NULL,67},
		case 68:
		{
			double tmp = io->prob_lambda;
			io->prob_lambda = atof(optarg);
			if(io->prob_lambda < 0.0) io->prob_lambda = tmp;
			break;
		}
		//		{"prob_lambda",       required_argument,NULL,68},
		case 69:
		{
			double tmp = io->prob_rr;
			io->prob_rr = atof(optarg);
			if(io->prob_rr < 0.0) io->prob_rr = tmp;
			break;
		}
		//		{"prob_rr",           required_argument,NULL,69},
		case 70:
		{
			double tmp = io->prob_rate_proportion;
			io->prob_rate_proportion = atof(optarg);
			if(io->prob_rate_proportion < 0.0) io->prob_rate_proportion = tmp;
			break;
		}
		//		{"prob_rate_proportion",required_argument,NULL,70},
		case 71:
		{
			double tmp = io->prob_topology;
			io->prob_topology = atof(optarg);
			if(io->prob_topology < 0.0) io->prob_topology = tmp;
			break;
		}
		//		{"prob_topology",     required_argument,NULL,71},
		case 72:
		{
			double tmp = io->prob_pinvar;
			io->prob_pinvar = atof(optarg);
			if(io->prob_pinvar < 0.0) io->prob_pinvar = tmp;
			break;
		}
		//		{"prob_pinvar",       required_argument,NULL,72},
		case 73:
		{
			double tmp = io->prob_pi;
			io->prob_pi = atof(optarg);
			if(io->prob_pi < 0.0) io->prob_pi = tmp;
			break;
		}
		//		{"prob_pi",           required_argument,NULL,73},
		case 74:
		{
			double tmp = io->prob_emig;
			io->prob_emig = atof(optarg);
			if(io->prob_emig < 0.0) io->prob_emig = tmp;
			break;
		}
		case 75:
		{
			double tmp = io->tau_start;
			io->tau_start = atof(optarg);
			if(io->tau_start < 0.0) io->tau_start = tmp;
			break;
		}
		case 76:
		{
			double tmp = io->tau_end;
			io->tau_end = atof(optarg);
			if(io->tau_end < 0.0) io->tau_end = tmp;
			break;
		}

		case 'u':case 15:
		{
			char *tmp;
			tmp = (char *)mCalloc(T_MAX_FILE, sizeof(char));
			if (strlen (optarg) > T_MAX_FILE -11)
			{
				char choix;
				strcpy (tmp, "\n. The file name'");
				strcat (tmp, optarg);
				strcat (tmp, "' is too long.\n");
				PhyML_Printf("%s",tmp);
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else if (! Filexists (optarg))
			{
				char choix;
				strcpy (tmp, "\n. The file '");
				strcat (tmp, optarg);
				strcat (tmp, "' doesn't exist.\n");
				PhyML_Printf("%s",tmp);
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else
			{
				strcpy(io->in_tree_file, optarg);
				io->in_tree = 2;
				io->fp_in_tree = Openfile(io->in_tree_file,0);
			}
			Free(tmp);
			break;
		}

		case 'v':case 13:
		{
			if ((strcmp (optarg, "e") == 0) ||
					(strcmp (optarg, "E") == 0) ||
					(strcmp (optarg, "estimated") == 0) ||
					(strcmp (optarg, "ESTIMATED") == 0)) {
				io->mod->s_opt->opt_num_param = 1;
				io->mod->s_opt->opt_pinvar    = 1;
				io->mod->invar                = 1;
			}
			else if ((atof(optarg) < 0.0) || (atof(optarg) > 1.0))
			{
				char choix;
				PhyML_Printf("\n. The proportion of invariable site must be a number between 0.0 and 1.0\n");
				PhyML_Printf("\n. Type any key to exit.");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			else
			{
				io->mod->pinvar = (m3ldbl)atof(optarg);
				if (io->mod->pinvar > 0.0+MDBL_MIN)
					io->mod->invar = 1;
				else
					io->mod->invar = 0;
				io->mod->s_opt->opt_pinvar = 0;
			}
			break;
		}
		case 'o':
		{
			if(!strcmp(optarg,"tlr"))
			{
				io->mod->s_opt->opt_topo = 1;
				io->mod->s_opt->opt_bl   = 1;
			}
			else if(!strcmp(optarg,"tl"))
			{
				io->mod->s_opt->opt_topo = 1;
				io->mod->s_opt->opt_bl   = 1;
			}
			else if(!strcmp(optarg,"t"))
			{
				Warn_And_Exit("\n. You can't optimize the topology without adjusting branch length too...\n");
			}
			else if(!strcmp(optarg,"lr"))
			{
				io->mod->s_opt->opt_topo = 0;
				io->mod->s_opt->opt_bl   = 1;
			}
			else if(!strcmp(optarg,"l"))
			{
				io->mod->s_opt->opt_topo = 0;
				io->mod->s_opt->opt_bl   = 1;
			}
			else if(!strcmp(optarg,"r"))
			{
				io->mod->s_opt->opt_topo = 0;
				io->mod->s_opt->opt_bl   = 0;
			}
			else if(!strcmp(optarg,"none") || !strcmp(optarg,"n"))
			{
				io->mod->s_opt->opt_topo = 0;
				io->mod->s_opt->opt_bl   = 0;
			}
			else
			{
				char choix;
				PhyML_Printf ("\n. The optimization parameter must be 'tlr' or 'tl' or 'lr' or 'l' or 'r' or ''.");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}
			break;
		}
//		case '?':
//		{
//			char choix;
//			if (isprint (optopt))
//				PhyML_Printf ("\n. Unknown option `-%c'.\n", optopt);
//			else
//				PhyML_Printf ("\n. Unknown option character `\\x%x'.\n", optopt);
//			PhyML_Printf("\n. Type any key to exit.\n");
//			if(!scanf("%c",&choix)) Exit("\n");
//			Exit("\n");
//			break;
//		}

		default:
			Usage();
		}
	}


#ifndef PHYML
	if((open_ps_file) || (io->m4_model == YES))
	{
		strcpy(io->out_ps_file,io->in_seq_file);
		strcat(io->out_ps_file, "_mc_tree.ps");
		io->fp_out_ps = Openfile(io->out_ps_file,1);
	}
#endif





	if((io->mod->s_opt->n_rand_starts)           &&
			(io->mod->s_opt->topo_search == NNI_MOVE) &&
			(io->mod->s_opt->random_input_tree))
	{
		Warn_And_Exit("\n. The random starting tree option is only compatible with SPR based search options.\n");
	}

	if ((io->mod->datatype == NT) && (io->mod->whichmodel > 10))
	{
		char choix;
		PhyML_Printf("\n. Err: model incompatible with the data type. Please use JC69, K80, F81, HKY, F84, TN93 or GTR\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Warn_And_Exit("\n");
	}
	else if ((io->mod->datatype == AA) && (io->mod->whichmodel < 11))
	{
		char choix;
		PhyML_Printf("\n. Err: model incompatible with the data type. Please use LG, Dayhoff, JTT, MtREV, WAG, DCMut, RtREV, CpREV, VT, Blosum62, MtMam, MtArt, HIVw or HIVb.\n");
		PhyML_Printf("\n. Type any key to exit.\n");
		if(!scanf("%c",&choix)) Exit("\n");
		Exit("\n");
	}

	if(io->m4_model == YES)
	{
#ifdef M4
		io->mod->ns *= io->mod->m4mod->n_h;
		io->mod->use_m4mod = 1;
		M4_Make_Complete(io->mod->m4mod->n_h,
				io->mod->m4mod->n_o,
				io->mod->m4mod);
#endif
	}
	else
	{
		io->mod->s_opt->opt_cov_delta      = 0;
		io->mod->s_opt->opt_cov_alpha      = 0;
		io->mod->s_opt->opt_cov_free_rates = 0;
	}

	if((io->mod->s_opt->opt_cov_free_rates) && (io->mod->s_opt->opt_cov_alpha))
	{
		io->mod->s_opt->opt_cov_free_rates = 0;
		io->mod->m4mod->use_cov_alpha      = 0;
		io->mod->m4mod->use_cov_free       = 1;
	}

	if(io->print_site_lnl)
	{
		strcpy(io->out_lk_file,io->in_seq_file);
		strcat(io->out_lk_file, "_phyml_lk");
		if(io->append_run_ID) { strcat(io->out_lk_file,"_"); strcat(io->out_lk_file,io->run_id_string); }
		strcat(io->out_lk_file, ".txt");
		io->fp_out_lk = Openfile(io->out_lk_file,1);
	}

	if(io->print_trace)
	{
		strcpy(io->out_trace_file,io->in_seq_file);
		strcat(io->out_trace_file,"_phyml_trace");
		if(io->append_run_ID) { strcat(io->out_trace_file,"_"); strcat(io->out_trace_file,io->run_id_string); }
		strcat(io->out_trace_file,".txt");
		io->fp_out_trace = Openfile(io->out_trace_file,1);
	}

	if(io->mod->s_opt->random_input_tree)
	{
		strcpy(io->out_trees_file,io->in_seq_file);
		strcat(io->out_trees_file,"_phyml_trees");
		if(io->append_run_ID) { strcat(io->out_trees_file,"_"); strcat(io->out_trees_file,io->run_id_string); }
		strcat(io->out_trees_file,".txt");
		io->fp_out_trees = Openfile(io->out_trees_file,1);
	}

	if((io->print_boot_trees) && (io->mod->bootstrap > 0))
	{
		strcpy(io->out_boot_tree_file,io->in_seq_file);
		strcat(io->out_boot_tree_file,"_phyml_boot_trees");
		if(io->append_run_ID) { strcat(io->out_boot_tree_file,"_"); strcat(io->out_boot_tree_file,io->run_id_string); }
		strcat(io->out_boot_tree_file,".txt");
		io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);

		strcpy(io->out_boot_stats_file,io->in_seq_file);
		strcat(io->out_boot_stats_file,"_phyml_boot_stats");
		if(io->append_run_ID) { strcat(io->out_boot_stats_file,"_"); strcat(io->out_boot_stats_file,io->run_id_string); }
		strcat(io->out_boot_stats_file,".txt");
		io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);
	}

	if(io->append_run_ID)
	{
		strcat(io->out_tree_file,"_");
		strcat(io->out_stats_file,"_");
		strcat(io->out_tree_file,io->run_id_string);
		strcat(io->out_stats_file,io->run_id_string);
	}
	strcat(io->out_tree_file,".txt");
	strcat(io->out_stats_file,".txt");

	io->fp_out_tree  = Openfile(io->out_tree_file,writemode);
	io->fp_out_stats = Openfile(io->out_stats_file,writemode);

	return;
}

/*********************************************************/
