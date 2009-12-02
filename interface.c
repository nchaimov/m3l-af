/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

 */

#include "utilities.h"
#include "options.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "interface.h"
#ifdef MG
#include "multigene.h"
#endif


void Launch_Interface(option *io)
{
	Launch_Interface_Input(io);
	io->ready_to_go = 0;
	do
	{
		switch(io->curr_interface)
		{
		case INTERFACE_DATA_TYPE :
		{
			Launch_Interface_Data_Type(io);
			break;
		}
		case INTERFACE_MULTIGENE :
		{
			Launch_Interface_Multigene(io);
			break;
		}
		case INTERFACE_MODEL :
		{
			Launch_Interface_Model(io);
			break;
		}
		case INTERFACE_TOPO_SEARCH :
		{
			Launch_Interface_Topo_Search(io);
			break;
		}
		case INTERFACE_BRANCH_SUPPORT :
		{
			Launch_Interface_Branch_Support(io);
			break;
		}
		case INTERFACE_MBL_MODEL :
		{
			Launch_Interface_MBL_Model(io);
			break;
		}
		default :
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Exit("");
			break;
		}
		}
	}while(!io->ready_to_go);


	if(io->in_tree == 2)
	{
		PhyML_Printf("\n. Enter the name of the input tree file > ");
		Getstring_Stdin(io->in_tree_file);
		io->fp_in_tree = Openfile(io->in_tree_file,0);
	}


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
		strcat(io->out_tree_file,".txt");
		strcat(io->out_stats_file,".txt");
	}

	io->fp_out_tree  = Openfile(io->out_tree_file,1);
	io->fp_out_stats = Openfile(io->out_stats_file,1);
}

/*********************************************************/

void Clear()
{
#ifdef WIN32
	system("cls");
#elif UNIX
	PhyML_Printf("\033[2J\033[H");
#endif
}

/*********************************************************/
void Launch_Interface_MBL_Model(option *io)
{
	char choix;
	char *n_branches;
	int i,n_trial;
	Clear();
	Print_Banner(stdout);

	PhyML_Printf("\n\n");

	PhyML_Printf("                        .................................                                              \n");
	PhyML_Printf("                         Menu : Mixed Branch Length Model                                               \n");
	PhyML_Printf("                       ...................................                                           \n");

	PhyML_Printf("\n\n");

	PhyML_Printf("                [+] "
			".................................... Next sub-menu\n");

	PhyML_Printf("                [-] "
			"................................ Previous sub-menu\n");

	PhyML_Printf("                [Y] "
			".............................. Launch the analysis\n");

	PhyML_Printf("\n");
	PhyML_Printf("                [N] "
			"............. Number of Branch lengths per edge:  "
			" %i \n", io->mod->n_l);
	if (io->mod->n_l > 1)
	{
		PhyML_Printf("                [P] "
				"Initial props in branch length set: [");
		For(i,io->mod->n_l){
			if(i+1 != io->mod->n_l){
				PhyML_Printf(" %lf,",io->mod->bl_props[i]);
			}else{
				PhyML_Printf(" %lf ",io->mod->bl_props[i]);
			}
		}
		PhyML_Printf("]\n");
		PhyML_Printf("                [F] "
				"...... Fixed Starting Proportions (Yes/No) "
				" %-15s \n",
				(io->fixed_props)?("Yes"):("No"));
	}

	PhyML_Printf("\n\n. Are these settings correct ? "
			"(type '+', '-', 'Y' or other letter for one to change)  ");


	if(!scanf("%c",&choix)) Exit("\n");
	if(choix != '\n') getchar(); /* \n */
	fflush(NULL);

	Uppercase(&choix);

	switch(choix)
	{
	/*     case '\n': */
	/*       { */
	/* 	io->curr_interface++; */
	/* 	break; */
	/*       } */
	case '+' :
	{
		io->curr_interface = (io->multigene)?(INTERFACE_MULTIGENE):(INTERFACE_MODEL);
		break;
	}
	case '-' :
	{	io->curr_interface = INTERFACE_DATA_TYPE;
	break;
	}
	case 'Y' :
	{
		io->ready_to_go = 1;
		break;
	}
	case 'N' :
	{
		Clear();
		Print_Banner(stdout);
		PhyML_Printf("\n\n");

		PhyML_Printf("                        .................................                                              \n");
		PhyML_Printf("                         Menu : Mixed Branch Length Model                                               \n");
		PhyML_Printf("                       ...................................                                           \n");

		PhyML_Printf("\n\n");


		PhyML_Printf("Number of Branch Length Sets > ");
		n_branches = (char *)mCalloc(10000,sizeof(char));
		Getstring_Stdin(n_branches);
		n_trial = 0;
		while((!atoi(n_branches)) || (atoi(n_branches) < 0))
		{
			if(++n_trial > 10) Exit("\n. Err : the number of branch length sets must be a positive integer.");
			PhyML_Printf("\n. The number of branch length sets must be a positive integer.");
			PhyML_Printf("\n. Enter a new value > ");
			Getstring_Stdin(n_branches);
		}
		io->mod->n_l = atoi(n_branches);
		if(io->mod->n_l == 1)
		{   //if one branch length category, then fix the branch length proportions (i.e. 100% for 0th branch length class)
			io->fixed_props = 1;
			io->mod->s_opt->opt_props = io->fixed_props;
		}
		else
		{
			io->fixed_props = 0;
			io->mod->s_opt->opt_props = 1;
			io->mod->s_opt->opt_five_branch = 0;
			//Default to Simulated Thermal Annealing if there is more than 1 branch length category
			// The user can change back to NNI, SPR, etc., but we *think* simulated annealing provides a
			// better optimization hueristic for complex likelihood landscapes.
			if(io->user_topo == 0)
			{
				io->mod->s_opt->topo_search = SIMULATED_THERMAL_ANNEALING;
			}
		}
		Free(n_branches);
		break;
	}
	case 'F':
	{
		if(io->mod->n_l == 1){
			io->fixed_props = 1;
			io->mod->s_opt->opt_props = io->fixed_props;
		} else {
			io->fixed_props = (io->fixed_props)?(0):(1);
			io->mod->s_opt->opt_props = 1;
		}
		break;
	}
	case 'P':
	{
		int i;
		int n_l = io->mod->n_l;

		For(i,n_l){
			Clear();
			Print_Banner(stdout);
			PhyML_Printf("\n\n");

			PhyML_Printf("                        .................................                                              \n");
			PhyML_Printf("                         Menu : Mixed Branch Length Model                                               \n");
			PhyML_Printf("                       ...................................                                           \n");

			PhyML_Printf("\n\n");
			PhyML_Printf("Old proportion for set %i: %lf \n",i,io->mod->bl_props[i]);
			PhyML_Printf("Enter a new proportion between 0.0 and 1.0 for set %i > ",i);
			n_branches = (char *)mCalloc(10000,sizeof(char));
			Getstring_Stdin(n_branches);
			n_trial = 0;
			while((!atof(n_branches)) || (atof(n_branches) <= 0.0) || (atof(n_branches) > 1.0))
			{
				if(++n_trial > 10) Exit("\n. Err : each site proportion must be a decimal greater than 0.0 and less than or equal to 1.0\n");
				PhyML_Printf("\n. The number of branch length sets must be a decimal value greater than 0.0 and less than or equal to 1.0\n");
				PhyML_Printf("\n. Enter a new proportion for branch length set number %i > ",i);
				Getstring_Stdin(n_branches);
			}
			io->mod->bl_props[i] = atof(n_branches);
			Free(n_branches);
		}
		io->user_props = 1;

		break;
	}

	default :
	{
		break;
	}
	}


}



/*********************************************************/

void Launch_Interface_Input(option *io)
{
	char choix;
	int n_trial;

	Clear();
	Print_Banner(stdout);


#ifdef EVOLVE

	char *n_data_sets;

	PhyML_Printf("\n\n");
	PhyML_Printf("\n. Enter the tree file name > "); fflush(NULL);
	Getstring_Stdin(io->in_tree_file);
	io->fp_in_tree = Openfile(io->in_tree_file,0);
	PhyML_Printf("\n");

	PhyML_Printf("\n. Enter the reference sequence file name > "); fflush(NULL);
	Getstring_Stdin(io->in_seq_file);
	io->fp_in_seq = Openfile(io->in_seq_file,0);
	PhyML_Printf("\n");

	PhyML_Printf("\n. Number of data sets > ");
	n_data_sets = (char *)mCalloc(10000,sizeof(char));
	Getstring_Stdin(n_data_sets);
	n_trial = 0;
	while((!atoi(n_data_sets)) || (atoi(n_data_sets) < 0))
	{
		if(++n_trial > 10) Exit("\n. Err : the number of sets must be a positive integer");
		PhyML_Printf("\n. The number of sets must be a positive integer");
		PhyML_Printf("\n. Enter a new value > ");
		Getstring_Stdin(n_data_sets);
	}
	io->n_data_set_asked = atoi(n_data_sets);
	Free(n_data_sets);

#elif OPTIMIZ

	PhyML_Printf("\n. Enter the tree file name > "); fflush(NULL);
	Getstring_Stdin(io->in_tree_file);
	io->fp_in_tree = Openfile(io->in_tree_file,0);
	PhyML_Printf("\n");

	PhyML_Printf("\n. Enter the reference sequence file name > "); fflush(NULL);
	Getstring_Stdin(io->in_seq_file);
	io->fp_in_seq = Openfile(io->in_seq_file,0);
	PhyML_Printf("\n");

#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)

	PhyML_Printf("\n. Enter the sequence file name > "); fflush(NULL);
	Getstring_Stdin(io->in_seq_file);
	io->fp_in_seq = Openfile(io->in_seq_file,0);

#endif


#if defined(PHYML) || defined(MG) || defined(PHYML_INSERT)

	strcpy(io->out_stats_file,io->in_seq_file);
	strcat(io->out_stats_file,"_phyml_stats.txt");

	strcpy(io->out_tree_file,io->in_seq_file);
	strcat(io->out_tree_file,"_phyml_tree.txt");

	strcpy(io->out_lk_file,io->in_seq_file);
	strcat(io->out_lk_file,"_phyml_lk.txt");


#endif


#ifdef WIN32
#ifdef EVOLVE
	if(Filexists("evolve_out.txt"));
#elif OPTIMIZ
	if(Filexists("optimiz_out.txt"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
		if(Filexists(io->out_stats_file))
#endif
#elif UNIX
#ifdef EVOLVE
			if(Filexists("evolve_out"));
#elif OPTIMIZ
	if(Filexists("optimiz_out"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
		if(Filexists(io->out_stats_file))
#endif
#endif
		{
			PhyML_Printf("\n");
#ifdef EVOLVE
			PhyML_Printf("\n. A file 'evolve_out' already exists");
#elif OPTIMIZ
			PhyML_Printf("\n. A file 'optimiz_out' already exists");
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
			PhyML_Printf("\n. A file '%s' already exists",io->out_stats_file);
#endif
			PhyML_Printf("\n. Do you want to Replace it or Append to it ? ");
			n_trial = 0;
			do
			{
				PhyML_Printf("\n. Please type R or A > ");
				if(!scanf("%c",&choix)) Exit("\n");
				if(choix == '\n') choix = 'r';
				else getchar();
				if(++n_trial>10) Exit("\n");
				Uppercase(&choix);
			}
			while((choix != 'R') && (choix != 'A'));
			if(choix == 'R') io->out_stats_file_open_mode = 1;
			else             io->out_stats_file_open_mode = 2;
		}

	io->fp_out_stats = Openfile(io->out_stats_file,io->out_stats_file_open_mode);

#ifdef WIN32
#ifdef EVOLVE
	if(Filexists("evolve_seq.txt"))
#elif OPTIMIZ
		if(Filexists("optimiz_tree.txt"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
			if(Filexists(io->out_tree_file))
#endif
#elif UNIX
#ifdef EVOLVE
				if(Filexists("evolve_seq"))
#elif OPTIMIZ
					if(Filexists("optimiz_tree"))
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
						if(Filexists(io->out_tree_file))
#endif
#endif
						{
							PhyML_Printf("\n");
#ifdef EVOLVE
							PhyML_Printf("\n. A file 'evolve_seq' already exists");
#elif OPTIMIZ
							PhyML_Printf("\n. A file 'optimiz_tree' already exists");
#elif defined(PHYML) || defined(MG) || defined(PHYML_INSERT)
							PhyML_Printf("\n. A file '%s' already exists",io->out_tree_file);
#endif
							PhyML_Printf("\n. Do you want to Replace it or Append to it ? ");
							n_trial = 0;
							do
							{
								PhyML_Printf("\n. Please type R or A > ");
								if(!scanf("%c",&choix)) Exit("\n");
								if(choix == '\n') choix = 'r';
								else getchar();
								Uppercase(&choix);
								if(++n_trial>10) Exit("\n");
							}
							while((choix != 'R') && (choix != 'A'));
							if(choix == 'R') io->out_tree_file_open_mode = 1;
							else             io->out_tree_file_open_mode = 2;
						}
	io->fp_out_tree = Openfile(io->out_tree_file,io->out_tree_file_open_mode);
}

/*********************************************************/

void Launch_Interface_Data_Type(option *io)
{
	char choix;
	char *s;
	char *buff;

	Clear();
	Print_Banner(stdout);
	if(io->config_multigene) Print_Data_Set_Number(io,stdout);

	s    = (char *)mCalloc(100,sizeof(char));
	buff = (char *)mCalloc(100,sizeof(char));


	PhyML_Printf("\n\n");

	PhyML_Printf("                                   ...................                                              \n");
	PhyML_Printf("                                    Menu : Input Data                                               \n");
	PhyML_Printf("                                .........................                                           \n");

	PhyML_Printf("\n\n");

	PhyML_Printf("                [+] "
			".................................... Next sub-menu\n");

	PhyML_Printf("                [-] "
			"................................ Previous sub-menu\n");

	PhyML_Printf("                [Y] "
			".............................. Launch the analysis\n");

	PhyML_Printf("\n");

	PhyML_Printf("                [D] "
			"............................... Data type (DNA/AA) "
			" %-15s \n",
			(io->mod->datatype)?("AA"):("DNA"));

	PhyML_Printf("                [I] "
			"...... Input sequences interleaved (or sequential) "
			" %-15s \n",
			(io->interleaved)?("interleaved"):("sequential"));

	strcpy(s,"");
	sprintf(s," (%d sets)",io->n_data_sets);
	strcpy(buff,(io->n_data_sets > 1)?("yes"):("no"));
	buff=strcat(buff,(io->n_data_sets > 1)?(s):("\0"));
	PhyML_Printf("                [M] "
			"....................... Analyze multiple data sets "
			" %-15s \n",buff);

	if(!io->append_run_ID) strcpy(s,"none");
	else strcpy(s,io->run_id_string);
	PhyML_Printf("                [R] "
			"........................................... Run ID "
			" %-15s \n",s);


	PhyML_Printf("\n\n. Are these settings correct ? "
			"(type '+', '-', 'Y' or other letter for one to change)  ");


	if(!scanf("%c",&choix)) Exit("\n");
	if(choix != '\n') getchar(); /* \n */
	fflush(NULL);

	Uppercase(&choix);

	switch(choix)
	{
	/*     case '\n': */
	/*       { */
	/* 	io->curr_interface++; */
	/* 	break; */
	/*       } */
	case 'R' :
	{
		io->append_run_ID = (io->append_run_ID)?(0):(1);
		PhyML_Printf("\n. Enter a run ID (any string of characters) > ");
		Getstring_Stdin(io->run_id_string);
		break;
	}


	case 'M' :
	{
		char *c;
		int n_trial;

		PhyML_Printf("\n. How many data sets > ");
		c = (char *)mCalloc(100,sizeof(char));
		Getstring_Stdin(c);
		n_trial = 0;
		while((!atoi(c)) || (atoi(c) < 0))
		{
			if(++n_trial > 10) Exit("\n. Err : The number of data sets must be a positive integer");
			PhyML_Printf("\n. The number of data sets must be a positive integer");
			PhyML_Printf("\n. Enter a new value > ");
			Getstring_Stdin(c);
		}
		io->n_data_sets = atoi(c);

#ifdef MG
		if(io->n_data_sets > 1)
		{
			io->multigene = 1;
		}
#endif

		if((io->mod->bootstrap > 1) && (io->n_data_sets > 1))
		{
			PhyML_Printf("\n. Bootstrap option is not allowed with multiple data sets !\n");
			PhyML_Printf("\n. Type any key to exit.");
			if(!scanf("%c",&choix)) Exit("\n");
			Exit("\n");
		}

		Free(c);
		break;
	}
	case 'I' :
	{
		if(io->interleaved)
			io->interleaved = 0;
		else io->interleaved = 1;
		break;
	}
	case 'D' :
	{
		if(io->mod->datatype == NT)
		{
			io->mod->datatype         = 1;
			io->mod->stepsize         = 1;
			io->mod->ns               = 20;
			io->mod->s_opt->opt_kappa = 0;
			io->mod->whichmodel       = LG;
			strcpy(io->mod->modelname,"LG");
		}
		else
		{
			io->mod->datatype         = 0;
			io->mod->stepsize         = 1;
			io->mod->ns               = 4;
			io->mod->whichmodel       = HKY85;
			strcpy(io->mod->modelname,"HKY85");
			strcpy(io->nt_or_cd,"nucleotides");
		}
		break;
	}
	case '-' :
	{
		io->curr_interface = (io->config_multigene)?(INTERFACE_MODEL):(INTERFACE_BRANCH_SUPPORT);
		break;
	}
	case '+' :
	{
		io->curr_interface = INTERFACE_MBL_MODEL;
		break;
	}
	case 'Y' :
	{
		io->ready_to_go = 1;
		break;
	}
	default :
	{
		break;
	}
	}

	Free(s);
	Free(buff);
}
/*********************************************************/

void Launch_Interface_Model(option *io)
{
	char choix;
	char *s;

	s = (char *)mCalloc(100,sizeof(char));

	Clear();
	Print_Banner(stdout);
	if(io->config_multigene) Print_Data_Set_Number(io,stdout);


	PhyML_Printf("\n\n");

	PhyML_Printf("                                ...........................                                      \n");
	PhyML_Printf("                                 Menu : Substitution Model                                       \n");
	PhyML_Printf("                             .................................                                   \n");

	PhyML_Printf("\n\n");

	PhyML_Printf("                [+] "
			".................................... Next sub-menu\n");

	PhyML_Printf("                [-] "
			"................................ Previous sub-menu\n");

	PhyML_Printf("                [Y] "
			".............................. Launch the analysis\n");

	PhyML_Printf("\n");

	if (io->mod->datatype == NT)
	{
		if(!strcmp(io->nt_or_cd,"nucleotides"))
		{
			PhyML_Printf("                [M] "
					"................. Model of nucleotide substitution "
					" %-15s \n", io->mod->modelname);

			if((io->mod->whichmodel == F81)   ||
					(io->mod->whichmodel == HKY85) ||
					(io->mod->whichmodel == F84)   ||
					(io->mod->whichmodel == TN93)  ||
					(io->mod->whichmodel == GTR)   ||
					(io->mod->whichmodel == CUSTOM))
			{
				/* 	      PhyML_Printf("                [F] " */
				/* 		     ".......... Base frequency estimates (empirical/ML) " */
				/* 		     " %-15s \n", */
				/* 		     (io->mod->s_opt->opt_state_freq)?("ML"):("empirical")); */
				/* 	    } */
				/* 	  else if(io->mod->whichmodel == CUSTOM) */
				/* 	    { */
				PhyML_Printf("                [F] "
						"................. Optimise equilibrium frequencies "
						" %-15s \n",
						(io->mod->s_opt->opt_state_freq)?("yes"):("no"));
			}


			if(io->mod->whichmodel == CUSTOM)
			{
				if(!io->mod->s_opt->opt_state_freq)
				{
					PhyML_Printf("                [E] "
							"......... Equilibrium frequencies (empirical/user) "
							" %-15s \n",
							(io->mod->s_opt->user_state_freq)?("user defined"):("empirical"));
				}
				PhyML_Printf("                [K] "
						"............................. Current custom model "
						" %-15s \n", io->mod->custom_mod_string);

				PhyML_Printf("                [O] "
						"................ Optimise relative rate parameters "
						" %-15s \n",(io->mod->s_opt->opt_rr)?("yes"):("no"));
			}
		}
	}
	else
	{
		PhyML_Printf("                [M] "
				"................ Model of amino-acids substitution "
				" %-15s \n", io->mod->modelname);

		PhyML_Printf("                [F] "
				". Amino acid frequencies (empirical/model defined) "
				" %-15s \n",
				(io->mod->s_opt->opt_state_freq)?("empirical"):("model"));
	}


	if ((io->mod->datatype    == NT)   &&
			((io->mod->whichmodel == K80)  ||
					(io->mod->whichmodel == HKY85)||
					(io->mod->whichmodel == F84)  ||
					(io->mod->whichmodel == TN93)))
	{
		strcpy(s,(io->mod->s_opt->opt_kappa)?("estimated"):("fixed"));
		(io->mod->s_opt->opt_kappa)?((char *)strcat(s,"")):((char *)strcat(s," (ts/tv = "));
		/*       (io->mod->s_opt->opt_kappa)?((char *)strcat(s,"")):((char *)sprintf(s+(int)strlen((char *)s),"%3.2f)",io->mod->kappa)); */
		if(io->mod->s_opt->opt_kappa)
		{
			strcat((char *)s,"");
		}
		else
		{
			sprintf((char *)(s+(int)strlen(s)),"%3.2f)",io->mod->kappa);
		}

		PhyML_Printf("                [T] "
				".................... Ts/tv ratio (fixed/estimated) "
				" %-15s \n",s);
	}


	(io->mod->s_opt->opt_pinvar)?(strcpy(s,"estimated")):(strcpy(s,"fixed"));
	(io->mod->s_opt->opt_pinvar)?((char *)strcat(s,"")):((char *)strcat(s," (p-invar = "));

	if(io->mod->s_opt->opt_pinvar)
	{
		strcat(s,"");
	}
	else
	{
		sprintf((char *)(s+strlen((char *)s)),"%3.2f)",io->mod->pinvar);
	}

	PhyML_Printf("                [V] "
			". Proportion of invariable sites (fixed/estimated)"
			"  %-15s \n",s);

	PhyML_Printf("                [R] "
			"....... One category of substitution rate (yes/no) "
			" %-15s \n",
			(io->mod->n_catg > 1)?("no"):("yes"));

	if(io->mod->n_catg > 1)
	{
		PhyML_Printf("                [C] "
				"........... Number of substitution rate categories "
				" %-15d \n",
				io->mod->n_catg);
	}


	if(io->mod->n_catg > 1)
	{
		strcpy(s,(io->mod->s_opt->opt_alpha)?("estimated"):("fixed"));
		(io->mod->s_opt->opt_alpha)?(strcat(s, "")):(strcat(s," (alpha = "));

		if(io->mod->s_opt->opt_alpha)
		{
			strcat(s, "");
		}
		else
		{
			sprintf(s+strlen(s),"%3.2f)",io->mod->alpha);
		}

		PhyML_Printf("                [A] "
				"... Gamma distribution parameter (fixed/estimated) "
				" %-15s \n",s);


		strcpy(s,(io->mod->gamma_median)?("median"):("mean"));

		PhyML_Printf("                [G] "
				".........'Middle' of each rate class (mean/median) "
				" %-15s \n",s);


	}

	PhyML_Printf("\n\n. Are these settings correct ? "
			"(type '+', '-', 'Y' or other letter for one to change)  ");


	if(!scanf("%c",&choix)) Exit("\n");
	if(choix != '\n') getchar(); /* \n */

	Uppercase(&choix);

	switch(choix)
	{
	/*     case '\n': */
	/*       { */
	/* 	io->curr_interface++; */
	/* 	break; */
	/*       } */

	case 'G' :
	{
		io->mod->gamma_median = (io->mod->gamma_median)?(0):(1);
		break;
	}

	case 'O' :
	{
		io->mod->s_opt->opt_rr = (io->mod->s_opt->opt_rr)?(0):(1);
		break;
	}

	case 'K' :
	{
		int i,j;
		char **rr_param,*rr;
		model *mod;
		int curr_param;
		int n_trial;

		if(io->mod->whichmodel == CUSTOM)
		{
			rr_param = (char **)mCalloc(6,sizeof(char *));
			For(i,6) rr_param[i] = (char *)mCalloc(10,sizeof(char));
			rr = (char *)mCalloc(50,sizeof(char));

			mod = io->mod;
			mod->s_opt->opt_rr = 1;

			n_trial = 0;
			do
			{
				PhyML_Printf("\n. Enter a new custom model > ");
				Getstring_Stdin(io->mod->custom_mod_string);
				if(strlen(io->mod->custom_mod_string) == 6)
				{
					For(i,6)
					{
						while(!isdigit((int)io->mod->custom_mod_string[i]))
						{
							if(++n_trial > 10) Exit("\n. Err : this string is not valid !\n");
							PhyML_Printf("\n. This string is not valid\n");
							PhyML_Printf("\n. Enter a new model > ");
							Getstring_Stdin(io->mod->custom_mod_string);
						}
					}
					if(i == 6) break;
				}
				else
				{
					PhyML_Printf("\n. The string should be of length 6\n");
					n_trial++;
				}
			}while(n_trial < 10);
			if(n_trial == 10) Exit("");

			Translate_Custom_Mod_String(io->mod);

			strcpy(rr_param[0],"A<->C");
			strcpy(rr_param[1],"A<->G");
			strcpy(rr_param[2],"A<->T");
			strcpy(rr_param[3],"C<->G");
			strcpy(rr_param[4],"C<->T");
			strcpy(rr_param[5],"G<->T");

			PhyML_Printf("\n. Set the relative rate values\n");
			curr_param = 0;
			For(i,mod->n_diff_rr)
			{
				sprintf(rr,"\n. [");
				For(j,6)
				{
					if(mod->rr_num[j] == i)
					{
						sprintf(rr+strlen(rr),"%s = ",rr_param[j]);
					}
				}
				sprintf(rr+strlen(rr)-3,"]");
				PhyML_Printf("%s",rr);

				PhyML_Printf("  (current=%.2f) > ",mod->rr_val[i]);

				Getstring_Stdin(rr);

				if(rr[0] != '\0')
				{
					n_trial = 0;
					while((atof(rr) < .0))
					{
						if(++n_trial > 10)
							Exit("\n. Err : the value of this parameter must be a positive number\n");
						PhyML_Printf("\n. The value of this parameter must be a positive number\n");
						PhyML_Printf("\n. Enter a new value > ");
						Getstring_Stdin(rr);
					}
					io->mod->rr_val[i] = (m3ldbl)atof(rr);
				}
			}

			For(i,5) Free(rr_param[i]);
			Free(rr_param);
			Free(rr);
		}
		break;
	}
	case 'E' :
	{
		int i;

		if(io->mod->whichmodel == CUSTOM)
		{
			io->mod->s_opt->user_state_freq = (io->mod->s_opt->user_state_freq)?(0):(1);

			if(io->mod->s_opt->user_state_freq)
			{
				if(!io->mod->s_opt->opt_state_freq)
				{
					char **bases;
					char *bs;
					m3ldbl sum;
					int n_trial;

					bases = (char **)mCalloc(4,sizeof(char *));
					For(i,4) bases[i] = (char *)mCalloc(50,sizeof(char));
					bs = (char *)mCalloc(100,sizeof(char));

					strcpy(bases[0],". f(A)> ");
					strcpy(bases[1],". f(C)> ");
					strcpy(bases[2],". f(G)> ");
					strcpy(bases[3],". f(T)> ");

					PhyML_Printf("\n. Set nucleotide frequencies \n");
					sum = .0;
					For(i,4)
					{
						PhyML_Printf("%s",bases[i]);
						Getstring_Stdin(bs);
						n_trial = 0;

						while((atof(bs) < .0001) || (bs[0] == '\0'))
						{
							if(++n_trial > 10)
								Exit("\n. Err : the value of this parameter must be a positive number\n");
							PhyML_Printf("\n. The value of this parameter must be a positive number\n");
							PhyML_Printf("\n. Enter a new value > ");
							Getstring_Stdin(bs);
						}
						io->mod->user_b_freq[i] = (m3ldbl)atof(bs);
						sum += io->mod->user_b_freq[i];
					}

					For(i,4) io->mod->user_b_freq[i] /= sum;

					if(sum > 1.0 || sum < 1.0)
					{
						PhyML_Printf("\n. The nucleotide frequencies have to be normalised in order to sum to 1.0.\n");
						PhyML_Printf("\n. The frequencies are now : f(A)=%f, f(C)=%f, f(G)=%f, f(T)=%f.\n",
								io->mod->user_b_freq[0],
								io->mod->user_b_freq[1],
								io->mod->user_b_freq[2],
								io->mod->user_b_freq[3]);
						PhyML_Printf("\n. Enter any key to continue.\n");
						if(!scanf("%c",bs)) Exit("\n");
					}

					For(i,4) Free(bases[i]);
					Free(bases);
					Free(bs);
				}
				else
				{
					Warn_And_Exit("\n. 'E' is not a valid option with these model settings.\n");
				}
			}
		}
		break;
	}
	case 'A' :
	{
		char answer;
		int n_trial;

		switch(io->mod->s_opt->opt_alpha)
		{
		case 0 :
		{
			PhyML_Printf("\n. Optimise alpha ? [Y/n] ");
			if(!scanf("%c",&answer)) Exit("\n");
			if(answer == '\n') answer = 'Y';
			else getchar();
			break;
		}
		case 1 :
		{
			PhyML_Printf("\n. Optimise alpha ? [N/y] ");
			if(!scanf("%c",&answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
			break;
		}
		default : Exit("\n");
		}

		n_trial = 0;
		while((answer != 'Y') && (answer != 'y') &&
				(answer != 'N') && (answer != 'n'))
		{
			if(++n_trial > 10) Exit("\n. Err : wrong answers !");
			PhyML_Printf("\n. Optimise alpha ? [N/y] ");
			if(!scanf("%c",&answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
		}

		switch(answer)
		{
		case 'Y' : case 'y' :
		{
			io->mod->s_opt->opt_alpha = 1;
			io->mod->s_opt->opt_num_param = 1;
			break;
		}
		case 'N' : case 'n' :
		{
			char *a;
			a = (char *)mCalloc(100,sizeof(char));
			io->mod->alpha = 10.0;
			io->mod->s_opt->opt_alpha = 0;
			PhyML_Printf("\n. Enter your value of alpha > ");
			Getstring_Stdin(a);
			n_trial = 0;
			while((!atof(a)) || (atof(a) < 1.E-10))
			{
				if(++n_trial > 10) Exit("\n. Err : alpha must be > 1.E-10\n");
				PhyML_Printf("\n. Alpha must be 1.E-10\n");
				PhyML_Printf("\n. Enter a new value > ");
				Getstring_Stdin(a);
			}
			io->mod->alpha = (m3ldbl)atof(a);
			Free(a);
			io->mod->s_opt->opt_alpha  = 0;
			break;
		}
		}
		break;
	}

	case 'C' :
	{
		char *c;
		int n_trial;

		PhyML_Printf("\n. Enter your number of categories > ");
		c = (char *)mCalloc(100,sizeof(char));
		Getstring_Stdin(c);
		n_trial = 0;
		while((!atoi(c)) || (atoi(c) < 0))
		{
			if(++n_trial > 10) Exit("\n. Err : the number of categories must be a positive integer\n");
			PhyML_Printf("\n. The number of categories must be a positive integer\n");
			PhyML_Printf("\n. Enter a new value > ");
			Getstring_Stdin(c);
		}
		io->mod->n_catg = atoi(c);
		Free(c);
		break;
	}

	case 'R' :
	{
		(io->mod->n_catg == 1)?(io->mod->n_catg = 4):(io->mod->n_catg = 1);
		break;
	}

	case 'V' :
	{
		char answer;
		int n_trial;

		switch(io->mod->s_opt->opt_pinvar)
		{
		case 0 :
		{
			PhyML_Printf("\n. Optimise p-invar ? [Y/n] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'Y';
			else getchar();
			break;
		}
		case 1 :
		{
			PhyML_Printf("\n. Optimise p-invar ? [N/y] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
			break;
		}
		default : Exit("\n");
		}

		n_trial = 0;
		while((answer != 'Y') && (answer != 'y') &&
				(answer != 'N') && (answer != 'n'))
		{
			if(++n_trial > 10) Exit("\n. Err : wrong answers !");
			PhyML_Printf("\n. Optimise p-invar ? [N/y] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
		}

		switch(answer)
		{
		case 'Y' : case 'y' :
		{
			io->mod->s_opt->opt_num_param = 1;
			io->mod->s_opt->opt_pinvar = 1;
			io->mod->pinvar = 0.2;
			io->mod->invar  = 1;
			break;
		}
		case 'N' : case 'n' :
		{
			char *p;
			p = (char *)mCalloc(100,sizeof(char));
			PhyML_Printf("\n. Enter your value of p-invar > ");
			Getstring_Stdin(p);
			n_trial = 0;
			while((atof(p) < 0.0) || (atof(p) > 1.0))
			{
				if(++n_trial > 10)
					Exit("\n. Err : the proportion of invariable sites must be a positive number between 0.0 and 1.0");
				PhyML_Printf("\n. The proportion must be a positive number between 0.0 and 1.0\n");
				PhyML_Printf("\n. Enter a new value > ");
				Getstring_Stdin(p);
			}
			io->mod->pinvar = (m3ldbl)atof(p);

			if(io->mod->pinvar > 0.0+MDBL_MIN) io->mod->invar = 1;
			else                             io->mod->invar = 0;

			Free(p);

			io->mod->s_opt->opt_pinvar = 0;
			break;
		}
		}
		break;
	}

	case 'T' :
	{
		char answer;
		int n_trial;

		if((io->mod->datatype   == AA)  ||
				(io->mod->whichmodel == JC69)||
				(io->mod->whichmodel == F81) ||
				(io->mod->whichmodel == GTR) ||
				(io->mod->whichmodel == CUSTOM))
		{
			PhyML_Printf("\n. 'K' is not a valid choice for this model\n");
			PhyML_Printf("\n. Type any key to exit.\n");
			if(!scanf("%c",&choix)) Exit("\n");
			Exit("\n");
		}

		switch(io->mod->s_opt->opt_kappa)
		{
		case 0 :
		{
			PhyML_Printf("\n. Optimise ts/tv ratio ? [Y/n] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'Y';
			else getchar();
			break;
		}
		case 1 :
		{
			PhyML_Printf("\n. Optimise ts/tv ratio ? [N/y] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
			break;
		}
		default : Exit("\n");
		}

		n_trial = 0;
		while((answer != 'Y') && (answer != 'y') &&
				(answer != 'N') && (answer != 'n'))
		{
			if(++n_trial > 10) Exit("\n. Err : wrong answers !");
			PhyML_Printf("\n. Optimise ts/tv ratio ? [N/y] ");
			if(!scanf("%c", &answer)) Exit("\n");
			if(answer == '\n') answer = 'N';
			else getchar();
		}

		switch(answer)
		{
		case 'Y' : case 'y' :
		{
			io->mod->kappa = 4.0;
			io->mod->s_opt->opt_num_param = 1;
			io->mod->s_opt->opt_kappa = 1;
			if(io->mod->whichmodel == TN93)
				io->mod->s_opt->opt_lambda = 1;
			break;
		}
		case 'N' : case 'n' :
		{
			char *t;
			t = (char *)mCalloc(100,sizeof(char));
			io->mod->s_opt->opt_kappa = 0;
			PhyML_Printf("\n. Enter your value of the ts/tv ratio > ");
			Getstring_Stdin(t);
			n_trial = 0;
			while((!atof(t)) || (atof(t) < .0))
			{
				if(++n_trial > 10) Exit("\n. Err : the ts/tv ratio must be a positive number\n");
				PhyML_Printf("\n. The ratio must be a positive number");
				PhyML_Printf("\n. Enter a new value > ");
				Getstring_Stdin(t);
			}
			io->mod->kappa = (m3ldbl)atof(t);
			io->mod->s_opt->opt_kappa  = 0;
			io->mod->s_opt->opt_lambda = 0;
			Free(t);
			break;
		}
		}
		break;
	}

	case 'F' :
	{
		if((io->mod->whichmodel == JC69) ||
				(io->mod->whichmodel == K80))
		{
			Warn_And_Exit("\n. 'F' is not a valid choice with these model settings.\n");
		}
		io->mod->s_opt->opt_state_freq = (io->mod->s_opt->opt_state_freq)?(0):(1);
		break;
	}

	case 'M' :
	{
		if(io->mod->datatype == NT)
		{
			if(!strcmp(io->nt_or_cd,"nucleotides"))
			{
				if(io->mod->whichmodel == JC69)
				{
					io->mod->whichmodel = K80;
					strcpy(io->mod->modelname,"K80");
				}
				else if(io->mod->whichmodel == K80)
				{
					io->mod->whichmodel = F81;
					strcpy(io->mod->modelname,"F81");
					io->mod->s_opt->opt_kappa  = 0;
					io->mod->s_opt->opt_lambda = 0;
				}
				else if(io->mod->whichmodel == F81)
				{
					io->mod->whichmodel = HKY85;
					strcpy(io->mod->modelname,"HKY85");
				}
				else if(io->mod->whichmodel == HKY85)
				{
					io->mod->whichmodel = F84;
					strcpy(io->mod->modelname,"F84");
				}
				else if(io->mod->whichmodel == F84)
				{
					io->mod->whichmodel = TN93;
					strcpy(io->mod->modelname,"TN93");
					if(io->mod->s_opt->opt_kappa) io->mod->s_opt->opt_lambda = 1;
				}
				else if(io->mod->whichmodel == TN93)
				{
					io->mod->whichmodel = GTR;
					strcpy(io->mod->modelname,"GTR");
					io->mod->s_opt->opt_kappa  = 0;
					io->mod->s_opt->opt_lambda = 0;
					io->mod->s_opt->opt_rr     = 1;
				}
				else if(io->mod->whichmodel == GTR)
				{
					io->mod->whichmodel = CUSTOM;
					strcpy(io->mod->modelname,"custom");
					io->mod->s_opt->opt_kappa  = 0;
					io->mod->s_opt->opt_lambda = 0;
				}

				else if(io->mod->whichmodel == CUSTOM)
				{
					io->mod->whichmodel = JC69;
					strcpy(io->mod->modelname,"JC69");
					io->mod->s_opt->opt_kappa  = 0;
					io->mod->s_opt->opt_lambda = 0;
				}
			}
		}
		else
		{
			if(io->mod->whichmodel == LG)
			{
				io->mod->whichmodel = WAG;
				strcpy(io->mod->modelname,"WAG");
			}
			else if(io->mod->whichmodel == WAG)
			{
				io->mod->whichmodel = DAYHOFF;
				strcpy(io->mod->modelname,"Dayhoff");
			}
			else if(io->mod->whichmodel == DAYHOFF)
			{
				io->mod->whichmodel = JTT;
				strcpy(io->mod->modelname,"JTT");
			}
			else if(io->mod->whichmodel == JTT)
			{
				io->mod->whichmodel = BLOSUM62;
				strcpy(io->mod->modelname,"Blossum62");
			}
			else if(io->mod->whichmodel == BLOSUM62)
			{
				io->mod->whichmodel = MTREV;
				strcpy(io->mod->modelname,"MtREV");
			}
			else if(io->mod->whichmodel == MTREV)
			{
				io->mod->whichmodel = RTREV;
				strcpy(io->mod->modelname,"RtREV");
			}
			else if(io->mod->whichmodel == RTREV)
			{
				io->mod->whichmodel = CPREV;
				strcpy(io->mod->modelname,"CpREV");
			}
			else if(io->mod->whichmodel == CPREV)
			{
				io->mod->whichmodel = DCMUT;
				strcpy(io->mod->modelname,"DCMut");
			}
			else if(io->mod->whichmodel == DCMUT)
			{
				io->mod->whichmodel = VT;
				strcpy(io->mod->modelname,"VT");
			}
			else if(io->mod->whichmodel == VT)
			{
				io->mod->whichmodel = MTMAM;
				strcpy(io->mod->modelname,"MtMam");
			}
			else if(io->mod->whichmodel == MTMAM)
			{
				io->mod->whichmodel = MTART;
				strcpy(io->mod->modelname,"MtART");
			}
			else if(io->mod->whichmodel == MTART)
			{
				io->mod->whichmodel = HIVW;
				strcpy(io->mod->modelname,"HIVw");
			}
			else if(io->mod->whichmodel == HIVW)
			{
				io->mod->whichmodel = HIVB;
				strcpy(io->mod->modelname,"HIVb");
			}
			else if(io->mod->whichmodel == HIVB)
			{
				io->mod->whichmodel = CUSTOMAA;
				strcpy(io->mod->modelname,"Read from file");
			}
			else if(io->mod->whichmodel == CUSTOMAA)
			{
				io->mod->whichmodel = LG;
				strcpy(io->mod->modelname,"LG");
			}
		}
		break;
	}
	case '-' :
	{
		io->curr_interface = INTERFACE_MBL_MODEL;
		break;
	}
	case '+' :
	{
		io->curr_interface = (io->config_multigene)?(INTERFACE_DATA_TYPE):(INTERFACE_TOPO_SEARCH);
		break;
	}
	case 'Y' :
	{
		io->ready_to_go = 1;
		break;
	}
	default :
	{
		break;
	}
	}

	if(io->mod->s_opt->opt_alpha  ||
			io->mod->s_opt->opt_kappa  ||
			io->mod->s_opt->opt_lambda ||
			io->mod->s_opt->opt_pinvar ||
			io->mod->s_opt->opt_rr) io->mod->s_opt->opt_num_param = 1;
	else                       io->mod->s_opt->opt_num_param = 0;





	Free(s);
}

/*********************************************************/

void Launch_Interface_Topo_Search(option *io)
{
	char choix;
	char *s;

	s = (char *)mCalloc(100,sizeof(char));

	Clear();
	Print_Banner(stdout);

	PhyML_Printf("\n\n");


	PhyML_Printf("                                   .......................                                     \n");
	PhyML_Printf("                                    Menu : Tree Searching                                        \n");
	PhyML_Printf("                                .............................                                  \n");


	PhyML_Printf("\n\n");

	PhyML_Printf("                [+] "
			".................................... Next sub-menu\n");

	PhyML_Printf("                [-] "
			"................................ Previous sub-menu\n");

	PhyML_Printf("                [Y] "
			".............................. Launch the analysis\n");

	PhyML_Printf("\n");

	PhyML_Printf("                [O] "
			"........................... Optimise tree topology "
			" %-15s \n",
			(io->mod->s_opt->opt_topo)?("yes"):("no"));


	if(io->mod->s_opt->opt_topo)
	{
		switch(io->in_tree)
		{
			case 0: { strcpy(s,"BioNJ");     break; }
			case 1: { strcpy(s,"parsimony"); break; }
			case 2: { strcpy(s,"user tree"); break; }
		}

		PhyML_Printf("                [U] "
				"........ Starting tree (BioNJ/parsimony/user tree) "
				" %-15s \n",s);
	}
	else
	{
		PhyML_Printf("                [U] "
				"..................... Input tree (BioNJ/user tree) "
				" %-15s \n",
				(!io->in_tree)?("BioNJ"):("user tree"));
	}

	if(io->mod->s_opt->opt_topo)
	{
		char *s;

		s = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
		if(io->mod->s_opt->topo_search == NNI_MOVE)
		{
			strcpy(s,"NNI moves (fast, approximate)\0");
		}
		else if(io->mod->s_opt->topo_search == SPR_MOVE)
		{
			strcpy(s,"SPR moves (slow, accurate)\0");
		}
		else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR)
		{
			strcpy(s,"Best of NNI and SPR \0");
		}else if(io->mod->s_opt->topo_search == SIMULATED_THERMAL_ANNEALING)
		{
			strcpy(s,"Simulated Thermal Annealing \0");
		}
		else if(io->mod->s_opt->topo_search == EMPIRICAL_BAYES)
		{
			strcpy(s,"Empirial Bayes MCMC \0");
		}
/*		else if(io->mod->s_opt->topo_search == SIMULATED_QUANTUM_ANNEALING)
		{
			strcpy(s,"Simulated Quantum Annealing \0");
		}
*/


		PhyML_Printf("                [S] "
				".................. Tree topology search operations "
				" %-15s \n",s);

		Free(s);

		if(io->mod->s_opt->topo_search != NNI_MOVE)
		{
			PhyML_Printf("                [R] "
					"........................ Add random starting trees "
					" %-15s \n",
					(io->mod->s_opt->random_input_tree)?("yes"):("no"));

			if(io->mod->s_opt->random_input_tree)
			{
				PhyML_Printf("                [N] "
						".................. Number of random starting trees "
						" %-15d \n",io->mod->s_opt->n_rand_starts);
			}
		}

		if(io->mod->s_opt->topo_search == EMPIRICAL_BAYES)
		{
			PhyML_Printf("                [G] "
					"........................ MCMC generations "
					" %-15d \n",io->eb_n_gens);
		}
	}
	else
	{
		PhyML_Printf("                [L] "
				".......................... Optimize branch lengths "
				" %-15s \n",
				(io->mod->s_opt->opt_bl)?("yes"):("no"));
	}

	PhyML_Printf("\n\n. Are these settings correct ? "
			"(type '+', '-', 'Y' or other letter for one to change)  ");


	if(!scanf("%c",&choix)) Exit("\n");
	if(choix != '\n') getchar(); /* \n */

	Free(s);

	Uppercase(&choix);


	switch(choix)
	{
	case '-' :
	{
		io->curr_interface = INTERFACE_MODEL;
		break;
	}
	case '+' :
	{
		io->curr_interface = INTERFACE_BRANCH_SUPPORT;
		break;
	}
	case 'Y' :
	{
		io->ready_to_go = 1;
		break;
	}
	case 'U' :
	{
		io->in_tree++;
		if(io->in_tree == 3) io->in_tree = 0;
		break;
	}
	case 'N' :
	{
		char *n;
		int n_trial;

		PhyML_Printf("\n. Enter your number of starting trees > ");
		n = (char *)mCalloc(100,sizeof(char));
		Getstring_Stdin(n);
		n_trial = 0;
		while(atoi(n) < 1)
		{
			if(++n_trial > 10) Exit("\n. Err : the number of starting trees must be a positive integer\n");
			PhyML_Printf("\n. The number of starting trees must be a positive integer\n");
			PhyML_Printf("\n. Enter a new value > ");
			Getstring_Stdin(n);
		}
		io->mod->s_opt->n_rand_starts = atoi(n);
		io->n_trees = 1;
		Free(n);
		break;
	}
	case 'O' :
	{
		io->mod->s_opt->opt_topo = (io->mod->s_opt->opt_topo)?(0):(1);
		break;
	}

	case 'L' :
	{
		if(!io->mod->s_opt->opt_topo)
		{
			io->mod->s_opt->opt_bl = (io->mod->s_opt->opt_bl)?(0):(1);
		}
		break;
	}

	case 'S' :
	{
		if(io->mod->s_opt->topo_search == NNI_MOVE)
		{
			io->mod->s_opt->topo_search         = SPR_MOVE;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
		else if(io->mod->s_opt->topo_search == SPR_MOVE)
		{
			io->mod->s_opt->topo_search         = BEST_OF_NNI_AND_SPR;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
		else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR)
		{
			io->mod->s_opt->topo_search         = SIMULATED_THERMAL_ANNEALING;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
		else if(io->mod->s_opt->topo_search == SIMULATED_THERMAL_ANNEALING)
		{
			io->mod->s_opt->topo_search         = EMPIRICAL_BAYES;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
		else if(io->mod->s_opt->topo_search == EMPIRICAL_BAYES){

			io->mod->s_opt->topo_search         = NNI_MOVE;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
/*		else if(io->mod->s_opt->topo_search == SIMULATED_QUANTUM_ANNEALING){

			io->mod->s_opt->topo_search         = EMPIRICAL_BAYES;
			io->mod->s_opt->n_rand_starts       = 1;
			io->user_topo 						= 1;
			io->mod->s_opt->random_input_tree   = 0;
			io->mod->s_opt->greedy              = 0;
		}
*/
		break;
	}
	case 'R' :
	{
		io->mod->s_opt->random_input_tree = (io->mod->s_opt->random_input_tree)?(0):(1);

		if(io->mod->s_opt->random_input_tree)
		{
			if(io->fp_in_tree) fclose(io->fp_in_tree);
			/* 	    io->in_tree                   = 0; */
			io->n_trees                   = 1;
			io->mod->s_opt->n_rand_starts = 5;

			strcpy(io->out_trees_file,io->in_seq_file);
			strcat(io->out_trees_file,"_phyml_trees.txt");
		}
		break;
	}
	case 'G' :
	{
		char *n;
		int n_trial;

		PhyML_Printf("\n. How many generations to run MCMC? > ");
		n = (char *)mCalloc(100,sizeof(char));
		Getstring_Stdin(n);
		n_trial = 0;
		while(atoi(n) < 1)
		{
			if(++n_trial > 10) Exit("\n. Err : the number of generations must be a positive integer\n");
			PhyML_Printf("\n. The number of generations must be a positive integer\n");
			PhyML_Printf("\n. Enter a new value > ");
			Getstring_Stdin(n);
		}
		io->eb_n_gens = atoi(n);
		Free(n);
		break;
	}
	default :
	{
		break;
	}
	}
}

/*********************************************************/

void Launch_Interface_Branch_Support(option *io)
{
	char choix;
	char *s;

	s = (char *)mCalloc(100,sizeof(char));

	Clear();
	Print_Banner(stdout);

	PhyML_Printf("\n\n");

	PhyML_Printf("                                  .......................                                          \n");
	PhyML_Printf("                                   Menu : Branch Support                                           \n");
	PhyML_Printf("                               .............................                                       \n");

	PhyML_Printf("\n\n");

	PhyML_Printf("                [+] "
			".................................... Next sub-menu\n");

	PhyML_Printf("                [-] "
			"................................ Previous sub-menu\n");

	PhyML_Printf("                [Y] "
			".............................. Launch the analysis\n");

	PhyML_Printf("\n");


	strcpy(s,(io->mod->bootstrap > 0)?("yes"):("no"));
	if(io->mod->bootstrap > 0) sprintf(s+strlen(s)," (%d replicate%s)",
			io->mod->bootstrap,
			(io->mod->bootstrap>1)?("s"):(""));

	/*   PhyML_Printf("                [+] " */
	PhyML_Printf("                [B] "
			"................ Non parametric bootstrap analysis "
			" %-15s \n",s);

	if(io->ratio_test == 0)
	{
		strcpy(s,"no");
	}
	else if(io->ratio_test == 1)
	{
		strcpy(s,"yes / aLRT statistics");
	}
	else if(io->ratio_test == 2)
	{
		strcpy(s,"yes / Chi2-based supports");
	}
	else if(io->ratio_test == 3)
	{
		strcpy(s,"yes / min of SH-like & Chi2-based supports");
	}
	else
	{
		strcpy(s,"yes / SH-like supports");
	}


	/*   PhyML_Printf("                [+] " */
	PhyML_Printf("                [A] "
			"................ Approximate likelihood ratio test "
			" %-15s \n",s);

	PhyML_Printf("\n. Are these settings correct ? "
			"(type '+', '-', 'Y' or other letter for one to change)  ");


	if(!scanf("%c",&choix)) Exit("\n");
	if(choix != '\n') getchar(); /* \n */

	Uppercase(&choix);

	Free(s);

	switch(choix)
	{
	/*     case '\n': */
	/*       { */
	/* 	io->curr_interface++; */
	/* 	break; */
	/*       } */

	case 'B' :
	{
		if(io->mod->bootstrap > 0) io->mod->bootstrap = 0;
		else
		{
			char *r;
			char answer;
			int n_trial;

			io->ratio_test = 0;

			if(io->n_data_sets > 1)
			{
				PhyML_Printf("\n. Bootstrap option is not allowed with multiple data sets.\n");
				PhyML_Printf("\n. Type any key to exit.\n");
				if(!scanf("%c",&choix)) Exit("\n");
				Exit("\n");
			}


			PhyML_Printf("\n. Number of replicates > ");
			r = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
			Getstring_Stdin(r);
			n_trial = 0;
			while((!atoi(r)) || (atoi(r) < 0))
			{
				if(++n_trial > 10) Exit("\n. Err : the number of replicates must be a positive integer\n");
				PhyML_Printf("\n. The number of replicates must be a positive integer");
				PhyML_Printf("\n. Enter a new value > ");
				Getstring_Stdin(r);
			}
			io->mod->bootstrap = atoi(r);

			PhyML_Printf("\n. Print bootstrap trees (and statistics) ? (%s) > ",
					(io->print_boot_trees)?("Y/n"):("y/N"));

			if(!scanf("%c",&answer)) Exit("\n");
			if(answer == '\n') answer = (io->print_boot_trees)?('Y'):('N');
			else getchar();

			switch(answer)
			{
			case 'Y' : case 'y' :
			{
				io->print_boot_trees = 1;

				strcpy(io->out_boot_tree_file,io->in_seq_file);
				strcat(io->out_boot_tree_file,"_phyml_boot_trees.txt");
				io->fp_out_boot_tree = Openfile(io->out_boot_tree_file,1);

				strcpy(io->out_boot_stats_file,io->in_seq_file);
				strcat(io->out_boot_stats_file,"_phyml_boot_stats.txt");
				io->fp_out_boot_stats = Openfile(io->out_boot_stats_file,1);

				break;
			}
			case 'N' : case 'n' :
			{
				io->print_boot_trees  = 0;
				io->fp_out_boot_tree  = NULL;
				io->fp_out_boot_stats = NULL;
				break;
			}
			}
			Free(r);
		}
		break;
	}
	case 'A' :
	{
		io->mod->bootstrap = 0;

		switch(io->ratio_test)
		{
		case 0 :
		{
			io->ratio_test = 1;
			break;
		}
		case 1 :
		{
			io->ratio_test = 2;
			break;
		}
		case 2 :
		{
			/* 	      io->ratio_test = 3; */
			io->ratio_test = 4;
			break;
		}
		case 3 :
		{
			io->ratio_test = 4;
			break;
		}
		case 4 :
		{
			io->ratio_test = 0;
			break;
		}
		}
		break;
	}

	case '-' :
	{
		io->curr_interface = INTERFACE_TOPO_SEARCH;
		break;
	}
	case '+' :
	{
		io->curr_interface = INTERFACE_DATA_TYPE;
		break;
	}
	case 'Y' :
	{
		io->ready_to_go = 1;
		break;
	}
	default : break;
	}
}

/*********************************************************/

void Launch_Interface_Multigene(option *io)
{

#ifdef MG

	if((io->n_data_sets > 1) && (io->multigene))
	{
		int set, n_trial;
		char choix;

		io->n_gt     = io->n_data_sets;
		io->st       = (superarbre *)Mg_Make_Superarbre_Light(io);
		io->st->n_gt = io->n_data_sets;

		For(set,io->n_gt)
		{
			io->st->optionlist[set] = Make_Input();
			Set_Defaults_Input(io->st->optionlist[set]);
			Set_Defaults_Model(io->st->optionlist[set]->mod);
			Set_Defaults_Optimiz(io->st->optionlist[set]->mod->s_opt);
			io->st->optionlist[set]->curr_gt = set;
			PhyML_Printf("\n. Enter the sequence file name [data set %2d] > ",set+1); fflush(NULL);
			Getstring_Stdin(io->st->optionlist[set]->in_seq_file);
			io->st->optionlist[set]->fp_in_seq = Openfile(io->st->optionlist[set]->in_seq_file,0);
		}


		PhyML_Printf("\n. Do you want to use your own initial tree ? [N/y] ");
		if(!scanf("%c", &choix)) Exit("\n");
		if(choix != '\n') getchar();

		n_trial = 0;
		while((choix != 'Y') && (choix != 'y') &&
				(choix != 'N') && (choix != 'n'))
		{
			if(++n_trial > 10) Exit("\n. Err : wrong answers !");
			PhyML_Printf("\n. Do you want to use your own initial tree ? [N/y] ? ");
			if(!scanf("%c", &choix)) Exit("\n");
			if(choix == '\n') choix = 'N';
			else getchar();
		}

		switch(choix)
		{
		case '\n' : break;
		case 'Y' : case 'y' :
		{
			io->in_tree = 1;
			PhyML_Printf("\n. Enter the name of the tree file > ");
			Getstring_Stdin(io->in_tree_file);
			io->fp_in_tree = Openfile(io->in_tree_file,0);
			break;
		}
		case 'N' : case 'n' :
		{
			io->in_tree = 0;
			break;
		}
		default : break;
		}

		io->curr_gt = 0;

		do
		{
			io->st->optionlist[io->curr_gt]->config_multigene = 1;
			do
			{
				switch(io->st->optionlist[io->curr_gt]->curr_interface)
				{
				case INTERFACE_DATA_TYPE :
				{
					Launch_Interface_Data_Type(io->st->optionlist[io->curr_gt]);
					break;
				}
				case INTERFACE_MODEL :
				{
					Launch_Interface_Model(io->st->optionlist[io->curr_gt]);
					break;
				}
				default :
				{
					PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
					Exit("");
					break;
				}
				}
			}while(!io->st->optionlist[io->curr_gt]->ready_to_go);
			io->curr_gt++;
		}while(io->curr_gt < io->n_gt);
	}
	io->ready_to_go = 1;
	Mg_Fill_Model_Partitions_Table(io->st);
#endif
}
