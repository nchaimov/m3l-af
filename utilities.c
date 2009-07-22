/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

 */

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "bionj.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"
#include "m4.h"
#include "mc.h"
#include "rates.h"
#include "numeric.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef MG
#include "mg.h"
#endif



/*********************************************************/

/* void Make_All_Edges_Light(node *a, node *d, int *curr_num_edge) */
/* { */
/*   int i; */

/*   Make_Edge_Light(a,d,*curr_num_edge); */
/*   (*curr_num_edge)++; */
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if(d->v[i] != a) */
/* 	    Make_All_Edges_Light(d,d->v[i],curr_num_edge); */
/* 	} */
/*     } */
/* } */

/*********************************************************/

void Make_All_Edges_Lk(node *a, node *d, arbre *tree)
{
	int i;

	For(i,3) if((a->v[i]) && (a->v[i] == d)) Make_Edge_Lk(a->b[i],tree);
	if(d->tax) return;
	else
	{
		For(i,3)
		{
			if(d->v[i] != a)
				Make_All_Edges_Lk(d,d->v[i],tree);
		}
	}
}

/*********************************************************/

arbre *Read_Tree(char *s_tree)
{
	char **subs;
	int i,n_ext,n_int,n_otu, n_l;
	arbre *tree;
	int degree;

	//JSJ: not good! need to modify so that commas are taxa or bl set deliniators
	n_l = 0;
	n_otu=0;
	//For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
	For(i,(int)strlen(s_tree)){
		if(s_tree[i] == '['){
			int tmp_nl = 0;
			for(;;){ //continue until close bracket
				i++;
				if(s_tree[i] == ','){
					tmp_nl++;
				} else if(s_tree[i] == ']'){
					if(n_l == 0){
						n_l = tmp_nl;
					}else if(n_l != tmp_nl){//we have a problem!!
						printf("JSJ: number of branch length sets not properly counted!\n");
						exit(1);
					}
					break; //done counting this round of branch_length sets
				}
			}
		}else if(s_tree[i] == ','){
			n_otu++;
		}
	}
	n_otu+=1;
	n_l +=1;


	//JSJ: Initialize memory and/or put in default values
	tree = (arbre *)Make_Tree(n_otu, n_l); //JSJ: Init_Tree called once here...
	Init_Tree(tree,tree->n_otu, tree->n_l); //JSJ: this is called twice... why??
	Make_All_Tree_Nodes(tree);
	Make_All_Tree_Edges(tree);
	Make_Tree_Path(tree); //just allocates memory
	Make_List_Of_Reachable_Tips(tree);

	//JSJ: not sure what this does...
	tree->noeud[n_otu]->num = n_otu;
	tree->noeud[n_otu]->tax = 0;
	//JSJ: added initialization of n_l
	tree->noeud[n_otu]->n_l = n_l;

	subs = Sub_Trees(s_tree,&degree); //degree is currently uninitialized
	Clean_Multifurcation(subs,degree,3);
	if(degree == 2) Unroot_Tree(subs);
	degree = 3;

	tree->has_branch_lengths = 0;
	tree->num_curr_branch_available = 0;
	n_int = n_ext = 0;
	For(i,degree) R_rtree(s_tree,subs[i],tree->noeud[n_otu],tree,&n_int,&n_ext);

	For(i,NODE_DEG_MAX) Free(subs[i]);
	Free(subs);
	return tree;
}

/*********************************************************/


/*********************************************************/
/* 'a' in node a stands for ancestor. 'd' stands for descendant */
void R_rtree(char *s_tree_a, char *s_tree_d, node *a, arbre *tree, int *n_int, int *n_ext)
{
	int i,j;
	node *d;
	int n_otu = tree->n_otu;

	if(strstr(s_tree_a," ")) Warn_And_Exit("\n Err : tree must not contain a ' ' character\n");

	if(s_tree_d[0] == '(')
	{
		char **subs;
		int degree;

		(*n_int)+=1;
		d      = tree->noeud[n_otu+*n_int];
		d->num = n_otu+*n_int;
		d->tax = 0;
		d->n_l = tree->n_l;

		Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
		Read_Branch_Lengths(s_tree_d,s_tree_a,tree);
		//JSJ: works now
		For(i,3)
		{
			if(!a->v[i])
			{
				//printf("JSJ: Entered utilities.c line 174\n");
				a->v[i]=d;
				For(j,tree->n_l){
					d->l[j][0]=tree->t_edges[tree->num_curr_branch_available]->l[j];
					a->l[j][i]=tree->t_edges[tree->num_curr_branch_available]->l[j];
					//printf("JSJ: The branch length here is: %f\n",tree->t_edges[tree->num_curr_branch_available]->l[j]);
				}
				break;
			}
		}
		d->v[0]=a;

		Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
		tree->num_curr_branch_available++;

		subs=Sub_Trees(s_tree_d,&degree);
		Clean_Multifurcation(subs,degree,2);
		R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
		R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
		For(i,NODE_DEG_MAX) Free(subs[i]);
		Free(subs);
	}

	else
	{
		int i,j;

		d      = tree->noeud[*n_ext];
		d->tax = 1;
		d->n_l = tree->n_l;

		Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
		Read_Branch_Lengths(s_tree_d,s_tree_a,tree);
		Read_Node_Name(d,s_tree_d,tree);

		For(i,3)
		{
			if(!a->v[i])
			{
				For(j,tree->n_l){
					//printf("JSJ: Entered utilities.c line 214\n");
					a->v[i]=d;
					d->l[j][0]=tree->t_edges[tree->num_curr_branch_available]->l[j];
					a->l[j][i]=tree->t_edges[tree->num_curr_branch_available]->l[j];
					//printf("JSJ: The branch length here is: %f\n",tree->t_edges[tree->num_curr_branch_available]->l[j]);
				}
				break;
			}
		}
		d->v[0]=a;

		Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
		tree->num_curr_branch_available++;

		d->num=*n_ext;
		(*n_ext)+=1;
	}
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, edge *b)
{
	char *sub_tp;
	char *p;
	int i,pos;

	sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	strcpy(sub_tp,s_d);
	strcat(sub_tp,"#");
	p = strstr(s_a,sub_tp);
	i = 0;
	b->n_labels = 0;
	if(p)
	{
		if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
		b->n_labels++;

		pos = 0;
		do
		{
			b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
			i++;
			pos++;
			if(p[i+strlen(s_d)+1] == '#')
			{
				b->labels[b->n_labels-1][pos] = '\0';
				b->n_labels++;
				if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
				i++;
				pos=0;
			}
		}
		while((p[i+strlen(s_d)+1] != ':') &&
				(p[i+strlen(s_d)+1] != ',') &&
				(p[i+strlen(s_d)+1] != '(') &&
				(p[i+strlen(s_d)+1] != '['));

		b->labels[b->n_labels-1][pos] = '\0';
	}

	if(p)
	{
		if(b->n_labels == 1)
			PhyML_Printf("\n. Found label '%s' on edge %3d.",b->labels[0],b->num);
		else
		{
			PhyML_Printf("\n. Found labels ");
			For(i,b->n_labels) PhyML_Printf("'%s' ",b->labels[i]);
			PhyML_Printf("on edge %3d.",b->num);
		}
	}

	Free(sub_tp);
}

/*********************************************************/
// JSJ: Modified so that it searches for the first open bracket as the beginning of the
// branch length set
void Read_Branch_Lengths(char *s_d, char *s_a, arbre *tree)
{
	char *sub_tp;
	//	char *p, *p2;
	char *p;
	char *brac;
	edge *b;
	int i,j;
	brac = strstr(s_a,"[");
	b = tree->t_edges[tree->num_curr_branch_available];
	b->n_l = tree->n_l; //not sure, but just in case...

	sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

	For(i,b->n_labels)
	{
		strcat(s_d,"#");
		strcat(s_d,b->labels[i]);
	}

	strcpy(sub_tp,s_d);
	if(brac){
		strcat(sub_tp,":[");
	}else{
		strcat(sub_tp,":");
	}
	//	printf("JSJ:Made it here!\n");
	//	//p = strstr(s_a,sub_tp);
	//printf("JSJ:The string about to be strstr: %s\n",s_a);
	//	//p2 = strstr((char *)p+(int)strlen(sub_tp),'[');
	p = strstr(s_a,sub_tp);
	if(p)
	{
		if(!brac){
			b->l[0] = atof((char *)p+(int)strlen(sub_tp));
		}else{
			//now we need to grab the string up until ',' and do it n_l times
			int glob = strlen(sub_tp);
			For(i,tree->n_l){
				char tmp[100];
				j = 0;
				glob++;//JSJ: increment glob past '[' and ','
				while ((p[glob] != ']') &&
						(p[glob] != ',')){
					tmp[j] = p[glob];
					j++;
					glob++;
				}
				tmp[j]='\0'; //JSJ: end the string
				b->l[i] = atof(tmp);
			}
		}
		//JSJ: Print all of the bls read in
		//For(i,b->n_l) printf("JSJ: Reading a branch from set %i of length %f\n",i,(double)b->l[i]);
		tree->has_branch_lengths = 1;
	}
	Free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(node *d, char *s_tree_d, arbre *tree)
{
	int i;

	if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
	{
		strcpy(d->name,s_tree_d);
	}
	else
	{
		i = 0;
		do
		{
			d->name[i] = s_tree_d[i];
			i++;
		}
		while(s_tree_d[i] != '#');
		d->name[i] = '\0';
	}
}
/*********************************************************/

void Unroot_Tree(char **subtrees)
{
	char **tmp_sub;
	int degree,i,j;

	PhyML_Printf("\n. Removing the root...\n");

	tmp_sub = Sub_Trees(subtrees[0],&degree);
	if(degree >= 2)
	{
		strcpy(subtrees[2],subtrees[1]);
		Clean_Multifurcation(tmp_sub,degree,2);
		For(j,2) strcpy(subtrees[j],tmp_sub[j]);
	}
	else
	{
		tmp_sub = Sub_Trees(subtrees[1],&degree);
		strcpy(subtrees[2],subtrees[0]);
		Clean_Multifurcation(tmp_sub,degree,2);
		For(j,2) strcpy(subtrees[j],tmp_sub[j]);
	}

	For(i,degree) Free(tmp_sub[i]);
	Free(tmp_sub);
}

/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

	if(current_deg <= end_deg) return;
	else
	{
		char *s_tmp;
		int i;

		s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

		strcat(s_tmp,"(\0");
		strcat(s_tmp,subtrees[0]);
		strcat(s_tmp,",\0");
		strcat(s_tmp,subtrees[1]);
		strcat(s_tmp,")\0");
		Free(subtrees[0]);
		subtrees[0] = s_tmp;

		for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

		Clean_Multifurcation(subtrees,current_deg-1,end_deg);
	}
}

/*********************************************************/
//modify so that it does not get tripped up over multiple branch lengths
char **Sub_Trees(char *tree, int *degree)
{
	char **subs;
	int posbeg,posend;
	int i;

	if(tree[0] != '(') {*degree = 1; return NULL;}

	subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

	For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));


	posbeg=posend=1;
	(*degree)=0;
	do
	{
		posbeg = posend;
		if(tree[posend] != '(')
		{
			while((tree[posend] != ',' ) &&
					(tree[posend] != ':' ) &&
					(tree[posend] != '#' ) &&
					(tree[posend] != ')' ))
			{
				//JSJ: skip over brackets...
				if (tree[posend] == '[') posend = Next_Brac(tree,posend);
				posend++ ;
			}
			posend -= 1;
		}
		else posend=Next_Par(tree,posend);

		while((tree[posend+1] != ',') &&
				(tree[posend+1] != ':') &&
				(tree[posend+1] != '#') &&
				(tree[posend+1] != ')')) {
			//JSJ: skip over brackets...
			if(tree[posend] == '[') posend = Next_Brac(tree,posend);
			posend++;

		}


		strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
		strcat(subs[(*degree)],"\0");

		posend += 1;
		while((tree[posend] != ',') &&
				(tree[posend] != ')')) {
			//JSJ: Skip over brackets
			if(tree[posend]== '[') posend = Next_Brac(tree,posend);
			posend++;

		}
		posend+=1;


		(*degree)++;
		if((*degree) == NODE_DEG_MAX)
		{
			For(i,(*degree))
			//JSJ: Changed print statement so it is sure to come out...
			//PhyML_Printf("\n. Subtree %d : %s\n",i+1,subs[i]);
			printf("\n. Subtree %d : %s\n",i+1,subs[i]);

			PhyML_Printf("\n. The degree of a node cannot be greater than %d\n",NODE_DEG_MAX);
			Warn_And_Exit("\n");
		}
	}
	while(tree[posend-1] != ')');

	return subs;
}


/*********************************************************/

int Next_Par(char *s, int pos)
{
	int curr;

	curr=pos+1;

	while(*(s+curr) != ')')
	{
		if(*(s+curr) == '(') curr=Next_Par(s,curr);
		curr++;
	}

	return curr;
}

//JSJ: Helper method to skip past brackets in tree string.
int Next_Brac(char *s, int pos)
{
	int curr;

	curr=pos+1;

	while(*(s+curr) != ']')
	{
		if(*(s+curr) == '['){
			printf("JSJ: nested brackets!!!! Bad!!!!");//curr=Next_Brac(s,curr);
			//this shouldn't ever happen!
			exit(1);
		}
		curr++;
	}

	return curr;
}

/*********************************************************/

void Print_Tree(FILE *fp, arbre *tree)
{
	char *s_tree;
	int i;

	s_tree = (char *)Write_Tree(tree);

	if(OUTPUT_TREE_FORMAT == 0) PhyML_Fprintf(fp,"%s\n",s_tree);
	else if(OUTPUT_TREE_FORMAT == 1)
	{
		PhyML_Fprintf(fp,"#NEXUS\n");
		PhyML_Fprintf(fp,"BEGIN TREES;\n");
		PhyML_Fprintf(fp,"\tTRANSLATE\n");
		For(i,tree->n_otu) PhyML_Fprintf(fp,"\t%3d\t%s,\n",i+1,tree->noeud[i]->name);
		PhyML_Fprintf(fp,"\tUTREE PAUP_1=\n");
		PhyML_Fprintf(fp,"%s\n",s_tree);
		PhyML_Fprintf(fp,"ENDBLOCK;");
	}
	Free(s_tree);
}

/*********************************************************/
void Print_Tree_Screen(arbre *tree){
	char *tree_s;
	tree_s = (char *)Write_Tree(tree);
	PhyML_Printf("%s\n",tree_s);
	Free(tree_s);
}

char *Write_Tree(arbre *tree)
{

	char *s;
	int i;

	s=(char *)mCalloc(T_MAX_LINE,sizeof(char));

	s[0]='(';

#ifdef PHYML
	tree->n_root = NULL;
	tree->e_root = NULL;
#endif


#ifdef MC
	if(tree->rates->bl_from_rt) RATES_Get_Br_Len(tree);
#endif

	if(!tree->n_root)
	{
		i = 0;
		while((!tree->noeud[tree->n_otu+i]->v[0]) ||
				(!tree->noeud[tree->n_otu+i]->v[1]) ||
				(!tree->noeud[tree->n_otu+i]->v[2])) i++;

		R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],s,tree);
		R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],s,tree);
		R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],s,tree);
	}
	else
	{
		R_wtree(tree->n_root,tree->n_root->v[0],s,tree);
		R_wtree(tree->n_root,tree->n_root->v[1],s,tree);
	}

	s[(int)strlen(s)-1]=')';
	s[(int)strlen(s)]=';';

	return s;
}

/*********************************************************/
// JSJ: this method will have to be dealt with...
void R_wtree(node *pere, node *fils, char *s_tree, arbre *tree)
{
	int i,p,j;

	p = -1;
	if(fils->tax)
	{
		if(OUTPUT_TREE_FORMAT == 0)
			strcat(s_tree,fils->name);
		else
			sprintf(s_tree+(int)strlen(s_tree),"%d",fils->num+1);
		//JSJ: temporary change of l
		if((fils->b) && (fils->b[0]) && (fils->b[0]->l[0] != -1))
		{
			if(tree->print_labels)
			{
				if(fils->b[0]->n_labels < 10)
					For(i,fils->b[0]->n_labels) sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[0]->labels[i]);
				else
					sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[0]->n_labels);
			}

			strcat(s_tree,":");

			if(pere != tree->n_root){
				if(tree->n_l == 1){
					sprintf(s_tree+(int)strlen(s_tree),"%.10f",fils->b[0]->l[0]);
				}else{
					sprintf(s_tree+(int)strlen(s_tree),"[");
					For(j,tree->n_l){
						if(j+1 == tree->n_l){
							sprintf(s_tree+(int)strlen(s_tree),"%.10f]",fils->b[0]->l[j]);
						}else{
							sprintf(s_tree+(int)strlen(s_tree),"%.10f,",fils->b[0]->l[j]);
						}
					}
				}
			}
			else
			{
				if(tree->n_root->v[0] == fils)
				{
					if(tree->n_l == 1){
						sprintf(s_tree+(int)strlen(s_tree),"%.10f",tree->n_root->l[0][0]);
					}else{
						sprintf(s_tree+(int)strlen(s_tree),"[");
						For(j,tree->n_l){
							if(j+1 == tree->n_l){
								sprintf(s_tree+(int)strlen(s_tree),"%.10f]",tree->n_root->l[j][0]);
							}else{
								sprintf(s_tree+(int)strlen(s_tree),"%.10f,",tree->n_root->l[j][0]);
							}
						}
					}
				}
				else
				{
					if(tree->n_l == 1){
						sprintf(s_tree+(int)strlen(s_tree),"%.10f",tree->n_root->l[0][1]);
					}else{
						sprintf(s_tree+(int)strlen(s_tree),"[");
						For(j,tree->n_l){
							if(j+1 == tree->n_l){
								sprintf(s_tree+(int)strlen(s_tree),"%.10f]",tree->n_root->l[j][1]);
							}else{
								sprintf(s_tree+(int)strlen(s_tree),"%.10f,",tree->n_root->l[j][1]);
							}
						}
					}
				}
			}
		}
		sprintf(s_tree+(int)strlen(s_tree),",");
	}
	else
	{
		s_tree[(int)strlen(s_tree)]='(';

		if(tree->n_root)
		{
			For(i,3)
			{
				if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
					R_wtree(fils,fils->v[i],s_tree,tree);
				else p=i;
			}
		}
		else
		{
			For(i,3)
			{
				if(fils->v[i] != pere)
					R_wtree(fils,fils->v[i],s_tree,tree);
				else p=i;
			}
		}

		s_tree[(int)strlen(s_tree)-1]=')';
		if((fils->b) && (fils->b[0]->l[0] != -1))
		{
			if(tree->print_boot_val)
				sprintf(s_tree+(int)strlen(s_tree),"%d",fils->b[p]->bip_score);
			else if(tree->print_alrt_val)
				sprintf(s_tree+(int)strlen(s_tree),"%.10f",fils->b[p]->ratio_test);

			if(tree->print_labels)
			{
				if(fils->b[p]->n_labels < 10)
					For(i,fils->b[p]->n_labels) sprintf(s_tree+(int)strlen(s_tree),"#%s",fils->b[p]->labels[i]);
				else
					sprintf(s_tree+(int)strlen(s_tree),"#%d_labels",fils->b[p]->n_labels);
			}

			strcat(s_tree,":");

			if(pere != tree->n_root){
				if(tree->n_l == 1){
					sprintf(s_tree+(int)strlen(s_tree),"%.10f",fils->b[p]->l[0]);
				}else{
					sprintf(s_tree+(int)strlen(s_tree),"[");
					For(j,tree->n_l){
						if(j+1 == tree->n_l){
							sprintf(s_tree+(int)strlen(s_tree),"%.10f]",fils->b[p]->l[j]);
						}else{
							sprintf(s_tree+(int)strlen(s_tree),"%.10f,",fils->b[p]->l[j]);
						}
					}
				}
			}
			else
			{
				if(tree->n_root->v[0] == fils)
				{
					if(tree->n_l == 1){
						sprintf(s_tree+(int)strlen(s_tree),"%.10f",tree->n_root->l[0][0]);
					}else{
						sprintf(s_tree+(int)strlen(s_tree),"[");
						For(j,tree->n_l){
							if(j+1 == tree->n_l){
								sprintf(s_tree+(int)strlen(s_tree),"%.10f]",tree->n_root->l[j][0]);
							}else{
								sprintf(s_tree+(int)strlen(s_tree),"%.10f,",tree->n_root->l[j][0]);
							}
						}
					}
				}
				else
				{
					if(tree->n_l == 1){
						sprintf(s_tree+(int)strlen(s_tree),"%.10f",tree->n_root->l[0][1]);
					}else{
						sprintf(s_tree+(int)strlen(s_tree),"[");
						For(j,tree->n_l){
							if(j+1 == tree->n_l){
								sprintf(s_tree+(int)strlen(s_tree),"%.10f]",tree->n_root->l[j][1]);
							}else{
								sprintf(s_tree+(int)strlen(s_tree),"%.10f,",tree->n_root->l[j][1]);
							}
						}
					}
				}
			}
		}
		strcat(s_tree,",");
	}
}

/*********************************************************/
//JSJ: No changes yet, but may have to change later to
// initialize the proportion of sites in each partition
void Init_Tree(arbre *tree, int n_otu, int n_l)
{
	int i;
	For(i,n_l) tree->props[i]       = 1.0/n_l; //JSJ: initialize with equal proportions of sites in each set
	tree->n_otu                     = n_otu;
	tree->n_l					    = n_l;
	tree->best_tree                 = NULL;
	tree->old_tree                  = NULL;
	tree->mat                       = NULL;
	tree->n_root                    = NULL;
	tree->e_root                    = NULL;
	tree->ps_tree                   = NULL;

	tree->depth_curr_path           = 0;
	tree->has_bip                   = 0;
	tree->n_moves                   = 0;
	tree->n_improvements            = 0;
	tree->number_of_lk_calls        = 0;
	tree->number_of_branch_lk_calls = 0;
	tree->bl_from_node_stamps       = 0;
	tree->lock_topo                 = 0;
	tree->ps_page_number            = 0;
	tree->init_lnL                  = UNLIKELY;
	tree->best_lnL                  = UNLIKELY;
	tree->c_lnL                     = UNLIKELY;
	tree->n_swap                    = 0;
	tree->best_pars                 = 1E+5;

	tree->n_pattern                 = -1;
	tree->prop_of_sites_to_consider = 1.;
	tree->n_root_pos                = -1.;
	tree->print_labels              = 1;

	tree->print_boot_val            = 0;
	tree->print_alrt_val            = 0;
	tree->num_curr_branch_available = 0;
	Normalize_Props(tree); //JSJ: Correct for floating point rounding error in props
}

/*********************************************************/
//JSJ: Function to normalize the site proportions of a tree to 1
//     This does a slight fudge to correct for floating point rounding error.
void Normalize_Props(arbre *tree){
	int i;
	m3ldbl sum = 0.;
	For(i,tree->n_l)sum           += tree->props[i]; //get the sum of props to normalize
	For(i,tree->n_l)tree->props[i] = tree->props[i] / sum; //normalize proportions to correct for float round error
}
//JSJ: Function to normalize the site proportions of an option struct to 1
//     This does a slight fudge to correct for floating point rounding error.
void Normalize_Props_IO(option *io){
	int i;
	m3ldbl sum = 0.;
	For(i,io->n_l)sum           += io->props[i]; //get the sum of props to normalize
	For(i,io->n_l)io->props[i] = io->props[i] / sum; //normalize proportions to correct for float round error
}

/*********************************************************/

void Make_New_Edge_Label(edge *b)
{
	int i;

	b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

	if(!b->labels)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	else
	{
		for(i=b->n_labels;i<b->n_labels+100;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
	}
}

/*********************************************************/
//JSJ: Just modified so that it also takes the n_l as an argument
edge *Make_Edge_Light(node *a, node *d, int num, int n_l)
{
	edge *b;
	int i;
	b = (edge *)mCalloc(1,sizeof(edge));
	b->n_l = n_l;

	Init_Edge_Light(b,num);

	if(a && b)
	{
		b->left = a;  b->rght = d;
		if(a->tax) {b->rght = a; b->left = d;} /* root */
		/* a tip is necessary on the right side of the edge */

		(b->left == a)?
				(Make_Edge_Dirs(b,a,d)):
					(Make_Edge_Dirs(b,d,a));
				For(i,n_l){
					b->l[i]                    = a->l[i][b->l_r];
					if(a->tax) b->l[i]         = a->l[i][b->r_l];
					if(b->l[i] < BL_MIN)  b->l[i] = BL_MIN;
					else if(b->l[i] > BL_MAX) b->l[i] = BL_MAX;
					b->l_old[i]                = b->l[i];
				}
	}
	else
	{
		b->left = NULL;
		b->rght = NULL;
	}

	return b;

}

/*********************************************************/

void Init_Edge_Light(edge *b, int num)
{
	int i;
	b->num                  = num;
	b->bip_score            = 0;
	b->dist_btw_edges       = .0;
	b->topo_dist_btw_edges  = 0;
	For(i,b->n_l) b->has_zero_br_len[i] = 0;
	b->is_p_lk_l_u2d        = 0;
	b->is_p_lk_r_u2d        = 0;
	b->n_jumps              = 0;

	b->p_lk_left            = NULL;
	b->p_lk_rght            = NULL;
	For(i,b->n_l)b->Pij_rr[i]= NULL;
}

/*********************************************************/

void Init_Node_Light(node *n, int num)
{
	n->list_of_reachable_tips = NULL;
	n->num                    = num;
	n->tax                    = -1;
	n->dist_to_root           = .0;
	n->common                 = 1;
}

/*********************************************************/

void Make_Edge_Dirs(edge *b, node *a, node *d)
{
	int i;

	if(a == b->rght)
	{
		PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	if(d == b->left)
	{
		PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	b->l_r = b->r_l = -1;
	For(i,3)
	{
		if((a->v[i]) && (a->v[i] == d))
		{
			b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
			a->b[i] = b;
		}
		if((d->v[i]) && (d->v[i] == a))
		{
			b->r_l  = i; /* we consider here that 'd' is on the right handside of 'b'*/
			d->b[i] = b;
		}
	}

	if(a->tax) {b->r_l = 0; For(i,3) if(d->v[i]==a) {b->l_r = i; break;}}

	b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
	For(i,3)
	{
		if(b->left->v[i] != b->rght)
		{
			if(b->l_v1 < 0) b->l_v1 = i;
			else            b->l_v2 = i;
		}

		if(b->rght->v[i] != b->left)
		{
			if(b->r_v1 < 0) b->r_v1 = i;
			else            b->r_v2 = i;
		}
	}
}

/*********************************************************/

void Make_Edge_Pars(edge *b, arbre *tree)
{
	/*   int site; */

	b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
	b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));

	b->ui_l = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
	b->ui_r = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));

	b->p_pars_l = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
	b->p_pars_r = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
}

/*********************************************************/

void Make_Edge_Lk(edge *b, arbre *tree)
{
	int i;
	For(i,tree->n_l){
		b->l_old[i] = b->l[i];
	}

	b->div_post_pred_left = (short int *)mCalloc((tree->mod->datatype == NT)?(4):(20),sizeof(short int));
	b->div_post_pred_rght = (short int *)mCalloc((tree->mod->datatype == NT)?(4):(20),sizeof(short int));
	//	b->Pij_rr             = (m3ldbl **)mCalloc(tree->n_l,sizeof(m3ldbl *));
	For(i,tree->n_l) b->Pij_rr[i] = (m3ldbl *)mCalloc(tree->mod->n_catg*tree->mod->ns*tree->mod->ns,sizeof(m3ldbl));

	b->scale_left = b->scale_rght = 0;

	if(!b->left->tax)
		b->sum_scale_f_left = (plkflt *)mCalloc(tree->data->crunch_len,sizeof(plkflt ));
	else
		b->sum_scale_f_left = NULL;

	if(!b->rght->tax)
		b->sum_scale_f_rght = (plkflt *)mCalloc(tree->data->crunch_len,sizeof(plkflt ));
	else
		b->sum_scale_f_rght = NULL;


	if((!b->left->tax) || (tree->mod->s_opt->greedy))
	{
		b->p_lk_left = (plkflt *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(plkflt));
		b->p_lk_tip_l = NULL;
	}
	else if(b->left->tax)
	{
		b->p_lk_left   = NULL;
		b->p_lk_tip_l  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int ));
	}

	if((!b->rght->tax) || (tree->mod->s_opt->greedy))
	{
		b->p_lk_rght = (plkflt *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(plkflt));
		b->p_lk_tip_r = NULL;
	}
	else if(b->rght->tax)
	{
		b->p_lk_rght = NULL;
		b->p_lk_tip_r  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int));
	}
}

/*********************************************************/

void Make_Edge_NNI(edge *b)
{
	b->nni    = Make_NNI(b->n_l);
	b->nni->b = b;
	b->nni->left = b->left;
	b->nni->rght = b->rght;
}

/*********************************************************/

nni *Make_NNI(int n_l)
{
	nni *a_nni;
	a_nni = (nni *)mCalloc(1,sizeof(nni ));
	Init_NNI(a_nni,n_l);
	return a_nni;
}

/*********************************************************/

void Init_NNI(nni *a_nni, int n_l)
{
	a_nni->left         = NULL;
	a_nni->rght         = NULL;
	a_nni->b            = NULL;
	a_nni->n_l			= n_l;
	a_nni->init_lk      = .0;
	a_nni->score        = +1.0;
	a_nni->swap_node_v1 = NULL;
	a_nni->swap_node_v2 = NULL;
	a_nni->swap_node_v3 = NULL;
	a_nni->swap_node_v4 = NULL;
	a_nni->lk0          = UNLIKELY;
	a_nni->lk1          = UNLIKELY;
	a_nni->lk2          = UNLIKELY;

	/**
	* JSJ: initialize arrays of branch lengths
	*/
	int i;
	For(i,n_l){
		a_nni->init_l[i] = -1.;
		a_nni->best_l[i] = -1.;
		a_nni->l0[i]     = -1.0;
		a_nni->l1[i]     = -1.0;
		a_nni->l2[i]     = -1.0;
	}
}

/*********************************************************/
////use this if the number of bl sets is not knonw, defaults to 1
//node *Make_Node_Light(int num)
//{
//  node *n;
//  n        = (node *)mCalloc(1,sizeof(node));
//  n->v     = (node **)mCalloc(3,sizeof(node *));
//  n->l     = (m3ldbl **)mCalloc(1,sizeof(m3ldbl*)); //JSJ: initialize into a single branch length set
//  n->l[0]  = (m3ldbl *)mCalloc(3,sizeof(m3ldbl)); //JSJ: if multiple, be sure to deallocate and reallocate
//  n->b     = (edge **)mCalloc(3,sizeof(edge *));
//  n->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
//  n->score = (m3ldbl *)mCalloc(3,sizeof(m3ldbl));
//  Init_Node_Light(n,num);
//  return n;
//}

/*********************************************************/
//JSJ: use this function if the number of bl sets is known
node *Make_Node_Light(int num, int num_bl_set)
{
//	int i;
	node *n;
	n        = (node *)mCalloc(1,sizeof(node));
	n->v     = (node **)mCalloc(3,sizeof(node *));
	n->b     = (edge **)mCalloc(3,sizeof(edge *));
	n->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	n->score = (m3ldbl *)mCalloc(3,sizeof(m3ldbl));
	n->n_l = num_bl_set;
	Init_Node_Light(n,num);
	return n;
}

/*********************************************************/


void Make_Node_Lk(node *n)
{
	/*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
	return;
}

/*********************************************************/

seq **Get_Seq(option *io,  int rw)
{
	seq **data;
	int i,j;
	char **buff;
	int n_unkn,n_removed,pos;
	int *remove;


	/*   rewind(fp_seq); */

	if(io->interleaved) data = Read_Seq_Interleaved(io->fp_in_seq,&(io->mod->n_otu));
	else                data = Read_Seq_Sequential(io->fp_in_seq,&(io->mod->n_otu));


	if(data)
	{
		buff = (char **)mCalloc(io->mod->n_otu,sizeof(char *));
		For(i,io->mod->n_otu) buff[i] = (char *)mCalloc(data[0]->len,sizeof(char));
		remove = (int *)mCalloc(data[0]->len,sizeof(int));

		n_removed = 0;

		For(i,data[0]->len)
		{
			For(j,io->mod->n_otu)
			{
				if((data[j]->state[i] == '?') || (data[j]->state[i] == '-')) data[j]->state[i] = 'X';
				if((io->mod->datatype == NT) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
				if(data[j]->state[i] == 'U') data[j]->state[i] = 'T';
			}

			n_unkn = 0;
			For(j,io->mod->n_otu) if(data[j]->state[i] == 'X') n_unkn++;

			if(n_unkn == io->mod->n_otu)
			{
				remove[i] = 1;
				n_removed++;
			}

			For(j,io->mod->n_otu) buff[j][i] = data[j]->state[i];
		}

		if(n_removed > 0)
		{
			if(io->mod->datatype == NT)
			{
				if(!io->quiet) PhyML_Printf("\n. %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n",n_removed);
			}
			else
			{
				if(!io->quiet) PhyML_Printf("\n. %d sites are made from completely undetermined states ('X', '-', '?')...\n",n_removed);
			}
		}

		pos = 0;
		For(i,data[0]->len)
		{
			/* 	  if(!remove[i]) */
			/* 	    { */
			For(j,io->mod->n_otu) data[j]->state[pos] = buff[j][i];
			pos++;
			/* 	    } */
		}

		For(i,io->mod->n_otu) data[i]->len = pos;
		For(i,io->mod->n_otu) Free(buff[i]);
		Free(buff);
		Free(remove);
	}
	return data;
}

/*********************************************************/

seq **Read_Seq_Sequential(FILE *in, int *n_otu)
{
	int i;
	char *line;
	int len,readok;
	seq **data;
	char c;
	char *format = (char *)mCalloc(T_MAX_NAME, sizeof(char));

	line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

	readok = len = 0;
	do
	{
		if(fscanf(in,"%s",line) == EOF)
		{
			Free(line); return NULL;
		}
		else
		{
			if(strcmp(line,"\n") && strcmp(line,"\n") && strcmp(line,"\t"))
			{
				*n_otu = atoi(line);
				data = (seq **)mCalloc(*n_otu,sizeof(seq *));
				if(*n_otu <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
				if(!fscanf(in,"%s",line)) Exit("\n");
				len = atoi(line);
				if(len <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
				else readok = 1;
			}
		}
	}while(!readok);


	/*   while((c=fgetc(in))!='\n'); */
	while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));

	For(i,*n_otu)
	{
		data[i] = (seq *)mCalloc(1,sizeof(seq));
		data[i]->len = 0;
		data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
		data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
		data[i]->is_ambigu = NULL;
		sprintf(format, "%%%ds", T_MAX_NAME);
		if(!fscanf(in, format, data[i]->name)) Exit("\n");

		while(data[i]->len < len) Read_One_Line_Seq(&data,i,in);

		if(data[i]->len != len)
		{
			PhyML_Printf("\n. Err: Problem with species %s's sequence (check the format)\n",
					data[i]->name);
			Warn_And_Exit("");
		}
	}

	Free(format);
	Free(line);
	return data;
}

/*********************************************************/

seq **Read_Seq_Interleaved(FILE *in, int *n_otu)
{
	int i,end,num_block;
	char *line;
	int len,readok;
	seq **data;
	char c;
	char *format;

	line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	format = (char *)mCalloc(T_MAX_NAME, sizeof(char));

	readok = len = 0;
	do
	{
		if(fscanf(in,"%s",line) == EOF)
		{
			Free(format);
			Free(line);
			return NULL;
		}
		else
		{
			if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
			{
				*n_otu = atoi(line);
				data = (seq **)mCalloc(*n_otu,sizeof(seq *));
				if(*n_otu <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
				if(!fscanf(in,"%s",line)) Exit("\n");
				len = atoi(line);
				if(len <= 0) Warn_And_Exit("\n. Problem with sequence format\n");
				else readok = 1;
			}
		}
	}while(!readok);


	while(((c=fgetc(in))!='\n') && (c != ' ') && (c != '\r') && (c != '\t'));


	end = 0;
	For(i,*n_otu)
	{
		data[i] = (seq *)mCalloc(1,sizeof(seq));
		data[i]->len = 0;
		data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
		data[i]->state = (char *)mCalloc(len+100,sizeof(char));
		data[i]->is_ambigu = NULL;
		sprintf(format, "%%%ds", T_MAX_NAME);
		/*       sprintf(format, "%%%ds", 10); */
		if(!fscanf(in, format, data[i]->name)) Exit("\n");
		if(!Read_One_Line_Seq(&data,i,in))
		{
			end = 1;
			if((i != *n_otu) && (i != *n_otu-1))
			{
				PhyML_Printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
				Warn_And_Exit("");
			}
			break;
		}
	}

	if(data[0]->len == len) end = 1;

	if(!end)
	{
		end = 0;

		num_block = 1;
		do
		{
			num_block++;

			/* interblock */
			if(!fgets(line,T_MAX_LINE,in)) break;

			if(line[0] != 13 && line[0] != 10)
			{
				PhyML_Printf("\n. One or more missing sequences in block %d\n",num_block-1);
				Warn_And_Exit("");
			}

			For(i,*n_otu)
			if(data[i]->len != len)
				break;

			if(i == *n_otu) break;


			For(i,*n_otu)
			{
				if(data[i]->len > len)
				{
					PhyML_Printf("\n. Observed length=%d expected length=%d\n",data[i]->len,len);
					PhyML_Printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
					Warn_And_Exit("");
				}
				else if(!Read_One_Line_Seq(&data,i,in))
				{
					end = 1;
					if((i != *n_otu) && (i != *n_otu-1))
					{
						PhyML_Printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
						Warn_And_Exit("");
					}
					break;
				}
			}
		}while(!end);
	}

	For(i,*n_otu)
	{
		if(data[i]->len != len)
		{
			PhyML_Printf("\n. Check sequence '%s' length...\n",data[i]->name);
			Warn_And_Exit("");
		}
	}

	Free(format);
	Free(line);
	return data;
}

/*********************************************************/

int Read_One_Line_Seq(seq ***data, int num_otu, FILE *in)
{
	char c;
	int nchar;

	nchar = 0;
	c=' ';
	while(1)
	{
		/*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */


		if((c == 13) || (c == 10))
		{
			/* 	  printf("[%d %d]\n",c,nchar); fflush(NULL); */
			if(!nchar)
			{
				c=(char)fgetc(in);
				continue;
			}
			else
			{
				/* 	      printf("break\n");  */
				break;
			}
		}
		else if(c == EOF)
		{
			/* 	  printf("EOL\n"); */
			break;
		}
		else if((c == ' ') || (c == '\t') || (c == 32))
		{
			/* 	  printf("[%d]",c); */
			c=(char)fgetc(in);
			continue;
		}

		nchar++;

		Uppercase(&c);


		if (strchr("ABCDEFGHIKLMNOPQRSTUVWXYZ?-.", c) == NULL)
		{
			PhyML_Printf("\n. Err: bad symbol: \"%c\" at position %d of species %s\n",
					c,(*data)[num_otu]->len,(*data)[num_otu]->name);
			Warn_And_Exit("");
		}

		if(c == '.')
		{
			c = (*data)[0]->state[(*data)[num_otu]->len];
			if(!num_otu) //num_otu == 0
				Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
		}
		(*data)[num_otu]->state[(*data)[num_otu]->len]=c;
		(*data)[num_otu]->len++;
		/*       printf("%c",c); */
		c = (char)fgetc(in);
	}

	if(c == EOF) return 0;
	else return 1;
}

/*********************************************************/

void Uppercase(char *ch)
{
	/* convert ch to upper case -- either ASCII or EBCDIC */
	*ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

/*********************************************************/

allseq *Compact_Seq(seq **data, option *io)
{
	allseq *alldata_tmp,*alldata;
	int i,j,k,site;
	int n_patt,which_patt,n_invar;
	char **sp_names;
	int n_otu, n_sites;
	pnode *proot;
	int compress;
	int n_ambigu,is_ambigu;

	n_otu        = io->mod->n_otu;
	n_patt       = 0;
	which_patt   = 0;

	sp_names = (char **)mCalloc(n_otu,sizeof(char *));
	For(i,n_otu)
	{
		sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
		strcpy(sp_names[i],data[i]->name);
	}

	alldata_tmp = Make_Cseq(n_otu,data[0]->len,data[0]->len,sp_names);
	proot       = (pnode *)Create_Pnode(T_MAX_ALPHABET);

	For(i,n_otu) Free(sp_names[i]);
	Free(sp_names);


	if(data[0]->len%io->mod->stepsize)
	{
		PhyML_Printf("\n. Sequence length is not a multiple of %d\n",io->mod->stepsize);
		Warn_And_Exit("");
	}

	compress = io->compress_seq;
	n_ambigu = 0;
	is_ambigu = 0;

	Fors(site,data[0]->len,io->mod->stepsize)
	{
		if(io->rm_ambigu)
		{
			is_ambigu = 0;
			For(j,n_otu)
			{
				if(Is_Ambigu(data[j]->state+site,io->mod->datatype,io->mod->stepsize)) break;
			}
			if(j != n_otu)
			{
				is_ambigu = 1;
				n_ambigu++;
			}
		}

		if(!is_ambigu)
		{
			if(compress)
			{
				which_patt = -1;
				Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,io,proot);
				if(which_patt == n_patt-1) /* New pattern found */
				{
					n_patt--;
					k=n_patt;
				}
				else
				{
					k = n_patt-10;
				}
			}
			else
			{
				PhyML_Printf("\n. WARNING: sequences are not compressed !");
				k = n_patt;
			}

			if(k == n_patt) /* add a new site pattern */
			{
				For(j,n_otu)
				Copy_One_State(data[j]->state+site,
						alldata_tmp->c_seq[j]->state+n_patt,
						io->mod->stepsize);


				For(i,n_otu)
				{
					For(j,n_otu)
					{
						if(!(Are_Compatible(alldata_tmp->c_seq[i]->state+n_patt,
								alldata_tmp->c_seq[j]->state+n_patt,
								io->mod->stepsize,
								io->mod->datatype))) break;
					}
					if(j != n_otu) break;
				}

				if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible -> the site is invariant */
				{
					For(j,n_otu)
					{
						alldata_tmp->invar[n_patt] = Assign_State(alldata_tmp->c_seq[j]->state+n_patt,
								io->mod->datatype,
								io->mod->stepsize);
						if(alldata_tmp->invar[n_patt] > -1.) break;
					}
				}
				else alldata_tmp->invar[n_patt] = -1;

				alldata_tmp->sitepatt[site] = n_patt;
				alldata_tmp->wght[n_patt]  += 1;
				n_patt                     += io->mod->stepsize;
			}
			else
			{
				alldata_tmp->sitepatt[site]    = which_patt;
				alldata_tmp->wght[which_patt] += 1;
			}
		}
	}

	data[0]->len -= n_ambigu;

	alldata_tmp->init_len                   = data[0]->len;
	alldata_tmp->crunch_len                 = n_patt;
	For(i,n_otu) alldata_tmp->c_seq[i]->len = n_patt;

	if(!io->quiet) PhyML_Printf("\n. %d patterns found. (out of a total of %d sites) \n",n_patt,data[0]->len);

	if((io->rm_ambigu) && (n_ambigu))
	{
		PhyML_Printf("\n. Removed %d columns of the alignment as the contain ambiguous characters (e.g., gaps) \n",n_ambigu);
	}

	n_invar=0;
	For(i,alldata_tmp->crunch_len) if(alldata_tmp->invar[i] > -1.) n_invar+=(int)alldata_tmp->wght[i];

	if(!io->quiet) PhyML_Printf("\n. %d sites without polymorphism (%.2f%c).\n",n_invar,100.*(m3ldbl)n_invar/data[0]->len,'%');

	alldata_tmp->obs_pinvar = (m3ldbl)n_invar/data[0]->len;

	n_sites = 0;
	For(i,alldata_tmp->crunch_len) n_sites += alldata_tmp->wght[i];
	if(n_sites != data[0]->len)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	if(io->mod->datatype == NT) Get_Base_Freqs(alldata_tmp);
	else                        Get_AA_Freqs(alldata_tmp);

	/*   PhyML_Fprintf(io->fp_out_stats,"\n. State frequencies: "); */
	/*   For(i,io->mod->ns) PhyML_Fprintf(io->fp_out_stats,"%f ",alldata_tmp->b_frq[i]); */
	/*   PhyML_Printf("\n"); */

	alldata = Copy_Cseq(alldata_tmp, alldata_tmp->crunch_len, io->mod->ns);

	Free_Cseq(alldata_tmp);
	Free_Prefix_Tree(proot,T_MAX_ALPHABET);

	return alldata;
}

/*********************************************************/

allseq *Compact_CSeq(allseq *data, model *mod)
{
	allseq *alldata;
	int i,j,k,site;
	int n_patt,which_patt;
	int n_otu;

	n_otu = data->n_otu;

	alldata         = (allseq *)mCalloc(1,sizeof(allseq));
	alldata->n_otu  = n_otu;
	alldata->c_seq  = (seq **)mCalloc(n_otu,sizeof(seq *));
	alldata->wght   = (int *)mCalloc(data->crunch_len,sizeof(int));
	alldata->b_frq  = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl));
	alldata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
	alldata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));

	alldata->crunch_len = alldata->init_len = -1;
	For(j,n_otu)
	{
		alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
		alldata->c_seq[j]->name      = (char *)mCalloc(T_MAX_NAME,sizeof(char));
		strcpy(alldata->c_seq[j]->name,data->c_seq[j]->name);
		alldata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
		alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
		alldata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
	}

	n_patt = which_patt =  0;

	Fors(site,data->crunch_len,mod->stepsize)
	{
		if(data->wght[site])
		{
			Fors(k,n_patt,mod->stepsize)
			{
				For(j,n_otu)
				{
					if(strncmp(alldata->c_seq[j]->state+k,
							data->c_seq[j]->state+site,
							mod->stepsize))
						break;
				}

				if(j == n_otu)
				{
					which_patt = k;
					break;
				}
			}

			/*       /\* TO DO *\/ */
			/*       k = n_patt; */

			if(k == n_patt)
			{
				For(j,n_otu) Copy_One_State(data->c_seq[j]->state+site,
						alldata->c_seq[j]->state+n_patt,
						mod->stepsize);

				For(i,n_otu)
				{
					For(j,n_otu)
					{
						if(!(Are_Compatible(alldata->c_seq[i]->state+n_patt,
								alldata->c_seq[j]->state+n_patt,
								mod->stepsize,
								mod->datatype))) break;
					}
					if(j != n_otu) break;
				}

				if((j == n_otu) && (i == n_otu))
				{
					For(j,n_otu)
					{
						alldata->invar[n_patt] = Assign_State(alldata->c_seq[j]->state+n_patt,
								mod->datatype,
								mod->stepsize);
						if(alldata->invar[n_patt] > -1.) break;
					}
				}
				else alldata->invar[n_patt] = -1;

				alldata->wght[n_patt] += data->wght[site];
				n_patt+=mod->stepsize;
			}
			else alldata->wght[which_patt] += data->wght[site];

			/*       Print_Site(alldata,k,n_otu,"\n",mod->stepsize); */
		}
	}

	alldata->init_len   = data->crunch_len;
	alldata->crunch_len = n_patt;
	For(i,n_otu) alldata->c_seq[i]->len = n_patt;

	(mod->datatype == NT)?
			(Get_Base_Freqs(alldata)):
				(Get_AA_Freqs(alldata));

			return alldata;
}

/*********************************************************/

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, seq **data, option *io, pnode *n)
{
	int ret_val;

	ret_val = -1;

	if(seqnum == io->mod->n_otu-1)
	{
		n->weight++;
		if(n->weight == 1)
		{
			n->num = *n_patt;
			(*n_patt) += 1;
		}
		(*patt_num) = n->num;
		return;
	}
	else
	{
		int next_state;

		next_state = -1;
		next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site,
				io->mod->datatype,
				io->mod->stepsize);

		if(!n->next[next_state]) n->next[next_state] = Create_Pnode(T_MAX_ALPHABET);
		Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,io,n->next[next_state]);
	}
}

/*********************************************************/

pnode *Create_Pnode(int size)
{
	pnode *n;
	int i;

	n = (pnode *)mCalloc(1,sizeof(pnode ));
	n->next = (pnode **)mCalloc(size,sizeof(pnode *));
	For(i,size) n->next[i] = NULL;
	n->weight = 0;
	n->num = -1;
	return n;
}
/*********************************************************/
/*********************************************************/

void Get_Base_Freqs(allseq *data)
{
	int i,j,k;
	m3ldbl A,C,G,T;
	m3ldbl fA,fC,fG,fT;
	int w;

	fA = fC = fG = fT = .25;

	For(k,8)
	{
		A = C = G = T = .0;
		For(i,data->n_otu)
		{
			For(j,data->crunch_len)
			{
				w = data->wght[j];
				if(w)
				{
					switch(data->c_seq[i]->state[j])
					{
					case 'A' : A+=w;
					break;
					case 'C' : C+=w;
					break;
					case 'G' : G+=w;
					break;
					case 'T' : T+=w;
					break;
					case 'U' : T+=w;
					break;
					case 'M' : C+=w*fC/(fC+fA); A+=w*fA/(fA+fC);
					break;
					case 'R' : G+=w*fG/(fA+fG); A+=w*fA/(fA+fG);
					break;
					case 'W' : T+=w*fT/(fA+fT); A+=w*fA/(fA+fT);
					break;
					case 'S' : C+=w*fC/(fC+fG); G+=w*fG/(fC+fG);
					break;
					case 'Y' : C+=w*fC/(fC+fT); T+=w*fT/(fT+fC);
					break;
					case 'K' : G+=w*fG/(fG+fT); T+=w*fT/(fT+fG);
					break;
					case 'B' : C+=w*fC/(fC+fG+fT); G+=w*fG/(fC+fG+fT); T+=w*fT/(fC+fG+fT);
					break;
					case 'D' : A+=w*fA/(fA+fG+fT); G+=w*fG/(fA+fG+fT); T+=w*fT/(fA+fG+fT);
					break;
					case 'H' : A+=w*fA/(fA+fC+fT); C+=w*fC/(fA+fC+fT); T+=w*fT/(fA+fC+fT);
					break;
					case 'V' : A+=w*fA/(fA+fC+fG); C+=w*fC/(fA+fC+fG); G+=w*fG/(fA+fC+fG);
					break;
					case 'N' : case 'X' : case '?' : case 'O' : case '-' :
						A+=w*fA; C+=w*fC; G+=w*fG; T+=w*fT; break;
					default : break;
					}
				}
			}
		}
		fA = A/(A+C+G+T);
		fC = C/(A+C+G+T);
		fG = G/(A+C+G+T);
		fT = T/(A+C+G+T);
	}

	data->b_frq[0] = fA;
	data->b_frq[1] = fC;
	data->b_frq[2] = fG;
	data->b_frq[3] = fT;
}

/*********************************************************/

void Get_AA_Freqs(allseq *data)
{
	int i,j,k;
	m3ldbl A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
	m3ldbl fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
	int w;
	m3ldbl sum;

	fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
			fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;

	For(k,8)
	{
		A = C = D = E = F = G = H = I = K = L =
				M = N = P = Q = R = S = T = V = W = Y = .0;

		For(i,data->n_otu)
		{
			For(j,data->crunch_len)
			{
				w = data->wght[j];
				if(w)
				{
					switch(data->c_seq[i]->state[j])
					{
					case 'A' : A+=w;		break;
					case 'C' : C+=w;		break;
					case 'D' : D+=w;		break;
					case 'E' : E+=w;		break;
					case 'F' : F+=w;		break;
					case 'G' : G+=w;		break;
					case 'H' : H+=w;		break;
					case 'I' : I+=w;		break;
					case 'K' : K+=w;		break;
					case 'L' : L+=w;		break;
					case 'M' : M+=w;		break;
					case 'N' : N+=w;		break;
					case 'P' : P+=w;		break;
					case 'Q' : Q+=w;		break;
					case 'R' : R+=w;		break;
					case 'S' : S+=w;		break;
					case 'T' : T+=w;		break;
					case 'V' : V+=w;		break;
					case 'W' : W+=w;		break;
					case 'Y' : Y+=w;		break;
					case 'Z' : Q+=w;		break;
					case 'X' : case '?' : case 'O' : case '-' :
						A+=w*fA;
						C+=w*fC;
						D+=w*fD;
						E+=w*fE;
						F+=w*fF;
						G+=w*fG;
						H+=w*fH;
						I+=w*fI;
						K+=w*fK;
						L+=w*fL;
						M+=w*fM;
						N+=w*fN;
						P+=w*fP;
						Q+=w*fQ;
						R+=w*fR;
						S+=w*fS;
						T+=w*fT;
						V+=w*fV;
						W+=w*fW;
						Y+=w*fY;
						break;
					default : break;
					}
				}
			}
		}
		sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
		fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
		fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
		fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
		fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
		fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
	}

	data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
	data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
	data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
	data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
	data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;
}

/*********************************************************/

arbre *Read_Tree_File(FILE *fp_input_tree)
{
	char *line;
	arbre *tree;
	int i;
	char c;

	line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

	do
	{
		c=fgetc(fp_input_tree);
	}
	while((c != '(') && (c != EOF));

	if(c==EOF)
	{
		Free(line);
		return NULL;
	}

	i=0;
	for(;;)
	{
		if((c == ' ') || (c == '\n'))
		{
			c=fgetc(fp_input_tree);
			if(c==EOF) break;
			else continue;
		}

		line[i]=c;
		i++;
		c=fgetc(fp_input_tree);
		if(c==EOF || c==';') break;
	}

	tree = Read_Tree(line);
	Free(line);
	return tree;
}

/*********************************************************/

void Connect_Edges_To_Nodes_Recur(node *a, node *d, arbre *tree)
{
	int i;

	Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
	tree->num_curr_branch_available += 1;

	if(d->tax) return;
	else For(i,3) if(d->v[i] != a) Connect_Edges_To_Nodes_Recur(d,d->v[i],tree);
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(node *a, node *d, edge *b, arbre *tree)
{
	int i,dir_a_d;
	dir_a_d = -1;
	For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}

	if(dir_a_d == -1)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}


	a->b[dir_a_d] = b;
	b->num        = tree->num_curr_branch_available;
	b->left       = a;
	b->rght       = d;
	if(a->tax) {b->rght = a; b->left = d;} /* root */
	/* a tip is necessary on the right hand side of the edge */

	(b->left == a)?
			(Make_Edge_Dirs(b,a,d)):
				(Make_Edge_Dirs(b,d,a));
			For(i,tree->n_l){ //JSJ: modified to add multiple lengths
				b->l[i]                    = a->l[i][b->l_r];
				if(a->tax) b->l[i]         = a->l[i][b->r_l];
				if(b->l[i] < BL_MIN)  b->l[i] = BL_MIN;
				else if(b->l[i] > BL_MAX) b->l[i] = BL_MAX;
				b->l_old[i]                = b->l[i];
				//printf("JSJ: Connecting node %s to node %s with an edge of length %f \n",a->name, d->name, b->l[i]);
			}
}

/*********************************************************/

void Update_Dirs(arbre *tree)
{
	int i;
	int buff;
	edge *b;

	b = NULL;
	buff = -1;
	For(i,2*tree->n_otu-3)
	{
		b = tree->t_edges[i];

		if((!b->left->tax) && (b->left->v[b->l_v1]->num < b->left->v[b->l_v2]->num))
		{
			buff    = b->l_v1;
			b->l_v1 = b->l_v2;
			b->l_v2 = buff;
		}
		if((!b->rght->tax) && (b->rght->v[b->r_v1]->num < b->rght->v[b->r_v2]->num))
		{
			buff    = b->r_v1;
			b->r_v1 = b->r_v2;
			b->r_v2 = buff;
		}
	}

}

/*********************************************************/

void Exit(char *message)
{
	fflush(NULL);
	PhyML_Fprintf(stderr,"%s",message);
	exit(1);
}

/*********************************************************/

void *mCalloc(int nb, size_t size)
{
	void *allocated;

	if((allocated = calloc((size_t)nb,(size_t)size)) != NULL)
	{
		return allocated;
	}
	else
		printf("Called mCalloc(nb= %i, size= %zu)\n",nb,size);
	Warn_And_Exit("\n. Err: low memory.\n mCalloc unsuccessful\n");

	return NULL;
}

/*********************************************************/

void *mRealloc(void *p,int nb, size_t size)
{
	if((p = realloc(p,(size_t)nb*size)) != NULL)
		return p;
	else
		printf("Called mRealloc(nb= %i, size= %zu)\n",nb,size);
	Warn_And_Exit("\n. Err: low memory\n mRealloc unsuccessful\n");

	return NULL;
}

/*********************************************************/

/* arbre *Make_Light_Tree_Struct(int n_otu) */
/* { */
/*   arbre *tree; */
/*   int i; */

/*   tree          = (arbre *)mCalloc(1,sizeof(arbre )); */
/*   tree->t_edges = (edge **)mCalloc(2*n_otu-3,sizeof(edge *)); */
/*   tree->noeud   = (node **)mCalloc(2*n_otu-2,sizeof(node *)); */
/*   tree->n_otu   = n_otu; */

/*   For(i,2*n_otu-3) */
/*     tree->t_edges[i] = Make_Edge_Light(NULL,NULL,i); */

/*   For(i,2*n_otu-2) */
/*     tree->noeud[i] = Make_Node_Light(i); */

/*   return tree; */
/* } */

/*********************************************************/

int Sort_Phydbl_Decrease(const void *a, const void *b)
{
	if((*(m3ldbl *)(a)) >= (*(m3ldbl *)(b))) return -1;
	else return 1;
}

/*********************************************************/
/* Sort in ascending order. Elements in B (if provided) are also re-ordered according to the ordering of A  */
void Qksort(m3ldbl *A, m3ldbl *B, int ilo, int ihi)
{
	m3ldbl pivot;	// pivot value for partitioning array
	int ulo, uhi;	// indices at ends of unpartitioned region
	int ieq;		// least index of array entry with value equal to pivot
	m3ldbl tempEntry;	// temporary entry used for swapping

	if (ilo >= ihi) {
		return;
	}
	// Select a pivot value.
	pivot = A[(ilo + ihi)/2];
	// Initialize ends of unpartitioned region and least index of entry
	// with value equal to pivot.
	ieq = ulo = ilo;
	uhi = ihi;
	// While the unpartitioned region is not empty, try to reduce its size.
	while (ulo <= uhi) {
		if (A[uhi] > pivot) {
			// Here, we can reduce the size of the unpartitioned region and
			// try again.
			uhi--;
		} else {
			// Here, A[uhi] <= pivot, so swap entries at indices ulo and
			// uhi.
			tempEntry = A[ulo];
			A[ulo]    = A[uhi];
			A[uhi]    = tempEntry;

			if(B)
			{
				tempEntry = B[ulo];
				B[ulo]    = B[uhi];
				B[uhi]    = tempEntry;
			}



			// After the swap, A[ulo] <= pivot.
			if (A[ulo] < pivot) {
				// Swap entries at indices ieq and ulo.
				tempEntry = A[ieq];
				A[ieq] = A[ulo];
				A[ulo] = tempEntry;


				if(B)
				{
					tempEntry = B[ieq];
					B[ieq] = B[ulo];
					B[ulo] = tempEntry;
				}


				// After the swap, A[ieq] < pivot, so we need to change
				// ieq.
				ieq++;
				// We also need to change ulo, but we also need to do
				// that when A[ulo] = pivot, so we do it after this if
				// statement.
			}
			// Once again, we can reduce the size of the unpartitioned
			// region and try again.
			ulo++;
		}
	}
	// Now, all entries from index ilo to ieq - 1 are less than the pivot
	// and all entries from index uhi to ihi + 1 are greater than the
	// pivot.  So we have two regions of the array that can be sorted
	// recursively to put all of the entries in order.
	Qksort(A, B, ilo, ieq - 1);
	Qksort(A, B, uhi + 1, ihi);
}

/********************************************************/

void Qksort_Matrix(m3ldbl **A, int col, int ilo, int ihi)
{
	m3ldbl pivot;	// pivot value for partitioning array
	int ulo, uhi;	// indices at ends of unpartitioned region
	int ieq;		// least index of array entry with value equal to pivot
	m3ldbl *tempEntry;	// temporary entry used for swapping

	tempEntry = NULL;

	if (ilo >= ihi) {
		return;
	}
	// Select a pivot value.
	pivot = A[(ilo + ihi)/2][col];
	// Initialize ends of unpartitioned region and least index of entry
	// with value equal to pivot.
	ieq = ulo = ilo;
	uhi = ihi;
	// While the unpartitioned region is not empty, try to reduce its size.
	while (ulo <= uhi) {
		if (A[uhi][col] > pivot) {
			// Here, we can reduce the size of the unpartitioned region and
			// try again.
			uhi--;
		} else {
			// Here, A[uhi] <= pivot, so swap entries at indices ulo and
			// uhi.
			tempEntry = A[ulo];
			A[ulo] = A[uhi];
			A[uhi] = tempEntry;
			// After the swap, A[ulo] <= pivot.
			if (A[ulo][col] < pivot) {
				// Swap entries at indices ieq and ulo.
				tempEntry = A[ieq];
				A[ieq] = A[ulo];
				A[ulo] = tempEntry;
				// After the swap, A[ieq] < pivot, so we need to change
				// ieq.
				ieq++;
				// We also need to change ulo, but we also need to do
				// that when A[ulo] = pivot, so we do it after this if
				// statement.
			}
			// Once again, we can reduce the size of the unpartitioned
			// region and try again.
			ulo++;
		}
	}
	// Now, all entries from index ilo to ieq - 1 are less than the pivot
	// and all entries from index uhi to ihi + 1 are greater than the
	// pivot.  So we have two regions of the array that can be sorted
	// recursively to put all of the entries in order.
	Qksort_Matrix(A, col, ilo, ieq - 1);
	Qksort_Matrix(A, col, uhi + 1, ihi);
}

/********************************************************/

void Print_Site(allseq *alldata, int num, int n_otu, char *sep, int stepsize)
{
	int i,j;
	For(i,n_otu)
	{
		PhyML_Printf("%s   ",alldata->c_seq[i]->name);
		For(j,stepsize)
		PhyML_Printf("%c",alldata->c_seq[i]->state[num+j]);
		PhyML_Printf("%s",sep);
	}
	PhyML_Fprintf(stderr,"%s",sep);
}

/*********************************************************/

void Print_Site_Lk(arbre *tree, FILE *fp)
{
	int site;
	int catg;
	char *s;
	m3ldbl postmean;

	if(!tree->io->print_site_lnl)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	if(!tree->io->print_trace)
	{
		s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

		PhyML_Fprintf(fp,"Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)\n");
		if(tree->mod->n_catg > 1 || tree->mod->invar)
			fprintf(fp,"P(D|M,rr[x]) is the probability of site D given the model M and the relative rate\nof evolution rr[x], where x is the class of rate to be considered.\nWe have P(D|M) = \\sum_x P(x) x P(D|M,rr[x]).\n");
		PhyML_Fprintf(fp,"\n\n");

		sprintf(s,"Site");
		PhyML_Fprintf(fp, "%-7s",s);

		sprintf(s,"P(D|M)");
		PhyML_Fprintf(fp,"%-16s",s);

		if(tree->mod->n_catg > 1)
		{
			For(catg,tree->mod->n_catg)
			{
				sprintf(s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->gamma_rr[catg]);
				PhyML_Fprintf(fp,"%-22s",s);
			}

			sprintf(s,"Posterior mean");
			PhyML_Fprintf(fp,"%-22s",s);
		}


		if(tree->mod->invar)
		{
			sprintf(s,"P(D|M,rr[0]=0)");
			PhyML_Fprintf(fp,"%-16s",s);
		}
		PhyML_Fprintf(fp,"\n");

		For(site,tree->data->init_len)
		{
			PhyML_Fprintf(fp,"%-7d",site+1);
			PhyML_Fprintf(fp,"%-16g",(m3ldbl)exp(tree->site_lk[tree->data->sitepatt[site]]));
			if(tree->mod->n_catg > 1)
			{
				For(catg,tree->mod->n_catg)
				fprintf(fp,"%-22g",(m3ldbl)exp(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]));

				postmean = .0;
				For(catg,tree->mod->n_catg)
				postmean +=
						tree->mod->gamma_rr[catg] *
						exp(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]) *
						tree->mod->gamma_r_proba[catg];
				postmean /= exp(tree->site_lk[tree->data->sitepatt[site]]);

				PhyML_Fprintf(fp,"%-22g",postmean);
			}
			if(tree->mod->invar)
			{
				if((m3ldbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5)
					fprintf(fp,"%-16g",tree->mod->pi[tree->data->invar[tree->data->sitepatt[site]]]);
				else
					fprintf(fp,"%-16g",0.0);
			}
			PhyML_Fprintf(fp,"\n");
		}
		Free(s);
	}
	else
	{
		For(site,tree->data->init_len)
		fprintf(fp,"%.2f\t",tree->site_lk[tree->data->sitepatt[site]]);
		PhyML_Fprintf(fp,"\n");
	}
}


/*********************************************************/

void Print_Seq(seq **data, int n_otu)
{
	int i,j;

	PhyML_Printf("%d\t%d\n",n_otu,data[0]->len);
	For(i,n_otu)
	{
		For(j,20)
		{
			if(j<(int)strlen(data[i]->name))
				putchar(data[i]->name[j]);
			else putchar(' ');
		}
		/*       PhyML_Printf("%10d  ",i); */
		For(j,data[i]->len)
		{
			PhyML_Printf("%c",data[i]->state[j]);
		}
		PhyML_Printf("\n");
	}
}

/*********************************************************/

void Print_CSeq(FILE *fp, allseq *alldata)
{
	int i,j,k;
	int n_otu;

	n_otu = alldata->n_otu;
	PhyML_Fprintf(fp,"%d\t%d\n",n_otu,alldata->init_len);
	For(i,n_otu)
	{
		For(j,50)
		{
			if(j<(int)strlen(alldata->c_seq[i]->name))
				fputc(alldata->c_seq[i]->name[j],fp);
			else fputc(' ',fp);
		}

		For(j,alldata->crunch_len)
		{
			For(k,alldata->wght[j])
			PhyML_Fprintf(fp,"%c",alldata->c_seq[i]->state[j]);
		}
		PhyML_Fprintf(fp,"\n");
	}
	PhyML_Fprintf(fp,"\n");

	/*   PhyML_Printf("\t"); */
	/*   For(j,alldata->crunch_len) */
	/*     PhyML_Printf("%.0f ",alldata->wght[j]); */
	/*   PhyML_Printf("\n"); */
}

/*********************************************************/

void Order_Tree_Seq(arbre *tree, seq **data)
{
	int i,j,n_otu;
	seq *buff;

	n_otu = tree->n_otu;

	For(i,n_otu)
	{
		For(j,n_otu)
		{
			if(!strcmp(tree->noeud[i]->name,data[j]->name))
				break;
		}
		buff = data[j];
		data[j] = data[i];
		data[i] = buff;
	}
}

/*********************************************************/

void Order_Tree_CSeq(arbre *tree, allseq *data)
{
	int i,j,n_otu_tree,n_otu_seq;
	seq *buff;


	n_otu_tree = tree->n_otu;
	n_otu_seq  = data->n_otu;


	if(n_otu_tree != n_otu_seq)
	{
		/*       PhyML_Printf("%d(tree) != %d(seq) \n",n_otu_tree,n_otu_seq); */
		Warn_And_Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");
	}

	For(i,MAX(n_otu_tree,n_otu_seq))
	{
		For(j,MIN(n_otu_tree,n_otu_seq))
		{
			if(!strcmp(tree->noeud[i]->name,data->c_seq[j]->name))
				break;
		}

		if(j==MIN(n_otu_tree,n_otu_seq))
		{
			PhyML_Printf("\n. Err: %s is not found in sequence data set\n",
					tree->noeud[i]->name);
			Warn_And_Exit("");
		}

		buff           = data->c_seq[j];
		data->c_seq[j] = data->c_seq[i];
		data->c_seq[i] = buff;
	}
}

/*********************************************************/

matrix *Make_Mat(int n_otu)
{
	matrix *mat;
	int i;

	mat = (matrix *)mCalloc(1,sizeof(matrix));

	mat->n_otu = n_otu;

	mat->P        = (m3ldbl **)mCalloc(n_otu,sizeof(m3ldbl *));
	mat->Q        = (m3ldbl **)mCalloc(n_otu,sizeof(m3ldbl *));
	mat->dist     = (m3ldbl **)mCalloc(n_otu,sizeof(m3ldbl *));
	mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
	mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
	mat->tip_node = (node **)mCalloc(n_otu,sizeof(node *));


	For(i,n_otu)
	{
		mat->P[i]    = (m3ldbl *)mCalloc(n_otu,sizeof(m3ldbl));
		mat->Q[i]    = (m3ldbl *)mCalloc(n_otu,sizeof(m3ldbl));
		mat->dist[i] = (m3ldbl *)mCalloc(n_otu,sizeof(m3ldbl));
		mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	}

	return mat;
}

/*********************************************************/

void Init_Mat(matrix *mat, allseq *data)
{
	int i;

	mat->n_otu = data->n_otu;
	mat->r = mat->n_otu;
	mat->curr_int = mat->n_otu;
	mat->method = 1;

	For(i,data->n_otu)
	{
		strcpy(mat->name[i],data->c_seq[i]->name);
		mat->on_off[i] = 1;
	}
}

/*********************************************************/

arbre *Make_Tree_From_Scratch(int n_otu, allseq *data, int n_l)
{
	arbre *tree;

	tree = Make_Tree(n_otu,n_l);
	Make_All_Tree_Nodes(tree);
	Make_All_Tree_Edges(tree);
	Make_Tree_Path(tree);
	Make_List_Of_Reachable_Tips(tree);
	if(data)
	{
		Copy_Tax_Names_To_Tip_Labels(tree,data);
		tree->data = data;
	}
	return tree;
}

/*********************************************************/

arbre *Make_Tree(int n_otu, int n_l)
{
	arbre *tree;
	int i;
	tree = (arbre *)mCalloc(1,sizeof(arbre ));
	tree->t_dir = (int **)mCalloc(2*n_otu-2,sizeof(int *)); //JSJ: not sure what t_dir does!!!
	For(i,2*n_otu-2) tree->t_dir[i] = (int *)mCalloc(2*n_otu-2,sizeof(int));
	Init_Tree(tree,n_otu, n_l);
	return tree;
}

/*********************************************************/

void Make_Tree_Path(arbre *tree)
{
	tree->curr_path = (node **)mCalloc(tree->n_otu,sizeof(node *));
}

/*********************************************************/

void Make_All_Tree_Nodes(arbre *tree)
{
	int i;

	tree->noeud          = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *));
	/*   tree->t_dead_nodes   = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *)); */

	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i] = (node *)Make_Node_Light(i, tree->n_l);
		if(i < tree->n_otu) tree->noeud[i]->tax = 1;
		else                tree->noeud[i]->tax = 0;
	}
}

/*********************************************************/

void Make_All_Tree_Edges(arbre *tree)
{
	int i;

	tree->t_edges      = (edge **)mCalloc(2*tree->n_otu-3,sizeof(edge *));
	/*   tree->t_dead_edges = (edge **)mCalloc(2*tree->n_otu-3,sizeof(edge *)); */

	For(i,2*tree->n_otu-3) tree->t_edges[i] = (edge *)Make_Edge_Light(NULL,NULL,i,tree->n_l);
}

/*********************************************************/

void Copy_Tax_Names_To_Tip_Labels(arbre *tree, allseq *data)
{
	int i;

	For(i,tree->n_otu)
	{
		strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
		tree->noeud[i]->tax = 1;
		tree->noeud[i]->num = i;
	}
}

/*********************************************************/

void Print_Dist(matrix *mat)
{
	int i,j;

	For(i,mat->n_otu)
	{
		PhyML_Printf("%s ",mat->name[i]);

		For(j,mat->n_otu)
		PhyML_Printf("%9.6f ",mat->dist[i][j]);
		PhyML_Printf("\n");
	}
}

/*********************************************************/

void Print_Node(node *a, node *d, arbre *tree)
{
	int i;
	int dir;
	dir = -1;
	For(i,3) if(a->v[i] == d) {dir = i; break;}
	PhyML_Printf("Node nums = %3d %3d  (dir=%d);",a->num,d->num,dir);
	PhyML_Printf("Node names = '%s' '%s' ; ",a->name,d->name);
	For(i,3) if(a->v[i] == d)
	{
		PhyML_Printf("Branch num = %3d (%d %d) %f",
				a->b[i]->num,a->b[i]->left->num,
				a->b[i]->rght->num,a->b[i]->l);
		if(a->b[i]->left->tax) PhyML_Printf(" WARNING LEFT->TAX!");
		break;
	}
	PhyML_Printf("\n");

	if(d->tax) return;
	else
		For(i,3)
		if(d->v[i] != a) Print_Node(d,d->v[i],tree);
}

/*********************************************************/

void Share_Lk_Struct(arbre *t_full, arbre *t_empt)
{
	int i,j,k,n_otu;
	edge *b_e,*b_f;
	node *n_e, *n_f;

	n_otu                   = t_full->n_otu;
	t_empt->n_root          = t_full->n_root;
	t_empt->e_root          = t_full->e_root;
	t_empt->c_lnL_sorted    = t_full->c_lnL_sorted;
	t_empt->log_site_lk_cat = t_full->log_site_lk_cat;
	t_empt->site_lk         = t_full->site_lk;
	t_empt->triplet_struct  = t_full->triplet_struct;
	t_empt->log_lks_aLRT    = t_full->log_lks_aLRT;

	For(i,2*n_otu-3)
	{
		b_f = t_full->t_edges[i];
		b_e = t_empt->t_edges[i];

		For(k,t_full->n_l){
			b_e->Pij_rr[k] = b_f->Pij_rr[k];
		}

		b_e->nni = b_f->nni;
	}


	for(i=n_otu;i<2*n_otu-2;i++)
	{
		n_f = t_full->noeud[i];
		n_e = t_empt->noeud[i];

		For(j,3)
		{
			if(n_f->b[j]->left == n_f)
			{
				if(n_e->b[j]->left == n_e)
				{
					n_e->b[j]->p_lk_left        = n_f->b[j]->p_lk_left;
					n_e->b[j]->sum_scale_f_left = n_f->b[j]->sum_scale_f_left;
					n_e->b[j]->p_lk_tip_l       = n_f->b[j]->p_lk_tip_l;
				}
				else
				{
					n_e->b[j]->p_lk_rght        = n_f->b[j]->p_lk_left;
					n_e->b[j]->sum_scale_f_rght = n_f->b[j]->sum_scale_f_left;
					n_e->b[j]->p_lk_tip_r       = n_f->b[j]->p_lk_tip_l;
				}
			}
			else
			{
				if(n_e->b[j]->rght == n_e)
				{
					n_e->b[j]->p_lk_rght        = n_f->b[j]->p_lk_rght;
					n_e->b[j]->sum_scale_f_rght = n_f->b[j]->sum_scale_f_rght;
					n_e->b[j]->p_lk_tip_r       = n_f->b[j]->p_lk_tip_r;
				}
				else
				{
					n_e->b[j]->p_lk_left        = n_f->b[j]->p_lk_rght;
					n_e->b[j]->sum_scale_f_left = n_f->b[j]->sum_scale_f_rght;
					n_e->b[j]->p_lk_tip_l       = n_f->b[j]->p_lk_tip_r;
				}
			}
		}
	}

	For(i,n_otu)
	{
		n_f = t_full->noeud[i];
		n_e = t_empt->noeud[i];

		if(n_f->b[0]->rght == n_f)
		{
			n_e->b[0]->p_lk_rght        = n_f->b[0]->p_lk_rght;
			n_e->b[0]->sum_scale_f_rght = n_f->b[0]->sum_scale_f_rght;
			n_e->b[0]->p_lk_tip_r       = n_f->b[0]->p_lk_tip_r;
		}
		else
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}
	}
}

/*********************************************************/

void Share_Spr_Struct(arbre *t_full, arbre *t_empt)
{
	t_empt->size_spr_list = t_full->size_spr_list;
	t_empt->spr_list      = t_full->spr_list;
	t_empt->best_spr      = t_full->best_spr;
}

/*********************************************************/

void Share_Pars_Struct(arbre *t_full, arbre *t_empt)
{
	int i;

	t_empt->site_pars = t_full->site_pars;
	t_empt->step_mat  = t_full->step_mat;

	For(i,2*t_full->n_otu-3)
	{
		t_empt->t_edges[i]->ui_l     = t_full->t_edges[i]->ui_l;
		t_empt->t_edges[i]->ui_r     = t_full->t_edges[i]->ui_r;

		t_empt->t_edges[i]->pars_l   = t_full->t_edges[i]->pars_l;
		t_empt->t_edges[i]->pars_r   = t_full->t_edges[i]->pars_r;

		t_empt->t_edges[i]->p_pars_l = t_full->t_edges[i]->p_pars_l;
		t_empt->t_edges[i]->p_pars_r = t_full->t_edges[i]->p_pars_r;
	}
}

/*********************************************************/

void Share_List_Of_Reachable_Tips_Struct(arbre *t_full, arbre *t_empt)
{
	int i;

	For(i,2*t_full->n_otu-2)
	{
		t_empt->noeud[i]->list_of_reachable_tips = t_full->noeud[i]->list_of_reachable_tips;
		t_empt->noeud[i]->n_of_reachable_tips    = t_full->noeud[i]->n_of_reachable_tips;
	}
}

/*********************************************************/

void Print_Mat(matrix *mat)
{
	int i,j;

	PhyML_Printf("%d",mat->n_otu);
	PhyML_Printf("\n");

	For(i,mat->n_otu)
	{
		For(j,13)
		{
			if(j>=(int)strlen(mat->name[i])) putchar(' ');
			else putchar(mat->name[i][j]);
		}

		For(j,mat->n_otu)
		{
			if(mat->dist[i][j] == -1)
				PhyML_Printf("   -     ");
			else
				PhyML_Printf("%7.8f  ",mat->dist[i][j]);
		}
		PhyML_Printf("\n");
	}
}

/*********************************************************/

int Sort_Edges_NNI_Score(arbre *tree, edge **sorted_edges, int n_elem)
{
	int i,j;
	edge *buff;

	For(i,n_elem-1)
	{
		for(j=i+1;j<n_elem;j++)
		{
			if(sorted_edges[j]->nni->score  < sorted_edges[i]->nni->score)
			{
				buff = sorted_edges[j];
				sorted_edges[j] = sorted_edges[i];
				sorted_edges[i] = buff;
			}
		}
	}
	return 1;
}

/*********************************************************/

int Sort_Edges_Depth(arbre *tree, edge **sorted_edges, int n_elem)
{
	int i,j;
	edge *buff;
	m3ldbl *depth,buff_depth;

	depth = (m3ldbl *)mCalloc(n_elem,sizeof(m3ldbl));

	For(i,n_elem)
	depth[i] =
			sorted_edges[i]->left->bip_size[sorted_edges[i]->l_r] *
			sorted_edges[i]->rght->bip_size[sorted_edges[i]->r_l] ;


	For(i,n_elem-1)
	{
		for(j=i+1;j<n_elem;j++)
		{
			if(depth[i] > depth[j])
			{
				buff = sorted_edges[i];
				sorted_edges[i] = sorted_edges[j];
				sorted_edges[j] = buff;

				buff_depth = depth[i];
				depth[i] = depth[j];
				depth[j] = buff_depth;
			}
		}
	}

	Free(depth);

	return 1;
}

/*********************************************************/
void NNI(arbre *tree, edge *b_fcus, int do_swap)
{
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
	node *v1,*v2,*v3,*v4;
	m3ldbl lk0, lk1, lk2;
	m3ldbl lk0_init, lk1_init, lk2_init;
	/*   m3ldbl lk_infa, lk_infb, lk_max; */
	m3ldbl lk_init;
	int i;

	m3ldbl bl_init[MAX_BL_SET];
	m3ldbl l0[MAX_BL_SET];
	m3ldbl l1[MAX_BL_SET];
	m3ldbl l2[MAX_BL_SET];
	m3ldbl l_infa[MAX_BL_SET];
	m3ldbl l_infb[MAX_BL_SET];
	m3ldbl l_max[MAX_BL_SET];



	lk_init                = tree->c_lnL;
	For(i,tree->n_l){
		bl_init[i]                = b_fcus->l[i];
		b_fcus->nni->init_l[i]    = b_fcus->l[i]; //JSJ: Fixed!
	}
	b_fcus->nni->init_lk  = tree->c_lnL;;

	b_fcus->nni->best_conf = 0;
	b_fcus->nni->score     = +1.0;

	lk0 = lk1 = lk2        = UNLIKELY;
	v1 = v2 = v3 = v4      = NULL;

	l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

	l_r                    = b_fcus->l_r;
	r_l                    = b_fcus->r_l;

	v1                     = b_fcus->left->v[b_fcus->l_v1];
	v2                     = b_fcus->left->v[b_fcus->l_v2];
	v3                     = b_fcus->rght->v[b_fcus->r_v1];
	v4                     = b_fcus->rght->v[b_fcus->r_v2];

	if(v1->num < v2->num)
	{

		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	if(v3->num < v4->num)
	{

		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	For(i,tree->n_l) l0[i] = l1[i] = l2[i] = -1.;


	/***********/
	Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
	tree->both_sides = 1;

	lk1_init = Update_Lk_At_Given_Edge(b_fcus,tree);
	//	For(i,tree->n_l){
	//		l_infa[i] = 10.*b_fcus->l[i]; //JSJ: Changed so it compiles
	//		l_max[i]  = b_fcus->l[i]; //JSJ: Changed so it compiles
	//		l_infb[i] = BL_MIN;
	//	}

	if(tree->mod->s_opt->fast_nni)
	{
		Fast_Br_Len(b_fcus,tree,1);
		lk1 = Lk_At_Given_Edge(b_fcus,tree);
	}
	else
	{
		For(i,tree->n_l){
		//	PhyML_Printf("Tree lnL (line 2983 of utilities.c): %lf \n",tree->c_lnL);
			lk1 = Br_Len_Brent_Iter(10.*b_fcus->l[i],b_fcus->l[i],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,
					tree->mod->s_opt->quickdirty,i);
		}
	}

	if(lk1 < lk1_init - tree->mod->s_opt->min_diff_lk_local)
	{
		PhyML_Printf("%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
		PhyML_Printf("%f -- %f \n",lk1_init,lk1);
		PhyML_Printf("\n. Err. in NNI (1)\n");
	}

	For(i,tree->n_l) l1[i]  = b_fcus->l[i];
	Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
	/***********/


	/***********/
	Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	For(i,tree->n_l) b_fcus->l[i] = bl_init[i];
	tree->both_sides = 1;

	lk2_init = Update_Lk_At_Given_Edge(b_fcus,tree);
	//	For(i,tree->n_l){
	//		l_infa[i] = 10.*b_fcus->l[i];
	//		l_max[i]  = b_fcus->l[i];
	//		l_infb[i] = BL_MIN;
	//	}

	if(tree->mod->s_opt->fast_nni)
	{
		Fast_Br_Len(b_fcus,tree,1);
		lk2 = Lk_At_Given_Edge(b_fcus,tree);
	}
	else
	{
		For(i,tree->n_l){
		//	PhyML_Printf("Tree lnL (line 3024 of utilities.c): %lf \n",tree->c_lnL);
			lk2 = Br_Len_Brent_Iter(10.*b_fcus->l[i],b_fcus->l[i],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,
					tree->mod->s_opt->quickdirty,i);
		}
	}

	if(lk2 < lk2_init - tree->mod->s_opt->min_diff_lk_local)
	{
		PhyML_Printf("%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
		PhyML_Printf("%f -- %f \n",lk2_init,lk2);
		PhyML_Printf("\n. Err. in NNI (2)\n");
	}

	For(i,tree->n_l) l2[i]  = b_fcus->l[i];
	Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
	/***********/



	/***********/
	For(i,tree->n_l) b_fcus->l[i] = bl_init[i];
	tree->both_sides = 1;

	lk0_init = Update_Lk_At_Given_Edge(b_fcus,tree);

	if(fabs(lk0_init - lk_init) > tree->mod->s_opt->min_diff_lk_local)
	{

		PhyML_Printf("\n. lk_init = %f; lk = %f diff = %f\n",
				lk_init,
				lk0_init,
				lk_init-lk0_init);
		PhyML_Printf("\n. Curr_lnL = %f\n",Return_Lk(tree));
		//Warn_And_Exit("\n. Err. in NNI (3)\n");
	}
	//	For(i,tree->n_l){
	//		l_infa[i] = 10.*b_fcus->l[i];
	//		l_max[i]  = b_fcus->l[i];
	//		l_infb[i] = BL_MIN;
	//	}

	if(tree->mod->s_opt->fast_nni)
	{
		Fast_Br_Len(b_fcus,tree,1);
		lk0 = Lk_At_Given_Edge(b_fcus,tree);
	}
	else
	{
		For(i,tree->n_l){
		//	PhyML_Printf("Tree lnL (line 3076 of utilities.c): %lf \n",tree->c_lnL);
			lk0 = Br_Len_Brent_Iter(10.*b_fcus->l[i],b_fcus->l[i],BL_MIN,
					tree->mod->s_opt->min_diff_lk_local,
					b_fcus,tree,
					tree->mod->s_opt->brent_it_max,
					tree->mod->s_opt->quickdirty,i);
		}
	}

	if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_local)
	{
		PhyML_Printf("\n\n%f %f %f %f\n",l_infa,l_max,l_infb,b_fcus->l);
		PhyML_Printf("%f -- %f \n",lk0_init,lk0);
		PhyML_Printf("\n. Err. in NNI (3)\n");

		//Warn_And_Exit("\n");
	}

	/***********/

	b_fcus->nni->lk0 = lk0;
	b_fcus->nni->lk1 = lk1;
	b_fcus->nni->lk2 = lk2;
	//JSJ: do over array...
	For(i,tree->n_l){
		l0[i]  = b_fcus->l[i];
		b_fcus->nni->l0[i]  = l0[i];
		b_fcus->nni->l1[i]  = l1[i];
		b_fcus->nni->l2[i]  = l2[i];
	}


	b_fcus->nni->score = lk0 - MAX(lk1,lk2);

	if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
			(b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
	{
		b_fcus->nni->score = .0;
		b_fcus->nni->lk1 = b_fcus->nni->lk0;
		b_fcus->nni->lk2 = b_fcus->nni->lk0;
	}

	if(lk0 > MAX(lk1,lk2))
	{
		b_fcus->nni->best_conf    = 0;
		For(i,tree->n_l) b_fcus->nni->best_l[i]       = l0[i]; //JSJ: fixed
		b_fcus->nni->swap_node_v1 = NULL;
		b_fcus->nni->swap_node_v2 = NULL;
		b_fcus->nni->swap_node_v3 = NULL;
		b_fcus->nni->swap_node_v4 = NULL;
	}
	else if(lk1 > MAX(lk0,lk2))
	{
		b_fcus->nni->best_conf    = 1;
		For(i,tree->n_l) b_fcus->nni->best_l[i]       = l1[i]; //JSJ: fixed
		b_fcus->nni->swap_node_v1 = v2;
		b_fcus->nni->swap_node_v2 = b_fcus->left;
		b_fcus->nni->swap_node_v3 = b_fcus->rght;
		b_fcus->nni->swap_node_v4 = v3;
	}
	else if(lk2 > MAX(lk0,lk1))
	{
		b_fcus->nni->best_conf    = 2;
		For(i,tree->n_l) b_fcus->nni->best_l[i]       = l2[i]; //JSJ: fixed
		b_fcus->nni->swap_node_v1 = v2;
		b_fcus->nni->swap_node_v2 = b_fcus->left;
		b_fcus->nni->swap_node_v3 = b_fcus->rght;
		b_fcus->nni->swap_node_v4 = v4;
	}
	else
	{
		b_fcus->nni->score        = +1.0;
		b_fcus->nni->best_conf    = 0;
		For(i,tree->n_l) b_fcus->nni->best_l[i]       = l0[i]; //JSJ: fixed
		b_fcus->nni->swap_node_v1 = NULL;
		b_fcus->nni->swap_node_v2 = NULL;
		b_fcus->nni->swap_node_v3 = NULL;
		b_fcus->nni->swap_node_v4 = NULL;
	}

	if((do_swap) && ((lk1 > lk0) || (lk2 > lk0)))
	{
		tree->n_swap++;
		PhyML_Printf("Swap edge %d -> %f\n",b_fcus->num,MAX(lk1,lk2));

		if(lk1 > lk2)
		{
			tree->best_lnL = lk1;
			Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
			For(i,tree->n_l) b_fcus->l[i] = l1[i];
			tree->both_sides = 1;
			Lk(tree);
		}
		else
		{
			tree->best_lnL = lk2;
			Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
			For(i,tree->n_l) b_fcus->l[i] = l2[i];
			tree->both_sides = 1;
			Lk(tree);
		}
	}
	else
	{
		For(i,tree->n_l) b_fcus->l[i] = bl_init[i];
		Update_PMat_At_Given_Edge(b_fcus,tree);
		tree->c_lnL = lk_init;
	}

}

/*********************************************************/

void NNI_Pars(arbre *tree, edge *b_fcus, int do_swap)
{
	int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
	node *v1,*v2,*v3,*v4;
	int pars0, pars1, pars2;
	int pars_init;

	pars_init              = tree->c_pars;
	b_fcus->nni->best_conf = 0;
	b_fcus->nni->score     = +1.0;

	pars0 = pars1 = pars2  = 0;
	v1 = v2 = v3 = v4      = NULL;

	l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

	l_r                    = b_fcus->l_r;
	r_l                    = b_fcus->r_l;

	v1                     = b_fcus->left->v[b_fcus->l_v1];
	v2                     = b_fcus->left->v[b_fcus->l_v2];
	v3                     = b_fcus->rght->v[b_fcus->r_v1];
	v4                     = b_fcus->rght->v[b_fcus->r_v2];

	if(v1->num < v2->num)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	if(v3->num < v4->num)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}


	/***********/
	Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
	tree->both_sides = 1;
	pars1 = Update_Pars_At_Given_Edge(b_fcus,tree);
	Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
	/***********/

	/***********/
	Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	tree->both_sides = 1;
	pars2 = Update_Pars_At_Given_Edge(b_fcus,tree);
	Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
	/***********/


	/***********/
	tree->both_sides = 1;
	pars0 = Update_Pars_At_Given_Edge(b_fcus,tree);

	if(pars0 != pars_init)
	{
		PhyML_Printf("\n. pars_init = %d; pars0 = %d\n",
				pars_init,
				pars0);
		Warn_And_Exit("\n. Err. in NNI (3)\n");
	}
	/***********/

	tree->c_pars = pars0;

	b_fcus->nni->score = MIN(pars1,pars2) - pars0;

	if(pars0 < MIN(pars1,pars2))
	{
		b_fcus->nni->best_conf    = 0;
		b_fcus->nni->swap_node_v1 = NULL;
		b_fcus->nni->swap_node_v2 = NULL;
		b_fcus->nni->swap_node_v3 = NULL;
		b_fcus->nni->swap_node_v4 = NULL;
	}
	else if(pars1 < MIN(pars0,pars2))
	{
		b_fcus->nni->best_conf    = 1;
		b_fcus->nni->swap_node_v1 = v2;
		b_fcus->nni->swap_node_v2 = b_fcus->left;
		b_fcus->nni->swap_node_v3 = b_fcus->rght;
		b_fcus->nni->swap_node_v4 = v3;
	}
	else if(pars2 > MIN(pars0,pars1))
	{
		b_fcus->nni->best_conf    = 2;
		b_fcus->nni->swap_node_v1 = v2;
		b_fcus->nni->swap_node_v2 = b_fcus->left;
		b_fcus->nni->swap_node_v3 = b_fcus->rght;
		b_fcus->nni->swap_node_v4 = v4;
	}
	else
	{
		b_fcus->nni->score        = +1.0;
		b_fcus->nni->swap_node_v1 = NULL;
		b_fcus->nni->swap_node_v2 = NULL;
		b_fcus->nni->swap_node_v3 = NULL;
		b_fcus->nni->swap_node_v4 = NULL;
	}
}

/*********************************************************/

void Swap(node *a, node *b, node *c, node *d, arbre *tree)
{
	int ab, ba, cd, dc, bc;
	int i;


	/* \             /d      \             /a
	 *  \           /         \           /
	 *   \b__...__c/    ->     \b__...__c/
	 *   /         \	       /		 \
	 *  /           \	      /		      \
	 * /a            \  	 /d            \
	 *
	 * nodes b and c are not necessarily on the same branch
	 */


#ifdef DEBUG
	if(!a || !b || !c || !d)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
#endif


	ab = ba = cd = dc = bc = -1;

	For(i,3) if(a->v[i] == b) { ab = i; break; }
	For(i,3) if(b->v[i] == a) { ba = i; break; }
	For(i,3) if(c->v[i] == d) { cd = i; break; }
	For(i,3) if(d->v[i] == c) { dc = i; break; }
	For(i,3) if(b->v[i] == c) { bc = i; break; }

#ifdef DEBUG
	if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
	{
		PhyML_Printf("\n. Nodes %d %d %d %d\n",a->num,b->num,c->num,d->num);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
#endif

	a->v[ab] = c;
	d->v[dc] = b;
	b->v[ba] = d;
	c->v[cd] = a;
	b->b[ba] = d->b[dc];
	c->b[cd] = a->b[ab];

	(a->b[ab]->left == b)?
			(a->b[ab]->left = c):
				(a->b[ab]->rght = c);

			(d->b[dc]->left == c)?
					(d->b[dc]->left = b):
						(d->b[dc]->rght = b);

					For(i,3)
					{
						if(a->b[ab]->left->v[i] == a->b[ab]->rght) a->b[ab]->l_r = i;
						if(a->b[ab]->rght->v[i] == a->b[ab]->left) a->b[ab]->r_l = i;
						if(d->b[dc]->left->v[i] == d->b[dc]->rght) d->b[dc]->l_r = i;
						if(d->b[dc]->rght->v[i] == d->b[dc]->left) d->b[dc]->r_l = i;
					}


					a->b[ab]->l_v1 = a->b[ab]->l_v2 =
							a->b[ab]->r_v1 = a->b[ab]->r_v2 =
									d->b[dc]->l_v1 = d->b[dc]->l_v2 =
											d->b[dc]->r_v1 = d->b[dc]->r_v2 = -1;


					For(i,3)
					{
						if(i != a->b[ab]->l_r)
						{
							if(a->b[ab]->l_v1 < 0) a->b[ab]->l_v1 = i;
							else a->b[ab]->l_v2 = i;
						}
						if(i != a->b[ab]->r_l)
						{
							if(a->b[ab]->r_v1 < 0) a->b[ab]->r_v1 = i;
							else a->b[ab]->r_v2 = i;
						}
						if(i != d->b[dc]->l_r)
						{
							if(d->b[dc]->l_v1 < 0) d->b[dc]->l_v1 = i;
							else d->b[dc]->l_v2 = i;
						}
						if(i != d->b[dc]->r_l)
						{
							if(d->b[dc]->r_v1 < 0) d->b[dc]->r_v1 = i;
							else d->b[dc]->r_v2 = i;
						}
					}
					Update_Dirs(tree);
}

/*********************************************************/

void Update_All_Partial_Lk(edge *b_fcus, arbre *tree)
{

	Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v1],
			b_fcus->left,
			b_fcus->left->v[b_fcus->l_v1],
			tree);

	Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v2],
			b_fcus->left,
			b_fcus->left->v[b_fcus->l_v2],
			tree);

	Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v1],
			b_fcus->rght,
			b_fcus->rght->v[b_fcus->r_v1],
			tree);

	Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v2],
			b_fcus->rght,
			b_fcus->rght->v[b_fcus->r_v2],
			tree);

	tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
}

/*********************************************************/

void Update_SubTree_Partial_Lk(edge *b_fcus, node *a, node *d, arbre *tree)
{
	int i;

	Update_P_Lk(tree,b_fcus,a);
	if(d->tax) return;
	else For(i,3) if(d->v[i] != a)
		Update_SubTree_Partial_Lk(d->b[i],d,d->v[i],tree);
}

/*********************************************************/

allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names)
{
	allseq *alldata;
	int j;

	alldata                        = (allseq *)mCalloc(1,sizeof(allseq));
	alldata->n_otu                 = n_otu;
	alldata->c_seq                 = (seq **)mCalloc(n_otu,sizeof(seq *));
	alldata->b_frq                 = (m3ldbl *)mCalloc(T_MAX_ALPHABET,sizeof(m3ldbl));
	alldata->wght                  = (int *)mCalloc(crunch_len,sizeof(int));
	alldata->ambigu                = (short int *)mCalloc(crunch_len,sizeof(short int));
	alldata->invar                 = (short int *)mCalloc(crunch_len,sizeof(short int));
	alldata->sitepatt              = (int *)mCalloc(  init_len,sizeof(int ));

	alldata->crunch_len = crunch_len;
	alldata->init_len   = init_len;
	alldata->obs_pinvar = .0;

	For(j,n_otu)
	{
		alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
		alldata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
		strcpy(alldata->c_seq[j]->name,sp_names[j]);
		alldata->c_seq[j]->state     = (char *)mCalloc(crunch_len,sizeof(char));
		alldata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
	}

	return alldata;
}

/*********************************************************/

arbrelist *Make_Treelist(int list_size)
{
	arbrelist *tlist;

	tlist = (arbrelist *)mCalloc(1,sizeof(arbrelist));
	tlist->list_size = list_size;
	tlist->tree = (arbre **)mCalloc(list_size,sizeof(arbre *));

	return tlist;
}


/*********************************************************/

void Copy_Seq_Names_To_Tip_Labels(arbre *tree, allseq *data)
{
	int i;
	For(i,tree->n_otu)
	{
		strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
	}
}

/*********************************************************/

allseq *Copy_Cseq(allseq *ori, int len, int ns)
{
	allseq *new;
	int i,j,n_otu;
	char **sp_names;

	n_otu = ori->n_otu;

	sp_names = (char **)mCalloc(n_otu,sizeof(char *));
	For(i,n_otu)
	{
		sp_names[i] = (char *)mCalloc(strlen(ori->c_seq[i]->name)+1,sizeof(char));
		strcpy(sp_names[i],ori->c_seq[i]->name);
	}

	new = Make_Cseq(n_otu, len, ori->init_len, sp_names);

	new->obs_pinvar = ori->obs_pinvar;

	For(i,ori->init_len) new->sitepatt[i] = ori->sitepatt[i];

	For(j,ori->crunch_len)
	{
		For(i,ori->n_otu)
		{
			new->c_seq[i]->state[j]     = ori->c_seq[i]->state[j];
			new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
		}

		new->wght[j]   = ori->wght[j];
		new->ambigu[j] = ori->ambigu[j];
		new->invar[j]  = ori->invar[j];
	}

	For(i,ori->n_otu)
	{
		new->c_seq[i]->len = ori->c_seq[i]->len;
		strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
	}

	new->init_len           = ori->init_len;
	new->clean_len          = ori->clean_len;
	new->crunch_len         = ori->crunch_len;
	For(i,ns) new->b_frq[i] = ori->b_frq[i];
	new->n_otu              = ori->n_otu;

	For(i,n_otu) Free(sp_names[i]);
	Free(sp_names);

	return new;
}

/*********************************************************/

optimiz *Alloc_Optimiz()
{
	optimiz *s_opt;
	s_opt = (optimiz *)mCalloc(1,sizeof(optimiz));
	return s_opt;
}

/*********************************************************/


int Filexists(char *filename)
{
	FILE *fp;
	fp =fopen(filename,"r");
	if (fp) {
		fclose(fp);
		return 1;
	} else
		return 0;
}

/*********************************************************/

FILE *Openfile(char *filename, int mode)
{
	/* mode = 0 -> read */
	/* mode = 1 -> write */
	/* mode = 2 -> append */

	FILE *fp;
	char *s;
	int open_test=0;

	/*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */

	/*   strcpy(s,filename); */

	s = filename;

	fp = NULL;

	switch(mode)
	{
	case 0 :
	{
		while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10)
		{
			PhyML_Printf("\n. Can't open file '%s', enter a new name : ",s);
			Getstring_Stdin(s);
		}
		break;
	}
	case 1 :
	{
		fp = (FILE *)fopen(s,"w");
		break;
	}
	case 2 :
	{
		fp = (FILE *)fopen(s,"a");
		break;
	}

	default : break;

	}

	/*   Free(s); */

	return fp;
}

/*********************************************************/

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *io, int n_data_set, int num_tree)
{
	char *s;
	div_t hour,min;


	/*   int i; */

	/*   For(i,2*tree->n_otu-3) fprintf(fp_out,"\n. * Edge %3d: %f",i,tree->t_edges[i]->l); */


	if((!n_data_set) || (!num_tree))
	{
		Print_Banner_Small(fp_out);
	}

	PhyML_Fprintf(fp_out,"\n\n. Sequence filename: \t\t\t%s", Basename(io->in_seq_file));
	PhyML_Fprintf(fp_out,"\n\n. Data set: \t\t\t\t#%d",n_data_set);

	if(io->mod->s_opt->random_input_tree)
		PhyML_Fprintf(fp_out,"\n\n. Random init tree: \t\t\t#%d",num_tree+1);
	else if(io->n_trees > 1)
		PhyML_Fprintf(fp_out,"\n\n. Starting tree number: \t\t\t#%d",num_tree+1);

	if(io->mod->s_opt->opt_topo)
	{
		if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Fprintf(fp_out,"\n\n. Tree topology search : \t\tNNIs");
		else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Fprintf(fp_out,"\n\n. Tree topology search : \t\tSPRs");
		else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Fprintf(fp_out,"\n\n. Tree topology search : \t\tBest of NNIs and SPRs");
	}
	else
	{
		PhyML_Fprintf(fp_out,"\n\n. Tree topology: \t\t\tfixed");
	}


	/* was after Sequence file ; moved here FLT */
	s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	if(io->in_tree == 2)
	{
		strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")");
	}
	else
	{
		if(!io->mod->s_opt->random_input_tree)
		{
			if(io->in_tree == 0)
				strcat(s,"BioNJ");
			if(io->in_tree == 1)
				strcat(s,"parsimony");
		}
		else
		{
			strcat(s,"random tree");
		}
	}

	PhyML_Fprintf(fp_out,"\n\n. Initial tree: \t\t\t%s",s);
	Free(s);

	(tree->mod->datatype == NT)?
			(fprintf(fp_out,"\n\n. Model of nucleotides substitution: \t%s",io->mod->modelname)):
				(fprintf(fp_out,"\n\n. Model of amino acids substitution: \t%s",io->mod->modelname));


			PhyML_Fprintf(fp_out,"\n\n. Number of taxa: \t\t\t%d",tree->n_otu);/*added FLT*/

			PhyML_Fprintf(fp_out,"\n\n. Log-likelihood: \t\t\t%.5f",tree->c_lnL);/*was last ; moved here FLT*/

			Unconstraint_Lk(tree);
			PhyML_Fprintf(fp_out,"\n\n. Unconstrained likelihood: \t\t%.5f",tree->unconstraint_lk);

			PhyML_Fprintf(fp_out,"\n\n. Parsimony: \t\t\t\t%d",tree->c_pars);

			PhyML_Fprintf(fp_out,"\n\n. Tree size: \t\t\t\t%.5f",tree->size);

			PhyML_Fprintf(fp_out,"\n\n. Discrete gamma model: \t\t%s",
					(tree->mod->n_catg>1)?("Yes"):("No"));
			if(tree->mod->n_catg > 1)
			{
				PhyML_Fprintf(fp_out,"\n  - Number of categories: \t\t%d",tree->mod->n_catg);
				PhyML_Fprintf(fp_out,"\n  - Gamma shape parameter: \t\t%.3f",tree->mod->alpha);
			}

			if(tree->mod->invar) PhyML_Fprintf(fp_out,"\n\n. Proportion of invariant: \t\t%.3f",tree->mod->pinvar);

			/*was before Discrete gamma model ; moved here FLT*/
			if((tree->mod->whichmodel == K80)   ||
					(tree->mod->whichmodel == HKY85) ||
					(tree->mod->whichmodel == F84))
				PhyML_Fprintf(fp_out,"\n\n. Transition/transversion ratio: \t%.3f",tree->mod->kappa);
			else if(tree->mod->whichmodel == TN93)
			{
				PhyML_Fprintf(fp_out,"\n\n. Transition/transversion ratio for purines: \t\t\t%.3f",
						tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
				PhyML_Fprintf(fp_out,"\n\n. Transition/transversion ratio for pyrimidines: \t\t\t%.3f",
						tree->mod->kappa*2./(1.+tree->mod->lambda));
			}

			if(tree->mod->datatype == NT)
			{
				PhyML_Fprintf(fp_out,"\n\n. Nucleotides frequencies:\n");
				PhyML_Fprintf(fp_out,"\n  - f(A)=%8.5f",tree->mod->pi[0]);
				PhyML_Fprintf(fp_out,"\n  - f(C)=%8.5f",tree->mod->pi[1]);
				PhyML_Fprintf(fp_out,"\n  - f(G)=%8.5f",tree->mod->pi[2]);
				PhyML_Fprintf(fp_out,"\n  - f(T)=%8.5f",tree->mod->pi[3]);
			}

			/*****************************************/
			if((tree->mod->whichmodel == GTR) ||
					(tree->mod->whichmodel == CUSTOM))
			{
				int i,j;

				Update_Qmat_GTR(tree->mod->rr,
						tree->mod->rr_val,
						tree->mod->rr_num,
						tree->mod->pi,
						tree->mod->qmat);

				PhyML_Fprintf(fp_out,"\n\n");
				PhyML_Fprintf(fp_out,". GTR relative rate parameters : \n\n");
				PhyML_Fprintf(fp_out,"  A <-> C   %8.5f\n",  tree->mod->rr[0]);
				PhyML_Fprintf(fp_out,"  A <-> G   %8.5f\n",  tree->mod->rr[1]);
				PhyML_Fprintf(fp_out,"  A <-> T   %8.5f\n",  tree->mod->rr[2]);
				PhyML_Fprintf(fp_out,"  C <-> G   %8.5f\n",  tree->mod->rr[3]);
				PhyML_Fprintf(fp_out,"  C <-> T   %8.5f\n",  tree->mod->rr[4]);
				PhyML_Fprintf(fp_out,"  G <-> T   %8.5f\n",tree->mod->rr[5]);


				PhyML_Fprintf(fp_out,"\n\n. Instantaneous rate matrix : \n");
				PhyML_Fprintf(fp_out,"\n  [A---------C---------G---------T------]\n");
				For(i,4)
				{
					PhyML_Fprintf(fp_out,"  ");
					For(j,4)
					PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
					PhyML_Fprintf(fp_out,"\n");
				}
				PhyML_Fprintf(fp_out,"\n");
			}
			/*****************************************/


			if(io->ratio_test == 1)
			{
				PhyML_Fprintf(fp_out,". aLRT statistics to test branches");
			}
			else if(io->ratio_test == 2)
			{
				PhyML_Fprintf(fp_out,". aLRT branch supports (cubic approximation, mixture of Chi2s distribution)");
			}


			hour = div(t_end-t_beg,3600);
			min  = div(t_end-t_beg,60  );

			min.quot -= hour.quot*60;

			PhyML_Fprintf(fp_out,"\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
			PhyML_Fprintf(fp_out,"\n\n. %d seconds\n",(int)(t_end-t_beg));

			PhyML_Fprintf(fp_out,"\n");
			PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
			PhyML_Fprintf(fp_out," Suggested citation:\n");
			PhyML_Fprintf(fp_out," S. Guindon & O. Gascuel\n");
			PhyML_Fprintf(fp_out," \"A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood\"\n");
			PhyML_Fprintf(fp_out," Systematic Biology. 2003. 52(5):696-704.\n");
			PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");



}

/*********************************************************/
/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *io, int n_data_set)
{
	char *s;
	/*div_t hour,min;*/

	if (n_data_set==1)
	{

		fprintf(fp_out,". Sequence file : [%s]\n\n", io->in_seq_file);

		(tree->mod->datatype == NT)?
				(fprintf(fp_out,". Model of nucleotides substitution : %s\n\n",io->mod->modelname)):
					(fprintf(fp_out,". Model of amino acids substitution : %s\n\n",io->mod->modelname));

				s = (char *)mCalloc(100,sizeof(char));

				switch(io->in_tree)
				{
				case 0: { strcpy(s,"BioNJ");     break; }
				case 1: { strcpy(s,"parsimony"); break; }
				case 2: { strcpy(s,"user tree (");
				strcat(s,io->in_tree_file);
				strcat(s,")");         break; }
				}

				fprintf(fp_out,". Initial tree : [%s]\n\n",s);

				Free(s);

				fprintf(fp_out,"\n");

				/*headline 1*/
				fprintf(fp_out, ". Data\t");

				fprintf(fp_out,"Nb of \t");

				fprintf(fp_out,"Likelihood\t");

				fprintf(fp_out, "Discrete   \t");

				if(tree->mod->n_catg > 1)
					PhyML_Fprintf(fp_out, "Number of \tGamma shape\t");

				fprintf(fp_out,"Proportion of\t");

				if(tree->mod->whichmodel <= 6)
					PhyML_Fprintf(fp_out,"Transition/ \t");

				fprintf(fp_out,"Nucleotides frequencies               \t");

				if((tree->mod->whichmodel == GTR) ||
						(tree->mod->whichmodel == CUSTOM))
					PhyML_Fprintf(fp_out,"Instantaneous rate matrix              \t");

				/*    PhyML_Fprintf(fp_out,"Time\t");*/

				fprintf(fp_out, "\n");


				/*headline 2*/
				fprintf(fp_out, "  set\t");

				fprintf(fp_out,"taxa\t");

				fprintf(fp_out,"loglk     \t");

				fprintf(fp_out, "gamma model\t");

				if(tree->mod->n_catg > 1)
					PhyML_Fprintf(fp_out, "categories\tparameter  \t");

				fprintf(fp_out,"invariant    \t");

				if(tree->mod->whichmodel <= 6)
					PhyML_Fprintf(fp_out,"transversion\t");

				fprintf(fp_out,"f(A)      f(C)      f(G)      f(T)    \t");

				if((tree->mod->whichmodel == GTR) ||
						(tree->mod->whichmodel == CUSTOM))
					PhyML_Fprintf(fp_out,"[A---------C---------G---------T------]\t");

				/*    PhyML_Fprintf(fp_out,"used\t");*/

				fprintf(fp_out, "\n");


				/*headline 3*/
				if(tree->mod->whichmodel == TN93)
				{
					PhyML_Fprintf(fp_out,"    \t      \t          \t           \t");
					if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"         \t         \t");
					PhyML_Fprintf(fp_out,"             \t");
					PhyML_Fprintf(fp_out,"purines pyrimid.\t");
					PhyML_Fprintf(fp_out, "\n");
				}

				PhyML_Fprintf(fp_out, "\n");
	}


	/*line items*/

	PhyML_Fprintf(fp_out,"  #%d\t",n_data_set);

	PhyML_Fprintf(fp_out,"%d   \t",tree->n_otu);

	PhyML_Fprintf(fp_out,"%.5f\t",tree->c_lnL);

	PhyML_Fprintf(fp_out,"%s        \t",
			(tree->mod->n_catg>1)?("Yes"):("No "));
	if(tree->mod->n_catg > 1)
	{
		PhyML_Fprintf(fp_out,"%d        \t",tree->mod->n_catg);
		PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->alpha);
	}

	/*if(tree->mod->invar)*/
	PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->pinvar);

	if(tree->mod->whichmodel <= 5)
	{
		PhyML_Fprintf(fp_out,"%.3f     \t",tree->mod->kappa);
	}
	else if(tree->mod->whichmodel == TN93)
	{
		PhyML_Fprintf(fp_out,"%.3f   ",
				tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
		PhyML_Fprintf(fp_out,"%.3f\t",
				tree->mod->kappa*2./(1.+tree->mod->lambda));
	}


	if(tree->mod->datatype == NT)
	{
		PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[0]);
		PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[1]);
		PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[2]);
		PhyML_Fprintf(fp_out,"%8.5f\t",tree->mod->pi[3]);
	}
	/*
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    PhyML_Fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
	 */

	/*****************************************/
	if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
	{
		int i,j;

		For(i,4)
		{
			if (i!=0) {
				/*format*/
				PhyML_Fprintf(fp_out,"      \t     \t          \t           \t");
				if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"          \t           \t");
				PhyML_Fprintf(fp_out,"             \t                                      \t");
			}
			For(j,4)
			PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
			if (i<3) PhyML_Fprintf(fp_out,"\n");
		}
	}
	/*****************************************/

	PhyML_Fprintf(fp_out, "\n\n");
}

/*********************************************************/

matrix *K80_dist(allseq *data, m3ldbl g_shape)
{
	int i,j,k;
	int diff;
	m3ldbl unc_len;
	matrix *mat;
	m3ldbl **len;

	len = (m3ldbl **)mCalloc(data->n_otu,sizeof(m3ldbl *));
	For(i,data->n_otu)
	len[i] = (m3ldbl *)mCalloc(data->n_otu,sizeof(m3ldbl));

	unc_len = .0;

	mat = Make_Mat(data->n_otu);
	Init_Mat(mat,data);

	For(i,data->c_seq[0]->len)
	{
		For(j,data->n_otu-1)
		{
			for(k=j+1;k<data->n_otu;k++)
			{
				if(((data->c_seq[j]->state[i] == 'A' || data->c_seq[j]->state[i] == 'G') &&
						(data->c_seq[k]->state[i] == 'C' || data->c_seq[k]->state[i] == 'T'))||
						((data->c_seq[j]->state[i] == 'C' || data->c_seq[j]->state[i] == 'T') &&
								(data->c_seq[k]->state[i] == 'A' || data->c_seq[k]->state[i] == 'G')))
				{
					diff++;
					mat->Q[j][k]+=data->wght[i];
					len[j][k]+=data->wght[i];
					len[k][j]=len[j][k];
				}

				else
					if(((data->c_seq[j]->state[i] == 'A' && data->c_seq[k]->state[i] == 'G') ||
							(data->c_seq[j]->state[i] == 'G' && data->c_seq[k]->state[i] == 'A'))||
							((data->c_seq[j]->state[i] == 'C' && data->c_seq[k]->state[i] == 'T') ||
									(data->c_seq[j]->state[i] == 'T' && data->c_seq[k]->state[i] == 'C')))
					{
						diff++;
						mat->P[j][k]+=data->wght[i];
						len[j][k]+=data->wght[i];
						len[k][j]=len[j][k];
					}
					else
						if((data->c_seq[j]->state[i] == 'A' ||
								data->c_seq[j]->state[i] == 'C' ||
								data->c_seq[j]->state[i] == 'G' ||
								data->c_seq[j]->state[i] == 'T')&&
								(data->c_seq[k]->state[i] == 'A' ||
										data->c_seq[k]->state[i] == 'C' ||
										data->c_seq[k]->state[i] == 'G' ||
										data->c_seq[k]->state[i] == 'T'))
						{
							len[j][k]+=data->wght[i];
							len[k][j]=len[j][k];
						}
			}
		}
	}


	For(i,data->n_otu-1)
	for(j=i+1;j<data->n_otu;j++)
	{
		if(len[i][j])
		{
			mat->P[i][j] /= len[i][j];
			mat->Q[i][j] /= len[i][j];
		}
		else
		{
			mat->P[i][j] = .5;
			mat->Q[i][j] = .5;
		}

		mat->P[j][i] = mat->P[i][j];
		mat->Q[j][i] = mat->Q[i][j];


		if((1-2*mat->P[i][j]-mat->Q[i][j] <= .0) || (1-2*mat->Q[i][j] <= .0))
		{
			mat->dist[i][j] = -1.;
			mat->dist[j][i] = -1.;
			continue;
		}

		mat->dist[i][j] = (g_shape/2)*
				(pow(1-2*mat->P[i][j]-mat->Q[i][j],-1./g_shape) +
						0.5*pow(1-2*mat->Q[i][j],-1./g_shape) - 1.5);


		if(mat->dist[i][j] > DIST_MAX)
		{
			mat->dist[i][j] = DIST_MAX;
		}
		mat->dist[j][i] = mat->dist[i][j];
	}

	For(i,data->n_otu) free(len[i]);
	free(len);
	return mat;
}

/*********************************************************/

matrix *JC69_Dist(allseq *data, model *mod)
{
	int site,i,j,k;
	m3ldbl unc_len;
	matrix *mat;
	m3ldbl **len;


	len = (m3ldbl **)mCalloc(data->n_otu,sizeof(m3ldbl *));
	For(i,data->n_otu)
	len[i] = (m3ldbl *)mCalloc(data->n_otu,sizeof(m3ldbl));

	unc_len = .0;

	mat = Make_Mat(data->n_otu);
	Init_Mat(mat,data);

	Fors(site,data->c_seq[0]->len,mod->stepsize)
	{
		For(j,data->n_otu-1)
		{
			for(k=j+1;k<data->n_otu;k++)
			{
				if((!Is_Ambigu(data->c_seq[j]->state+site,mod->datatype,mod->stepsize)) &&
						(!Is_Ambigu(data->c_seq[k]->state+site,mod->datatype,mod->stepsize)))
				{
					len[j][k]+=data->wght[site];
					len[k][j]=len[j][k];
					if(strncmp(data->c_seq[j]->state+site,data->c_seq[k]->state+site,mod->stepsize))
						mat->P[j][k]+=data->wght[site];
				}
			}
		}
	}


	For(i,data->n_otu-1)
	for(j=i+1;j<data->n_otu;j++)
	{
		if(len[i][j])
		{
			mat->P[i][j] /= len[i][j];
		}
		else
		{
			mat->P[i][j] = 1.;
		}

		mat->P[j][i] = mat->P[i][j];

		if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0)
		{
			mat->dist[i][j] = DIST_MAX;
		}
		else
			mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(m3ldbl)log(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);


		/* 	PhyML_Printf("\n. Incorrect JC distances"); */
		/* 	mat->dist[i][j] = len[i][j]; */


		if(mat->dist[i][j] > DIST_MAX)
		{
			mat->dist[i][j] = DIST_MAX;
		}

		mat->dist[j][i] = mat->dist[i][j];
	}

	For(i,data->n_otu) free(len[i]);
	free(len);

	return mat;
}

/*********************************************************/

matrix *Hamming_Dist(allseq *data, model *mod)
{
	int i,j,k;
	m3ldbl unc_len;
	matrix *mat;
	m3ldbl **len;


	len = (m3ldbl **)mCalloc(data->n_otu,sizeof(m3ldbl *));
	For(i,data->n_otu)
	len[i] = (m3ldbl *)mCalloc(data->n_otu,sizeof(m3ldbl));

	unc_len = .0;

	mat = Make_Mat(data->n_otu);
	Init_Mat(mat,data);

	For(i,data->c_seq[0]->len)
	{
		For(j,data->n_otu-1)
		{
			for(k=j+1;k<data->n_otu;k++)
			{
				if((!Is_Ambigu(data->c_seq[j]->state+i,mod->datatype,mod->stepsize)) &&
						(!Is_Ambigu(data->c_seq[k]->state+i,mod->datatype,mod->stepsize)))
				{
					len[j][k]+=data->wght[i];
					len[k][j]=len[j][k];
					if(data->c_seq[j]->state[i] != data->c_seq[k]->state[i])
						mat->P[j][k]+=data->wght[i];
				}
			}
		}
	}


	For(i,data->n_otu-1)
	for(j=i+1;j<data->n_otu;j++)
	{
		if(len[i][j])
		{
			mat->P[i][j] /= len[i][j];
		}
		else
		{
			mat->P[i][j] = 1.;
		}

		mat->P[j][i] = mat->P[i][j];

		mat->dist[i][j] = mat->P[i][j];


		if(mat->dist[i][j] > DIST_MAX)
		{
			mat->dist[i][j] = DIST_MAX;
		}
		mat->dist[j][i] = mat->dist[i][j];
	}

	For(i,data->n_otu) free(len[i]);
	free(len);

	return mat;
}

/*********************************************************/
/* Test if the given site pattern is invariant. Does not handle ambiguities */

int Is_Invar(int patt_num, int stepsize, int datatype, allseq *data)
{
	int i, j;

	For(i,data->n_otu)
	{
		For(j,data->n_otu)
		{
			if(!(Are_Compatible(data->c_seq[i]->state+patt_num,
					data->c_seq[j]->state+patt_num,
					stepsize,
					datatype)))
			{
				break;
			}
		}
		if(j != data->n_otu) break;
	}

	if(i == data->n_otu) return 1;
	else                 return 0;
}


/*********************************************************/

int Is_Ambigu(char *state, int datatype, int stepsize)
{
	int val,i;

	val = -1;
	if(datatype == NT)
	{
		For(i,stepsize)
		{
			switch(state[i])
			{
			case 'A' : case 'C' : case 'G' : case 'T' : case 'U' : { val=0; break; }
			default : { val=1; break; }
			}
			if(val == 1) break;
			/*       if(strchr("MRWSYKBDHVNXO?-.",(int)state)) return 1; */
		}
	}
	else
	{
		switch(state[0])
		{
		case 'X' : case '?' : case '-' : case '.' : {val=1; break; }
		default : { val=0; break; }
		}
		/*       if(strchr("X?-.",(int)state)) return 1; */
	}

	return val;
}

/*********************************************************/

void Check_Ambiguities(allseq *data, int datatype, int stepsize)
{
	int i,j;

	Fors(j,data->crunch_len,stepsize)
	{
		For(i,data->n_otu)
		{
			data->ambigu[j]              = 0;
			data->c_seq[i]->is_ambigu[j] = 0;
		}

		For(i,data->n_otu)
		{
			if(Is_Ambigu(data->c_seq[i]->state+j,
					datatype,
					stepsize))
			{
				data->ambigu[j]              = 1;
				data->c_seq[i]->is_ambigu[j] = 1;
			}
		}
	}
}

/*********************************************************/

int Get_State_From_Ui(int ui, int datatype)
{
	if(datatype == NT)
	{
		switch(ui)
		{
		case 1 : {return 0; break;}
		case 2 : {return 1; break;}
		case 4 : {return 2; break;}
		case 8 : {return 3; break;}
		default :
		{
			PhyML_Printf("\n. ui=%d",ui);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
			break;
		}
		}
	}
	else if(datatype == AA)
	{
		switch(ui)
		{
		case 1 :      {return 0;  break;}
		case 2 :      {return 1;  break;}
		case 4 :      {return 2;  break;}
		case 8 :      {return 3;  break;}
		case 16 :     {return 4;  break;}
		case 32 :     {return 5;  break;}
		case 64 :     {return 6;  break;}
		case 128 :    {return 7;  break;}
		case 256 :    {return 8;  break;}
		case 512 :    {return 9;  break;}
		case 1024 :   {return 10; break;}
		case 2048 :   {return 11; break;}
		case 4096 :   {return 12; break;}
		case 8192 :   {return 13; break;}
		case 16384 :  {return 14; break;}
		case 32768 :  {return 15; break;}
		case 65536 :  {return 16; break;}
		case 131072 : {return 17; break;}
		case 262144 : {return 18; break;}
		case 524288 : {return 19; break;}
		default :
		{
			PhyML_Printf("\n. ui=%d",ui);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}
		}
	}
	else
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
	return -1;
}

/*********************************************************/

int Assign_State(char *c, int datatype, int stepsize)
{
	int state[3];
	int i;

	state[0] = state[1] = state[2] = -1;
	if(datatype == NT)
	{
		For(i,stepsize)
		{
			switch(c[i])
			{
			case 'A' : {state[i]=0;  break;}
			case 'C' : {state[i]=1;  break;}
			case 'G' : {state[i]=2;  break;}
			case 'T' : {state[i]=3;  break;}
			case 'U' : {state[i]=3;  break;}
			default  : {state[i]=-1; break;}
			}
		}
		return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
	}
	else
	{
		switch(c[0])
		{
		case 'A' : {state[0]=0 ; break;}
		case 'R' : {state[0]=1 ; break;}
		case 'N' : {state[0]=2 ; break;}
		case 'D' : {state[0]=3 ; break;}
		case 'C' : {state[0]=4 ; break;}
		case 'Q' : {state[0]=5 ; break;}
		case 'E' : {state[0]=6 ; break;}
		case 'G' : {state[0]=7 ; break;}
		case 'H' : {state[0]=8 ; break;}
		case 'I' : {state[0]=9 ; break;}
		case 'L' : {state[0]=10; break;}
		case 'K' : {state[0]=11; break;}
		case 'M' : {state[0]=12; break;}
		case 'F' : {state[0]=13; break;}
		case 'P' : {state[0]=14; break;}
		case 'S' : {state[0]=15; break;}
		case 'T' : {state[0]=16; break;}
		case 'W' : {state[0]=17; break;}
		case 'Y' : {state[0]=18; break;}
		case 'V' : {state[0]=19; break;}

		case 'B' : {state[0] = 2; break;}
		case 'Z' : {state[0] = 5; break;}
		default  : {state[0]=-1;  break;}
		}
		return state[0];
	}
	return -1;
}

/*********************************************************/

char Reciproc_Assign_State(int i_state, int datatype)
{

	if(datatype == NT)
	{
		i_state = i_state%4;
		switch(i_state)
		{
		case 0 :   {return 'A';  break;}
		case 1 :   {return 'C';  break;}
		case 2 :   {return 'G';  break;}
		case 3 :   {return 'T';  break;}
		default  :
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
			break;
		}
		}
	}
	else
	{
		i_state = i_state%20;
		switch(i_state)
		{
		case 0  : {return 'A' ; break;}
		case 1  : {return 'R' ; break;}
		case 2  : {return 'N' ; break;}
		case 3  : {return 'D' ; break;}
		case 4  : {return 'C' ; break;}
		case 5  : {return 'Q' ; break;}
		case 6  : {return 'E' ; break;}
		case 7  : {return 'G' ; break;}
		case 8  : {return 'H' ; break;}
		case 9  : {return 'I' ; break;}
		case 10 : {return 'L';  break;}
		case 11 : {return 'K';  break;}
		case 12 : {return 'M';  break;}
		case 13 : {return 'F';  break;}
		case 14 : {return 'P';  break;}
		case 15 : {return 'S';  break;}
		case 16 : {return 'T';  break;}
		case 17 : {return 'W';  break;}
		case 18 : {return 'Y';  break;}
		case 19 : {return 'V';  break;}
		default  :
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
			break;
		}
		}
	}
	return -1;
}

/*********************************************************/

int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize)
{
	int state[3];
	int i;

	state[0] = state[1] = state[2] = -1;
	if(datatype == NT)
	{
		For(i,stepsize)
		{
			switch(c[i])
			{
			case 'A' : {state[i]= 0;  break;}
			case 'C' : {state[i]= 1;  break;}
			case 'G' : {state[i]= 2;  break;}
			case 'T' : {state[i]= 3;  break;}
			case 'U' : {state[i]= 3;  break;}
			case 'M' : {state[i]= 4;  break;}
			case 'R' : {state[i]= 5;  break;}
			case 'W' : {state[i]= 6;  break;}
			case 'S' : {state[i]= 7;  break;}
			case 'Y' : {state[i]= 8;  break;}
			case 'K' : {state[i]= 9;  break;}
			case 'B' : {state[i]=10;  break;}
			case 'D' : {state[i]=11;  break;}
			case 'H' : {state[i]=12;  break;}
			case 'V' : {state[i]=13;  break;}
			case 'N' : case 'X' : case '?' : case 'O' : case '-' : {state[i]=14;  break;}
			default :
			{
				PhyML_Printf("\n. Unknown character state : %c\n",state[i]);
				Warn_And_Exit("\n. Init failed (check the data type)\n");
				break;
			}
			}
			return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
		}
	}
	else
	{
		switch(c[0])
		{
		case 'A' : {state[0]= 0; break;}
		case 'R' : {state[0]= 1; break;}
		case 'N' : {state[0]= 2; break;}
		case 'D' : {state[0]= 3; break;}
		case 'C' : {state[0]= 4; break;}
		case 'Q' : {state[0]= 5; break;}
		case 'E' : {state[0]= 6; break;}
		case 'G' : {state[0]= 7; break;}
		case 'H' : {state[0]= 8; break;}
		case 'I' : {state[0]= 9; break;}
		case 'L' : {state[0]=10; break;}
		case 'K' : {state[0]=11; break;}
		case 'M' : {state[0]=12; break;}
		case 'F' : {state[0]=13; break;}
		case 'P' : {state[0]=14; break;}
		case 'S' : {state[0]=15; break;}
		case 'T' : {state[0]=16; break;}
		case 'W' : {state[0]=17; break;}
		case 'Y' : {state[0]=18; break;}
		case 'V' : {state[0]=19; break;}
		case 'B' : {state[0]= 2; break;}
		case 'Z' : {state[0]= 5; break;}
		case 'X' : case '?' : case '-' : {state[0]=20; break;}
		default  :
		{
			PhyML_Printf("\n. Unknown character state : %c\n",state[0]);
			Warn_And_Exit("\n. Init failed (check the data type)\n");
			break;
		}
		}
		return state[0];
	}
	return -1;
}

/*********************************************************/

void Clean_Tree_Connections(arbre *tree)
{

	int i;
	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i]->v[0] = NULL;
		tree->noeud[i]->v[1] = NULL;
		tree->noeud[i]->v[2] = NULL;
		tree->noeud[i]->b[0] = NULL;
		tree->noeud[i]->b[1] = NULL;
		tree->noeud[i]->b[2] = NULL;
	}
}

/*********************************************************/

void Bootstrap(arbre *tree)
{
	int *site_num, n_site;
	int replicate,j,k;
	int position,init_len;
	allseq *boot_data;
	arbre *boot_tree;
	model *boot_mod;
	matrix *boot_mat;
	char *s;
	/*   m3ldbl rf; */

	tree->print_boot_val = 1;
	tree->print_alrt_val = 0;
	boot_tree            = NULL;

	site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));

	Alloc_Bip(tree);
	Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

	n_site = 0;
	For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
	{
		site_num[n_site] = j;
		n_site++;
	}

	boot_data = Copy_Cseq(tree->data, tree->data->crunch_len, tree->mod->ns);

	PhyML_Printf("\n\n. Non parametric bootstrap analysis \n\n");
	PhyML_Printf("  [");


	For(replicate,tree->mod->bootstrap)
	{
		For(j,boot_data->crunch_len) boot_data->wght[j] = 0;

		init_len = 0;
		For(j,boot_data->init_len)
		{
			position = Rand_Int(0,(int)(tree->data->init_len-1.0));
			boot_data->wght[site_num[position]] += 1;
			init_len++;
		}

		if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");

		init_len = 0;
		For(j,boot_data->crunch_len) init_len += boot_data->wght[j];

		if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");

		(tree->mod->datatype == NT)?
				(Get_Base_Freqs(boot_data)):
					(Get_AA_Freqs(boot_data));

				if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

				boot_mod = Copy_Model(tree->mod);
				Init_Model(boot_data,boot_mod);

				if(tree->io->in_tree == 2)
				{
					rewind(tree->io->fp_in_tree);
					boot_tree = Read_Tree_File(tree->io->fp_in_tree);
				}
				else
				{
					boot_mat = ML_Dist(boot_data,boot_mod);
					boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data, tree->n_l);
					Fill_Missing_Dist(boot_mat);
					Bionj(boot_mat);
					boot_tree = boot_mat->tree;
					boot_tree->mat = boot_mat;
				}

				boot_tree->mod                = boot_mod;
				boot_tree->io                 = tree->io;
				boot_tree->data               = boot_data;
				boot_tree->both_sides         = 1;
				boot_tree->mod->s_opt->print  = 0;
				boot_tree->n_pattern          = boot_tree->data->crunch_len/
						boot_tree->mod->stepsize;
				boot_tree->io->print_site_lnl = 0;
				boot_tree->io->print_trace    = 0;

				if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
				Order_Tree_CSeq(boot_tree,boot_data);
				Share_Lk_Struct(tree,boot_tree);
				Share_Spr_Struct(tree,boot_tree);
				Share_Pars_Struct(tree,boot_tree);
				Fill_Dir_Table(boot_tree);
				Update_Dirs(boot_tree);

				if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
				else                         Init_P_Lk_Tips_Int(boot_tree);
				Init_Ui_Tips(boot_tree);
				Init_P_Pars_Tips(boot_tree);
				Br_Len_Not_Involving_Invar(boot_tree);


				if(boot_tree->mod->s_opt->opt_topo)
				{
					if(boot_tree->mod->s_opt->topo_search == NNI_MOVE)
					{
						Simu_Loop(boot_tree);
					}
					else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
							(boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
					{
						Speed_Spr_Loop(boot_tree);
					}
				}
				else
				{
					if(boot_tree->mod->s_opt->opt_num_param || boot_tree->mod->s_opt->opt_bl)
						Round_Optimize(boot_tree,boot_tree->data,ROUND_MAX);
					else
						Lk(boot_tree);
				}


				Alloc_Bip(boot_tree);

				Get_Bip(boot_tree->noeud[0],
						boot_tree->noeud[0]->v[0],
						boot_tree);

				Compare_Bip(tree,boot_tree);

				Br_Len_Involving_Invar(boot_tree);

				if(tree->io->print_boot_trees)
				{
					s = Write_Tree(boot_tree);
					PhyML_Fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
					Free(s);
					Print_Fp_Out_Lines(tree->io->fp_out_boot_stats,0,0,boot_tree,tree->io,replicate+1);
				}


				/*       rf = .0; */
				/*       For(j,2*tree->n_otu-3)  */
				/* 	rf += tree->t_edges[j]->bip_score; */


				PhyML_Printf(".");
#ifndef QUIET
				fflush(stdout);
#endif
				if(!((replicate+1)%20))
				{
					PhyML_Printf("] %4d/%4d\n  ",replicate+1,tree->mod->bootstrap);
					if(replicate != tree->mod->bootstrap-1) PhyML_Printf("[");
				}

				if(boot_tree->mat) Free_Mat(boot_tree->mat);
				Free_Tree(boot_tree);
				Free_Model(boot_mod);
	}

	if(((replicate)%20)) PhyML_Printf("] %4d/%4d\n ",replicate,tree->mod->bootstrap);

	tree->lock_topo = 1; /* Topology should not be modified afterwards */

	if(tree->io->print_boot_trees)
	{
		fclose(tree->io->fp_out_boot_tree);
		fclose(tree->io->fp_out_boot_stats);
	}

	Free_Cseq(boot_data);
	Free(site_num);
}

/*********************************************************/
//JSJ: temp update to the length
void Br_Len_Involving_Invar(arbre *tree)
{
	int i,k;
	For(k,tree->n_l){
		For(i,2*tree->n_otu-3) tree->t_edges[i]->l[k] *= (1.0-tree->mod->pinvar);
	}
}

/*********************************************************/
//JSJ: Another temp update
void Br_Len_Not_Involving_Invar(arbre *tree)
{
	int i,k;
	For(k,tree->n_l){
		For(i,2*tree->n_otu-3) tree->t_edges[i]->l[k] /= (1.0-tree->mod->pinvar);
	}
}

/*********************************************************/

void Getstring_Stdin(char *file_name)
{
	if(!fgets(file_name,T_MAX_LINE,stdin)) Exit("");
	if (strchr(file_name, '\n') != NULL)
		*strchr(file_name, '\n') = '\0';
}

/*********************************************************/

void Print_Freq(arbre *tree)
{

	switch(tree->mod->datatype)
	{
	case NT:
	{
		PhyML_Printf("A : %f\n",tree->mod->pi[0]);
		PhyML_Printf("C : %f\n",tree->mod->pi[1]);
		PhyML_Printf("G : %f\n",tree->mod->pi[2]);
		PhyML_Printf("T : %f\n",tree->mod->pi[3]);

		PhyML_Printf("U : %f\n",tree->mod->pi[4]);
		PhyML_Printf("M : %f\n",tree->mod->pi[5]);
		PhyML_Printf("R : %f\n",tree->mod->pi[6]);
		PhyML_Printf("W : %f\n",tree->mod->pi[7]);
		PhyML_Printf("S : %f\n",tree->mod->pi[8]);
		PhyML_Printf("Y : %f\n",tree->mod->pi[9]);
		PhyML_Printf("K : %f\n",tree->mod->pi[10]);
		PhyML_Printf("B : %f\n",tree->mod->pi[11]);
		PhyML_Printf("D : %f\n",tree->mod->pi[12]);
		PhyML_Printf("H : %f\n",tree->mod->pi[13]);
		PhyML_Printf("V : %f\n",tree->mod->pi[14]);
		PhyML_Printf("N : %f\n",tree->mod->pi[15]);
		break;
	}
	case AA:
	{
		PhyML_Printf("A : %f\n",tree->mod->pi[0]);
		PhyML_Printf("R : %f\n",tree->mod->pi[1]);
		PhyML_Printf("N : %f\n",tree->mod->pi[2]);
		PhyML_Printf("D : %f\n",tree->mod->pi[3]);
		PhyML_Printf("C : %f\n",tree->mod->pi[4]);
		PhyML_Printf("Q : %f\n",tree->mod->pi[5]);
		PhyML_Printf("E : %f\n",tree->mod->pi[6]);
		PhyML_Printf("G : %f\n",tree->mod->pi[7]);
		PhyML_Printf("H : %f\n",tree->mod->pi[8]);
		PhyML_Printf("I : %f\n",tree->mod->pi[9]);
		PhyML_Printf("L : %f\n",tree->mod->pi[10]);
		PhyML_Printf("K : %f\n",tree->mod->pi[11]);
		PhyML_Printf("M : %f\n",tree->mod->pi[12]);
		PhyML_Printf("F : %f\n",tree->mod->pi[13]);
		PhyML_Printf("P : %f\n",tree->mod->pi[14]);
		PhyML_Printf("S : %f\n",tree->mod->pi[15]);
		PhyML_Printf("T : %f\n",tree->mod->pi[16]);
		PhyML_Printf("W : %f\n",tree->mod->pi[17]);
		PhyML_Printf("Y : %f\n",tree->mod->pi[18]);
		PhyML_Printf("V : %f\n",tree->mod->pi[19]);

		PhyML_Printf("N : %f\n",tree->mod->pi[20]);
		break;
	}
	default : {break;}
	}
}

/*********************************************************/

m3ldbl Num_Derivatives_One_Param(m3ldbl (*func)(arbre *tree), arbre *tree,
		m3ldbl f0, m3ldbl *param, m3ldbl stepsize,
		m3ldbl *err, int precise)
{
	int i,j;
	m3ldbl errt,fac,hh,**a,ans;
	int n_iter;
	a = (m3ldbl **)mCalloc(11,sizeof(m3ldbl *));
	For(i,11) a[i] = (m3ldbl *)mCalloc(11,sizeof(m3ldbl));


	n_iter = 10; /* */

	ans  = .0;

	if (stepsize == 0.0) Warn_And_Exit("\n. h must be nonzero in Dfridr.");

	hh=stepsize;

	if(!precise)
	{

		*param   = *param+hh;
		a[0][0]  = (*func)(tree);
		a[0][0]  -= f0;
		a[0][0]  /= hh;
		*param   = *param-hh;

		ans =  a[0][0];
	}
	else
	{
		*param   = *param+hh;
		a[0][0]  = (*func)(tree);
		/*   *param   = *param-2*hh; */
		/*   a[0][0] -= (*func)(tree); */
		/*   a[0][0] /= (2.0*hh); */
		/*   *param   = *param+hh; */
		a[0][0]  -= f0;
		a[0][0]  /= hh;
		*param   = *param-hh;


		*err=1e30;
		for(i=1;i<n_iter;i++)
		{
			hh /= 1.4;

			/*       *param   = *param+hh; */
			/*       a[0][i]  = (*func)(tree); */
			/*       *param   = *param-2*hh; */
			/*       a[0][i] -= (*func)(tree); */
			/*       a[0][i] /= (2.0*hh); */
			/*       *param   = *param+hh; */


			*param   = *param+hh;
			a[0][i]  = (*func)(tree);
			/*   *param   = *param-2*hh; */
			/*   a[0][i] -= (*func)(tree); */
			/*   a[0][i] /= (2.0*hh); */
			/*   *param   = *param+hh; */
			a[0][i]  -= f0;
			a[0][i]  /= hh;
			*param   = *param-hh;


			fac=1.4*1.4;
			for (j=1;j<=i;j++)
			{
				a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
				fac=1.4*1.4*fac;

				errt=MAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));

				if (errt <= *err)
				{
					*err=errt;
					ans=a[j][i];
				}
			}

			if(fabs(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
		}
	}
	For(i,11) Free(a[i]);
	Free(a);

	return ans;
}

/*********************************************************/

void Num_Derivative_Several_Param(arbre *tree, m3ldbl *param, int n_param, m3ldbl stepsize,
		m3ldbl (*func)(arbre *tree), m3ldbl *derivatives)
{
	int i;
	m3ldbl err,f0;

	f0 = (*func)(tree);

	For(i,n_param)
	{
		derivatives[i] = Num_Derivatives_One_Param(func,
				tree,
				f0,
				param+i,
				stepsize,
				&err,
				0
		);
	}
}

/*********************************************************/

int Compare_Two_States(char *state1, char *state2, int state_size)
{
	/* 1 the two states are identical */
	/* 0 the two states are different */
	int i;

	For(i,state_size) if(state1[i] != state2[i]) break;

	return (i==state_size)?(1):(0);
}

/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
	int i;
	For(i,state_size) to[i] = from[i];
}

/*********************************************************/

model *Make_Model_Basic()
{
	model *mod;

	mod                     = (model *)mCalloc(1,sizeof(model));

	mod->modelname          = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	mod->custom_mod_string  = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
	mod->user_b_freq        = (m3ldbl *)mCalloc(T_MAX_OPTION,sizeof(m3ldbl));

	mod->rr                 = (m3ldbl *)mCalloc(6,sizeof(m3ldbl));
	mod->rr_val             = (m3ldbl *)mCalloc(6,sizeof(m3ldbl));
	mod->rr_num             = (int *)mCalloc(6,sizeof(int *));
	mod->n_rr_per_cat       = (int *)mCalloc(6,sizeof(int));
	mod->s_opt              = (optimiz *)Alloc_Optimiz();

	return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod)
{

	mod->pi             = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl));
	mod->gamma_r_proba  = (m3ldbl *)mCalloc(mod->n_catg,sizeof(m3ldbl));
	mod->gamma_rr       = (m3ldbl *)mCalloc(mod->n_catg,sizeof(m3ldbl));
	mod->pi_unscaled    = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl));

	mod->Pij_rr   = (double *)mCalloc(mod->n_catg*mod->ns*mod->ns,sizeof(double));

	mod->qmat      = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
	mod->qmat_buff = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
	mod->eigen     = (eigen *)Make_Eigen_Struct(mod);

	if(mod->n_rr_branch)
	{
		mod->rr_branch   = (m3ldbl *)mCalloc(mod->n_rr_branch,sizeof(m3ldbl));
		mod->p_rr_branch = (m3ldbl *)mCalloc(mod->n_rr_branch,sizeof(m3ldbl));
	}
}

/*********************************************************/

void Copy_Dist(m3ldbl **cpy, m3ldbl **orig, int n)
{
	int i,j;
	For(i,n) For(j,n) cpy[i][j] = orig[i][j];
}

/*********************************************************/

model *Copy_Model(model *ori)
{
	model *cpy;

	cpy = Make_Model_Basic();

	Copy_Optimiz(ori->s_opt,cpy->s_opt);

	cpy->ns      = ori->ns;
	cpy->n_catg  = ori->n_catg;

	Make_Model_Complete(cpy);

	Record_Model(ori,cpy);

#ifdef M4
	if(ori->m4mod) cpy->m4mod = M4_Copy_M4_Model(ori, ori->m4mod);
#endif

	return cpy;
}

/*********************************************************/

void Record_Model(model *ori, model *cpy)
{
	int i;

	cpy->ns           = ori->ns;
	cpy->n_catg       = ori->n_catg;
	cpy->datatype     = ori->datatype;
	cpy->n_otu        = ori->n_otu;
	cpy->alpha_old    = ori->alpha_old;
	cpy->kappa_old    = ori->alpha_old;
	cpy->lambda_old   = ori->lambda_old;
	cpy->pinvar_old   = ori->pinvar_old;
	cpy->whichmodel   = ori->whichmodel;
	cpy->seq_len      = ori->seq_len;
	cpy->update_eigen = ori->update_eigen;
	cpy->kappa        = ori->kappa;
	cpy->alpha        = ori->alpha;
	cpy->lambda       = ori->lambda;
	cpy->bootstrap    = ori->bootstrap;
	cpy->invar        = ori->invar;
	cpy->pinvar       = ori->pinvar;
	cpy->stepsize     = ori->stepsize;
	cpy->n_diff_rr    = ori->n_diff_rr;

	For(i,6) cpy->rr_num[i] = ori->rr_num[i];

	For(i,6)
	{
		cpy->rr_val[i]  = ori->rr_val[i];
		cpy->rr[i] = cpy->rr[i];
	}

	For(i,cpy->ns)
	{
		cpy->pi[i]          = ori->pi[i];
		cpy->pi_unscaled[i] = ori->pi_unscaled[i];
		cpy->user_b_freq[i] = ori->user_b_freq[i];
	}

	For(i,cpy->ns*cpy->ns) cpy->qmat[i] = ori->qmat[i];

	For(i,cpy->n_catg)
	{
		cpy->gamma_r_proba[i] = ori->gamma_r_proba[i];
		cpy->gamma_rr[i]      = ori->gamma_rr[i];
	}

#ifndef PHYML
	cpy->use_m4mod = ori->use_m4mod;
#endif

	cpy->eigen->size = ori->eigen->size;
	For(i,2*ori->ns)       cpy->eigen->space[i]       = ori->eigen->space[i];
	For(i,2*ori->ns)       cpy->eigen->space_int[i]   = ori->eigen->space_int[i];
	For(i,ori->ns)         cpy->eigen->e_val[i]       = ori->eigen->e_val[i];
	For(i,ori->ns)         cpy->eigen->e_val_im[i]    = ori->eigen->e_val_im[i];
	For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
	For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
	For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect_im[i] = ori->eigen->r_e_vect_im[i];
	For(i,ori->ns*ori->ns) cpy->eigen->l_e_vect[i]    = ori->eigen->l_e_vect[i];
	For(i,ori->ns*ori->ns) cpy->eigen->q[i]           = ori->eigen->q[i];
}

/*********************************************************/

option *Make_Input()
{
	option* io               = (option *)mCalloc(1,sizeof(option));
	io->mod                  = (model *)Make_Model_Basic();

	io->in_seq_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->in_tree_file         = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_tree_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_trees_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_boot_tree_file   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_boot_stats_file  = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_stats_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_lk_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_ps_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->out_trace_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->nt_or_cd             = (char *)mCalloc(T_MAX_FILE,sizeof(char));
	io->run_id_string        = (char *)mCalloc(T_MAX_OPTION,sizeof(char));

	return io;
}

/*********************************************************/
/*
 * JSJ: takes the number of branch lengths per set from io
 * 	and reallocates memory to the props, and resets the default
 * 	value of proportions (equal props).
 *
 */
void Update_Default_Props(option *io){
	int n_l = io->n_l;
	int i;
	For(i,n_l){
		io->props[i] = 1.0/n_l;
	}
	Normalize_Props_IO(io);
}

void Set_Defaults_Input(option* io)
{
	io->fp_in_seq                  = NULL;
	io->fp_in_tree                 = NULL;
	io->fp_out_tree                = NULL;
	io->fp_out_trees               = NULL;
	io->fp_out_boot_tree           = NULL;
	io->fp_out_boot_stats          = NULL;
	io->fp_out_stats               = NULL;
	io->n_l 					   = 1; //JSJ: initialize to 1 branch length set
	io->props[0]				   = 1.0; //JSJ: initialize to all sites falling under the one set
	io->fixed_props                = 1; //JSJ: don't optimize props on single branch length
	io->user_props				   = 0; //JSJ: assume user won't supply proportions...
	io->user_topo				   = 0; //JSJ: assume user won't supply search algo
	io->tree                       = NULL;
	io->mod->datatype              = 0;
	strcpy(io->nt_or_cd,"nucleotides");
	io->n_data_sets                = 1;
	io->interleaved                = 1;
	io->in_tree                    = 0;
	io->out_tree_file_open_mode    = 1;
	io->out_stats_file_open_mode   = 1;
	io->seq_len                    = -1;
	io->n_data_set_asked           = -1;
	io->print_boot_trees           = 1;
	io->n_gt                       = 1;
	io->ratio_test		           = 4;
	io->multigene                  = 0;
	io->config_multigene           = 0;
	io->curr_interface             = 0;
	io->r_seed                     = -1;
	io->collapse_boot              = 0;
	io->random_boot_seq_order      = 1;
	io->print_trace                = 0;
	io->print_site_lnl             = 0;
	io->m4_model                   = NO;
	io->rm_ambigu                  = 0;
	io->compress_seq               = 1;
	io->append_run_ID              = 0;
	io->quiet                      = 0;
}

/*********************************************************/

void Set_Defaults_Model(model *mod)
{
	int i;

	strcpy(mod->modelname,"HKY85");
	strcpy(mod->custom_mod_string,"000000");
	mod->whichmodel              = HKY85;
	mod->n_catg                  = 4;
	mod->kappa                   = 4.0;
	mod->alpha                   = 1.0;
	mod->lambda                  = 1.0;
	mod->bootstrap               = 0;
	mod->invar                   = 0;
	mod->pinvar                  = 0.0;
	mod->stepsize                = 1;
	mod->ns                      = 4;
	mod->n_diff_rr               = 0;
	For(i,6) mod->rr_val[i]      = 1.0;
	For(i,4) mod->user_b_freq[i] = -1.;
	mod->m4mod                   = NULL;
	mod->use_m4mod               = 0;
	mod->n_rr_branch             = 0;
	mod->rr_branch_alpha         = 0.1;
	mod->gamma_median            = 0;
}

/*********************************************************/

void Reset_Model_Parameters(model *mod)
{
	int i;

	PhyML_Printf("\n. Resetting model parameters... \n");

	mod->kappa                   = 4.0;
	mod->alpha                   = 1.0;
	mod->lambda                  = 1.0;
	mod->pinvar                  = 0.0;

	For(i,6)
	{
		mod->rr[i]     = 1.0;
		mod->rr_val[i] = 1.0;
	}
}

/*********************************************************/

void Set_Defaults_Optimiz(optimiz *s_opt)
{
	s_opt->print                = 1;
	s_opt->last_opt             = 1;
	s_opt->opt_alpha            = 1;
	s_opt->opt_kappa            = 1;
	s_opt->opt_bl               = 1;
	s_opt->opt_props			= 0; //JSJ: defaults to 0 because optimizing 1 proportion (default) makes no sense...
	s_opt->opt_lambda           = 0;
	s_opt->opt_pinvar           = 0;
	s_opt->opt_num_param        = 0;
	s_opt->opt_cov_delta        = 0;
	s_opt->opt_cov_alpha        = 0;
	s_opt->opt_cov_free_rates   = 0;
	s_opt->opt_rr               = 0;
	s_opt->init_lk              = UNLIKELY;
	s_opt->n_it_max             = 1000;
	s_opt->opt_topo             = 1;
	s_opt->topo_search          = NNI_MOVE;
	s_opt->random_input_tree    = 0;
	s_opt->n_rand_starts        = 5;
	s_opt->brent_it_max         = 500;
	s_opt->steph_spr            = 1;
	s_opt->user_state_freq      = 0;
	s_opt->min_diff_lk_local    = MIN_DIFF_LK_LOCAL; //JSJ: try playing around with these values
	s_opt->min_diff_lk_global   = MIN_DIFF_LK_GLOBAL;
	s_opt->min_diff_lk_move     = MIN_DIFF_LK_MOVE;
	s_opt->p_moves_to_examine   = 0.15;
	s_opt->fast_nni             = 0;
	s_opt->greedy               = 0;
	s_opt->general_pars         = 0;
	s_opt->tree_size_mult       = 1;
	s_opt->opt_five_branch      = 1;
	s_opt->pars_thresh          = 5;
	s_opt->hybrid_thresh        = 0;
	s_opt->quickdirty           = 0;
	s_opt->spr_pars             = 1;
	s_opt->spr_lnL              = 0;
	s_opt->min_depth_path       = 0;
	s_opt->max_depth_path       = 20;
	s_opt->deepest_path         = 20;
	s_opt->max_delta_lnL_spr    = 50.;

	s_opt->wim_n_rgrft          = -1;
	s_opt->wim_n_globl          = -1;
	s_opt->wim_max_dist         = -1;
	s_opt->wim_n_optim          = -1;
	s_opt->wim_n_best           = -1;
	s_opt->wim_inside_opt       =  0;
}

/*********************************************************/

void Copy_Optimiz(optimiz *ori, optimiz *cpy)
{
	cpy->print                =   ori->print                ;
	cpy->last_opt             =   ori->last_opt             ;
	cpy->opt_alpha            =   ori->opt_alpha            ;
	cpy->opt_kappa            =   ori->opt_kappa            ;
	cpy->opt_bl               =   ori->opt_bl               ;
	cpy->opt_props			  =   ori->opt_props            ;
	cpy->opt_lambda           =   ori->opt_lambda           ;
	cpy->opt_pinvar           =   ori->opt_pinvar           ;
	cpy->opt_num_param        =   ori->opt_num_param        ;
	cpy->opt_cov_delta        =   ori->opt_cov_delta        ;
	cpy->opt_cov_alpha        =   ori->opt_cov_alpha        ;
	cpy->opt_cov_free_rates   =   ori->opt_cov_free_rates   ;
	cpy->opt_rr               =   ori->opt_rr               ;
	cpy->init_lk              =   ori->init_lk              ;
	cpy->n_it_max             =   ori->n_it_max             ;
	cpy->opt_topo             =   ori->opt_topo             ;
	cpy->topo_search          =   ori->topo_search          ;
	cpy->random_input_tree    =   ori->random_input_tree    ;
	cpy->n_rand_starts        =   ori->n_rand_starts        ;
	cpy->brent_it_max         =   ori->brent_it_max         ;
	cpy->steph_spr            =   ori->steph_spr            ;
	cpy->user_state_freq      =   ori->user_state_freq      ;
	cpy->min_diff_lk_local    =   ori->min_diff_lk_local    ;
	cpy->min_diff_lk_global   =   ori->min_diff_lk_global   ;
	cpy->min_diff_lk_move     =   ori->min_diff_lk_move     ;
	cpy->p_moves_to_examine   =   ori->p_moves_to_examine   ;
	cpy->fast_nni             =   ori->fast_nni             ;
	cpy->greedy               =   ori->greedy               ;
	cpy->general_pars         =   ori->general_pars         ;
	cpy->tree_size_mult       =   ori->tree_size_mult       ;
	cpy->opt_five_branch      =   ori->opt_five_branch      ;
	cpy->pars_thresh          =   ori->pars_thresh          ;
	cpy->hybrid_thresh        =   ori->hybrid_thresh        ;
	cpy->quickdirty           =   ori->quickdirty           ;
	cpy->spr_pars             =   ori->spr_pars             ;
	cpy->spr_lnL              =   ori->spr_lnL              ;
	cpy->min_depth_path       =   ori->min_depth_path       ;
	cpy->max_depth_path       =   ori->max_depth_path       ;
	cpy->deepest_path         =   ori->deepest_path         ;
	cpy->max_delta_lnL_spr    =   ori->max_delta_lnL_spr    ;

	cpy->wim_n_rgrft          =   ori->wim_n_rgrft          ;
	cpy->wim_n_globl          =   ori->wim_n_globl          ;
	cpy->wim_max_dist         =   ori->wim_max_dist         ;
	cpy->wim_n_optim          =   ori->wim_n_optim          ;
	cpy->wim_n_best           =   ori->wim_n_best           ;
	cpy->wim_inside_opt       =   ori->wim_inside_opt       ;
}

/*********************************************************/

void Get_Bip(node *a, node *d, arbre *tree)
{
	int i,j;

	if(d->tax)
	{
		if(d->common)
		{
			d->bip_node[0][0] = d;
			d->bip_size[0]    = 1;
			strcpy(d->bip_name[0][0],d->name);

			For(i,3)
			{
				if(a->v[i] == d)
				{
					a->bip_size[i] = 0;
					For(j,tree->n_otu)
					{
						if(strcmp(tree->noeud[j]->name,d->name))
						{
							a->bip_node[i][a->bip_size[i]] = d;
							strcpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name);
							a->bip_size[i]++;
						}
					}
					qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);
					break;
				}
			}
		}
		return;
	}
	else
	{
		int k;
		int d_a;

		d_a = -1;

		For(i,3)
		{
			if(d->v[i] != a) Get_Bip(d,d->v[i],tree);
			else d_a = i;
		}

		d->bip_size[d_a] = 0;
		For(i,3)
		if(d->v[i] != a)
		{
			For(j,3)
			{
				if(d->v[i]->v[j] == d)
				{
					For(k,d->v[i]->bip_size[j])
					{
						d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
						strcpy(d->bip_name[d_a][d->bip_size[d_a]],d->v[i]->bip_node[j][k]->name);
						d->bip_size[d_a]++;
					}
					break;
				}
			}
		}

		qsort(d->bip_name[d_a],d->bip_size[d_a],sizeof(char *),Sort_String);

		For(i,3)
		if(a->v[i] == d)
		{
			a->bip_size[i] = 0;
			For(j,tree->n_otu)
			{
				For(k,d->bip_size[d_a])
				{
					if(d->bip_node[d_a][k] == tree->noeud[j])
						break;
				}

				if((k == d->bip_size[d_a]) && (tree->noeud[j]->common))
					/* 		if(k == d->bip_size[d_a]) */
				{
					a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
					strcpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name);
					a->bip_size[i]++;
				}
			}

			qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);

			/* 	    if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a]) */
			/* 	      { */
			/* 		PhyML_Printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]); */
			/* 		Warn_And_Exit("\n. Problem in counting bipartitions \n"); */
			/* 	      } */
			break;
		}
	}
}

/*********************************************************/

void Alloc_Bip(arbre *tree)
{
	int i,j,k;

	tree->has_bip = 1;

	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
		tree->noeud[i]->bip_node = (node ***)mCalloc(3,sizeof(node **));
		tree->noeud[i]->bip_name = (char ***)mCalloc(3,sizeof(char **));
		For(j,3)
		{
			tree->noeud[i]->bip_node[j] =
					(node **)mCalloc(tree->n_otu,sizeof(node *));

			tree->noeud[i]->bip_name[j] =
					(char **)mCalloc(tree->n_otu,sizeof(char *));

			For(k,tree->n_otu)
			tree->noeud[i]->bip_name[j][k] =
					(char *)mCalloc(T_MAX_NAME,sizeof(char ));
		}
	}
}

/*********************************************************/

int Sort_Phydbl_Increase(const void *a, const void *b)
{
	if((*(m3ldbl *)(a)) <= (*(m3ldbl *)(b))) return -1;
	else return 1;
}

/*********************************************************/

int Sort_String(const void *a, const void *b)
{
	return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

/*********************************************************/

//m3ldbl Compare_Bip_On_Existing_Edges(m3ldbl thresh_len, arbre *tree1, arbre *tree2)
//{
//	int i,j,k,m;
//	edge *b1,*b2;
//	char **bip1,**bip2;
//	int bip_size,n_edges1,n_edges2;
//	m3ldbl rf;
//	int bip_size1, bip_size2;
//
//	n_edges1 = 0;
//	For(i,2*tree1->n_otu-3)
//	{
//		if((!tree1->t_edges[i]->left->tax) &&
//				(!tree1->t_edges[i]->rght->tax) &&
//				(tree1->t_edges[i]->l[0] > thresh_len))
//		{
//			n_edges1++;
//		}
//	}
//	n_edges2 = 0;
//	For(i,2*tree2->n_otu-3)
//	{
//		if((!tree2->t_edges[i]->left->tax) &&
//				(!tree2->t_edges[i]->rght->tax) &&
//				(tree2->t_edges[i]->l[0] > thresh_len))
//		{
//			n_edges2++;
//		}
//	}
//
//
//	rf = 0.0;
//	For(i,2*tree1->n_otu-3)
//	{
//		b1 = tree1->t_edges[i];
//		bip_size1 = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
//
//		if((bip_size1 > 1) && (b1->l[0] > thresh_len))
//		{
//			For(j,2*tree2->n_otu-3)
//			{
//				b2 = tree2->t_edges[j];
//				bip_size2 = MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]);
//
//				if((bip_size2 > 1) && (b2->l[0] > thresh_len))
//				{
//					if(bip_size1 == bip_size2)
//					{
//						bip_size = bip_size1;
//
//						if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
//						{
//							if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
//							{
//								bip1 = b1->left->bip_name[b1->l_r];
//							}
//							else
//							{
//								bip1 = b1->rght->bip_name[b1->r_l];
//							}
//						}
//						else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
//						{
//							bip1 = b1->left->bip_name[b1->l_r];
//						}
//						else
//						{
//							bip1 = b1->rght->bip_name[b1->r_l];
//						}
//
//						if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
//						{
//							if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
//							{
//								bip2 = b2->left->bip_name[b2->l_r];
//							}
//							else
//							{
//								bip2 = b2->rght->bip_name[b2->r_l];
//							}
//						}
//						else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
//						{
//							bip2 = b2->left->bip_name[b2->l_r];
//						}
//						else
//						{
//							bip2 = b2->rght->bip_name[b2->r_l];
//						}
//
//						if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");
//
//
//						For(k,bip_size)
//						{
//							if(strcmp(bip1[k],bip2[k])) break;
//						}
//
//						if(k == bip_size)
//						{
//							b2->bip_score++;
//							b1->bip_score++;
//							rf+=1.0;
//							break;
//						}
//					}
//				}
//			}
//		}
//	}
//
//	PhyML_Printf("\n. rf= %f n_edges1=%d n_edges2=%d",rf,n_edges1,n_edges2);
//	rf /= MIN(n_edges1,n_edges2);
//	return 1-rf;
//}

/*********************************************************/

void Compare_Bip(arbre *tree1, arbre *tree2)
{
	int i,j,k;
	edge *b1,*b2;
	char **bip1,**bip2;
	int bip_size1, bip_size2, bip_size;



	For(i,2*tree1->n_otu-3)
	{
		b1 = tree1->t_edges[i];
		bip_size1 = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);

		if(bip_size1 > 1)
		{
			For(j,2*tree2->n_otu-3)
			{
				b2 = tree2->t_edges[j];
				bip_size2 = MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]);

				if(bip_size2 > 1)
				{
					if(bip_size1 == bip_size2)
					{
						bip_size = bip_size1;

						if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
						{
							if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
							{
								bip1 = b1->left->bip_name[b1->l_r];
							}
							else
							{
								bip1 = b1->rght->bip_name[b1->r_l];
							}
						}
						else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
						{
							bip1 = b1->left->bip_name[b1->l_r];
						}
						else
						{
							bip1 = b1->rght->bip_name[b1->r_l];
						}


						if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
						{
							if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
							{
								bip2 = b2->left->bip_name[b2->l_r];
							}
							else
							{
								bip2 = b2->rght->bip_name[b2->r_l];
							}
						}
						else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
						{
							bip2 = b2->left->bip_name[b2->l_r];
						}
						else
						{
							bip2 = b2->rght->bip_name[b2->r_l];
						}

						if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");


						For(k,bip_size)
						{
							if(strcmp(bip1[k],bip2[k])) break;
						}

						if(k == bip_size)
						{
							b1->bip_score++;
							b2->bip_score++;
							break;
						}
					}
				}
			}
		}
	}
}

/*********************************************************/

void Test_Multiple_Data_Set_Format(option *io)
{
	char *line;

	line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

	io->n_trees = 0;

	while(fgets(line,T_MAX_LINE,io->fp_in_tree)) if(strstr(line,";")) io->n_trees++;

	Free(line);

	if((io->mod->bootstrap > 1) && (io->n_trees > 1))
		Warn_And_Exit("\n. Bootstrap option is not allowed with multiple input trees !\n");

	rewind(io->fp_in_tree);

	return;
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype)
{
	int i,j;
	char a,b;


	if(datatype == NT)
	{
		For(i,stepsize)
		{
			a = statea[i];
			For(j,stepsize)
			{
				b = stateb[j];

				switch(a)
				{
				case 'A':
				{
					switch(b)
					{
					case 'A' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'G':
				{
					switch(b)
					{
					case 'G' :
					case 'R' :
					case 'S' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'C':
				{
					switch(b)
					{
					case 'C' :
					case 'M' :
					case 'S' :
					case 'Y' :
					case 'B' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'T':
				{
					switch(b)
					{
					case 'T' :
					case 'W' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'X' :
					{b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'M' :
				{
					switch(b)
					{
					case 'M' :
					case 'A' :
					case 'C' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' :
					{b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'R' :
				{
					switch(b)
					{
					case 'R' :
					case 'A' :
					case 'G' :
					case 'M' :
					case 'W' :
					case 'S' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}

				case 'W' :
				{
					switch(b)
					{
					case 'W' :
					case 'A' :
					case 'T' :
					case 'M' :
					case 'R' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}

				case 'S' :
				{
					switch(b)
					{
					case 'S' :
					case 'C' :
					case 'G' :
					case 'M' :
					case 'R' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}

				case 'Y' :
				{
					switch(b)
					{
					case 'Y' :
					case 'C' :
					case 'T' :
					case 'M' :
					case 'W' :
					case 'S' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}

				case 'K' :
				{
					switch(b)
					{
					case 'K' :
					case 'G' :
					case 'T' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'B' :
				{
					switch(b)
					{
					case 'B' :
					case 'C' :
					case 'G' :
					case 'T' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'K' :
					case 'D' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'D' :
				{
					switch(b)
					{
					case 'D' :
					case 'A' :
					case 'G' :
					case 'T' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'H' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'H' :
				{
					switch(b)
					{
					case 'H' :
					case 'A' :
					case 'C' :
					case 'T' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'V' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'V' :
				{
					switch(b)
					{
					case 'V' :
					case 'A' :
					case 'C' :
					case 'G' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'X' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				case 'X' :
				{
					switch(b)
					{
					case 'X' :
					case 'A' :
					case 'C' :
					case 'G' :
					case 'T' :
					case 'M' :
					case 'R' :
					case 'W' :
					case 'S' :
					case 'Y' :
					case 'K' :
					case 'B' :
					case 'D' :
					case 'H' :
					case 'V' : {b=b; break;}
					default : return 0;
					}
					break;
				}
				default :
				{
					PhyML_Printf("\n. Err. in Are_Compatible\n");
					PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
					PhyML_Printf("  correspond to existing nucleotides.\n");
					Warn_And_Exit("\n");
					return 0;
				}
				}
			}
		}
	}
	else
	{
		a = statea[0]; b = stateb[0];
		switch(a)
		{
		case 'A' :
		{
			switch(b)
			{
			case 'A' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'R' :
		{
			switch(b)
			{
			case 'R' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'N' :
		{
			switch(b)
			{
			case 'N' :
			case 'B' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'B' :
		{
			switch(b)
			{
			case 'N' :
			case 'B' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'D' :
		{
			switch(b)
			{
			case 'D' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'C' :
		{
			switch(b)
			{
			case 'C' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'Q' :
		{
			switch(b)
			{
			case 'Q' :
			case 'Z' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'Z' :
		{
			switch(b)
			{
			case 'Q' :
			case 'Z' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'E' :
		{
			switch(b)
			{
			case 'E' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'G' :
		{
			switch(b)
			{
			case 'G' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'H' :
		{
			switch(b)
			{
			case 'H' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'I' :
		{
			switch(b)
			{
			case 'I' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'L' :
		{
			switch(b)
			{
			case 'L' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'K' :
		{
			switch(b)
			{
			case 'K' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'M' :
		{
			switch(b)
			{
			case 'M' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'F' :
		{
			switch(b)
			{
			case 'F' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'P' :
		{
			switch(b)
			{
			case 'P' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'S' :
		{
			switch(b)
			{
			case 'S' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'T' :
		{
			switch(b)
			{
			case 'T' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'W' :
		{
			switch(b)
			{
			case 'W' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'Y' :
		{
			switch(b)
			{
			case 'Y' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'V' :
		{
			switch(b)
			{
			case 'V' :
			case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		case 'X' :
		{
			switch(b)
			{
			case 'A':case 'R':case 'N' :case 'B' :case 'D' :
			case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
			case 'H':case 'I':case 'L' :case 'K' :case 'M' :
			case 'F':case 'P':case 'S' :case 'T' :case 'W' :
			case 'Y':case 'V': case 'X' : {b=b; break;}
			default : return 0;
			}
			break;
		}
		default :
		{
			PhyML_Printf("\n. Err. in Are_Compatible\n");
			PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
			PhyML_Printf("  correspond to existing amino-acids.\n");
			Warn_And_Exit("\n");
			return 0;
		}
		}
	}
	return 1;
}

/*********************************************************/

void Hide_Ambiguities(allseq *data)
{
	int i;
	For(i,data->crunch_len) if(data->ambigu[i]) data->wght[i] = 0;
}

/*********************************************************/

void Copy_Tree(arbre *ori, arbre *cpy)
{
	int i,j,m;
	cpy->n_l = ori->n_l;

	//JSJ: copy rate proportions
	For(i,ori->n_l){
		cpy->props[i] = ori->props[i];
	}

	For(i,2*ori->n_otu-2)
	{
		For(j,3)
		{
			if(ori->noeud[i]->v[j])
			{
				cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
				For(m,ori->n_l) cpy->noeud[i]->l[m][j] = ori->noeud[i]->l[m][j];//JSJ: deep copy
				cpy->noeud[i]->n_l = ori->n_l;
				cpy->noeud[i]->b[j] = cpy->t_edges[ori->noeud[i]->b[j]->num];

			}
			else
			{
				cpy->noeud[i]->v[j] = NULL;
				cpy->noeud[i]->b[j] = NULL;
			}
		}
	}

	For(i,2*ori->n_otu-3)
	{
		For(j,ori->n_l){
			cpy->t_edges[i]->l[j]    = ori->t_edges[i]->l[j];//JSJ: deep copy
			cpy->t_edges[i]->best_l[j] = ori->t_edges[i]->best_l[j];
			cpy->t_edges[i]->l_old[j] = ori->t_edges[i]->l_old[j];
			cpy->t_edges[i]->has_zero_br_len[j] = ori->t_edges[i]->has_zero_br_len[j];
		}
		cpy->t_edges[i]->left = cpy->noeud[ori->t_edges[i]->left->num];
		cpy->t_edges[i]->rght = cpy->noeud[ori->t_edges[i]->rght->num];
		cpy->t_edges[i]->l_v1 = ori->t_edges[i]->l_v1;
		cpy->t_edges[i]->l_v2 = ori->t_edges[i]->l_v2;
		cpy->t_edges[i]->r_v1 = ori->t_edges[i]->r_v1;
		cpy->t_edges[i]->r_v2 = ori->t_edges[i]->r_v2;
		cpy->t_edges[i]->l_r  = ori->t_edges[i]->l_r;
		cpy->t_edges[i]->r_l  = ori->t_edges[i]->r_l;
	}

	For(i,ori->n_otu)
	{
		cpy->noeud[i]->tax = 1;
		strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
	}

	cpy->num_curr_branch_available = 0;
	/*   Connect_Edges_To_Nodes_Recur(cpy->noeud[0],cpy->noeud[0]->v[0],cpy); */
	/*   Update_Dirs(cpy); */
}

/*********************************************************/

void Prune_Subtree(node *a, node *d, edge **target, edge **residual, arbre *tree)
{
	node *v1, *v2;
	edge *b1, *b2;
	int dir_v1, dir_v2;
	int i;
	/*   plkflt ***buff_p_lk; */
	plkflt *buff_p_lk;
	plkflt *buff_scale;
	int *buff_p_pars, *buff_pars;
	unsigned int *buff_ui;
	short int *buff_p_lk_tip;

	if(a->tax)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	dir_v1 = dir_v2 = -1;
	For(i,3)
	{
		if(a->v[i] != d)
		{
			if(dir_v1 < 0) dir_v1 = i;
			else           dir_v2 = i;
		}
	}

	if(a->v[dir_v1]->num < a->v[dir_v2]->num)
	{
		v1 = a->v[dir_v1];
		v2 = a->v[dir_v2];
		b1 = a->b[dir_v1];
		b2 = a->b[dir_v2];
	}
	else
	{
		v1 = a->v[dir_v2];
		v2 = a->v[dir_v1];
		b1 = a->b[dir_v2];
		b2 = a->b[dir_v1];
	}

	if(v1->tax && v2->tax) PhyML_Printf("\n. Pruning is meaningless here.\n");

	a->v[dir_v1] = NULL;
	a->v[dir_v2] = NULL;
	a->b[dir_v1] = NULL;
	a->b[dir_v2] = NULL;

	if(v1 == b1->left)
	{
		b1->rght = v2;

		if(v2 == b2->left)
		{
			buff_p_lk            = b1->p_lk_rght;
			b1->p_lk_rght        = b2->p_lk_left;
			b2->p_lk_left        = buff_p_lk;

			buff_p_lk_tip        = b1->p_lk_tip_r;
			b1->p_lk_tip_r       = b2->p_lk_tip_l;
			b2->p_lk_tip_l       = buff_p_lk_tip;

			buff_scale           = b1->sum_scale_f_rght;
			b1->sum_scale_f_rght = b2->sum_scale_f_left;
			b2->sum_scale_f_left = buff_scale;

			buff_pars            = b1->pars_r;
			b1->pars_r           = b2->pars_l;
			b2->pars_l           = buff_pars;

			buff_ui              = b1->ui_r;
			b1->ui_r             = b2->ui_l;
			b2->ui_l             = buff_ui;

			buff_p_pars          = b1->p_pars_r;
			b1->p_pars_r         = b2->p_pars_l;
			b2->p_pars_l         = buff_p_pars;
		}
		else
		{
			buff_p_lk            = b1->p_lk_rght; /* b1->p_lk_rght = NULL if b1->rght->tax */
			b1->p_lk_rght        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
			b2->p_lk_rght        = buff_p_lk;

			buff_p_lk_tip        = b1->p_lk_tip_r;
			b1->p_lk_tip_r       = b2->p_lk_tip_r;
			b2->p_lk_tip_r       = buff_p_lk_tip;

			buff_scale           = b1->sum_scale_f_rght;
			b1->sum_scale_f_rght = b2->sum_scale_f_rght;
			b2->sum_scale_f_rght = buff_scale;

			buff_pars            = b1->pars_r;
			b1->pars_r           = b2->pars_r;
			b2->pars_r           = buff_pars;

			buff_ui              = b1->ui_r;
			b1->ui_r             = b2->ui_r;
			b2->ui_r             = buff_ui;

			buff_p_pars          = b1->p_pars_r;
			b1->p_pars_r         = b2->p_pars_r;
			b2->p_pars_r         = buff_p_pars;
		}
	}
	else
	{
		b1->left = v2;

		if(v2 == b2->left)
		{
			buff_p_lk            = b1->p_lk_left;
			b1->p_lk_left        = b2->p_lk_left;
			b2->p_lk_left        = buff_p_lk;

			buff_p_lk_tip        = b1->p_lk_tip_l;
			b1->p_lk_tip_l       = b2->p_lk_tip_l;
			b2->p_lk_tip_l       = buff_p_lk_tip;

			buff_scale           = b1->sum_scale_f_left;
			b1->sum_scale_f_left = b2->sum_scale_f_left;
			b2->sum_scale_f_left = buff_scale;

			buff_pars            = b1->pars_l;
			b1->pars_l           = b2->pars_l;
			b2->pars_l           = buff_pars;

			buff_ui              = b1->ui_l;
			b1->ui_l             = b2->ui_l;
			b2->ui_l             = buff_ui;

			buff_p_pars          = b1->p_pars_l;
			b1->p_pars_l         = b2->p_pars_l;
			b2->p_pars_l         = buff_p_pars;
		}
		else
		{
			buff_p_lk            = b1->p_lk_left;
			b1->p_lk_left        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
			b2->p_lk_rght        = buff_p_lk;

			buff_p_lk_tip        = b1->p_lk_tip_l;
			b1->p_lk_tip_l       = b2->p_lk_tip_r;
			b2->p_lk_tip_r       = buff_p_lk_tip;

			buff_scale           = b1->sum_scale_f_left;
			b1->sum_scale_f_left = b2->sum_scale_f_rght;
			b2->sum_scale_f_rght = buff_scale;

			buff_pars            = b1->pars_l;
			b1->pars_l           = b2->pars_r;
			b2->pars_r           = buff_pars;

			buff_ui              = b1->ui_l;
			b1->ui_l             = b2->ui_r;
			b2->ui_r             = buff_ui;

			buff_p_pars          = b1->p_pars_l;
			b1->p_pars_l         = b2->p_pars_r;
			b2->p_pars_r         = buff_p_pars;
		}
	}

	For(i,3)
	if(v2->v[i] == a)
	{
		v2->v[i] = v1;
		v2->b[i] = b1;
		break;
	}

#ifdef DEBUG
	if(i == 3)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
#endif

	For(i,3)
	if(v1->v[i] == a)
	{
		v1->v[i] = v2;
		break;
	}

#ifdef DEBUG
	if(i == 3)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}
#endif
	//JSJ: do for array of bls
	For(i,tree->n_l){
		b1->l[i] += b2->l[i];
	}

	(v1 == b1->left)?
			(Make_Edge_Dirs(b1,v1,v2)):
				(Make_Edge_Dirs(b1,v2,v1));

			if(target)   (*target)   = b1;
			if(residual) (*residual) = b2;


#ifdef DEBUG
			if(b1->left->tax)
			{
				printf("\n. b1->left->num = %d",b1->left->num);
				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
				Warn_And_Exit("");
			}
#endif


}

/*********************************************************/

void Graft_Subtree(edge *target, node *link, edge *residual, arbre *tree)
{
	node *v1, *v2;
	int i, dir_v1, dir_v2;
	plkflt *buff_p_lk;
	plkflt *buff_scale;
	int *buff_p_pars, *buff_pars;
	short int *buff_p_lk_tip;
	unsigned int *buff_ui;
	edge *b_up;

	dir_v1 = dir_v2 = -1;
	b_up = NULL;
	For(i,3)
	{
		if(!link->v[i])
		{
			if(dir_v1 < 0) dir_v1 = i;
			else           dir_v2 = i;
		}
		else b_up = link->b[i];
	}

	if(target->left->num < target->rght->num)
	{
		v1                           = target->left;
		v2                           = target->rght;

		buff_p_lk                    = residual->p_lk_rght;
		residual->p_lk_rght          = target->p_lk_rght;
		target->p_lk_rght            = buff_p_lk;

		buff_p_lk_tip                = residual->p_lk_tip_r;
		residual->p_lk_tip_r         = target->p_lk_tip_r;
		target->p_lk_tip_r           = buff_p_lk_tip;

		buff_scale                   = residual->sum_scale_f_rght;
		residual->sum_scale_f_rght   = target->sum_scale_f_rght;
		target->sum_scale_f_rght     = buff_scale;

		buff_pars                    = residual->pars_r;
		residual->pars_r             = target->pars_r;
		target->pars_r               = buff_pars;

		buff_ui                      = residual->ui_r;
		residual->ui_r               = target->ui_r;
		target->ui_r                 = buff_ui;

		buff_p_pars                  = residual->p_pars_r;
		residual->p_pars_r           = target->p_pars_r;
		target->p_pars_r             = buff_p_pars;
	}
	else
	{
		v1                           = target->rght;
		v2                           = target->left;

		buff_p_lk                    = residual->p_lk_rght;
		residual->p_lk_rght          = target->p_lk_left;
		target->p_lk_left            = buff_p_lk;

		buff_p_lk_tip                = residual->p_lk_tip_r;
		residual->p_lk_tip_r         = target->p_lk_tip_l;
		target->p_lk_tip_l           = buff_p_lk_tip;

		buff_scale                   = residual->sum_scale_f_rght;
		residual->sum_scale_f_rght   = target->sum_scale_f_left;
		target->sum_scale_f_left     = buff_scale;

		buff_pars                    = residual->pars_r;
		residual->pars_r             = target->pars_l;
		target->pars_l               = buff_pars;

		buff_ui                      = residual->ui_r;
		residual->ui_r               = target->ui_l;
		target->ui_l                 = buff_ui;

		buff_p_pars                  = residual->p_pars_r;
		residual->p_pars_r           = target->p_pars_l;
		target->p_pars_l             = buff_p_pars;
	}

	For(i,3)
	if(v2->b[i] == target)
	{
		v2->v[i] = link;
		v2->b[i] = residual;
		break;
	}

	link->v[dir_v2] = v2;
	link->b[dir_v2] = residual;

	residual->left  = link;
	residual->rght  = v2;

	(v1 == target->left)?(target->rght = link):(target->left = link);

	link->v[dir_v1] = v1;
	link->b[dir_v1] = target;

	For(i,3)
	if(v1->v[i] == v2)
	{
		v1->v[i] = link;
		break;
	}
	//JSJ: operate on array of edge lengths
	For(i,tree->n_l){
		target->l[i] /= 2.;
		residual->l[i] = target->l[i];
	}

	Make_Edge_Dirs(target,target->left,target->rght);
	Make_Edge_Dirs(residual,residual->left,residual->rght);
	Make_Edge_Dirs(b_up,b_up->left,b_up->rght);
}

/*********************************************************/

void Reassign_Node_Nums(node *a, node *d, int *curr_ext_node, int *curr_int_node, arbre *tree)
{
	node *buff;
	int i;

	if(a->tax)
	{
		buff = tree->noeud[*curr_ext_node];
		tree->noeud[*curr_ext_node] = a;
		tree->noeud[a->num] = buff;
		buff->num = a->num;
		a->num = *curr_ext_node;
		(*curr_ext_node)++;
	}

	if(d->tax)
	{
		buff = tree->noeud[*curr_ext_node];
		tree->noeud[*curr_ext_node] = d;
		tree->noeud[d->num] = buff;
		buff->num = d->num;
		d->num = *curr_ext_node;
		(*curr_ext_node)++;
		return;
	}
	else
	{
		buff = tree->noeud[*curr_int_node];
		tree->noeud[*curr_int_node] = d;
		tree->noeud[d->num] = buff;
		buff->num = d->num;
		d->num = *curr_int_node;
		(*curr_int_node)++;
	}

	For(i,3)
	{
		if(d->v[i] != a)
			Reassign_Node_Nums(d,d->v[i],curr_ext_node,curr_int_node,tree);
	}
}

/*********************************************************/

void Reassign_Edge_Nums(node *a, node *d, int *curr_br, arbre *tree)
{
	edge *buff;
	int i,j;

	For(i,3)
	if(a->v[i] == d)
	{
		buff = tree->t_edges[*curr_br];
		For(j,2*N_MAX_OTU-3) if(tree->t_edges[j] == a->b[i]) break;
		if(j == 2*N_MAX_OTU-3)
		{
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}
		tree->t_edges[*curr_br] = a->b[i];
		tree->t_edges[j] = buff;
		a->b[i]->num = *curr_br;
		(*curr_br)++;
		break;
	}

	if(d->tax) return;
	else
	{
		For(i,3)
		if(d->v[i] != a)
			Reassign_Edge_Nums(d,d->v[i],curr_br,tree);
	}
}

/*********************************************************/

void Make_List_Of_Reachable_Tips(arbre *tree)
{
	int i,j;

	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i]->list_of_reachable_tips = (node ***)mCalloc(3,sizeof(node **));
		tree->noeud[i]->n_of_reachable_tips    = (int *)mCalloc(3,sizeof(int));
		For(j,3)
		tree->noeud[i]->list_of_reachable_tips[j] = (node **)mCalloc(tree->n_otu,sizeof(node *));
	}
}

/*********************************************************/

void Get_List_Of_Reachable_Tips(arbre *tree)
{
	int i,j;

	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i]->n_of_reachable_tips[0] = 0;
		tree->noeud[i]->n_of_reachable_tips[1] = 0;
		tree->noeud[i]->n_of_reachable_tips[2] = 0;
		For(j,tree->n_otu)
		{
			tree->noeud[i]->list_of_reachable_tips[0][j] = NULL;
			tree->noeud[i]->list_of_reachable_tips[1][j] = NULL;
			tree->noeud[i]->list_of_reachable_tips[2][j] = NULL;
		}
	}

	Get_List_Of_Reachable_Tips_Post(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree);
	Get_List_Of_Reachable_Tips_Pre(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree);
}

/*********************************************************/

void Get_List_Of_Reachable_Tips_Post(node *a, node *d, arbre *tree)
{
	int i,j,k,cpt;

	if(d->tax)
	{
		For(i,3)
		if(a->v[i] == d)
		{
			a->list_of_reachable_tips[i][0] = d;
			a->n_of_reachable_tips[i]       = 1;
			break;
		}
		return;
	}
	else
	{
		For(i,3)
		if(d->v[i] != a)
			Get_List_Of_Reachable_Tips_Post(d,d->v[i],tree);

		For(i,3)
		{
			if(a->v[i] == d)
			{
				a->n_of_reachable_tips[i] = 0;
				cpt                       = 0;
				For(j,3)
				{
					if(d->v[j] != a)
					{
						For(k,d->n_of_reachable_tips[j])
						{
							a->list_of_reachable_tips[i][cpt] = d->list_of_reachable_tips[j][k];
							a->n_of_reachable_tips[i]++;
							cpt++;
						}
					}
				}
				break;
			}
		}
	}
}

/*********************************************************/

void Get_List_Of_Reachable_Tips_Pre(node *a, node *d, arbre *tree)
{
	int i,j,k,cpt;

	For(i,3)
	{
		if(d->v[i] == a)
		{
			if(a->tax)
			{
				d->list_of_reachable_tips[i][0] = a;
				d->n_of_reachable_tips[i]       = 1;
			}
			else
			{
				d->n_of_reachable_tips[i] = 0;
				cpt = 0;
				For(j,3)
				{
					if(a->v[j] != d)
					{
						For(k,a->n_of_reachable_tips[j])
						{
							d->list_of_reachable_tips[i][cpt] = a->list_of_reachable_tips[j][k];
							d->n_of_reachable_tips[i]++;
							cpt++;
						}
					}
				}
			}
			break;
		}
	}

	if(d->tax) return;
	else
	{
		For(i,3)
		if(d->v[i] != a)
			Get_List_Of_Reachable_Tips_Pre(d,d->v[i],tree);

	}
}

/*********************************************************/

int Compare_List_Of_Reachable_Tips(node **list1, int size_list1, node **list2, int size_list2)
{
	int i,j,n_matches;

	n_matches = 0;
	For(i,size_list1)
	{
		For(j,size_list2)
		{
			if(list1[i] == list2[j])
			{
				n_matches++;
			}
		}
	}
	return n_matches;
}

/*********************************************************/

void Find_Mutual_Direction(node *n1, node *n2, int *dir_n1_to_n2, int *dir_n2_to_n1)
{
	int scores[3][3];
	int n_zero_line, n_zero_col;
	int i,j;

	For(i,3) For(j,3) scores[i][j] = 0;

	For(i,3)
	{
		For(j,3)
		{
			scores[i][j] = Compare_List_Of_Reachable_Tips(n1->list_of_reachable_tips[i],
					n1->n_of_reachable_tips[i],
					n2->list_of_reachable_tips[j],
					n2->n_of_reachable_tips[j]);
		}
	}

	For(i,3)
	{
		n_zero_line = 0;
		For(j,3)
		{
			if(!scores[i][j]) n_zero_line++;
		}
		if(n_zero_line != 2) {*dir_n1_to_n2 = i; break;}
	}


	For(i,3)
	{
		n_zero_col = 0;
		For(j,3)
		{
			if(!scores[j][i]) n_zero_col++;
		}
		if(n_zero_col != 2) {*dir_n2_to_n1 = i; break;}
	}

}

/*********************************************************/

void Fill_Dir_Table(arbre *tree)
{
	int i,j,k,l;
	int found;

	Get_List_Of_Reachable_Tips(tree);

	For(i,tree->n_otu) For(j,2*tree->n_otu-2) tree->t_dir[i][j] = 0;

	for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
		For(j,tree->n_otu)
		{
		found = 0;
		For(k,3)
		{
			For(l,tree->noeud[i]->n_of_reachable_tips[k])
			{
				if(tree->noeud[i]->list_of_reachable_tips[k][l] == tree->noeud[j])
				{
					found = 1;
					tree->t_dir[i][j] = k;
					break;
				}
			}
			if(found) break;
		}
		}

	for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
		for(j=i;j<2*tree->n_otu-2;j++)
		{
			Find_Mutual_Direction(tree->noeud[i],tree->noeud[j],
					&(tree->t_dir[i][j]),
					&(tree->t_dir[j][i]));
		}
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

int Get_Subtree_Size(node *a, node *d)
{
	int size,i;

	if(d->tax) return 1;
	else
	{
		size = 0;
		For(i,3)
		if(d->v[i] != a)
			size += Get_Subtree_Size(d,d->v[i]);
	}
	return size;
}

/*********************************************************/
/**
* JSJ: This method is in need of attention. This is what we can use to
* modify branch lengths.
* @param b
* @param tree
* @param approx
*/
void Fast_Br_Len(edge *b, arbre *tree, int approx)
{
	m3ldbl sum;
	m3ldbl *prob, *F;
	int i, j, k, m, site;
	m3ldbl v_rght;
	int dim1,dim2,dim3;
	//m3ldbl eps_bl,old_l,new_l;
	int n_iter;

	n_iter = 0;
	dim1   = tree->mod->ns * tree->mod->n_catg;
	dim2   = tree->mod->ns ;
	dim3   = tree->mod->ns * tree->mod->ns;
	//eps_bl = BL_MIN;

	F    = tree->triplet_struct->F_bc;
	prob = tree->triplet_struct->F_cd;

	Update_PMat_At_Given_Edge(b,tree);

	For(i,dim1*dim2) F[i] = .0;

	For(site,tree->n_pattern)
	{
		/* Joint probabilities of the states at the two ends of the edge */
		v_rght = -1.;
		For(i,tree->mod->ns)
		{
			For(j,tree->mod->ns)
			{
				For(k,tree->mod->n_catg)
				{ //JSJ: temp fix of Pij_rr
					v_rght = (b->rght->tax)?((m3ldbl)(b->p_lk_tip_r[site*dim2+j])):(b->p_lk_rght[site*dim1+k*dim2+j]);
					For(m,tree->n_l){ //JSJ: made the summed probability include a mixed model
						prob[dim3*k+dim2*i+j]              +=
								tree->mod->gamma_r_proba[k]      *
								tree->mod->pi[i]                 *
								b->Pij_rr[m][k*dim3+i*dim2+j]    *
								b->p_lk_left[site*dim1+k*dim2+i] *
								v_rght							 *
								tree->props[m];
					}
				}
			}
		}

		/* Scaling */
		sum = .0;
		For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) sum += prob[dim3*k+dim2*i+j];
		For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) prob[dim3*k+dim2*i+j] /= sum;

		/* Expected number of each pair of states */
		For(i,tree->mod->ns) For(j,tree->mod->ns) For(k,tree->mod->n_catg)
		F[dim3*k+dim2*i+j] += tree->data->wght[site] * prob[dim3*k+dim2*i+j];
	}
	//JSJ: More temp fixes...
	//old_l = b->l[0];
	/**
	* Figure out how to fix the following function so that it works with a set of bls!
	*/
	Opt_Dist_F(b->l,F,tree->mod,tree); //JSJ: decide how to get array of bls to this fxn, one at a time, or simult
	//new_l = b->l[0];
	n_iter++;

	For(m,tree->n_l){
		if(b->l[m] < BL_MIN)      b->l[m] = BL_MIN;
		else if(b->l[m] > BL_MAX) b->l[m] = BL_MAX;
	}

	if(!approx){//JSJ: restored brent temporarily
		m3ldbl min,max;
		//		m3ldbl min[MAX_BL_SET];
		//		m3ldbl max[MAX_BL_SET];
		//
		//		For(m,tree->n_l){
		//			min[m] = b->l[m];
		//			max[m] = b->l[m];
		//			min[m] *= 0.02;
		//			max[m] *= 50.0;
		//		}
		For(m,tree->n_l){
			min = max = b->l[m];
			min *= 0.02;
			max *= 50.0;
			Br_Len_Brent_Iter(min,b->l[m],max,
					tree->mod->s_opt->min_diff_lk_local,
					b,tree,
					tree->mod->s_opt->brent_it_max,
					tree->mod->s_opt->quickdirty,m);
		}
	}
	else{
		Lk_At_Given_Edge(b,tree);
	}
}


/*********************************************************/

eigen *Make_Eigen_Struct(model *mod)
{
	eigen *eig;


	eig              = (eigen *)mCalloc(1,sizeof(eigen));
	eig->size        = mod->ns;
	eig->space       = (double *)mCalloc(2*mod->ns,sizeof(double));
	eig->space_int   = (int *)mCalloc(2*mod->ns,sizeof(int));
	eig->e_val       = (double *)mCalloc(mod->ns,sizeof(double));
	eig->e_val_im    = (double *)mCalloc(mod->ns,sizeof(double));
	eig->r_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
	eig->r_e_vect_im = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
	eig->l_e_vect    = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));
	eig->q           = (double *)mCalloc(mod->ns*mod->ns,sizeof(double));

	return eig;
}

/*********************************************************/

triplet *Make_Triplet_Struct(model *mod)
{
	int i,j,k;
	triplet *triplet_struct;

	triplet_struct                  = (triplet *)mCalloc(1,sizeof(triplet));
	triplet_struct->size            = mod->ns;
	triplet_struct->pi_bc           = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl ));
	triplet_struct->pi_cd           = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl ));
	triplet_struct->pi_bd           = (m3ldbl *)mCalloc(mod->ns,sizeof(m3ldbl ));
	triplet_struct->F_bc            = (m3ldbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(m3ldbl));
	triplet_struct->F_cd            = (m3ldbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(m3ldbl));
	triplet_struct->F_bd            = (m3ldbl *)mCalloc(mod->ns*mod->ns,sizeof(m3ldbl));
	triplet_struct->core            = (m3ldbl ****)mCalloc(mod->n_catg,sizeof(m3ldbl ***));
	triplet_struct->p_one_site      = (m3ldbl ***)mCalloc(mod->ns,sizeof(m3ldbl **));
	triplet_struct->sum_p_one_site  = (m3ldbl ***)mCalloc(mod->ns,sizeof(m3ldbl **));
	triplet_struct->eigen_struct    = (eigen *)Make_Eigen_Struct(mod);
	triplet_struct->mod             = mod;

	For(k,mod->n_catg)
	{
		triplet_struct->core[k]                = (m3ldbl ***)mCalloc(mod->ns,sizeof(m3ldbl **));
		For(i,mod->ns)
		{
			triplet_struct->core[k][i]         = (m3ldbl **)mCalloc(mod->ns,sizeof(m3ldbl *));
			For(j,mod->ns)
			triplet_struct->core[k][i][j]    = (m3ldbl  *)mCalloc(mod->ns,sizeof(m3ldbl ));
		}
	}

	For(i,mod->ns)
	{
		triplet_struct->p_one_site[i]          = (m3ldbl **)mCalloc(mod->ns,sizeof(m3ldbl *));
		For(j,mod->ns)
		triplet_struct->p_one_site[i][j]     = (m3ldbl  *)mCalloc(mod->ns,sizeof(m3ldbl ));
	}

	For(i,mod->ns)
	{
		triplet_struct->sum_p_one_site[i]      = (m3ldbl **)mCalloc(mod->ns,sizeof(m3ldbl *));
		For(j,mod->ns)
		triplet_struct->sum_p_one_site[i][j] = (m3ldbl  *)mCalloc(mod->ns,sizeof(m3ldbl ));
	}
	return triplet_struct;

}

/*********************************************************/

m3ldbl Triple_Dist(node *a, arbre *tree, int approx)
{
	if(a->tax) return UNLIKELY;
	else
	{
		Update_PMat_At_Given_Edge(a->b[1],tree);
		Update_PMat_At_Given_Edge(a->b[2],tree);

		Update_P_Lk(tree,a->b[0],a);
		Fast_Br_Len(a->b[0],tree,approx);
		/*       Br_Len_Brent (BL_MAX, a->b[0]->l,BL_MIN, 1.e-10,a->b[0],tree,50,0); */


		Update_P_Lk(tree,a->b[1],a);
		Fast_Br_Len(a->b[1],tree,approx);
		/*       Br_Len_Brent (BL_MAX, a->b[1]->l,BL_MIN, 1.e-10,a->b[1],tree,50,0); */


		Update_P_Lk(tree,a->b[2],a);
		Fast_Br_Len(a->b[2],tree,approx);
		/*       Br_Len_Brent (BL_MAX, a->b[2]->l,BL_MIN, 1.e-10,a->b[2],tree,50,0); */

		Update_P_Lk(tree,a->b[1],a);
		Update_P_Lk(tree,a->b[0],a);
	}

	return tree->c_lnL;

}


/*********************************************************/

void Make_Symmetric(m3ldbl **F, int size)
{
	int i,j;

	For(i,size)
	{
		for(j=i+1;j<size;j++)
		{
			(*F)[size*i+j] = ((*F)[size*i+j] + (*F)[size*j+i])/2.;
			(*F)[size*j+i] = (*F)[size*i+j];
		}
	}
}

/*********************************************************/

void Round_Down_Freq_Patt(m3ldbl **F, arbre *tree)
{
	int i,j;

	For(i,tree->mod->ns)
	{
		For(j,tree->mod->ns)
		{
			(*F)[tree->mod->ns*i+j] = rint((*F)[tree->mod->ns*i+j]);
		}
	}
}

/*********************************************************/

m3ldbl Get_Sum_Of_Cells(m3ldbl *F, arbre *tree)
{
	int i,j;
	m3ldbl sum = .0;

	For(i,tree->mod->ns)
	For(j,tree->mod->ns)
	sum += F[tree->mod->ns*i+j];

	return sum;
}


/*********************************************************/

void Divide_Cells(m3ldbl **F, m3ldbl div, arbre *tree)
{
	int i,j;

	For(i,tree->mod->ns)
	For(j,tree->mod->ns)
	(*F)[tree->mod->ns*i+j] /= div;
}

/*********************************************************/

void Divide_Mat_By_Vect(m3ldbl **F, m3ldbl *vect, int size)
{
	int i,j;
	For(i,size)
	For(j,size)
	(*F)[size*i+j] = (*F)[size*i+j] / vect[j];
}

/*********************************************************/

void Multiply_Mat_By_Vect(m3ldbl **F, m3ldbl *vect, int size)
{
	int i,j;
	For(i,size)
	For(j,size)
	(*F)[size*i+j] = (*F)[size*i+j] * vect[j];
}

/*********************************************************/

void Found_In_Subtree(node *a, node *d, node *target, int *match, arbre *tree)
{
	if(d->tax) return;
	else
	{
		int i;
		if(d == target) *match =  1;
		For(i,3)
		{
			if(d->v[i] != a)
				Found_In_Subtree(d,d->v[i],target,match,tree);
		}
	}
}

/*********************************************************/

void Get_List_Of_Target_Edges(node *a, node *d, edge **list, int *list_size, arbre *tree)
{
	int i;

	For(i,3)
	{
		if(a->v[i] && a->v[i] == d)
		{
			list[*list_size] = a->b[i];
			(*list_size)++;
		}
	}

	if(d->tax) return;
	else
	{
		For(i,3)
		{
			if(d->v[i] != a)
				Get_List_Of_Target_Edges(d,d->v[i],list,list_size,tree);
		}
	}
}

/*********************************************************/

void Fix_All(arbre *tree)
{
	int i,j;

	tree->mod->pinvar_old = tree->mod->pinvar;
	tree->mod->alpha_old  = tree->mod->alpha;
	tree->mod->kappa_old  = tree->mod->kappa;
	tree->mod->lambda_old = tree->mod->lambda;
	For(j,tree->n_l){
		for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
		{
			tree->noeud[i]->b[0]->l_old[j] = tree->noeud[i]->b[0]->l[j];
			tree->noeud[i]->b[1]->l_old[j] = tree->noeud[i]->b[1]->l[j];
			tree->noeud[i]->b[2]->l_old[j] = tree->noeud[i]->b[2]->l[j];
		}
	}
}

/*********************************************************/
//JSJ: fixed to record array of bls
void Record_Br_Len(m3ldbl **where, arbre *tree)
{
	int i,j;

	if(!where)
	{
		For(j,tree->n_l){
			For(i,2*tree->n_otu-3){
				tree->t_edges[i]->l_old[j] = tree->t_edges[i]->l[j];
			}
		}
	}
	else
	{
		For(j,tree->n_l){
			For(i,2*tree->n_otu-3){
				where[j][i] = tree->t_edges[i]->l[j];
			}
		}
	}
}

/*********************************************************/
//JSJ: changed so that an array of lengths are coppied
void Restore_Br_Len(m3ldbl **from, arbre *tree)
{
	int i,j;

	if(!from)
	{
		For(i,2*tree->n_otu-3){
			For(j,tree->n_l){
				tree->t_edges[i]->l[j] = tree->t_edges[i]->l_old[j];
			}
		}
	}
	else
	{
		For(i,2*tree->n_otu-3){
			For(j,tree->n_l) {
				tree->t_edges[i]->l[j] = from[j][i];
			}
		}
	}
}

/*********************************************************/
//JSJ: probably will eventually have to change
void Get_Dist_Btw_Edges(node *a, node *d, arbre *tree)
{
	int i;
	edge *b_fcus;

	b_fcus = NULL;
	For(i,3) if(a->v[i] == d) {b_fcus = a->b[i]; break;}

	if(d->tax) return;
	else
	{
		For(i,3)
		if(d->v[i] != a)
		{ //JSJ: Another temporary fix for compilation
			d->b[i]->topo_dist_btw_edges = b_fcus->topo_dist_btw_edges + 1;
			d->b[i]->dist_btw_edges      = b_fcus->dist_btw_edges + d->b[i]->l[0] / 2.;
			Get_Dist_Btw_Edges(d,d->v[i],tree);
		}
	}


}

/*********************************************************/
//JSJ: fixed to detects Polytomies on an array of bls
void Detect_Polytomies(edge *b, m3ldbl l_thresh, arbre *tree)
{
	int i;
	For(i,tree->n_l){
		if((b->l[i] < l_thresh) && (!b->left->tax) && (!b->rght->tax))
		{
			b->l[i]               = 0.0;
			b->has_zero_br_len[i] = 1;
		}
		else b->has_zero_br_len[i] = 0;
	}
}

/*********************************************************/
//JSJ: temporary fix for list of nodes in Polytomy
void Get_List_Of_Nodes_In_Polytomy(node *a, node *d, node ***list, int *size_list)
{
	if(d->tax) return;
	else
	{
		int i;

		For(i,3)
		{
			if(d->v[i] != a)
			{
				if(!d->b[i]->has_zero_br_len[0])
				{
					(*list)[*size_list] = d->v[i];
					(*size_list)++;
				}

				if(d->b[i]->has_zero_br_len[0])
					Get_List_Of_Nodes_In_Polytomy(d,d->v[i],list,size_list);
			}
		}
	}

}


/*********************************************************/

void Check_Path(node *a, node *d, node *target, arbre *tree)
{
	PhyML_Printf("path---------\n");
	if(d==target) return;
	else Check_Path(d,d->v[tree->t_dir[d->num][target->num]],target,tree);
}


/*********************************************************/

void Connect_Two_Nodes(node *a, node *d)
{
	a->v[0] = d;
	d->v[0] = a;
}

/*********************************************************/

void Get_List_Of_Adjacent_Targets(node *a, node *d, node ***node_list, edge ***edge_list, int *list_size)
{
	int i;

	For(i,3)
	if(a->v[i] == d)
	{
		(*node_list)[*list_size] = a;
		(*edge_list)[*list_size] = a->b[i];
		(*list_size)++;
	}
	if(d->tax) return;
	else
		For(i,3)
		if(d->v[i] != a) Get_List_Of_Adjacent_Targets(d,d->v[i],node_list,edge_list,list_size);
}

/*********************************************************/

void Sort_List_Of_Adjacent_Targets(edge ***list, int list_size)
{
	edge *buff_edge;
	int i,j;

	buff_edge = NULL;

	For(i,list_size-1)
	{
		for(j=i+1;j<list_size;j++)
			if((*list)[j]->topo_dist_btw_edges < (*list)[i]->topo_dist_btw_edges)
			{
				buff_edge = (*list)[j];
				(*list)[j] = (*list)[i];
				(*list)[i] = buff_edge;
			}
	}
}

/*********************************************************/


node *Common_Nodes_Btw_Two_Edges(edge *a, edge *b)
{
	if(a->left == b->left)      return b->left;
	else if(a->left == b->rght) return b->rght;
	else if(a->rght == b->left) return b->left;
	else if(a->rght == b->rght) return b->rght;

	PhyML_Printf("\n. First edge = %d (%d %d); Second edge = %d (%d %d)\n",
			a->num,a->left->num,a->rght->num,
			b->num,b->left->num,b->rght->num);
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");

	return NULL;
}

/*********************************************************/

int KH_Test(m3ldbl *site_lk_M1, m3ldbl *site_lk_M2, arbre *tree)
{
	m3ldbl *delta,mean,sd,obs_stat,threshold;
	int i;


	delta = (m3ldbl *)mCalloc(tree->data->init_len,sizeof(m3ldbl));

	threshold = .0;
	mean = .0;
	obs_stat = .0;
	For(i,tree->n_pattern)
	{
		delta[i] = site_lk_M1[i] - site_lk_M2[i];
		mean += ((int)tree->data->wght[i])*delta[i];
	}

	obs_stat = mean;

	mean /= tree->data->init_len;

	For(i,tree->data->init_len) delta[i] -= mean;

	sd = .0;
	For(i,tree->data->init_len) sd += pow(delta[i],2);
	sd /= (m3ldbl)(tree->data->init_len-1.);

	/*   threshold = tree->dnorm_thresh*sqrt(sd*tree->data->init_len); */


	/*   PhyML_Printf("\nObs stat = %f Threshold = %f\n",obs_stat,threshold); */
	Free(delta);

	if(obs_stat > threshold) return 1;
	else                     return 0;
}

/*********************************************************/

/*********************************************************/

void Random_Tree(arbre *tree)
{
	int *is_available,*list_of_nodes;
	int i,node_num,step,n_available;

	PhyML_Printf("\n. Randomising the tree...\n");

	is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
	list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));

	For(i,tree->n_otu) is_available[i]  = 1;
	For(i,tree->n_otu) list_of_nodes[i] = i;

	step = 0;
	do
	{
		/*       node_num = (int)rint(rand()/(m3ldbl)(RAND_MAX+1.0)*(tree->n_otu-1-step)); */
		node_num = Rand_Int(0,tree->n_otu-1-step);
		node_num = list_of_nodes[node_num];
		is_available[node_num] = 0;
		For(i,tree->n_otu) list_of_nodes[i] = -1;
		n_available = 0;
		For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

		tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
		tree->noeud[tree->n_otu+step]->v[1] = tree->noeud[node_num];

		/*       node_num = (int)rint(rand()/(m3ldbl)(RAND_MAX+1.0)*(tree->n_otu-2-step)); */
		node_num = Rand_Int(0,tree->n_otu-2-step);
		node_num = list_of_nodes[node_num];
		is_available[node_num] = 0;
		For(i,tree->n_otu) list_of_nodes[i] = -1;
		n_available = 0;
		For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

		tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
		tree->noeud[tree->n_otu+step]->v[2] = tree->noeud[node_num];

		is_available[tree->n_otu+step] = 1;
		For(i,tree->n_otu) list_of_nodes[i] = -1;
		n_available = 0;
		For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;

		step++;
	}while(step < tree->n_otu-2);

	tree->noeud[list_of_nodes[0]]->v[0] = tree->noeud[list_of_nodes[1]];
	tree->noeud[list_of_nodes[1]]->v[0] = tree->noeud[list_of_nodes[0]];

	tree->num_curr_branch_available = 0;
	Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
	/*   Print_Node(tree->noeud[0],tree->noeud[0]->v[0],tree); */

	Fill_Dir_Table(tree);
	Update_Dirs(tree);

	Free(is_available);
	Free(list_of_nodes);
}

/*********************************************************/


void Random_NNI(int n_moves, arbre *tree)
{
	int i,j;
	edge *b;
	node *n1,*n2,*n_target;

	n1 = n2 = NULL;
	b = NULL;
	For(i,n_moves)
	{
		n_target  = tree->noeud[tree->n_otu + (int)((m3ldbl)rand()/RAND_MAX * (2*tree->n_otu-3-tree->n_otu))];
		For(j,3) if(!n_target->v[j]->tax) {b = n_target->b[j]; break;}


		For(j,3) if(b->left->v[j] != b->rght) {n1 = b->left->v[j]; break;}
		For(j,3) if(b->rght->v[j] != b->left) {n2 = b->rght->v[j]; break;}

		Swap(n1,b->left,b->rght,n2,tree);
	}
}

/*********************************************************/


void Print_Settings(option *io)
{
	int answer,i;
	char *s;

	s = (char *)mCalloc(100,sizeof(char));

	PhyML_Printf("\n\n\n");
	PhyML_Printf("\n\n");

	PhyML_Printf("                                 ..........................                                      \n");
	PhyML_Printf(" ooooooooooooooooooooooooooooo        CURRENT SETTINGS        ooooooooooooooooooooooooooooooooooo\n");
	PhyML_Printf("                                 ..........................                                      \n");

	PhyML_Printf("\n                . Sequence filename : \t\t\t\t %s", Basename(io->in_seq_file));
	PhyML_Printf("\n                . Data type :             \t\t\t %s", (io->mod->datatype ? "aa" : "dna"));
	PhyML_Printf("\n                . Sequence format : \t\t\t\t %s", io->interleaved ? "interleaved" : "sequential");
	PhyML_Printf("\n                . Number of data sets : \t\t\t %d", io->n_data_sets);

	PhyML_Printf("\n                . Nb of bootstrapped data sets : \t\t %d", io->mod->bootstrap);

	if (io->mod->bootstrap > 0)
		PhyML_Printf("\n                . Compute approximate likelihood ratio test : \t no");
	else
	{
		if(io->ratio_test == 1)
			PhyML_Printf("\n                . Compute approximate likelihood ratio test : \t yes (aLRT statistics)");
		else if(io->ratio_test == 2)
			PhyML_Printf("\n                . Compute approximate likelihood ratio test : \t yes (Chi2-based parametric branch supports)");
		else if(io->ratio_test == 3)
			PhyML_Printf("\n                . Compute approximate likelihood ratio test : \t yes (Minimum of SH-like and Chi2-based branch supports)");
		else if(io->ratio_test == 4)
			PhyML_Printf("\n                . Compute approximate likelihood ratio test : \t yes (SH-like branch supports)");
	}

	PhyML_Printf("\n                . Model name : \t\t\t\t\t %s", io->mod->modelname);

	if (io->mod->datatype == NT)
	{
		if ((io->mod->whichmodel == K80)  ||
				(io->mod->whichmodel == HKY85)||
				(io->mod->whichmodel == F84)  ||
				(io->mod->whichmodel == TN93))
		{
			if (io->mod->s_opt->opt_kappa)
				PhyML_Printf("\n                . Ts/tv ratio : \t\t\t\t estimated");
			else
				PhyML_Printf("\n                . Ts/tv ratio : \t\t\t\t %f", io->mod->kappa);
		}
	}

	if (io->mod->s_opt->opt_pinvar)
		PhyML_Printf("\n                . Proportion of invariable sites :\t\t estimated");
	else
		PhyML_Printf("\n                . Proportion of invariable sites :\t\t %f", io->mod->pinvar);


	PhyML_Printf("\n                . Number of subst. rate categs : \t\t %d", io->mod->n_catg);
	if(io->mod->s_opt->opt_alpha)
		PhyML_Printf("\n                . Gamma distribution parameter : \t\t estimated");
	else
		PhyML_Printf("\n                . Gamma distribution parameter : \t\t %f", io->mod->alpha);

	if(io->mod->n_catg > 1)
		PhyML_Printf("\n                . 'Middle' of each rate class  : \t\t %s",(io->mod->gamma_median)?("median"):("mean"));


	if(io->mod->datatype == AA)
		PhyML_Printf("\n                . Amino acid equilibrium frequencies : \t\t %s", (io->mod->s_opt->opt_state_freq) ? ("empirical"):("model"));
	else if(io->mod->datatype == NT)
		if((io->mod->whichmodel != JC69) &&
				(io->mod->whichmodel != K80)  &&
				(io->mod->whichmodel != F81))
			PhyML_Printf("\n                . Nucleotide equilibrium frequencies : \t\t %s", (io->mod->s_opt->opt_state_freq) ? ("ML"):("empirical"));


	PhyML_Printf("\n                . Optimise tree topology : \t\t\t %s", (io->mod->s_opt->opt_topo) ? "yes" : "no");

	switch(io->in_tree)
	{
	case 0: { strcpy(s,"BioNJ");     break; }
	case 1: { strcpy(s,"parsimony"); break; }
	case 2: { strcpy(s,"user tree (");
	strcat(s,Basename(io->in_tree_file));
	strcat(s,")");         break; }
	}

	if(io->mod->s_opt->opt_topo)
	{
		if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n                . Tree topology search : \t\t\t NNIs");
		else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n                . Tree topology search : \t\t\t SPRs");
		else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Printf("\n                . Tree topology search : \t\t\t Best of NNIs and SPRs");



		PhyML_Printf("\n                . Starting tree : \t\t\t\t %s",s);

		PhyML_Printf("\n                . Add random input tree : \t\t\t %s", (io->mod->s_opt->random_input_tree) ? "yes" : "no");
		if(io->mod->s_opt->random_input_tree)
			PhyML_Printf("\n                . Number of random starting trees : \t\t %d", io->mod->s_opt->n_rand_starts);
	}
	else
		if(!io->mod->s_opt->random_input_tree)
			PhyML_Printf("\n                . Evaluted tree : \t\t\t\t file \"%s\"",s);

	PhyML_Printf("\n                . Optimise branch lengths : \t\t\t %s", (io->mod->s_opt->opt_bl) ? "yes" : "no");

	answer = 0;
	if(io->mod->s_opt->opt_alpha  ||
			io->mod->s_opt->opt_kappa  ||
			io->mod->s_opt->opt_lambda ||
			io->mod->s_opt->opt_pinvar ||
			io->mod->s_opt->opt_rr) answer = 1;

	PhyML_Printf("\n                . Optimise substitution model parameters : \t %s", (answer) ? "yes" : "no");

	PhyML_Printf("\n                . Number of branch length categories: \t %i", io->n_l);
	PhyML_Printf("\n                . Proportion of sites in each branch length category: [");
	For(i,io->n_l){
		if(i+1 == io->n_l){
			PhyML_Printf(" %lf ]",(double)io->props[i]);
		}else{
			PhyML_Printf(" %lf,",(double)io->props[i]);
		}
	}
	PhyML_Printf("\n                . Optimize proportion of sites in each branch length category: \t %s", (io->fixed_props == 0) ? ("Yes"):("No"));

	PhyML_Printf("\n                . Run ID : \t\t\t\t\t %s", (io->append_run_ID) ? (io->run_id_string) : ("none"));


	PhyML_Printf("\n\n oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

	PhyML_Printf("\n\n");
	fflush(NULL);

	Free(s);
}

/*********************************************************/

void Fill_Missing_Dist(matrix *mat)
{
	int i,j;
	For(i,mat->n_otu)
	{
		for(j=i+1;j<mat->n_otu;j++)
		{
			if(i != j)
			{
				if(mat->dist[i][j] < .0)
				{
					Fill_Missing_Dist_XY(i,j,mat);
					mat->dist[j][i] = mat->dist[i][j];
				}
			}
		}
	}
}

/*********************************************************/

void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{

	int i,j;
	m3ldbl *local_mins,**S1S2;
	int cpt;
	int pos_best_estimate;
	m3ldbl min_crit, curr_crit;

	local_mins = (m3ldbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(m3ldbl ));
	S1S2       = (m3ldbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(m3ldbl *));
	For(i,mat->n_otu*mat->n_otu) S1S2[i] = (m3ldbl *)mCalloc(2,sizeof(m3ldbl));

	cpt = 0;
	For(i,mat->n_otu)
	{
		if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
		{
			For(j,mat->n_otu)
			{
				if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
				{
					if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
					{
						S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
						S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
						cpt++;
					}
				}
			}
		}
	}

	Qksort_Matrix(S1S2,0,0,cpt-1);

	local_mins[0] = S1S2[0][1];
	for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(m3ldbl)(i+1);

	pos_best_estimate = 0;
	min_crit = curr_crit = MDBL_MAX;

	For(i,cpt-1)
	{
		if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
		{
			curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
			if(curr_crit < min_crit)
			{
				min_crit = curr_crit;
				pos_best_estimate = i;
			}
		}
	}

	mat->dist[x][y] = local_mins[pos_best_estimate];
	mat->dist[y][x] = mat->dist[x][y];

	For(i,mat->n_otu*mat->n_otu) Free(S1S2[i]);
	Free(S1S2);
	Free(local_mins);
}

/*********************************************************/

m3ldbl Least_Square_Missing_Dist_XY(int x, int y, m3ldbl dxy, matrix *mat)
{
	int i,j;
	m3ldbl fit;

	fit = .0;
	For(i,mat->n_otu)
	{
		if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
		{
			For(j,mat->n_otu)
			{
				if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
				{
					if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
					{
						if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
						{
							fit += pow((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
						}
						else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
						{
							fit += pow(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
						}
						else
						{
							fit += pow(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
						}
					}
				}
			}
		}
	}
	return fit;
}

/*********************************************************/

void Print_Banner(FILE *fp)
{
	PhyML_Fprintf(fp,"\n");
	PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp,"                                   ---  PhyML %s  ---                                             \n",VERSION);
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp,"    A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood    \n");
	PhyML_Fprintf(fp,"                            Stephane Guindon & Olivier Gascuel                                      \n");
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp,"                           http://www.atgc-montpellier.fr/phyml                                          \n");
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/

void Print_Banner_Small(FILE *fp)
{
	PhyML_Fprintf(fp,"\n");
	PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	PhyML_Fprintf(fp,"                                    ---  PhyML %s  ---                                             \n",VERSION);
	PhyML_Fprintf(fp,"                            http://www.atgc-montpellier.fr/phyml                                          \n");
	PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
	PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/


void Print_Data_Set_Number(option *io, FILE *fp)
{
	PhyML_Fprintf(fp,"\n");
	PhyML_Fprintf(fp,"                                                                                                  \n");
	PhyML_Fprintf(fp,"                                 [ Data set number %3d ]                                           \n",io->curr_gt+1);
	PhyML_Fprintf(fp,"                                                                                                  \n");
}
/*********************************************************/


void Check_Memory_Amount(arbre *tree)
{
	/* Rough estimate of the amount of memory that has to be used */

	int nbytes;
	model *mod;

	mod = tree->mod;

	nbytes = 0;


	/* Pmat JSJ: an array of them...*/
	nbytes += (2*mod->n_otu-3) * tree->n_l * mod->n_catg * mod->ns * mod->ns * sizeof(m3ldbl);

	/* Partial Lk */
	nbytes += ((2*mod->n_otu-3) * 2 - tree->n_otu) * tree->n_pattern * mod->n_catg * mod->ns * sizeof(m3ldbl);

	/* Scaling factors */
	nbytes += ((2*mod->n_otu-3) * 2 - tree->n_otu) * tree->n_pattern * sizeof(m3ldbl);

	/* Partial Pars */
	nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->ns * sizeof(short int);
	nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->ns * sizeof(int);
	nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(int);
	nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(int);

	if(((m3ldbl)nbytes/(1.E+06)) > 256.)
	{
		char answer;
		PhyML_Printf("\n. WARNING: this analysis requires at least %.0fMo of memory space.\n",(m3ldbl)nbytes/(1.E+06));
#ifndef BATCH
		if (! tree->io->quiet) {
			PhyML_Printf("\n. Do you really want to continue ? [Y/n] ");
			if(scanf("%c", &answer))
			{
				if(answer == '\n') answer = 'Y';
				else if(answer == 'n' || answer == 'N') Warn_And_Exit("\n\n");
				else getchar();
			}
			else
			{
				Warn_And_Exit("\n\n");
			}
		}
#endif
	}
	else if(((m3ldbl)nbytes/(1.E+06)) > 100.)
	{
		if(!tree->io->quiet) PhyML_Printf("\n. WARNING: this analysis will use at least %.0fMo of memory space...\n",(m3ldbl)nbytes/(1.E+06));
	}
	else
	{
		if(!tree->io->quiet) PhyML_Printf("\n. This analysis requires at least %.0fMo of memory space.\n",(m3ldbl)nbytes/(1.E+06));
	}
}

/*********************************************************/

int Get_State_From_P_Lk(m3ldbl *p_lk, int pos, arbre *tree)
{
	int i;
	For(i,tree->mod->ns) if(p_lk[pos+i] > .0) return i;
	return -1;
}

/*********************************************************/

int Get_State_From_P_Pars(short int *p_pars, int pos, arbre *tree)
{
	int i;
	For(i,tree->mod->ns) if(p_pars[pos+i] > .0) return i;
	return -1;
}

/*********************************************************/

void Print_Lk(arbre *tree, char *string)
{
	time(&(tree->t_current));
	PhyML_Printf("\n. (%5d sec) [%15.4f] %s",(int)(tree->t_current-tree->t_beg),tree->c_lnL,string);
#ifndef QUIET
	fflush(NULL);
#endif
}

/*********************************************************/

void Print_Pars(arbre *tree)
{
	time(&(tree->t_current));
	PhyML_Printf("\n. (%5d sec) [%5d]",(int)(tree->t_current-tree->t_beg),tree->c_pars);
#ifndef QUIET
	fflush(NULL);
#endif
}

/*********************************************************/

void Print_Lk_And_Pars(arbre *tree)
{
	time(&(tree->t_current));

	PhyML_Printf("\n. (%5d sec) [%15.4f] [%5d]",
			(int)(tree->t_current-tree->t_beg),
			tree->c_lnL,tree->c_pars);
#ifndef QUIET
	fflush(NULL);
#endif
}

/*********************************************************/

void Check_Dirs(arbre *tree)
{
	int i;

	For(i,2*tree->n_otu-3)
	{
		if(!tree->t_edges[i]->left->tax)
		{
			if(tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num <
					tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num)
			{
				PhyML_Printf("\n. Edge %d ; v1=%d v2=%d",
						tree->t_edges[i]->num,
						tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num,
						tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num);
				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
				Warn_And_Exit("");
			}
		}

		if(!tree->t_edges[i]->rght->tax)
		{
			if(tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num <
					tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num)
			{
				PhyML_Printf("\n. Edge %d ; v3=%d v4=%d",
						tree->t_edges[i]->num,
						tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num,
						tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num);
				PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
				Warn_And_Exit("");
			}
		}
	}
}

/*********************************************************/

void Warn_And_Exit(char *s)
{
	PhyML_Fprintf(stdout,"%s",s);
	fflush(NULL);
#ifndef BATCH
	//  if (! tree->io->quiet) {
	char c;
	PhyML_Fprintf(stdout,"\n. Type any key to exit.\n");
	if(!fscanf(stdin,"%c",&c)) Exit("");
	//  }
#endif
	Exit("\n");
}

/*********************************************************/

void Read_Qmat(double *daa, m3ldbl *pi, FILE *fp)
{
	int i,j;
	m3ldbl sum;

	for(i=1;i<20;i++)
	{
		For(j,19)
		{
			if(!fscanf(fp,"%lf",&(daa[i*20+j]))) Exit("\n");
			daa[j*20+i] = daa[i*20+j];
			if(j == i-1) break;
		}
	}

	For(i,20) { if(!fscanf(fp,"%lf",pi+i)) Exit("\n");}
	sum = .0;
	For(i,20) sum += pi[i];
	if(fabs(sum - 1.) > 1.E-06)
	{
		PhyML_Printf("\n. Scaling amino-acid frequencies...\n");
		For(i,20) pi[i] /= sum;
	}
}

/*********************************************************/

void Print_Qmat_AA(double *daa, m3ldbl *pi)
{
	int i,j,cpt;

	cpt = 0;
	For(i,20)
	{
		for(j=0;j<i;j++)
		{
			PhyML_Printf("daa[%2d*20+%2d] = %10f;  ",i,j,daa[i*20+j]);
			cpt++;
			if(!(cpt%4)) PhyML_Printf("\n");
		}
	}

	PhyML_Printf("\n\n");
	PhyML_Printf("for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];\n\n");
	For(i,20) PhyML_Printf("pi[%d] = %f; ",i,pi[i]);
	PhyML_Printf("\n");
	PhyML_Printf("Ala\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet\tPhe\tPro\tSer\tThr\tTrp\tTyr\tVal\n");
}


/*********************************************************/

void Randomize_Sequence_Order(allseq *data)
{
	int i,exchange_with;
	m3ldbl buff_dbl;
	char *buff_name,*buff_seq;
	short int *buff_ambigu;

	exchange_with = -1;
	For(i,data->n_otu)
	{
		buff_dbl  = rand();
		buff_dbl /= (RAND_MAX+1.);
		buff_dbl *= data->n_otu;
		exchange_with = (int)floor(buff_dbl);

		buff_name = data->c_seq[i]->name;
		data->c_seq[i]->name = data->c_seq[exchange_with]->name;
		data->c_seq[exchange_with]->name = buff_name;

		buff_seq = data->c_seq[i]->state;
		data->c_seq[i]->state = data->c_seq[exchange_with]->state;
		data->c_seq[exchange_with]->state = buff_seq;

		buff_ambigu = data->c_seq[i]->is_ambigu;
		data->c_seq[i]->is_ambigu = data->c_seq[exchange_with]->is_ambigu;
		data->c_seq[exchange_with]->is_ambigu = buff_ambigu;
	}
}

/*********************************************************/
//JSJ: the first one was n_root->l[0] or n_root->[1] before I modified
void Update_Root_Pos(arbre *tree)
{
	int i;
	For(i, tree->n_l){
		if(tree->n_root_pos > -1.0)
		{
			tree->n_root->l[i][0] = tree->e_root->l[i] * tree->n_root_pos;
			tree->n_root->l[i][1] = tree->e_root->l[i] * (1.-tree->n_root_pos);
		}
		else
		{
			tree->n_root->l[i][0] = tree->e_root->l[i] / 2.;
			tree->n_root->l[i][1] = tree->e_root->l[i] / 2.;
		}
	}
}

/*********************************************************/

void Add_Root(edge *target, arbre *tree)
{
	int i;
	PhyML_Printf("\n. Add root on edge %d left = %d right = %d",target->num,target->left->num,target->rght->num); fflush(NULL);
	tree->e_root = target;

	/* Create the root node if it does not exist yet */
	if((!tree->n_root) || (tree->n_root->num != 2*tree->n_otu-2))
	{
		tree->n_root = (node *)Make_Node_Light(2*tree->n_otu-2, tree->n_l);
	}

	tree->n_root->tax = 0;

	/* Set the position of the root */
	tree->n_root->v[0] = tree->e_root->left;
	tree->n_root->v[1] = tree->e_root->rght;

	tree->n_root->b[0] = tree->e_root;
	tree->n_root->b[1] = tree->e_root;
	For(i,tree->n_l){ //JSJ: iterate over set and do below for each member of set
		if(tree->n_root_pos > -1.0)
		{
			if(tree->n_root_pos < 1.E-6 &&  tree->n_root_pos > -1.E-6)
				printf("\n. WARNING: you put the root at a weird position...");

			/*       tree->n_root->l[0] = tree->e_root->l * (tree->n_root_pos/(1.+tree->n_root_pos)); */
			/*       tree->n_root->l[1] = tree->e_root->l - tree->n_root->l[0]; */
			//JSJ: simmilar compilation fix
			tree->n_root->l[i][0] = tree->e_root->l[i] * tree->n_root_pos;
			tree->n_root->l[i][1] = tree->e_root->l[i] * (1. - tree->n_root_pos);
		}
		else
		{
			tree->n_root->l[i][0] = tree->e_root->l[i] / 2.;
			tree->n_root->l[i][1] = tree->e_root->l[i] / 2.;
			tree->n_root_pos = 0.5;
		}
	}

	Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
	Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
	tree->n_root->anc = NULL;
}

/*********************************************************/

void Update_Ancestors(node *a, node *d, arbre *tree)
{
	d->anc = a;
	if(d->tax) return;
	else
	{
		int i;
		For(i,3)
		if((d->v[i] != d->anc) && (d->b[i] != tree->e_root))
			Update_Ancestors(d,d->v[i],tree);
	}
}

/*********************************************************/
/* Generate a random unrooted tree with 'n_otu' OTUs */
#ifdef MC
arbre *Generate_Random_Tree_From_Scratch(int n_otu, int rooted)
{
	arbre *tree;
	int *connected,*nonconnected,*available_nodes;
	int i,n_connected,n_nonconnected,n1,n2,new_n,n_internal,n_external,n_available;
	node *root,*curr_n,**internal_nodes, **external_nodes;
	m3ldbl *t,*tmp;

	tree = (arbre *)Make_Tree(n_otu);
	Init_Tree(tree,tree->n_otu);
	Make_All_Tree_Nodes(tree);
	Make_All_Tree_Edges(tree);
	Make_Tree_Path(tree);
	Make_List_Of_Reachable_Tips(tree);
	tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
	RATES_Init_Rate_Struct(tree->rates,tree->n_otu);

	For(i,2*tree->n_otu-2)
	{
		tree->noeud[i]->v[1] = NULL;
		tree->noeud[i]->v[2] = NULL;
	}

	root = (node *)Make_Node_Light(2*tree->n_otu-2);

	connected       = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
	nonconnected    = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
	available_nodes = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
	internal_nodes  = (node **)mCalloc(tree->n_otu-2,sizeof(node *));
	external_nodes  = (node **)mCalloc(tree->n_otu,  sizeof(node *));
	t               = (m3ldbl *)mCalloc(tree->n_otu-1,sizeof(m3ldbl ));
	tmp             = (m3ldbl *)mCalloc(2*tree->n_otu-2,sizeof(m3ldbl ));

	n_nonconnected = 2*n_otu-2;

	For(i,2*tree->n_otu-2) nonconnected[i] = i;

	available_nodes[0] = 2*n_otu-2;



	/* Node times are generated according to a Birth-death process.
     Formulae are as described by Yang and Rannala (1997) */
	m3ldbl    phi;
	m3ldbl    rho; /* sampling intensity */
	m3ldbl     mu; /* birth rate */
	m3ldbl lambda; /* death rate */
	m3ldbl      u; /* random U[0,1] */
	m3ldbl expval;

	/* rho = 1.0 and mu = 0.0 correspond to the Yule process */

	lambda = 6.7;
	mu     = 2.5;
	rho    = 9./150.;

	expval = exp(mu-lambda);
	phi = (rho*lambda*(expval-1.) + (mu-lambda)*expval)/(expval-1.); /* Equation 16 */

	For(i,tree->n_otu-1)
	{
		u = rand();
		u /= RAND_MAX;

		if(fabs(lambda - mu) > 1.E-4)
			t[i] = (log(phi-u*rho*lambda) - log(phi-u*rho*lambda + u*(lambda-mu)))/(mu-lambda); /* Equation 15 */
		else
			t[i] = u / (1.+lambda*rho*(1-u)); /* Equation 17 */
	}

	Qksort(t,NULL,0,tree->n_otu-2); /* Node times ordering in ascending order */

	For(i,tree->n_otu-1) tmp[i] =  t[tree->n_otu-2-i];
	For(i,tree->n_otu-1) t[i]   = -tmp[i];

	/* Rescale node times such that the time at the root node is -100 */
	for(i=1;i<tree->n_otu-1;i++)
	{
		t[i] /= -t[0];
		t[i] *= 1.E+02;
	}
	t[0] = -1.E+02;

	n_available = 1;
	curr_n = root;
	n_connected = 0;
	do
	{
		n1 = Rand_Int(0,n_nonconnected-1);
		n1 = nonconnected[n1];
		connected[n1] = 1;

		n_nonconnected = 0;
		For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}

		n2 = Rand_Int(0,n_nonconnected-1);
		n2 = nonconnected[n2];
		connected[n2] = 1;

		n_nonconnected = 0;
		For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}

		curr_n->v[1] = tree->noeud[n1];
		curr_n->v[2] = tree->noeud[n2];
		tree->noeud[n1]->v[0] = curr_n;
		tree->noeud[n2]->v[0] = curr_n;

		tree->rates->nd_t[curr_n->num] = t[n_connected/2];

		available_nodes[n_available] = tree->noeud[n1]->num;
		For(i,n_available)
		if(available_nodes[i] == curr_n->num)
		{
			available_nodes[i] = tree->noeud[n2]->num;
			break;
		}
		n_available++;

		new_n = Rand_Int(0,n_available-1);
		curr_n = tree->noeud[available_nodes[new_n]];

		n_connected+=2;

	}while(n_connected < 2*tree->n_otu-2);

	For(i,2*tree->n_otu-2) tmp[i] = tree->rates->nd_t[i];

	/* Unroot the tree */
	root->v[1]->v[0] = root->v[2];
	root->v[2]->v[0] = root->v[1];


	n_internal = n_external = 0;
	For(i,2*tree->n_otu-2)
	{
		if(tree->noeud[i]->v[1]) internal_nodes[n_internal++] = tree->noeud[i];
		else                     external_nodes[n_external++] = tree->noeud[i];
	}


	n_internal = n_external = 0;
	For(i,2*tree->n_otu-2)
	{
		if(i < tree->n_otu)
		{
			tree->noeud[i]      = external_nodes[n_external++];
			tree->noeud[i]->tax = 1;
		}
		else
		{
			tree->rates->nd_t[i] = tmp[internal_nodes[n_internal]->num];
			tree->noeud[i]        = internal_nodes[n_internal++];
			tree->noeud[i]->tax   = 0;
		}

		tree->noeud[i]->num = i;
	}

	For(i,tree->n_otu) tree->rates->nd_t[i] = 0.0;

	For(i,tree->n_otu)
	{
		strcpy(tree->noeud[i]->name,"x");
		sprintf(tree->noeud[i]->name+1,"%d",i);
	}


	tree->num_curr_branch_available = 0;
	Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
	Fill_Dir_Table(tree);
	Update_Dirs(tree);

	/* Add root */
	if(rooted)
	{
		For(i,2*tree->n_otu-3)
		{
			if(((tree->t_edges[i]->left == root->v[1]) || (tree->t_edges[i]->rght == root->v[1])) &&
					((tree->t_edges[i]->left == root->v[2]) || (tree->t_edges[i]->rght == root->v[2])))
			{
				Add_Root(tree->t_edges[i],tree);
				break;
			}
		}
	}
	/* Or not... */
	else
	{
		Free_Node(root);
	}

	tree->rates->use_rates = 0;
	MC_Bl_From_T(tree);

	if(rooted)
	{
		tree->n_root->l[0] = (tree->rates->nd_t[tree->n_root->v[0]->num] - t[0]);
		tree->n_root->l[1] = (tree->rates->nd_t[tree->n_root->v[1]->num] - t[0]);
		tree->e_root->l    = tree->n_root->l[0] + tree->n_root->l[1];
	}

	RATES_Random_Branch_Lengths(tree);

	Free(available_nodes);
	Free(connected);
	Free(nonconnected);
	Free(external_nodes);
	Free(internal_nodes);
	Free(t);
	Free(tmp);


	return tree;
}
#endif
/*********************************************************/

void Random_Lineage_Rates(node *a, node *d, edge *b, m3ldbl stick_prob, m3ldbl *rates, int curr_rate, int n_rates, arbre *tree)
{
	m3ldbl uni;
	int new_rate;
	int i,k;


	if(b)
	{
		uni  = rand();
		uni /= RAND_MAX;

		if(uni > stick_prob) /* Randomly pick a new rate */
		{
			uni  = rand();
			uni /= RAND_MAX;
			uni = (m3ldbl)(uni * (n_rates-1));
			if(uni-(int)(uni) > 0.5-MDBL_MAX) new_rate = (int)(uni)+1;
			else new_rate = (int)(uni);
		}
		else
		{
			new_rate = curr_rate;
		}
		For(k,tree->n_l){
			For(i,3)
			if(a->v[i] == d)
			{//JSJ: compilation fix
				a->b[i]->l[k] *= rates[new_rate];
				break;
			}
		}

		For(i,3)
		if(a->v[i] == d)
		{
			if(!(a->b[i]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(a->b[i]);
			if(rates[new_rate] > 1.0)      strcpy(a->b[i]->labels[a->b[i]->n_labels],"FAST");
			else if(rates[new_rate] < 1.0) strcpy(a->b[i]->labels[a->b[i]->n_labels],"SLOW");
			else                           strcpy(a->b[i]->labels[a->b[i]->n_labels],"MEDIUM");
			a->b[i]->n_labels++;
			break;
		}
		curr_rate = new_rate;
	}

	if(d->tax) return;
	else
	{
		For(i,3)
		if((d->v[i] != a) && (d->b[i] != tree->e_root))
			Random_Lineage_Rates(d,d->v[i],d->b[i],stick_prob,rates,curr_rate,n_rates,tree);
	}
}

/*********************************************************/

edge *Find_Edge_With_Label(char *label, arbre *tree)
{
	int i,j;

	For(i,2*tree->n_otu-3)
	{
		For(j,tree->t_edges[i]->n_labels)
		{
			if(!strcmp(tree->t_edges[i]->labels[j],label)) return tree->t_edges[i];
		}
	}
	return NULL;
}

/*********************************************************/

void Print_Square_Matrix_Generic(int n, m3ldbl *mat)
{
	int i,j;

	PhyML_Printf("\n");
	For(i,n)
	{
		For(j,n)
		{
			PhyML_Printf("%.3f ",mat[i*n+j]);
		}
		PhyML_Printf("\n");
	}
	PhyML_Printf("\n");
}

/*********************************************************/

void Evolve(allseq *data, model *mod, arbre *tree)
{
	int root_state, root_rate_class;
	int site,i;

	if(mod->use_m4mod) tree->print_labels = 1;

	/* Get the change probability matrices */
	Set_Model_Parameters(mod);
	For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);


	For(site,data->init_len)
	{
		root_state = root_rate_class = -1;

		/* Pick the root nucleotide/aa */
		root_state = Pick_State(mod->ns,mod->pi);
		data->c_seq[0]->state[site] = Reciproc_Assign_State(root_state,mod->datatype);

		/* Pick the rate class */
		root_rate_class = Pick_State(mod->n_catg,mod->gamma_r_proba);

		/* tree->noeud[0] is considered as the root node */
		Evolve_Recur(tree->noeud[0],
				tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],
				root_state,
				root_rate_class,
				site,
				data,
				mod,
				tree);

		/*       PhyML_Printf("%s\n",Write_Tree(tree)); */

		data->wght[site] = 1;
	}
	data->crunch_len = data->init_len;
}

/*********************************************************/

int Pick_State(int n, m3ldbl *prob)
{
	int pos;
	m3ldbl uni;

	do
	{
		pos  = rand();
		pos  = (pos % n);
		uni  = (m3ldbl)rand();
		uni /= (m3ldbl)RAND_MAX;
		if(uni < prob[pos]) break;
	}
	while(1);

	return (int)pos;
}

/*********************************************************/

void Evolve_Recur(node *a, node *d, edge *b, int a_state, int r_class, int site_num, allseq *gen_data, model *mod, arbre *tree)
{
	int d_state;
	int dim1,dim2;

	dim1 = tree->mod->ns * tree->mod->ns;
	dim2 = tree->mod->ns;
	//JSJ: temp fix of Pij_rr
	d_state = Pick_State(mod->ns,b->Pij_rr[0]+r_class*dim1+a_state*dim2);

	/*   PhyML_Printf("\n>> %c (%d,%d)",Reciproc_Assign_State(d_state,mod->datatype),d_state,(int)d_state/mod->m4mod->n_o); */

	if(mod->use_m4mod)
	{
		m3ldbl rrate; /* relative rate of substitutions */

		rrate = mod->m4mod->multipl[(int)d_state/mod->m4mod->n_o];
		if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
		if(rrate > 1.0) strcpy(b->labels[b->n_labels],"FASTER");
		else strcpy(b->labels[b->n_labels],"SLOWER");
		b->n_labels++;
	}

	if(d->tax)
	{
		gen_data->c_seq[d->num]->state[site_num] = Reciproc_Assign_State(d_state,mod->datatype);
		return;
	}
	else
	{
		int i;
		For(i,3)
		if(d->v[i] != a)
			Evolve_Recur(d,d->v[i],d->b[i],
					d_state,r_class,site_num,gen_data,
					mod,tree);
	}
}

/*********************************************************/

void Site_Diversity(arbre *tree)
{
	int i,j,k,ns;
	int *div,sum;

	ns = (tree->mod->datatype == NT)?(4):(20);

	div = (int *)mCalloc(ns,sizeof(int));

	Site_Diversity_Post(tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);
	Site_Diversity_Pre (tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);

	For(i,2*tree->n_otu-3)
	{
		For(j,ns)
		{
			tree->t_edges[i]->div_post_pred_left[j] = 0;
			tree->t_edges[i]->div_post_pred_rght[j] = 0;
		}
	}

	For(i,tree->n_pattern)
	{
		For(j,2*tree->n_otu-3)
		{
			Binary_Decomposition(tree->t_edges[j]->ui_l[i],div,ns);
			sum = 0;
			For(k,ns) sum += div[k];
			tree->t_edges[j]->div_post_pred_left[sum-1] += tree->data->wght[i];

			Binary_Decomposition(tree->t_edges[j]->ui_r[i],div,ns);
			sum = 0;
			For(k,ns) sum += div[k];
			tree->t_edges[j]->div_post_pred_rght[sum-1] += tree->data->wght[i];
		}
	}

	/*   For(j,2*tree->n_otu-3) */
	/*     { */
	/*       PhyML_Printf("\n. Edge %4d   div_left = %4d %4d %4d %4d -- div_rght = %4d %4d %4d %4d", */
	/* 	     j, */
	/* 	     tree->t_edges[j]->div_post_pred_left[0], */
	/* 	     tree->t_edges[j]->div_post_pred_left[1], */
	/* 	     tree->t_edges[j]->div_post_pred_left[2], */
	/* 	     tree->t_edges[j]->div_post_pred_left[3], */
	/* 	     tree->t_edges[j]->div_post_pred_rght[0], */
	/* 	     tree->t_edges[j]->div_post_pred_rght[1], */
	/* 	     tree->t_edges[j]->div_post_pred_rght[2], */
	/* 	     tree->t_edges[j]->div_post_pred_rght[3]); */
	/*     } */

	Free(div);
}

/*********************************************************/

void Site_Diversity_Post(node *a, node *d, edge *b, arbre *tree)
{
	if(d->tax) return;
	else
	{
		int i;

		For(i,3)
		if(d->v[i] != a)
			Site_Diversity_Post(d,d->v[i],d->b[i],tree);

		Subtree_Union(d,b,tree);
	}
}

/*********************************************************/

void Site_Diversity_Pre(node *a, node *d, edge *b, arbre *tree)
{
	if(d->tax) return;
	else
	{
		int i;

		For(i,3)
		if(d->v[i] != a)
		{
			Subtree_Union(d,d->b[i],tree);
			Site_Diversity_Pre(d,d->v[i],d->b[i],tree);
		}
	}
}

/*********************************************************/

void Subtree_Union(node *n, edge *b_fcus, arbre *tree)
{
	/*
           |
	   |<- b_cus
	   |
	   n
          / \
       	 /   \
       	/     \
	 */

	int site;
	unsigned int *ui, *ui_v1, *ui_v2;

	ui = ui_v1 = ui_v2 = NULL;

	if(n == b_fcus->left)
	{
		ui = b_fcus->ui_l;

		ui_v1 =
				(n == n->b[b_fcus->l_v1]->left)?
						(n->b[b_fcus->l_v1]->ui_r):
							(n->b[b_fcus->l_v1]->ui_l);

						ui_v2 =
								(n == n->b[b_fcus->l_v2]->left)?
										(n->b[b_fcus->l_v2]->ui_r):
											(n->b[b_fcus->l_v2]->ui_l);
	}
	else
	{
		ui = b_fcus->ui_r;

		ui_v1 =
				(n == n->b[b_fcus->r_v1]->left)?
						(n->b[b_fcus->r_v1]->ui_r):
							(n->b[b_fcus->r_v1]->ui_l);

						ui_v2 =
								(n == n->b[b_fcus->r_v2]->left)?
										(n->b[b_fcus->r_v2]->ui_r):
											(n->b[b_fcus->r_v2]->ui_l);
	}

	For(site,tree->n_pattern) ui[site] = ui_v1[site] | ui_v2[site];

}

/*********************************************************/

void Binary_Decomposition(int value, int *bit_vect, int size)
{
	int i,cumul;

	For(i,size) bit_vect[i] = 0;

	cumul = 0;
	for(i=size-1;i>=0;i--)
	{
		if(value - cumul < (int)pow(2,i))
		{
			bit_vect[i] = 0;
		}
		else
		{
			bit_vect[i] = 1;
			cumul += (int)pow(2,i);
		}
	}
}

/*********************************************************/

void Print_Diversity_Header(FILE *fp, arbre *tree)
{
	/*   PhyML_Fprintf(fp,"edge side mean\n");  */
	PhyML_Fprintf(fp,"edge side diversity count\n");
}

/*********************************************************/

void Print_Diversity(FILE *fp, arbre *tree)
{
	int ns;

	ns = (tree->mod->datatype == NT)?(4):(20);

	Print_Diversity_Pre(tree->noeud[0],
			tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],
			fp,
			tree);

	/*       mean_div_left = .0; */
	/*       For(k,ns)  */
	/* 	{ */
	/* 	  mean_div_left += (k+1) * tree->t_edges[j]->div_post_pred_left[k]; */
	/* 	} */
	/*       mean_div_rght = .0; */
	/*       For(k,ns) mean_div_rght += (k+1) * tree->t_edges[j]->div_post_pred_rght[k]; */

	/*       mean_div_left /= (m3ldbl)tree->data->init_len; */
	/*       mean_div_rght /= (m3ldbl)tree->data->init_len; */

	/*       PhyML_Fprintf(fp,"%4d 0 %f\n",j,mean_div_left); */
	/*       PhyML_Fprintf(fp,"%4d 1 %f\n",j,mean_div_rght); */


	/*       mean_div_left = .0; */
	/*       For(k,ns) mean_div_left += tree->t_edges[j]->div_post_pred_left[k]; */

	/*       mean_div_rght = .0; */
	/*       For(k,ns)  */
	/* 	{ */
	/* 	  mean_div_rght += tree->t_edges[j]->div_post_pred_rght[k]; */
	/* 	} */

	/*       if((mean_div_left != tree->data->init_len) || (mean_div_rght != tree->data->init_len)) */
	/* 	{ */
	/* 	  PhyML_Printf("\n. mean_div_left = %f mean_div_rght = %f init_len = %d", */
	/* 		 mean_div_left,mean_div_rght,tree->data->init_len); */
	/* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
	/* 	  Warn_And_Exit(""); */
	/* 	} */
}

/*********************************************************/

void Print_Diversity_Pre(node *a, node *d, edge *b, FILE *fp, arbre *tree)
{
	int k,ns;


	if(d->tax) return;
	else
	{
		ns = (tree->mod->datatype == NT)?(4):(20);
		if(d == b->left) For(k,ns) PhyML_Fprintf(fp,"%4d 0 %2d %4d\n",b->num,k,b->div_post_pred_left[k]);
		else             For(k,ns) PhyML_Fprintf(fp,"%4d 1 %2d %4d\n",b->num,k,b->div_post_pred_rght[k]);

		For(k,3) if(d->v[k] != a) Print_Diversity_Pre(d,d->v[k],d->b[k],fp,tree);
	}

}

/*********************************************************/
/* Estimation of density using kernel smoothing.
- where : point where I want to estimate the density,
- x : data vector,
- sample_size :  number of data points in x
 */
m3ldbl Univariate_Kernel_Density_Estimate(m3ldbl where, m3ldbl *x, int sample_size)
{
	m3ldbl sd,h;
	m3ldbl density,sqrt2pi,cons;
	int i;

	sqrt2pi = 2.506628;

	sd = sqrt(Var(x,sample_size));
	h = 1.06 * sd * pow(sample_size,-1./5.); /* Quick and dirty way to set the bandwidth */

	cons = (1./sample_size) * (1./h) * (1./sqrt2pi);

	density = .0;
	For(i,sample_size) density += exp(-0.5 * pow((x[i] - where)/h,2));
	density *= cons;

	return density;
}

/*********************************************************/

/* Estimation of a multivariate density using kernel smoothing.

- where : vector where I want to estimate the density,
- x : data matrix, i.e., sample of vectors,
- sample_size : number of vectors,
- vect_size : vector length.

See "Multivariate Density Estimation" by David Scott. pp 150.
 */
m3ldbl Multivariate_Kernel_Density_Estimate(m3ldbl *where, m3ldbl **x, int sample_size, int vect_size)
{
	m3ldbl sd,*h,cons,density,tmp;
	m3ldbl _2pi;
	int i,j;

	h = (m3ldbl *)mCalloc(vect_size,sizeof(m3ldbl));

	_2pi = 6.283185;

	For(i,vect_size)
	{
		sd = sqrt(Var(x[i],sample_size));
		/*       h[i] = pow(4./(vect_size+2.),1./(vect_size+4)) * sd * pow(sample_size,-1./(vect_size+4)); */
		h[i] = sd * pow(sample_size,-1./(vect_size+4));
		/*       PhyML_Printf("\n. sd = %f, h[i] = %f",sd,h[i]); */
	}

	cons = sample_size;
	For(i,vect_size) cons *= h[i];
	cons *= pow(_2pi,vect_size/2.);
	cons = 1./cons;

	density = .0;
	For(i,sample_size)
	{
		tmp = 1.0;
		For(j,vect_size)
		{
			tmp *= exp(-0.5 * pow((x[j][i] - where[j])/h[j],2));
		}
		density += tmp;
	}

	density *= cons;

	Free(h);

	return density;
}

/*********************************************************/

m3ldbl Var(m3ldbl *x, int n)
{
	m3ldbl mean, sum2;
	int i;

	mean = Mean(x,n);

	sum2 = .0;
	For(i,n) sum2 += x[i] * x[i];

	return (1./n) * (sum2 - n * pow(mean,2));
}

/*********************************************************/

m3ldbl Mean(m3ldbl *x, int n)
{
	int i;
	m3ldbl sum;

	sum = .0;

	For(i,n) sum += x[i];

	return sum / n;
}

/*********************************************************/

void Best_Of_NNI_And_SPR(arbre *tree)
{
	int i;
	int n_l = tree->n_l;
	if(tree->mod->s_opt->random_input_tree) Speed_Spr_Loop(tree); /* Don't do simultaneous NNIs if starting tree is random */
	else
	{
		arbre *ori_tree,*best_tree;
		model *ori_mod,*best_mod;
		m3ldbl **ori_bl,**best_bl;
		m3ldbl best_lnL,ori_lnL,nni_lnL,spr_lnL;
		ori_bl = (m3ldbl **)mCalloc(tree->n_l,sizeof(m3ldbl *));
		best_bl = (m3ldbl **)mCalloc(tree->n_l,sizeof(m3ldbl *));
		For(i,tree->n_l){ // VHS: I think John added this loop?
			ori_bl[i] = (m3ldbl *)mCalloc(2*tree->n_otu-3,sizeof(m3ldbl));
			best_bl[i] = (m3ldbl *)mCalloc(2*tree->n_otu-3,sizeof(m3ldbl));
		}

		ori_mod   = Copy_Model(tree->mod);
		best_mod  = Copy_Model(tree->mod);

		ori_tree = Make_Tree(tree->n_otu, tree->n_l);
		Init_Tree(ori_tree,tree->n_otu, tree->n_l);
		Make_All_Tree_Nodes(ori_tree);
		Make_All_Tree_Edges(ori_tree);

		best_tree = Make_Tree(tree->n_otu, tree->n_l);
		Init_Tree(best_tree,tree->n_otu, tree->n_l);
		Make_All_Tree_Nodes(best_tree);
		Make_All_Tree_Edges(best_tree);

		Copy_Tree(tree,ori_tree);
		Record_Br_Len(ori_bl,tree);

		best_lnL = UNLIKELY;
		Lk(tree);
		ori_lnL = tree->c_lnL; /* Record likelihood of the starting tree */

		Simu_Loop(tree); /* Perform simultaneous NNIs */
		best_lnL = tree->c_lnL; /* Record the likelihood */
		nni_lnL = tree->c_lnL;
		Copy_Tree(tree,best_tree); /* Record the tree topology and branch lengths */
		Record_Br_Len(best_bl,tree);
		Restore_Br_Len(best_bl,best_tree);
		Record_Model(tree->mod,best_mod);

		Copy_Tree(ori_tree,tree); /* Back to the original tree topology */
		Restore_Br_Len(ori_bl,tree); /* Back to the original branch lengths */
		Record_Model(ori_mod,tree->mod); /* Back to the original model */

		/* Make sure the tree is in its original form */
		Lk(tree);
		if(fabs(tree->c_lnL - ori_lnL) > tree->mod->s_opt->min_diff_lk_global)
		{
			PhyML_Printf("\n. ori_lnL = %f, c_lnL = %f",ori_lnL,tree->c_lnL);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}

		Speed_Spr_Loop(tree);
		spr_lnL = tree->c_lnL;
		if(tree->c_lnL > best_lnL)
		{
			best_lnL = tree->c_lnL;
			Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
			Record_Br_Len(best_bl,tree);
			Restore_Br_Len(best_bl,best_tree);
			Record_Model(tree->mod,best_mod);
		}

		Copy_Tree(best_tree,tree);
		Restore_Br_Len(best_bl,tree);
		Record_Model(best_mod,tree->mod);

		/* Make sure the current tree has the best topology, branch lengths and model parameters */
		Lk(tree);
		if(fabs(tree->c_lnL - best_lnL) > tree->mod->s_opt->min_diff_lk_global)
		{
			PhyML_Printf("\n. best_lnL = %f, c_lnL = %f",best_lnL,tree->c_lnL);
			PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		}

		if(tree->mod->s_opt->print)
		{
			PhyML_Printf("\n\n. Log likelihood obtained after NNI moves : %f",nni_lnL);
			PhyML_Printf("\n. Log likelihood obtained after SPR moves : %f",spr_lnL);
		}
		For(i,n_l){

			if(ori_bl[i]) Free(ori_bl[i]);
			if(best_bl[i]) Free(best_bl[i]);
		}
		if(ori_bl) Free(ori_bl);
		if(best_bl) Free(best_bl);
		//JSJ: currently problems in Free_Tree... causes bus error...
		if(ori_tree) Free_Tree(ori_tree);
		if(best_tree) Free_Tree(best_tree);

		if(ori_mod) Free_Model(ori_mod);
		if(best_mod) Free_Model(best_mod);

	}
}

/*********************************************************/

/* Polynomial interpolation. Adapted from "Numerical Recipes in C".
Press, Flannery, Teukolsky, Vetterling, 1988.
 */
int Polint(m3ldbl *xa, m3ldbl *ya, int n, m3ldbl x, m3ldbl *y, m3ldbl *dy)
{
	int i,m,ns=1;
	m3ldbl den,dif,dift,ho,hp,w;
	m3ldbl *c,*d;

	dif=fabs(x-xa[1]);

	c = (m3ldbl *)mCalloc(n,sizeof(m3ldbl));
	d = (m3ldbl *)mCalloc(n,sizeof(m3ldbl));

	for(i=1;i<=n;i++)
	{
		if((dift=fabs(x-xa[i])) < dif)
		{
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}

	*y=ya[ns--];

	for (m=1;m<n;m++)
	{
		for (i=1;i<=n-m;i++)
		{
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
			{
				/* 	       Rprintf("\n. Error in routine POLINT.\n"); */
				Exit("\n. Error in routine POLINT.\n");
				return(-1);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}

	Free(d);
	Free(c);
	return(0);
}

/*********************************************************/

void JF(arbre *tree)
{
	//printing loglk for each site, to compute SH-like tests */
	m3ldbl sum=0.0;
	PhyML_Printf("\n\nSITES LKS:\n");
	int n_patterns = (int)floor(tree->n_pattern*tree->prop_of_sites_to_consider);
	int site=0;
	For(site,n_patterns) {
		int wei=0;
		For(wei,tree->data->wght[site]) {
			PhyML_Printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]);
			sum+=tree->c_lnL_sorted[site] / tree->data->wght[site];
		}
	}

	PhyML_Printf("\n\nsum=%f\n\n",sum);
	int i=0;
	For(i,2*tree->n_otu-3)
	{
		if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax))
		{
			PhyML_Printf("%3d %f %f %f\n",
					tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l);
		}
	}


	/*   //printing loglk for each site, to compute SH-like tests */
	/*   m3ldbl sum=0.0; */
	/*   PhyML_Printf("\n\nSITES LKS:\n"); */
	/*   int n_patterns = (int)floor(tree->n_pattern*tree->prop_of_sites_to_consider); */
	/*   int site=0; */
	/*   For(site,n_patterns) { */
	/*     int wei=0; */
	/*     For(wei,tree->data->wght[site]) { */
	/*       PhyML_Printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]); */
	/*       sum+=tree->c_lnL_sorted[site] / tree->data->wght[site]; */
	/*     } */
	/*   } */

	/*   PhyML_Printf("\n\nsum=%f\n\n",sum); */

	/*   int i=0; */
	/*   For(i,2*tree->n_otu-3) */
	/*     { */
	/*       if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax)) */
	/* 	{ */
	/* 	  PhyML_Printf("%3d %f %f %f\n", */
	/* 		 tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l); */
	/* 	} */
	/*     } */
}

/*********************************************************/

arbre *Dist_And_BioNJ(allseq *alldata, model *mod, option *io)
{
	arbre *tree;
	matrix *mat;
//	int i;


	if(!io->quiet) PhyML_Printf("\n. Computing pairwise distances...\n");

	mat = ML_Dist(alldata,mod);
	Fill_Missing_Dist(mat);

	if(!io->quiet) PhyML_Printf("\n. Building BioNJ tree...\n");
	// JSJ: !!!!!!! need to get n_l from io
	// for now just use a temp n_l
	int n_l = 2;
	mat->tree = Make_Tree_From_Scratch(alldata->n_otu,alldata,n_l);
	Bionj(mat);
	tree      = mat->tree;
	tree->mat = mat;

	Fix_Tree_From_IO(tree,io);

	return tree;
}

/*********************************************************/

void Fix_Tree_From_IO(arbre *tree, option *io){
	int i,j,k;
	int num_node = 2*tree->n_otu-2; //node
	int num_edge = 2*tree->n_otu-3; //edge
	//	int num_path = 2*tree->n_otu; //path
	//JSJ: now copy the branch lengths from this tree into all members of the set
	tree->n_l = io->n_l;
	tree->props[io->n_l - 1] = io->props[io->n_l - 1];//since it doesn't get to last element, fill in here
	for(i = 1; i < tree->n_l; i++){
		tree->props[i-1] = io->props[i-1];
		//		For(j,num_path){
		//			tree->curr_path[j]->n_l = io->n_l;
		//			For(k,3){
		//				tree->curr_path[j]->l[i][k] = tree->curr_path[j]->l[0][k];
		//				tree->curr_path[j]->b[k]->best_l[i] = tree->curr_path[j]->b[k]->best_l[0];
		//				tree->curr_path[j]->b[k]->l[i] = tree->curr_path[j]->b[k]->l[0];
		//				tree->curr_path[j]->b[k]->l_old[i] =  tree->curr_path[j]->b[k]->l_old[0];
		//				tree->curr_path[j]->b[k]->has_zero_br_len[i] =  tree->curr_path[j]->b[k]->has_zero_br_len[0];
		//				tree->curr_path[j]->b[k]->n_l = io->n_l;
		//			}
		//		}

		For(j,num_node){
			tree->noeud[j]->n_l = io->n_l;
			For(k,3){
				tree->noeud[j]->l[i][k] = tree->noeud[j]->l[0][k];
				//				tree->noeud[j]->b[k]->best_l[i] = tree->noeud[j]->b[k]->best_l[0];
				//				tree->noeud[j]->b[k]->n_l = io->n_l;
				//				tree->noeud[j]->b[k]->has_zero_br_len[i] = tree->noeud[j]->b[k]->has_zero_br_len[0];
				//				tree->noeud[j]->b[k]->l[i] = tree->noeud[j]->b[k]->l[0];
				//				tree->noeud[j]->b[k]->l_old[i] = tree->noeud[j]->b[k]->l_old[0];
			}
		}

		For(j,num_edge){
			tree->t_edges[j]->l[i] = tree->t_edges[j]->l[0];
			tree->t_edges[j]->l_old[i] = tree->t_edges[j]->l_old[0];
			tree->t_edges[j]->best_l[i] = tree->t_edges[j]->best_l[0];
			tree->t_edges[j]->n_l = io->n_l;
			tree->t_edges[j]->has_zero_br_len[i] = tree->t_edges[j]->has_zero_br_len[0];
		}
	}

}



void Add_BioNJ_Branch_Lengths(arbre *tree, allseq *alldata, model *mod, option *io)
{
	matrix *mat;

	PhyML_Printf("\n. Computing branch length estimates...\n");

	Order_Tree_CSeq(tree,alldata);
	mat = ML_Dist(alldata,mod);
	mat->tree = tree;
	mat->method = 0;
	Bionj_Br_Length(mat);

	Fix_Tree_From_IO(tree, io); //JSJ: coppy over branch lengths

	Free_Mat(mat);
}

/*********************************************************/

arbre *Read_User_Tree(allseq *alldata, model *mod, option *io)
{
	arbre *tree;


	PhyML_Printf("\n. Reading tree...\n"); fflush(NULL);
	if(io->n_trees == 1) rewind(io->fp_in_tree);
	tree = Read_Tree_File(io->fp_in_tree);
	if(!tree) Exit("\n. Input tree not found...\n");
	/* Add branch lengths if necessary */
	if(!tree->has_branch_lengths) Add_BioNJ_Branch_Lengths(tree,alldata,mod,io);

	return tree;
}

/*********************************************************/

void Print_Time_Info(time_t t_beg, time_t t_end)
{
	div_t hour,min;

	hour = div(t_end-t_beg,3600);
	min  = div(t_end-t_beg,60  );
	min.quot -= hour.quot*60;

	PhyML_Printf("\n\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
	PhyML_Printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/

char *Bootstrap_From_String(char *s_tree, allseq *alldata, model *mod, option *io)
{
	arbre *tree;

	tree = Read_Tree(s_tree);

	if(!tree)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	tree->mod         = mod;
	tree->io          = io;
	tree->data        = alldata;
	tree->both_sides  = 1;
	tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

	Order_Tree_CSeq(tree,alldata);
	if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
	Fill_Dir_Table(tree);
	Update_Dirs(tree);
	Make_Tree_4_Pars(tree,alldata,alldata->init_len);
	Make_Tree_4_Lk(tree,alldata,alldata->init_len);
	tree->triplet_struct = Make_Triplet_Struct(mod);
	Br_Len_Not_Involving_Invar(tree);
	Make_Spr_List(tree);
	Make_Best_Spr(tree);

#ifdef MPI
	Bootstrap_MPI(tree);
#else
	Bootstrap(tree);
#endif

	Free(s_tree);
	s_tree = Write_Tree(tree);

	Free_Spr_List(tree);
	Free_One_Spr(tree->best_spr);
	Free_Triplet(tree->triplet_struct);
	Free_Tree_Pars(tree);
	Free_Tree_Lk(tree);
	Free_Tree(tree);

	return s_tree;
}

/*********************************************************/

char *aLRT_From_String(char *s_tree, allseq *alldata, model *mod, option *io)
{
	arbre *tree;

	tree = Read_Tree(s_tree);

	if(!tree)
	{
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("");
	}

	tree->mod         = mod;
	tree->io          = io;
	tree->data        = alldata;
	tree->both_sides  = 1;
	tree->n_pattern   = tree->data->crunch_len/tree->mod->stepsize;

	Order_Tree_CSeq(tree,alldata);
	if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
	Fill_Dir_Table(tree);
	Update_Dirs(tree);
	Make_Tree_4_Pars(tree,alldata,alldata->init_len);
	Make_Tree_4_Lk(tree,alldata,alldata->init_len);
	tree->triplet_struct = Make_Triplet_Struct(mod);
	Br_Len_Not_Involving_Invar(tree);
	Make_Spr_List(tree);
	Make_Best_Spr(tree);

	tree->both_sides = 1;
	Lk(tree);

	aLRT(tree);

	Free(s_tree);
	s_tree = Write_Tree(tree);

	Free_Spr_List(tree);
	Free_One_Spr(tree->best_spr);
	Free_Triplet(tree->triplet_struct);
	Free_Tree_Pars(tree);
	Free_Tree_Lk(tree);
	Free_Tree(tree);

	return s_tree;
}

/*********************************************************/

void Prepare_Tree_For_Lk(arbre *tree)
{
	Order_Tree_CSeq(tree,tree->data);
	if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
	Fill_Dir_Table(tree);
	Update_Dirs(tree);
	Make_Tree_4_Pars(tree,tree->data,tree->data->init_len);
	Make_Tree_4_Lk(tree,tree->data,tree->data->init_len);
	tree->triplet_struct = Make_Triplet_Struct(tree->mod);
	Br_Len_Not_Involving_Invar(tree);
	Make_Spr_List(tree);
	Make_Best_Spr(tree);
}

/*********************************************************/

void PhyML_Printf(char *format, ...)
{
	va_list ptr;

#ifdef MPI
	if(Global_myRank == 0)
	{
		va_start (ptr, format);
		vprintf (format, ptr);
		va_end(ptr);
	}
#else
	va_start (ptr, format);
	vprintf (format, ptr);
	va_end(ptr);
#endif

	fflush (NULL);
}

/*********************************************************/

void PhyML_Fprintf(FILE *fp, char *format, ...)
{
	va_list ptr;

#ifdef MPI
	if(Global_myRank == 0)
	{
		va_start (ptr, format);
		vfprintf (fp,format, ptr);
		va_end(ptr);
	}
#else
	va_start (ptr, format);
	vfprintf (fp,format, ptr);
	va_end(ptr);
#endif

	fflush (NULL);
}

/*********************************************************/

void Find_Common_Tips(arbre *tree1, arbre *tree2)
{
	int i,j;

	For(i,tree1->n_otu) tree1->noeud[i]->common = 0;
	For(i,tree2->n_otu) tree2->noeud[i]->common = 0;

	For(i,tree1->n_otu)
	{
		For(j,tree2->n_otu)
		{
			if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
			{
				tree1->noeud[i]->common = 1;
				tree2->noeud[j]->common = 1;
				break;
			}
		}
	}
}

/*********************************************************/
//JSJ: returns the summation of the weighted branch lengths from branch length sets
m3ldbl Get_Tree_Size(arbre *tree)
{
	int i,k;
	m3ldbl tree_size;

	tree_size = 0.0;
	//JSJ: temporary compilation fix
	For(k,tree->n_l){
		For(i,2*tree->n_otu-3) tree_size += (tree->t_edges[i]->l[k] * tree->props[k]);
	}

	tree->size = tree_size;

	return tree_size;

}

/*********************************************************/

/* check whther target_bip can be found within tree. WARNING: target_bip names
   must be sorted in alphabetical order */
int Find_Bipartition(char **target_bip, int bip_size, arbre *tree)
{
	int score,i,j;
	edge *b;

	score = 0;

	For(i,2*tree->n_otu-3)
	{
		b = tree->t_edges[i];

		/*       printf("\n. %d %d",b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

		if(b->left->bip_size[b->l_r] == bip_size)
		{
			For(j,bip_size)
			{
				/* 	      printf("%s %s\n",b->left->bip_name[b->l_r][j]); */
				if(strcmp(b->left->bip_name[b->l_r][j],target_bip[j])) break;
			}
			if(j == bip_size) score++;
		}

		if(b->rght->bip_size[b->r_l] == bip_size)
		{
			For(j,bip_size)
			{
				/* 	      printf("%s %s\n",b->rght->bip_name[b->r_l][j]); */
				if(strcmp(b->rght->bip_name[b->r_l][j],target_bip[j])) break;
			}
			if(j == bip_size) score++;
		}
	}

	return score;
}

/*********************************************************/

//JSJ: compilation fix, returns weighted distance to root
void Dist_To_Root_Pre(node *a, node *d, edge *b, arbre *tree)
{
	int i,k;

	if(b) {
		For(k,tree->n_l){
			d->dist_to_root = a->dist_to_root + (b->l[k] * tree->props[k]);
		}
	}

	if(d->tax) return;
	else
	{
		For(i,3)
		if((d->v[i] != a) && (d->b[i] != tree->e_root))
			Dist_To_Root_Pre(d,d->v[i],d->b[i],tree);
	}
}

/*********************************************************/
//JSJ: compilation fix, returns weighted distance to root
void Dist_To_Root(node *n_root, arbre *tree)
{
	int k;
	n_root->v[0]->dist_to_root = tree->n_root->l[0][0];
	n_root->v[1]->dist_to_root = tree->n_root->l[0][1];
	for(k=1;k<tree->n_l;k++){
		n_root->v[0]->dist_to_root += (tree->n_root->l[k][0] * tree->props[k]);
		n_root->v[1]->dist_to_root += (tree->n_root->l[k][1] * tree->props[k]);
	}
	Dist_To_Root_Pre(n_root,n_root->v[0],NULL,tree);
	Dist_To_Root_Pre(n_root,n_root->v[1],NULL,tree);
}

/*********************************************************/
/* 'Borrowed' fromn libgen */
char *Basename(char *path)
{
	char *p;

	if( path == NULL || *path == '\0' ) return ".";

	p = path + strlen(path) - 1;

	while( *p == '/' )
	{
		if( p == path ) return path;
		*p-- = '\0';
	}

	while( p >= path && *p != '/' ) p--;

	return p + 1;
}

/*********************************************************/

m3ldbl *Matrix_Mult(m3ldbl *A, m3ldbl *B, int nra, int nca, int nrb, int ncb)
{
	int i,j,k;
	m3ldbl *C;

	C = (m3ldbl *)mCalloc(nra*ncb,sizeof(m3ldbl));

	if(nca != nrb)
	{
		PhyML_Printf("\n. Matrices dimensions don't match.");
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("\n");
	}

	For(i,nra)
	For(j,ncb)
	For(k,nca)
	C[i*ncb+j] += A[i*nca+k] * B[k*ncb+j];

	return C;
}

/*********************************************************/

m3ldbl *Matrix_Transpose(m3ldbl *A, int dim)
{
	m3ldbl *tA,buff;
	int i,j;

	tA = (m3ldbl *)mCalloc(dim*dim,sizeof(m3ldbl));

	For(i,dim*dim) tA[i]=A[i];

	For(i,dim) for(j=i+1;j<dim;j++)
	{
		buff        = tA[i*dim+j];
		tA[i*dim+j] = tA[j*dim+i];
		tA[j*dim+i]  = buff;
	}

	return tA;
}

/*********************************************************/
