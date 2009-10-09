#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "m4.h"
#include "mc.h"

#ifdef MC
#include "rates.h"
#endif

#ifdef COMPRESS_SUBALIGNMENTS
#include "compress.h"

/*
 * Compute_Red_Arrays performs a post-order traversal of the tree and compresses subalignments
 * based on the topology.  (See emails between me and Bryan for a longer description of this
 * process).
 *
 * At the end of this method, d->red will contain information about compressable sites for node d.
 */
void Compute_Red_Arrays(node *a, node *d, arbre *tree)
{
	//PhyML_Printf("entered Post_Order_Foo(node %d, node %d, ...)\n", a->num, d->num);

	int i,site,state;
	if(d->tax)
	{
		//PhyML_Printf("Post_Order_Foo: node %d is a taxon\n", d->num);
		int *state_firstsite; // key = a state number, value = the first site at which that state appears
		state_firstsite = (int *)mCalloc(tree->mod->ns + 1,sizeof(int)); //VHS: I changed this to +1 so as t consider "-" and other ambiguities.
		for(state = 0; state < tree->mod->ns+1; state++)
		{
			//PhyML_Printf("Post_Order_Foo: state_firstsite[%d] = -1\n", state);
			state_firstsite[state] = -1;
		}

		for(site = 0; site < tree->n_pattern; site++)
		{
			//PhyML_Printf("considering site %d\n", site);
			char this_state = tree->data->c_seq[ d->num ]->state[ site ];
			int this_state_id = Assign_State_With_Ambiguity(&this_state, tree->mod->datatype, 1);
			if (state_firstsite[ this_state_id ] == -1)
			{
				//PhyML_Printf("this_state_id = %d -> site = %d\n", this_state_id, site);
				state_firstsite[ this_state_id ] = site;
			}
			else
			{
				//PhyML_Printf("Post_Order_Foo: for taxon %d, site %d shares the same state as site %d\n", d->num, site, state_firstsite[ this_state_id ]);
				d->red[ site ] = state_firstsite[ this_state_id ];
			}
		}
		return;
	}
	else
	{
		node *child1 = 0;
		node *child2 = 0;

		For(i,3)
		{
			if(d->v[i] != a)
			{
				if (child1 == 0)
				{
					//PhyML_Printf("pointing child1 to node %d\n", d->v[i]->num);
					child1 = d->v[i];
				}
				else if (child2 == 0)
				{
					//PhyML_Printf("pointing child2 to node %d\n", d->v[i]->num);
					child2 = d->v[i];
				}

				//PhyML_Printf("Post_Order_Foo: traversing from node %d to node %d\n", d->num, d->v[i]->num);
				Compute_Red_Arrays(d,d->v[i],tree);
			}
		}

		//PhyML_Printf("Post_Order_Foo: node %d has children %d and %d\n", d->num, child1->num, child2->num);

		for(site = 0; site < tree->n_pattern; site++)
		{
			if (child1->red[site] != -1)
			{
				if (child1->red[site] == child2->red[site])
				{
					//PhyML_Printf("Post_Order_Foo: compressing site %d for node d=%d\n", site, d->num);
					d->red[site] = child1->red[site];
				}
			}
		}
	}
}

void Post_Order_Lk_Red(node *a, node *d, arbre *tree, int site)
{
	int i,dir;

	dir = -1; // dir will be given the index of the edge (in d->b) connecting node d to node a.

	if(d->tax) // if 'd' is terminal, then there is no partial likelihood to compute.
	{	return;
	}
	else if (d->red[site] == -1) // if d->red[site] == -1, then the partial likelihood
								 // at 'site' needs to be recomputed, and we therefore
								 // need to traverse to subtrees.
	{
		For(i,3)
		{
			if(d->v[i] != a)
				Post_Order_Lk_Red(d,d->v[i],tree,site);
			else dir = i;
		}
	}
	else{ // we need to determine which edge (in d->b) connects d to a:
		//PhyML_Printf(" . debug: Post_Order_Lk_Red: d = %d, skipping traverse to lower nodes.\n", d->num);
		For(i,3)
		{
			if (d->v[i] == a)
			{
				dir = i;
				break; // skip the traverse
			}
		}

	}
	Update_P_Lk_Red(tree, d->b[dir], d, site);
}

void Pre_Order_Lk_Red(node *a, node *d, arbre *tree, int site)
{
	int i;

	if(d->tax) // if 'd' is terminal, then there is no partial likelihood to compute.
	{	return;
	}
	else if (d->red[site] == -1) // if d->red[site] == -1, then the partial likelihood
								 // at 'site' needs to be recomputed, and we therefore
								 // need to traverse to subtrees.
	{
		For(i,3)
		{
			if(d->v[i] != a)
			{
				Update_P_Lk_Red(tree, d->b[i], d, site);
				Pre_Order_Lk_Red(d,d->v[i],tree,site);
			}
		}
	}
}

/*
 * POSTCONDITION: b->p_lk will be updated.
 */
void Update_P_Lk_Red(arbre *tree, edge *b, node *d, int site)
{
	/*
				   |
				   |<- b
				   |
				   d
				  / \
		 dir1 -> /   \ <- dir2
				/     \
			  n_v1    n_v2

	 */

	//PhyML_Printf("entered Update_P_Lk_Red(..., edge %d, node %d)\n", b->num, d->num);
	//PhyML_Printf("left = %d, right = %d\n", b->left->num, b->rght->num);

	node *n_v1, *n_v2;
	m3ldbl p1_lk1,p2_lk2;
	plkflt *p_lk,*p_lk_v1,*p_lk_v2;
	double *Pij1, *Pij2;
	plkflt max_p_lk;
	plkflt *sum_scale, *sum_scale_v1, *sum_scale_v2;
	plkflt scale_v1, scale_v2;
	int i,j,k;
	int catg;
	int dir1,dir2;
	int n_patterns;

	/**
	* ambiguity_check_v1 is used to store the site specific ambiguity of the neighbor
	* node in direction dir1. Ambiguity check v2 does the same for dir2
	*/
	int ambiguity_check_v1,ambiguity_check_v2;
	int state_v1,state_v2;
	int dima, dimb, dimc, dimaa, dimd;
	dima = tree->mod->ns;
	dimb = dima * tree->mod->n_l;
	dimc = dimb * tree->mod->n_catg;
	dimaa = dima * dima; // = ns^2
	dimd = dimb * dima; // = number of BL sets * ns^2
	state_v1 = state_v2 = -1;
	ambiguity_check_v1 = ambiguity_check_v2 = -1;
	scale_v1 = scale_v2 = 0.0;
	p1_lk1 = p2_lk2 = .0;

	if(d->tax)
	{
		PhyML_Printf("\n. node %d is a leaf...",d->num);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Warn_And_Exit("\n");
	}

	n_patterns = tree->n_pattern;

	dir1=dir2=-1;
	For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i); //VHS: here we set dir1 and dir2 to point to the two
															// edges (of three possible edges) which aren't
															// the edge pointed to by edge *b.

	if((dir1 == -1) || (dir2 == -1))
	{
		PhyML_Printf("\n. d = %d",d->num);
		PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num);
		PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num);
		PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2);
		PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	}

	n_v1 = d->v[dir1];
	n_v2 = d->v[dir2];

	/*
	 * VHS: In the remainder of this method, we work on computing p_lk
	 * as the product of p_lk_v1 and p_lk_v2. However, first we must establish
	 * where p_lk_v1 and p_lk_v2 point within our tree.
	 *
	 * 'p_lk' becomes a shortcut pointer to b->p_lk_left OR b->p_lk_rght.
	 */
	if(d == b->left)
	{
		p_lk = b->p_lk_left;
		sum_scale = b->sum_scale_f_left;
	}
	else
	{
		p_lk = b->p_lk_rght;
		sum_scale = b->sum_scale_f_rght;
	}

	// VHS: p_lk_v1 is the partial likelihood for the subtree attached to branch dir1
	if(d == d->b[dir1]->left)
	{
		p_lk_v1 = d->b[dir1]->p_lk_rght;
		sum_scale_v1 = d->b[dir1]->sum_scale_f_rght;
	}
	else
	{
		p_lk_v1 = d->b[dir1]->p_lk_left;
		sum_scale_v1 = d->b[dir1]->sum_scale_f_left;
	}

	// VHS: p_lk_v2 is the partial likelihood for the subtree attached to branch dir2
	if(d == d->b[dir2]->left)
	{
		p_lk_v2 = d->b[dir2]->p_lk_rght;
		sum_scale_v2 = d->b[dir2]->sum_scale_f_rght;
	}
	else
	{
		p_lk_v2 = d->b[dir2]->p_lk_left;
		sum_scale_v2 = d->b[dir2]->sum_scale_f_left;
	}

	Pij1 = d->b[dir1]->Pij_rr;
	Pij2 = d->b[dir2]->Pij_rr;


#ifdef USE_OPENMP
	int chunk = n_patterns/omp_get_num_procs();
	//	int chunk = n_patterns/2;
	//printf("Chunk size: %i\n",chunk);
#pragma omp parallel for\
		default(shared) private(k,catg,i,j,site,scale_v1,scale_v2,\
				max_p_lk,state_v1,state_v2,ambiguity_check_v1,\
				ambiguity_check_v2,p1_lk1,p2_lk2)\
				schedule(static,chunk)
#endif

	/**
	* JSJ: If sum_scale_v1 was assigned a non-null value in the above if/else cascade,
	* the scale_v1 value is assigned as the site specific sum_scale_v1 from above.
	*/
	scale_v1 = (sum_scale_v1)?(sum_scale_v1[site]):(0.0);
	scale_v2 = (sum_scale_v2)?(sum_scale_v2[site]):(0.0);


	// Muy importante!
	sum_scale[site] = scale_v1 + scale_v2;

	max_p_lk = -MDBL_MAX;
	state_v1 = state_v2 = -1; //just ints
	ambiguity_check_v1 = ambiguity_check_v2 = -1;


	if (d->red[site] != -1) // i.e. this site is redundant, use the partial likelihood (that we already computed) for the previous site.
	{
		int red_site = d->red[site];

		For(catg,tree->mod->n_catg)
		{
			for(k = 0; k < tree->mod->n_l; k++)
			{
				For(i,tree->mod->ns)
				{

					//PhyML_Printf(" compress line 302: plk edge %d site %d bl %d catg %d state %d = %f\n", b->num, site, k, catg, i, p_lk[red_site*dimc + catg*dimb + k*dima + i]);
					p_lk[site*dimc + catg*dimb + k*dima + i] = p_lk[red_site*dimc + catg*dimb + k*dima + i];

					if( p_lk[site*dimc + catg*dimb + k*dima + i] > max_p_lk )
					{	max_p_lk = p_lk[site*dimc + catg*dimb + k*dima + i];
					}
				}
			}
		}

		// Scale all the likelihoods
		if((max_p_lk < LIM_SCALE_VAL) || (max_p_lk > (1./LIM_SCALE_VAL)))
		{
			For(catg,tree->mod->n_catg)
			{
				for(k = 0; k < tree->mod->n_l; k++)
				{
					For(i,tree->mod->ns)
					{
						p_lk[site*dimc + catg*dimb + k*dima + i] /= max_p_lk;
					}
				}
			}
			sum_scale[site] += (plkflt)log(max_p_lk);
		}

		//PhyML_Printf(" . debug compress.c 343: edge->num=%d, sum_scale[%d]=%f\n", b->num, site, sum_scale[site]);
		//PhyML_Printf(" . debug compress.c 343: max_p_lk = %f\n", max_p_lk);

		return; // We return, because we don't need to recompute partial likelihoods for this node.
	}


	/**
	* JSJ: If the model is set to greedy, and the node is terminal,
	*  then we check if the site is ambiguous
	* 	otherwise if it is not ambiguous we default to true
	*/
	if(!tree->mod->s_opt->greedy)
	{
		if(n_v1->tax)
		{
			ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
			if(!ambiguity_check_v1)
			{
				state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dima,tree);
			}
		}

		if(n_v2->tax)
		{
			ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
			if(!ambiguity_check_v2)
			{
				state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dima,tree);
			}
		}
	}

	if(tree->mod->use_m4mod)
	{
		ambiguity_check_v1 = 1;
		ambiguity_check_v2 = 1;
	}

	For(i,tree->mod->ns)
	{
		For(catg,tree->mod->n_catg)
		{
			for(k = 0; k < tree->mod->n_l; k++)
			{
				p1_lk1 = .0;

				// if n_v1 is terminal, then use it's known state...
				if((n_v1->tax) && (!tree->mod->s_opt->greedy))
				{
					if(!ambiguity_check_v1)
					{
						p1_lk1 = Pij1[catg*dimd + k*dimaa + i*dima + state_v1];
					}
					else
					{
						For(j,tree->mod->ns)
						{
							p1_lk1 += Pij1[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)n_v1->b[0]->p_lk_tip_r[site*dima+j];
						}
					}
				}
				else // otherwise, n_v1 is non-terminal, so we need to consider all possible states...
				{
					For(j,tree->mod->ns)
					{
						p1_lk1 += Pij1[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)p_lk_v1[site*dimc + catg*dimb + k*dima + j];   //[site*dim1+catg*dim2+j];
					}
				}

				p2_lk2 = .0;

				if((n_v2->tax) && (!tree->mod->s_opt->greedy))
				{
					if(!ambiguity_check_v2)
					{
						p2_lk2 = Pij2[catg*dimd + k*dimaa + i*dima + state_v2];
					}
					else
					{
						For(j,tree->mod->ns)
						{
							p2_lk2 += Pij2[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)n_v2->b[0]->p_lk_tip_r[site*dima+j];
						}
					}
				}
				else
				{
					For(j,tree->mod->ns)
					{
						p2_lk2 += Pij2[catg*dimd + k*dimaa + i*dima + j] * (m3ldbl)p_lk_v2[site*dimc + catg*dimb + k*dima + j];
					}
				}

				// Muy importante!
				p_lk[site*dimc + catg*dimb + k*dima + i] = (plkflt)(p1_lk1 * p2_lk2);

				if( p_lk[site*dimc + catg*dimb + k*dima + i] > max_p_lk )
				{
					max_p_lk = p_lk[site*dimc + catg*dimb + k*dima + i];
				}

			}
		}
	}


	// Scale all the likelihoods
	if((max_p_lk < LIM_SCALE_VAL) || (max_p_lk > (1./LIM_SCALE_VAL)))
	{
		For(catg,tree->mod->n_catg)
		{
			for(k = 0; k < tree->mod->n_l; k++)
			{
				For(i,tree->mod->ns)
				{
					p_lk[site*dimc + catg*dimb + k*dima + i] /= max_p_lk;
				}
			}
		}
		sum_scale[site] += (plkflt)log(max_p_lk);
	}

	//PhyML_Printf(" . debug compress.c 464: edge->num=%d, sum_scale[%d]=%f\n", b->num, site, sum_scale[site]);
}

/*********************************************************/
// This method assumes that tree->n_pattern is initialized to a real value (not -1).
//
// This method allocates memory for red arrays AND fills all the arrays with values -1.
void Make_All_Nodes_Red(arbre *tree)
{
	PhyML_Printf(". Allocating memory for red. arrays. . .\n");

	int i;
	int k;

	For(i,2*tree->n_otu-2)
	{
		// VHS: we allocate memory space for the redundant site array here, instead of inside
		// Make_Node_Light, because we need to know how many patterns are in the alignment.
		// An optimization (to-do) would be to leave ->red size 0, thus saving memory space,
		// and allocate more memory everytime we want to add information to ->red.
		//PhyML_Printf("Make_Node_Red: calling mCalloc(%d, %d)\n", tree->n_pattern, sizeof(int));
		tree->noeud[i]->red = (int *)mCalloc( tree->n_pattern, sizeof(int) );
		For(k, tree->n_pattern)
		{
			tree->noeud[i]->red[k] = -1; // we initialize all values to -1.
		}
	}
}

//
// This method assumes that the memory for red arrays has already been allocated.
// This method fills all (pre-allocated) red arrays with the value -1.
//
void Init_All_Nodes_Red(arbre *tree)
{
	//PhyML_Printf(". Resetting red arrays\n");

	int i;
	int k;

	For(i,2*tree->n_otu-2)
	{
		For(k, tree->n_pattern)
		{
			tree->noeud[i]->red[k] = -1; // we initialize all values to -1.
		}
	}
}

#endif // end of SUBALIGNMENT_COMPRESSION
