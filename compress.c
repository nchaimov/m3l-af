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
 * This method fills red arrays using a first pass post-order recursive tree traversal.
 * It's a so-called "firt pass" traversal because not *all* red arrays will be filled;
 * Every edge will have either it's red_left or red_right filled, but not both.
 * The other red array can be filled during a second pass, if necessary.  However,
 * I don't think a second pass is necessary because
 *
 *
 * b is the edge on which we'll fill the red. arrays
 * from is the node connected to b, indicating the direction from which we are traversing
 * tree is. . . the tree, duh.
 */
void Fill_Red_Arrays_First(edge *b, node* from, arbre *tree)
{

	node *down; // down points to the node at the opposite side of 'b' from 'from'.
	int* b_red; // b_red points to either b->red_left or b->red_right, whichever is on the down side of b.
	if (b->left->num == from->num)
	{	down = b->rght;
		b_red = b->red_right;
	}
	else
	{	down = b->left;
		b_red = b->red_left;
	}


	// 1. if 'down' is a terminal taxa, then fill down->red
	if (down->tax)
	{
		int *state_firstsite; // key = a state number, value = the first site at which that state appears
		state_firstsite = (int *)mCalloc(tree->mod->ns + 1,sizeof(int)); //VHS: I changed this to +1 so as to consider "-" and other ambiguities.
		int state;
		for(state = 0; state < tree->mod->ns+1; state++)
		{	state_firstsite[state] = -1;
		}

		int site;
		for(site = 0; site < tree->n_pattern; site++)
		{
			//PhyML_Printf("considering site %d\n", site);
			char this_state = tree->data->c_seq[ down->num ]->state[ site ];
			int this_state_id = Assign_State_With_Ambiguity(&this_state, tree->mod->datatype, 1);
			if (state_firstsite[ this_state_id ] == -1)
			{
				//PhyML_Printf("this_state_id = %d -> site = %d\n", this_state_id, site);
				state_firstsite[ this_state_id ] = site;
			}
			else
			{
				//PhyML_Printf("Post_Order_Foo: for taxon %d, site %d shares the same state as site %d\n", d->num, site, state_firstsite[ this_state_id ]);
				b_red[ site ] = state_firstsite[ this_state_id ];
			}
		}
		Free( state_firstsite );
	}
	else
	{
		// 2a. Otherwise, recur down the tree

		// 'down' has three branches.  We need to determine which two branches != b.
		// b2 and b3 will point to those two branches.
		edge *b2 = NULL;
		edge *b3 = NULL;
		int i;
		For(i,3)
		{
			if ( (down->b[i] != b) && (b2 == NULL) )
			{
				b2 = down->b[i];
			}
			else if ( (down->b[i] != b) && (b3 == NULL) )
			{
				b3 = down->b[i];
			}
		}
		Fill_Red_Arrays_First(b2, down, tree);
		Fill_Red_Arrays_First(b3, down, tree);

		int *b2_red; // points to the red array that we just filled in b2.
		int *b3_red;
		if (b2->left->num == down->num)
		{
			b2_red = b2->red_right;
		}
		else if (b2->rght->num == down->num)
		{
			b2_red = b2->red_left;
		}
		if (b3->left->num == down->num)
		{
			b3_red = b3->red_right;
		}
		else if (b3->rght->num)
		{
			b3_red = b3->red_left;
		}

		// 2b. And then fill down->red using the downstream red arrays (from the previous recursion).
		int site;
		for(site = 0; site < tree->n_pattern; site++)
		{
			if (b2_red[site] == b3_red[site])
			{
				b_red[site] = b2_red[site];
			}
		}
	}

}



/*********************************************************/
// This method assumes that tree->n_pattern is initialized to a real value (not -1).
//
// This method allocates memory for red arrays AND fills all the arrays with values -1.
void Make_All_Edges_Red(arbre *tree)
{
	PhyML_Printf(". Allocating memory for red. arrays. . .\n");

	int i;
	int k;

	For(i, 2*tree->n_otu-3) /* For every branch */
	{
		// VHS: we allocate memory space for the redundant site array here, instead of inside
		// Make_Node_Light, because we need to know how many patterns are in the alignment.
		// An optimization (to-do) would be to leave ->red size 0, thus saving memory space,
		// and allocate more memory everytime we want to add information to ->red.
		//PhyML_Printf("Make_Node_Red: calling mCalloc(%d, %d)\n", tree->n_pattern, sizeof(int));
		tree->t_edges[i]->red_left = (int *)mCalloc( tree->n_pattern, sizeof(int) );
		tree->t_edges[i]->red_right = (int *)mCalloc( tree->n_pattern, sizeof(int) );
	}
}

//
// This method assumes that the memory for red arrays has already been allocated.
// This method fills all (pre-allocated) red arrays with the value -1.
//
void Init_All_Edges_Red(arbre *tree)
{
	PhyML_Printf(" . debug compress.c 158: initializing red arrays\n");

	int i;
	int k;

	For(i, 2*tree->n_otu-3) /* For every branch */
	{
		For(k, tree->n_pattern)
		{
			tree->t_edges[i]->red_left[k] = -1;
			tree->t_edges[i]->red_right[k] = -1;
		}
	}

	Fill_Red_Arrays_First(tree->noeud[0]->b[0], tree->noeud[0], tree);

	/*
	For(i, 2*tree->n_otu-3) // For every branch
	{
		PhyML_Printf("[%d] ", i);
		For(k, tree->n_pattern)
		{
			PhyML_Printf(" [site %d:]", k);
			PhyML_Printf(" %d", tree->t_edges[i]->red_left[k]);
			PhyML_Printf(" %d", tree->t_edges[i]->red_right[k]);
		}
		PhyML_Printf("\n");
	}
	*/
}

#endif // end of SUBALIGNMENT_COMPRESSION
