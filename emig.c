#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"
#include "models.h"
#include "mc.h"
#include "rates.h"
#include "annealing.h"
#include "numeric.h"
#include "spr.h"

//Input: a directed acyclic graph with n leaves (i.e. a phylogeny).  Each leaf node corresponds to a unique taxon; each internal node corresponds to an ancestral taxon.  The length (or weight) of each edge indicates the phylogenetic distance between two taxa.
	//
	//Output: a new phylogeny, which will have different branch lengths than the input phylogeny and potentially a different branching pattern than the input tree.

void Emig_Swap(arbre *tree, edge *ea, node *v_ea, edge *eb, node *v_eb)
{
	// The swap operation should exchange ea's pointer to node v_ea with
	// eb's pointer to node ev_b.
	//
	// In other words: swap the branches, but maintain branch lengths and
	// the topology of attached clades.

	// Each edge has a left and right pointer to nodes. each node has a pointer
	// to the three or one edges and the three or one neighbor nodes as well.













	//JSJ: lets implement this using the existing Swap algorithm


	/*
	 * Swap(node *a, node *b, node *c, node *d, arbre *tree);
	 *
	 * \             /d      \             /a
	 *  \           /         \           /
	 *   \b__...__c/    ->     \b__...__c/
	 *   /         \	       /		 \
	 *  /           \	      /		      \
	 * /a            \  	 /d            \
	 *
	 * nodes b and c are not necessarily on the same branch, I think b == c is ok...
	 */

	node *a,*b,*c,*d; //four node pointers to use in swap.

	a = v_ea;
	d = v_eb;

	if(ea->left->num == a->num){
		b = ea->rght;
	}else{
		b = ea->left;
	}

	if(eb->left->num == d->num){
		c = eb->rght;
	}else{
		c = eb->left;
	}
//	if(b->num == c->num){ //this is expected, now it is safe to swap
		PhyML_Printf("Call to swap at line %d, in file %s\n",__LINE__,__FILE__);
		Swap(a,b,c,d,tree);
//	}else{
//		PhyML_Printf("Bad call to swap at line %d, in file %s\n",__LINE__,__FILE__);
//	}


}

//
// At the end of this method, th object 'tree' will be modified.
//
void Migrate_One_Edge(arbre *tree,annealing *ann)
{
	int ran,i,set,tmp;
	double r;
	int n_edges = (tree->n_otu * 2 - 3);
	edge *e1 = NULL; // VHS: all these NULL assignments are added to avoid warnings during compilation.
	edge *e2 = NULL;
	edge *e3 = NULL;
	edge *e4 = NULL;
	edge *e5 = NULL;
	node *v1,*v2;
	set = 0;
	//
	//	1. 	Select a maximum perturbation distance.
	//			Call this distance d_max.
	//			(Perhaps select d_max from a Gaussian distribution.)
	double d_max = gsl_ran_gaussian(ann->rng,ann->emig_sigma);
	//	2. 	d_total = 0.0
	double d_total = 0.0;
	//	3. 	while d_total < d_max:
	while(d_total < d_max){
	//	4. 	  	Select a random edge.  Call this edge e1.
		ran = gsl_rng_uniform_int(ann->rng,n_edges);
		e1 = tree->t_edges[ran];
	//	5.	  	Select a non-leaf vertex attached to e1.
	//		  		Call this vertex v1.
		//randomly check left or rght node first...
		r = gsl_rng_uniform(ann->rng);
		if(r > 0.5){
			if(!e1->left->tax){
				v1 = e1->left;
			}else{
				v1 = e1->rght;
			}
		}else{
			if(!e1->rght->tax){
				v1 = e1->rght;
			}else{
				v1 = e1->left;
			}
		}
	//	6.	  	Select two edges attached to v1.
	//		  		Call these edges e2 and e3.
		tmp = -1;
		for(i=0;i<3;i++){
			if(v1->b[i] == e1) tmp = i;
		}
		switch(tmp){
		case 0:
		{
			r = gsl_rng_uniform(ann->rng);
			if(r > 0.5){
				e2 = v1->b[1];
				e3 = v1->b[2];
			}else{
				e2 = v1->b[2];
				e3 = v1->b[1];
			}
			break;
		}
		case 1:
		{
			r = gsl_rng_uniform(ann->rng);
			if(r > 0.5){
				e2 = v1->b[0];
				e3 = v1->b[2];
			}else{
				e2 = v1->b[2];
				e3 = v1->b[0];
			}
			break;
		}
		case 2:
		{
			r = gsl_rng_uniform(ann->rng);
			if(r > 0.5){
				e2 = v1->b[1];
				e3 = v1->b[0];
			}else{
				e2 = v1->b[0];
				e3 = v1->b[1];
			}
			break;
		}
		default:
		{ //this is bad!
			PhyML_Printf("In default case at line %d, in file %s\n",__LINE__,__FILE__);
			break;
		}
		}//end switch(tmp)

	//		  // Migrate e1 along e3
	//
	//	7. 	  if length(e3) >= d_max - d_total:
		//since e2 and e3 are random, this move is also random
		if(e3->l[set] >= (d_max - d_total)){
	//	8.	     d_total = d_max
			d_total = d_max;
	//	9.	     length(e3) -= d_max - d_total
			e3->l[set] -= (d_max - d_total);
	//	10.	     length(e2) += d_max - d_total
			e2->l[set] += (d_max - d_total);
		}
	//	11.	  else if length(e3) < d_max - d_total:
		else if(e3->l[set] < (d_max - d_total)){
	//	12.	     d_total += length(e3)
			d_total += e3->l[set];
	//	13.	     length(e2) += length(e3)
			e2->l[set] += e3->l[set];
	//	14.	     length(e3) = 0
			e3->l[set] = 0.0;
	//	15.	     Select the other vertex attached to e3, which is not v1.
	//					Call this vertex v2.
			if(e3->left == v1){
				v2 = e3->rght;
			}else{
				v2 = e3->left;
			}
	//
	//		     // Deal with the case where e3 is terminal
	//
	//	16.	     if is_leaf(v2):
			if(v2->tax){
	//	17.	     	if length(e1) >= d_max - d_total:
				if(e1->l[set] >= (d_max - d_total)){
	//	18.		   		length(e1) -= d_max - d_total
					e1->l[set] -= (d_max - d_total);
	//	19.		   		length(e2) += d_max - d_total
					e2->l[set] += (d_max - d_total);
	//	20.		   		d_total = d_max
					d_total = d_max;
				}

	//	21.		 else if length(e1) < d_max - d_total:
				else if(e1->l[set] < (d_max - d_total)){
	//	22.		   	d_total += length(e1)
					d_total += e1->l[set];
	//	24.		   	length(e2) += length(e1)
					e2->l[set] += e1->l[set];
	//	23.		   	length(e1) = 0
					e1->l[set] = 0.0;
	//
				}
			}
	//		     // Deal with the case where e3 is internal
	//
	//	24.	     else:
			else{
	//	25.			Select two edges attached to v2.
	//				Call these edges e4 and e5
	//				e4 and e5 should not be e3; e4 and e5 should not be e1
				tmp = -1;
				for(i=0;i<3;i++){
					if(v2->b[i] == e3) tmp = i;
				}
				switch(tmp){
				case 0:
				{
					r = gsl_rng_uniform(ann->rng);
					if(r > 0.5){
						e4 = v2->b[1];
						e5 = v2->b[2];
					}else{
						e4 = v2->b[2];
						e5 = v2->b[1];
					}
					break;
				}
				case 1:
				{
					r = gsl_rng_uniform(ann->rng);
					if(r > 0.5){
						e4 = v2->b[0];
						e5 = v2->b[2];
					}else{
						e4 = v2->b[2];
						e5 = v2->b[0];
					}
					break;
				}
				case 2:
				{
					r = gsl_rng_uniform(ann->rng);
					if(r > 0.5){
						e4 = v2->b[1];
						e5 = v2->b[0];
					}else{
						e4 = v2->b[0];
						e5 = v2->b[1];
					}
					break;
				}
				default:
				{ //this is bad!
					PhyML_Printf("In default case at line %d, in file %s\n",__LINE__,__FILE__);
					break;
				}
				}//end switch(tmp)
	//	26.			Swap e1 with e4.
	//	27.			if length(e5) >= d_max - d_total:
				if(e5->l[set] >= (d_max - d_total)){
	//	28.		   		length(e5) -= d_max - d_total
					e5->l[set] -= (d_max - d_total);
					//if(e5->l[set] < 0.0) e5->l[set] = 0.0;
	//	29.		   		length(e3) += d_max - d_total
					e3->l[set] += (d_max - d_total);
					//if(e3->l[set] < 0.0) e3->l[set] = 0.0;
	//	30.		   		d_total = d_max
					d_total = d_max;
				}
	//	31.			else if length(e5) < d_max - d_total:
				else if(e5->l[set] < (d_max - d_total)){
	//	32.		   		d_total += length(e5)
					d_total += e5->l[set];
	//	33.		   		length(e5) = 0
					e5->l[set] = 0.0;
				}
				Emig_Swap(tree,e1,v1,e4,v2);
			}
		}
	}
}
