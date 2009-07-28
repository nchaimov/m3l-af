	//Input: a directed acyclic graph with n leaves (i.e. a phylogeny).  Each leaf node corresponds to a unique taxon; each internal node corresponds to an ancestral taxon.  The length (or weight) of each edge indicates the phylogenetic distance between two taxa.
	//
	//Output: a new phylogeny, which will have different branch lengths than the input phylogeny and potentially a different branching pattern than the input tree.

void Emig_swap(arbre *tree, edge *a, node *v_a, edge *b, node *v_b)
{
	// The swap operation should exchange a's pointer to node v_a with
	// b's pointer to node v_b.
	//
	// In other words: swap the branches, but maintain branch lengths and
	// the topology of attached clades.

}

//
// At the end of this method, th object 'tree' will be modified.
//
void Migrate_one_edge(arbre *tree)
{
	//
	//	1. 	Select a maximum perturbation distance.
	//			Call this distance d_max.
	//			(Perhaps select d_max from a Gaussian distribution.)
	//	2. 	d_total = 0.0
	//	3. 	while d_total < d_max:
	//	4. 	  	Select a random edge.  Call this edge e1.
	//	5.	  	Select a non-leaf vertex attached to e1.
	//		  		Call this vertex v1.
	//	6.	  	Select two edges attached to v1.
	//		  		Call these edges e2 and e3.
	//
	//		  // Migrate e1 along e3
	//
	//	7. 	  if length(e3) >= d_max - d_total:
	//	8.	     d_total = d_max
	//	9.	     length(e3) -= d_max - d_total
	//	10.	     length(e2) += d_max - d_total
	//	11.	  else if length(e3) < d_max - d_total:
	//	12.	     d_total += length(e3)
	//	13.	     length(e2) += length(e3)
	//	14.	     length(e3) = 0
	//	15.	     Select the other vertex attached to e3, which is not v1.
	//					Call this vertex v2.
	//
	//		     // Deal with the case where e3 is terminal
	//
	//	16.	     if is_leaf(v2):
	//	17.	     	if length(e1) >= d_max - d_total:
	//	18.		   		length(e1) -= d_max - d_total
	//	19.		   		length(e2) += d_max - d_total
	//	20.		   		d_total = d_max
	//	21.		 else if length(e1) < d_max - d_total:
	//	22.		   	d_total += length(e1)
	//	24.		   	length(e2) += length(e1)
	//	23.		   	length(e1) = 0
	//
	//		     // Deal with the case where e3 is internal
	//
	//	24.	     else:
	//	25.			Select two edges attached to v2.
	//				Call these edges e4 and e5
	//				e4 and e5 should not be e3; e4 and e5 should not be e1
	//	26.			Swap e1 with e4.
	//	27.			if length(e5) >= d_max - d_total:
	//	28.		   		length(e5) -= d_max - d_total
	//	29.		   		length(e3) += d_max - d_total
	//	30.		   		d_total = d_max
	//	31.			else if length(e5) < d_max - d_total:
	//	32.		   		d_total += length(e5)
	//	33.		   		length(e5) = 0
}
