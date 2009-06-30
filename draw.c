/*

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "utilities.h"
#include "draw.h"


/*********************************************************/

void DR_Get_Tree_Coord(arbre *tree)
{
  DR_Init_Tdraw_Struct(tree->ps_tree);
  DR_Get_Tree_Box_Width(tree->ps_tree,tree);
  if(!tree->n_root) 
    {
      PhyML_Printf("\n. Adding root before rendering the tree.");
      Add_Root(tree->t_edges[0],tree);
    }
  else Update_Root_Pos(tree);
  Dist_To_Root(tree->n_root,tree);
  tree->ps_tree->max_dist_to_root = DR_Get_Max_Dist_To_Root(tree);
  DR_Get_X_Coord(tree->ps_tree,tree);
  DR_Get_Y_Coord(tree->ps_tree,tree);
}

/*********************************************************/

void DR_Print_Postscript_Header(int n_pages, FILE *fp)
{
  if(!fp)
    {
      PhyML_Printf("\n. Failed to open the postscript file.");
      PhyML_Printf("\n. Did you forget the '--ps' option ?.");
      Warn_And_Exit("\n");
    }

  PhyML_Fprintf(fp,"%%!PS-Adobe-3.0\n");
  PhyML_Fprintf(fp,"%%%%DocumentFonts: Times-Roman Times-Roman\n");
  PhyML_Fprintf(fp,"%%%%Creator: Stephane Guindon\n");
  PhyML_Fprintf(fp,"%%%%Title: tree\n");
  PhyML_Fprintf(fp,"%%%%BeginFeature: *PageSize\n"); 
  PhyML_Fprintf(fp,"a4\n");
  PhyML_Fprintf(fp,"%%%%EndFeature\n");
  PhyML_Fprintf(fp,"%%%%EndComments\n");
  PhyML_Fprintf(fp,"%%%%Pages: %d\n",n_pages);

  PhyML_Fprintf(fp,"/lt {lineto} bind def\n");
  PhyML_Fprintf(fp,"/mt {moveto} bind def\n");
  PhyML_Fprintf(fp,"/sc {setrgbcolor} bind def\n");

  PhyML_Fprintf(fp,"/clipbox\n");
  PhyML_Fprintf(fp,"{\n");
  PhyML_Fprintf(fp,"newpath\n");
  PhyML_Fprintf(fp,"40 40 moveto\n");
  PhyML_Fprintf(fp,"560 40 lineto\n");
  PhyML_Fprintf(fp,"560 820 lineto\n");
  PhyML_Fprintf(fp,"40 820 lineto\n");
  PhyML_Fprintf(fp,"40 40 lineto\n");
  PhyML_Fprintf(fp,"closepath\n");
  PhyML_Fprintf(fp,"clip\n");
  PhyML_Fprintf(fp,"} bind def\n");
  
  PhyML_Fprintf(fp,"/Times-Roman findfont\n");
  PhyML_Fprintf(fp,"12 scalefont\n");
  PhyML_Fprintf(fp,"setfont\n");


}

/*********************************************************/

void DR_Print_Postscript_EOF(FILE *fp)
{
  PhyML_Fprintf(fp,"%%%%Trailer\n");
  PhyML_Fprintf(fp,"%%%%EOF\n");
}

/*********************************************************/

void DR_Print_Tree_Postscript(int page_num, FILE *fp, arbre *tree)
{
  int i;
  tdraw *draw;
  node *n_root;
  
  printf("coucou\n");

  draw = tree->ps_tree;
  DR_Get_Tree_Coord(tree);
  n_root = tree->n_root;


  PhyML_Fprintf(fp,"%%%%Page: %d %d\n",page_num,page_num); 
  PhyML_Fprintf(fp,"clipbox\n");
  PhyML_Fprintf(fp,"stroke\n");
  PhyML_Fprintf(fp,"50 50 translate\n");
  PhyML_Fprintf(fp,"newpath\n");

/*   if(b_root->prob_sel_regime <= 0.1) */
/*     PhyML_Fprintf(fp,".0 .0 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.1 && b_root->prob_sel_regime <= 0.2) */
/*     PhyML_Fprintf(fp,".0 .5 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.2 && b_root->prob_sel_regime <= 0.3) */
/*     PhyML_Fprintf(fp,".0 1. 1. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.3 && b_root->prob_sel_regime <= 0.4) */
/*     PhyML_Fprintf(fp,".0 1. .5 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.4 && b_root->prob_sel_regime <= 0.5) */
/*     PhyML_Fprintf(fp,".0 1. .0 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.5 && b_root->prob_sel_regime <= 0.6) */
/*     PhyML_Fprintf(fp,".5 1. .0 sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.6 && b_root->prob_sel_regime <= 0.7) */
/*     PhyML_Fprintf(fp,"1. 1. 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.7 && b_root->prob_sel_regime <= 0.8) */
/*     PhyML_Fprintf(fp,"1. .5 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.8 && b_root->prob_sel_regime <= 0.9) */
/*     PhyML_Fprintf(fp,"1. 0. 0. sc\n"); */
/*   else if(b_root->prob_sel_regime > 0.9) */
/*     PhyML_Fprintf(fp,"1. .0 .0 sc\n"); */

  PhyML_Fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[0]->num],draw->ycoord[n_root->v[0]->num]);
  PhyML_Fprintf(fp,"%d %d lt\n",0,draw->ycoord[n_root->v[0]->num]);
  PhyML_Fprintf(fp,"%d %d lt\n",0,draw->ycoord[n_root->v[1]->num]);
  PhyML_Fprintf(fp,"%d %d lt\n",draw->xcoord[n_root->v[1]->num],draw->ycoord[n_root->v[1]->num]);
  PhyML_Fprintf(fp,"stroke\n");


  PhyML_Fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[0]->num],draw->ycoord[n_root->v[0]->num]);
  if(n_root->v[0]->tax) PhyML_Fprintf(fp,"(%s) show\n",n_root->v[0]->name);
  else
    {
      For(i,3)
	if((n_root->v[0]->v[i]) && (n_root->v[0]->v[i] != n_root->v[1]))
	  DR_Print_Tree_Postscript_Pre(n_root->v[0],n_root->v[0]->v[i],fp,draw,tree);
    }

  PhyML_Fprintf(fp,"%d %d mt\n",draw->xcoord[n_root->v[1]->num],draw->ycoord[n_root->v[1]->num]);

  if(n_root->v[1]->tax) PhyML_Fprintf(fp,"(%s) show\n",n_root->v[1]->name);
  else
    {
      For(i,3)
	if((n_root->v[1]->v[i]) && (n_root->v[1]->v[i] != n_root->v[0]))
	  DR_Print_Tree_Postscript_Pre(n_root->v[1],n_root->v[1]->v[i],fp,draw,tree);
    }

  PhyML_Fprintf(fp,"closepath\n");
  PhyML_Fprintf(fp,"stroke\n");
  PhyML_Fprintf(fp,"showpage\n");
}

/*********************************************************/

void DR_Print_Tree_Postscript_Pre(node *a, node *d, FILE *fp, tdraw *w, arbre *tree)
{
  int i;

  PhyML_Fprintf(fp,"gsave\n");
  
  For(i,3)
    if(a->v[i] == d)
      {
/* 	if(a->b[i]->prob_sel_regime <= 0.1) */
/* 	  PhyML_Fprintf(fp,".0 .0 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.1 && a->b[i]->prob_sel_regime <= 0.2) */
/* 	  PhyML_Fprintf(fp,".0 .5 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.2 && a->b[i]->prob_sel_regime <= 0.3) */
/* 	  PhyML_Fprintf(fp,".0 1. 1. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.3 && a->b[i]->prob_sel_regime <= 0.4) */
/* 	  PhyML_Fprintf(fp,".0 1. .5 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.4 && a->b[i]->prob_sel_regime <= 0.5) */
/* 	  PhyML_Fprintf(fp,".0 1. .0 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.5 && a->b[i]->prob_sel_regime <= 0.6) */
/* 	  PhyML_Fprintf(fp,".5 1. .0 sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.6 && a->b[i]->prob_sel_regime <= 0.7) */
/* 	  PhyML_Fprintf(fp,"1. 1. 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.7 && a->b[i]->prob_sel_regime <= 0.8) */
/* 	  PhyML_Fprintf(fp,"1. .5 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.8 && a->b[i]->prob_sel_regime <= 0.9) */
/* 	  PhyML_Fprintf(fp,"1. 0. 0. sc\n"); */
/* 	else if(a->b[i]->prob_sel_regime > 0.9) */
/* 	  PhyML_Fprintf(fp,"1. .0 .0 sc\n"); */
	break;
      }

  PhyML_Fprintf(fp,"%d %d mt\n",w->xcoord[a->num],w->ycoord[a->num]);
  PhyML_Fprintf(fp,"%d %d lt\n",w->xcoord[a->num],w->ycoord[d->num]);
  PhyML_Fprintf(fp,"%d %d lt\n",w->xcoord[d->num],w->ycoord[d->num]);

  if(d->tax) 
    {
      PhyML_Fprintf(fp,"(%s) show \n",d->name);
      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"grestore\n");
      return;
    }
  else
    {
      PhyML_Fprintf(fp,"stroke\n");
      PhyML_Fprintf(fp,"grestore\n");
      For(i,3)
	if(d->v[i] != a) DR_Print_Tree_Postscript_Pre(d,d->v[i],fp,w,tree);
    }


  return;
}

/*********************************************************/


void DR_Get_X_Coord_Pre(node *a, node *d, edge *b, tdraw *w, arbre *tree)
{
  int i;

  if(b) w->xcoord[d->num] =  d->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  DR_Get_X_Coord_Pre(d,d->v[i],d->b[i],w,tree);
    }
}

/*********************************************************/

void DR_Get_X_Coord(tdraw *w, arbre *tree)
{
  w->xcoord[tree->n_root->v[0]->num] = tree->n_root->v[0]->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;
  w->xcoord[tree->n_root->v[1]->num] = tree->n_root->v[1]->dist_to_root * (double)w->tree_box_width/w->max_dist_to_root;
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[0],NULL,w,tree);
  DR_Get_X_Coord_Pre(tree->n_root,tree->n_root->v[1],NULL,w,tree);
}


/*********************************************************/

void DR_Get_Y_Coord(tdraw *w, arbre *tree)
{
  int next_y_slot;
  next_y_slot = 0;
  DR_Get_Y_Coord_Post(tree->e_root->left,tree->e_root->rght,tree->e_root,&next_y_slot,w,tree);
  DR_Get_Y_Coord_Post(tree->e_root->rght,tree->e_root->left,tree->e_root,&next_y_slot,w,tree);
}

/*********************************************************/

void DR_Get_Y_Coord_Post(node *a, node *d, edge *b, int *next_y_slot, tdraw *w, arbre *tree)
{
  int i;

  if(d->tax) 
    {
      w->ycoord[d->num] = *next_y_slot + (int)(w->page_height / (2.*tree->n_otu));
      (*next_y_slot) += (int)(w->page_height / (tree->n_otu));
    }
  else
    {
      int d1, d2;

      d1 = d2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      DR_Get_Y_Coord_Post(d,d->v[i],d->b[i],next_y_slot,w,tree);
	      if(d1<0) d1 = i;
	      else     d2 = i;
	    }
	}
      w->ycoord[d->num] = (w->ycoord[d->v[d1]->num] + w->ycoord[d->v[d2]->num])/2.; 
    }
}

/*********************************************************/

tdraw *DR_Make_Tdraw_Struct(arbre *tree)
{
  tdraw *w;

  w = (tdraw *)mCalloc(1,sizeof(tdraw));
  w->xcoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  w->ycoord = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));

  return w;
}

/*********************************************************/

void DR_Init_Tdraw_Struct(tdraw *w)
{
  w->page_width  = 510;
  w->page_height = 770;
}

/*********************************************************/

void DR_Get_Tree_Box_Width(tdraw *w, arbre *tree)
{
  int i;
  int max_name_len, curr_len;

  max_name_len = curr_len = 0;
  For(i,tree->n_otu)
    {
      curr_len = (int)strlen(tree->noeud[i]->name);
      if(curr_len > max_name_len) max_name_len = curr_len;
    }

  w->tree_box_width = w->page_width - max_name_len * 8.66667;
}

/*********************************************************/

double DR_Get_Max_Dist_To_Root(arbre *tree)
{
  double mx;
  int i;

  mx = .0;
  For(i,tree->n_otu)
    {
      if(tree->noeud[i]->dist_to_root > mx)
	{
	  mx = tree->noeud[i]->dist_to_root;
	}
    }

  return mx;
}

/*********************************************************/
