/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef INTERFACE_H
#define INTERFACE_H

void Launch_Interface(option *input);
void Clear();
void Launch_Interface_Input(option *input);
void Launch_Interface_Data_Type(option *input);
void Launch_Interface_Model(option *input);
void Launch_Interface_Topo_Search(option *input);
void Launch_Interface_Branch_Support(option *input);
void Launch_Interface_Multigene(option *input);

#endif
