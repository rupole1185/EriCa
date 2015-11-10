#include "../MESH/mesh.h"

void boundary( void);

double ghst2fld(int ighst, int ivar);

// Function to move from global to ghost and vice-versa
int glb2ghst( int icell);
int ghst2gbl( int ighost);

//Tools
int isitwall( int ighst);

void monitors( void);
