#include "mesh.h"
#include "meshset.h"

//Set Face Pointers -----------------------------------
int* setintface(int iface)
{
   return &faceintdata[iface*nfaceint];
}

double* setdblface(int iface)
{
   return &facedbldata[iface*nfacedbl];
}

//Set Cell Pointers----------------------------------
int* setintcell(int icell)
{
   return &cellintdata[icell*ncellint];
}

double* setdblcell(int icell)
{
   return &celldbldata[icell*ncelldbl];
}

//Set Ghost Pointers----------------------------------
int* setintghst(int ighost)
{
   return &ghstintdata[ighost*nghstint];
}

double* setdblghst(int ighost)
{
   return &ghstdbldata[ighost*nghstdbl];
}

//Set variables Pointers--------------------------------
double* setvar(int icell)
{
   return &var[icell*nvar];
}

//Set dependent Pointers
double* setdep(int icell)
{
   return &dep[icell*ndep];
}

// Other pointers for generality ----------------------
double* setsol(int icell, int nnvar, double *pntr)
{
   return &pntr[icell*nnvar];
}
