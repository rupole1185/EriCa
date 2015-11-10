#include "eqntns.h"

// Distance between two points (vector form) -------------
void pntsdist( double xxa[ndir], double xxb[ndir], double dx[ndir])
{
   int idir;

   for (idir=0;idir<ndir;idir++)
      dx[idir] = xxa[idir] - xxb[idir];
}

//Set facevar -----------------------------------------
double facevar(double alfa, double var_a, double var_b)
{
   return blend( alfa, var_a, var_b);
}
