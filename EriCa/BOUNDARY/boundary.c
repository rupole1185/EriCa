#include "boundary.h"

void boundary( void)
{
   int ighst, *GhstintPnt, icell, ivar, idir;
   double *GhstdblPnt, dx[ndir], alfa, *phi, *GhstCelldblPnt;
   double  phi_m[nvar], *CelldblPnt, *Ghstphi, dist;
   double *Ghstdep, *dep;

   for (ighst=0;ighst<nghst;ighst++) {
      GhstdblPnt = setdblghst(ighst);
      GhstintPnt = setintghst(ighst);
      GhstCelldblPnt = setdblcell( ghst2gbl( ighst));
      Ghstphi        = setvar( ghst2gbl(ighst));
      Ghstdep        = setdep( ghst2gbl(ighst));

      icell = GhstintPnt[GhstFidx];
      phi   = setvar(icell);
      dep   = setdep(icell);
      CelldblPnt = setdblcell( icell);

      alfa   = GhstdblPnt[GhstAlfa];

      // Mirror point for the FLD cell ---------
      for (ivar=0;ivar<nvar;ivar++)
         phi_m[ivar] = phi[ivar] + scalarp( ndir, &dep[ PGrd+ndir*ivar], \
                                                  &GhstdblPnt[GhstDx]);

      // Calculating Wall Variable -------------
      for (ivar=0;ivar<nvar;ivar++) {
         switch (GhstintPnt[GhstType+ivar]) {
         case -1:
         case 0:
            GhstdblPnt[GhstWallPhi+ivar] = GhstdblPnt[GhstBnd+ivar];
            break;
         case 1:
            GhstdblPnt[GhstWallPhi+ivar] = phi_m[ivar] - \
                       GhstdblPnt[GhstBnd+ivar] * GhstdblPnt[GhstDxAlfa];
            break;
         default:
            FatalError("Wrong BND condition selected!");
            break;
         }
      }

      // Application of Dirichlet BC -----------
      for (ivar=0;ivar<nvar;ivar++)
         Ghstphi[ivar] = GhstdblPnt[GhstWallPhi+ivar] / (1.0 - alfa) + phi_m[ivar] * alfa / (alfa - 1.0);

      // Flux Calculation ----------------------
      for (ivar=0;ivar<nvar;ivar++)     // Wall normal convention: positive if towards fluid
         GhstdblPnt[GhstWallDPDn+ivar] = (phi_m[ivar] - Ghstphi[ivar]) / GhstdblPnt[GhstDist];
   }
}

//Subroutine to know the relation between the ghost cell and the fluid one
double ghst2fld(int ighst, int ivar)
{
   int *GhstintPnt;
   double ghst2fld;

#if _DBG == 10
   if (ighst >= nghst)
      FatalError("Error in Ghst2Fld");
#endif

   GhstintPnt = setintghst(ighst);
   ghst2fld  = 0.0;
  
   switch (GhstintPnt[GhstType+ivar]) {
   case -1:
   case  0:
      ghst2fld =  0.0;
      break;
   case  1:
      ghst2fld = 1.0;
      break;
   default:
      FatalError("Wrong BND int Ghst2Fld!");
      break;
   }

   return ghst2fld;
}

//  Function that returns local ghost idx
int gbl2ghst( int icell)
{
   return icell - ncell;
}

int ghst2gbl( int ighst)
{
   return ighst + ncell;
}

// Is it wall function
int isitwall( int ighst)
{
   int *GhstintPnt;

   GhstintPnt = setintghst(ighst);

   if (GhstintPnt[GhstType+Uvar] == -1)
      return 1;
   else
      return 0;
}
