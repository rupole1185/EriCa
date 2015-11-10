#include "eqntns.h"

//Treatment of "pressure" equation fluxes into convection
void presflx(double phi[nvar], double F_i[nvar][nvar])
{
   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IPISO:
   case ICOUPLED:
   case IFROZENP:
      simple_conf( phi, F_i);
      break;
   case IDBIG:
      dnsbsd_conf( phi, F_i);
      break;
   }
}

//Time advancing matrix
void timematrix(double phi[nvar], double dep[ndep], double F_i[nvar][nvar])
{
   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IPISO:
   case ICOUPLED:
   case IFROZENP:
      simple_tmmt( phi, dep, F_i);
      break;
   case IDBIG:
      dnsbsd_tmmt( phi, dep, F_i);
      break;
   }
}

//Operations to perform before solving the system
void pv_adjustment(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   int icell;
   double *phi;

   // Copying VAR to SOL
   phi = setvar( 0); 
   memcpy( sol, phi, ncell*nvar* sizeof(double));

   // Modifications due to alghorithms 
   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IPISO:
   case IFROZENP:
      for (icell=0;icell<ncell;icell++) 
         sol[icell*nvar+Pvar] = 0.0;
      break;
   case ICOUPLED:
      coupled_press( Mat, Rhs);
      break;
   case IDBIG:
      break;
   }
}

//Time update subroutines
void update(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IFROZENP:
      simple_upd( Mat, Rhs, sol);
      break;
   case IPISO:
      piso_upd( Mat, Rhs, sol);
   case ICOUPLED:
      coupled_upd( Mat, Rhs, sol);
      break;
   case IDBIG:
      dnsbsd_upd( Mat, Rhs, sol);
      break;
   }
}

//Return Matrix Pattern
void pattern( int pattern[nvar][nvar])
{
   int ivar, iivar;

   for (ivar=0;ivar<nvar;ivar++)
      for (iivar=0;iivar<nvar;iivar++)
         pattern[ivar][iivar] = 1;

   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IPISO:
   case IFROZENP:
      for (ivar=0;ivar<nvar;ivar++)
         pattern[ivar][ivar] = 0;
      break;
   case ICOUPLED:
      for (ivar=0;ivar<Uvar+ndir;ivar++)
         for (iivar=0;iivar<Uvar+ndir;iivar++)
            pattern[ivar][iivar] = 0;

      for (ivar=Uvar+ndir;ivar<nvar;ivar++)
         pattern[ivar][ivar] = 0;

      break;
   case IDBIG:
      for (ivar=0;ivar<=Tvar;ivar++)
         for (iivar=0;iivar<=Tvar;iivar++)
            pattern[ivar][iivar] = 0;

      for (ivar=Tvar+1;ivar<nvar;ivar++)
         pattern[ivar][ivar] = 0;

      break;
   }
}
