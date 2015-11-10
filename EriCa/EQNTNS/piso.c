#include "eqntns.h"

void piso_upd(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   int icell, ivar, istart, idir;
   double *phi, *Vprime;
/*
   // Velocity coefficients -------------
   calcdd( Mat);

   // Pressure equation resolution ------
   simple_peqn( &p_sol, sol, NULL, setVELsol);

   // Velocity correction ---------------
   simple_uprm( &uprime, p_sol, NULL);

   for (icell=0;icell<ncell;icell++) {
      Vprime = setVELndir( icell, uprime);
      phi    = setvar(icell);

      for (idir=0;idir<ndir;idir++)
         Vprime[idir] += phi[Uvar+idir];
   }

   // Pressure p'' ----------------------
   simple_peqn( &pp_sol, uprime, sol, setVELndir);

   // Velocity v** ----------------------
   simple_uprm( &upprime, pp_sol);

   // Variables update ------------------
   for (icell=0;icell<ncell;icell++)  {
      istart = icell * nvar;
      phi    = setvar(icell);
  
      //Pressure update P = P + p'
      phi[Pvar] = phi[Pvar] + p_sol[icell] + pp_sol[icell];

      //Velocity update U = u* + u'
      for (idir=0;idir<ndir;idir++)
         phi[Uvar+idir] = uprime[icell*ndir+idir] + \
                          upprime[icell*ndir+idir];

      //Other variables are just the solution
      for (ivar=Uvar+ndir;ivar<nvar;ivar++)
         //si puo' usare um memcpy?
         phi[ivar] = sol[istart+ivar];
   }

   // Memory release --------------------
   free(uprime);
   free(upprime);
   free(p_sol);
   free(pp_sol);
   free(dd); */
}
