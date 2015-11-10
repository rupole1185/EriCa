#include "eqntns.h"

// Density Equation convective fluxes
void dnsbsd_conf(double phi[nvar], double F_i[nvar][nvar])
{
   phi[Pvar]       = 1.0;
}

void dnsbsd_tmmt( double phi[nvar], double dep[ndep], double F_i[nvar][nvar])
{
   int ivar, idir, iidir;
   
   // density jacobian ----------------
   F_i[Pvar][Pvar] = drdp( phi);
   F_i[Pvar][Tvar] = drdt( phi);

//printf(" %f %f\n", F_i[Pvar][Pvar], F_i[Pvar][Tvar]);

   for (idir=0;idir<ndir;idir++) {
   // momentum jacobian ---------------
      F_i[Uvar+idir][Pvar]      = drdp( phi) * phi[Uvar+idir];
      F_i[Uvar+idir][Uvar+idir] = dep[Rdep];
      F_i[Uvar+idir][Tvar]      = drdt( phi) * phi[Uvar+idir];

   // enthalpy jacobian/1 -------------
      F_i[Tvar][Uvar+idir] = dep[Rdep] * phi[Uvar+idir];
   }

   // enthalpy jacobian/2 -------------
   F_i[Tvar][Pvar] = drdp( phi) - ( 1.0 * dep[Rdep] * dhdp( phi));
   F_i[Tvar][Tvar] = drdt( phi) * enthalpy( phi) + dep[Rdep] * dhdt( phi);

   for (ivar=Tvar+1;ivar<nvar;ivar++)
      F_i[ivar][ivar] = dep[Rdep];
}

void dnsbsd_upd(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   int icell, istart, ivar;
   double *phi, *phi_new;

   for (icell=0;icell<ncell;icell++) {
      phi     = setvar( icell);
      phi_new = setsol( icell, nvar, sol);

      for (ivar=0;ivar<nvar;ivar++)
         phi[ivar] = phi_new[ivar];
   }
}
