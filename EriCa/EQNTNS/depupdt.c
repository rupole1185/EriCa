#include "eqntns.h"

void depupdt(void )
{
   calcrho();
   calccfl();

   if (MTdep != -1)
      calcmt();

#if _DBG == 10
   calcdbgdep();
#endif
}

// Cells density update ---------
void calcrho(void )
{
   int icell;
   double *phi, *dep;

   switch (eos) {
   case ISIMPLE: 
   case ISIMPLEC: 
   case IPISO:
   case ICOUPLED: 
   case IFROZENP:
      break;
   case IDBIG:
      for (icell=0;icell<ntot;icell++) {
         phi = setvar(icell);
         dep = setdep(icell);

         dep[Rdep] = phi[Pvar] / (R * phi[Tvar]);
      }
      break;
   default:
      FatalError("Wrong Equation of state");
      break;
   }
}

// density and enthalpy defs ----
double drdp( double phi[nvar])
{
   switch (eos) {
   case IDBIG:
      return ( 1.0 / (R * phi[Tvar]));
      break;
   default:
      FatalError("No drdp mode defined");
      break;
   }
}

double drdt( double phi[nvar])
{
   switch (eos) {
   case IDBIG:
      return (- 1.0 * ( phi[Pvar] / (R * phi[Tvar] * phi[Tvar])));
      break;
   default:
      FatalError("No drdt mode defined");
      break;
   }
}

double enthalpy( double phi[nvar])
{
   double hh;

   switch (eos) {
   case IDBIG:
      hh = Cp * phi[Tvar] + 0.5 * scalarp( ndir, &phi[Uvar], &phi[Uvar]);
      break;
   default:
      FatalError("No enthalpy mode defined");
      break;
   }

   return hh;
}

double dhdp( double phi[nvar])
{
   switch (eos) {
   case IDBIG:
      return 0.0;
      break;
   default:
      FatalError("No dhdp mode defined");
      break;
   }
}

double dhdt( double phi[nvar])
{
   switch (eos) {
   case IDBIG:
      return ( Cp);
      break;
   default:
      FatalError("No dhdt mode defined");
      break;
   }
}

void calccfl(void )
{
   int icell, idir;
   double *phi, *dep, *CelldblPnt, Vabs;

   if (cfl > 0.0) { // TimeStep calculation for STEADY flows
      timst = 1.0E10;

      for (icell=0;icell<ncell;icell++) {
         phi         = setvar(icell);
         dep         = setdep(icell);
         CelldblPnt  = setdblcell(icell);

         Vabs        = sqrt(scalarp( ndir, &phi[Uvar], &phi[Uvar]));

         if (Vabs < 0.01)
            Vabs = 0.01;

         timst = fmin( timst, cfl / Vabs * pow(CelldblPnt[CellVol], 1.0/ndir));
      }
   }

   for (icell=0;icell<ncell;icell++) {
      phi        = setvar(icell);
      dep        = setdep(icell);
      CelldblPnt = setdblcell(icell);

      Vabs      = sqrt(scalarp( ndir, &phi[Uvar], &phi[Uvar]));
      dep[CFLdep] = Vabs * timst / pow(CelldblPnt[CellVol],1.0/ndir);
      dep[TIMEdep]= timst;
   }
}

// Turbulent viscosity ----------
void calcmt(void )
{
   int icell, i, j;
   double *phi, *dep, omhat, om2, Sij[ndir][ndir], div, WW[ndir][ndir];
   double xi;

   switch (turb) {
   case SA:
      for (icell=0;icell<ntot;icell++) {
         phi = setvar(icell);
         dep = setdep(icell);

         dep[MTdep] = dep[Rdep] * phi[MTvar] * fv1( phi, dep);
      }
      break;
   case KW06:
      for (icell=0;icell<ntot;icell++) {
         phi = setvar(icell);
         dep = setdep(icell);

         strain( dep, Sij);
         div = divergence( dep);

         for (i=0;i<ndir;i++)
            Sij[i][i] -= div / 3.0;

         om2 = 2.0 * matsum( Sij);

         omhat = fmax( phi[OMvar], Clim * sqrt(om2 / Cmu));

         dep[MTdep] = dep[Rdep] * phi[Kvar] / omhat;
      }
      break;
   case KWSST:
      for (icell=0;icell<ntot;icell++) {
         phi = setvar(icell);
         dep = setdep(icell);

         rotation( dep, WW);

         dep[MTdep] = dep[Rdep] * a1 * phi[Kvar];
         dep[MTdep]/= fmax(a1 * phi[OMvar], \
                           sqrt(2.0*matsum(WW)) * F2( phi, dep));
      }
      break;
   }
}

// Velocity divergence ----------
double divergence( double dep[ndep])
{
   int idir;
   double div = 0.0;

   for (idir=0;idir<ndir;idir++)
      div += dep[UGrd+(ndir*idir)+idir];

   return div;
}

#if _DBG == 10
void calcdbgdep(void )
{
   int icell, idir;
   double *dep, *phi, tij[ndir][ndir];

   for (icell=0;icell<ncell;icell++) {
      dep       = setdep(icell);
      phi       = setvar(icell);

      // Divergence ----------------
      dep[Divdep] = divergence( dep);

      // Vorticity -----------------
      double Omij[ndir][ndir];
      rotation( dep, Omij);
      dep[Vrtdep]= 2.0 * matsum( Omij);

      // Pk ------------------------
      if (Kvar!=-1) {
         stress( phi, dep, tij);
         dep[Pkdep]  = prodk( dep, tij);
      }
   }
}
#endif
