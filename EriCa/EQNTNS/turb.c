#include "turb.h"

// GENERAL TURBULENCE FUNCTIONALITIES ###############################
// Definition of constants dependent from primary other constants ...
void init_turb_const( void)
{
   gamma1 = beta1 / Cmu - sigw1 * kappa * kappa /sqrt(Cmu);
   gamma2 = beta2 / Cmu - sigw2 * kappa * kappa /sqrt(Cmu);
   cw1 = cb1 /  kappa * kappa + ( 1.0 + cb2) / sigma;
}

// Stress tensor Tij (turbulent flow)
void stress( double phi[nvar], double dep[ndep], double tij[ndir][ndir])
{
   int idir, iidir;
   double Sij[ndir][ndir], div;

   switch (turb) {
   case KWSST:
   case KW06:
      strain( dep, Sij);
      div = divergence( dep);

      for (idir=0;idir<ndir;idir++) {
         for (iidir=0;iidir<ndir;iidir++) {
            tij[idir][iidir] = 2.0 * dep[MTdep] * Sij[idir][iidir];
         }

         tij[idir][idir] -= 2.0 / 3.0 * ( dep[MTdep] * div +\
                                          dep[Rdep] * phi[Kvar]);
      }
      break;
   default:
      FatalError("No TAUIJ method defined for this turbulence");
      break;
   }
}

// Velocity strain tensor
void strain( double dep[ndep], double Sij[ndir][ndir])
{
   int idir, iidir;

   for (idir=0;idir<ndir;idir++) {
      for (iidir=0;iidir<ndir;iidir++) {
         Sij[idir][iidir] = 0.5 * (dep[UGrd+idir*ndir+iidir] + \
                                   dep[UGrd+iidir*ndir+idir]);
      }
   }
}

// Rotation strain tensor (or vorticity)
void rotation( double dep[ndep], double Omij[ndir][ndir])
{
   int idir, iidir;

   for (idir=0;idir<ndir;idir++) {
      for (iidir=0;iidir<ndir;iidir++) {
         Omij[idir][iidir] = 0.5 * (dep[UGrd+idir*ndir+iidir] - \
                                    dep[UGrd+iidir*ndir+idir]);
      }
   }
}

// TKE production term
double prodk( double dep[ndep], double tij[ndir][ndir])
{
   int idir, iidir;
   double Pk = 0.0;

   for (idir=0;idir<ndir;idir++) {
   for (iidir=0;iidir<ndir;iidir++)
      Pk += tij[idir][iidir] * dep[UGrd+idir*ndir+iidir];
   }

   return Pk;
}

double matsum( double mat[ndir][ndir])
{
   int idir, iidir;
   double ss = 0.0;

   for (idir=0;idir<ndir;idir++) {
   for (iidir=0;iidir<ndir;iidir++)
      ss += mat[idir][iidir] * mat[idir][iidir];
   }

   return ss;
}

// KW SST peculiar ones #############################################
double F1( double phi[nvar], double dep[ndep])
{
   double arg1;

   arg1 = fmin( fmax( sqrt(phi[Kvar]) / (Cmu * phi[OMvar] * dep[Ddep]), \
                      500.0 * dep[MUdep] / ( dep[Ddep] * dep[Ddep] * phi[OMvar] )),
                4.0 * dep[Rdep] * sigw2 * phi[Kvar] / ( CDkw( phi, dep) * dep[Ddep] * dep[Ddep]) );

   return tanh( arg1 * arg1 * arg1 * arg1);
}

double CDkw( double phi[nvar], double dep[ndep])
{
   double CD;

   CD = fmax( 2.0 * dep[Rdep] * sigw2 / phi[OMvar] * \
               scalarp( ndir, &dep[KGrd], &dep[OMGrd]), 1.E-20);

   return CD;
}

double F2( double phi[nvar], double dep[ndep])
{
   double arg2;

   arg2 = fmax( 2.0 * sqrt(phi[Kvar]) / (Cmu * phi[OMvar] * dep[Ddep]), \
                500.0 * dep[MUdep] / ( dep[Ddep] * dep[Ddep] * phi[OMvar] ));

   return tanh( arg2 * arg2 );
}

double blend( double alfa, double vara, double varb)
{
   return (alfa * vara + (1.0 - alfa) * varb);
}

// SA peculiar ones #################################################
double fw( double phi[nvar], double dep[ndep])
{
   double r = fmin( phi[MTvar] / ( Shat( phi, dep) * kappa * kappa * dep[Ddep] * dep[Ddep]), 10.0);

   double g = r + cw2 * ( r*r*r*r*r*r - r);

   double cw3_6 = cw3*cw3*cw3*cw3*cw3*cw3;

   return ( g * pow( (1.0+cw3_6)/(g*g*g*g*g*g + cw3_6), 1.0/6.0));
}

double Shat( double phi[nvar], double dep[ndep])
{
   double Omij[ndir][ndir];
   rotation( dep, Omij);

   double Omega = sqrt( 2.0 * matsum( Omij));

   double out = Omega + phi[MTvar] / (kappa*kappa*dep[Ddep]*dep[Ddep]) * fv2( phi, dep);

   return fmax( out, 0.3 * Omega);
}

double fv1( double phi[nvar], double dep[ndep])
{
   double xi = phi[MTvar] / dep[MUdep];

   return ( xi * xi * xi / ( xi * xi * xi + cv1 * cv1 * cv1));
}

double fv2( double phi[nvar], double dep[ndep])
{
   double xi = phi[MTvar] / dep[MUdep];

   return ( 1.0 - xi / ( 1.0 + xi * fv1( phi, dep)));
}

double ft2( double phi[nvar], double dep[ndep])
{
   double xi = phi[MTvar] / dep[MUdep];

   return ( ct3 * exp( - ct4 * xi * xi));
}
