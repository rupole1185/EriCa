#include "eqntns.h"

//Source term depending on cell
void sources(struct CSRMat *Mat, struct RHS *Rhs)
{
   int icell, ivar, idir;
   double *CelldblPnt, *phi, *dep, Src[nvar];
   double tau[ndir][ndir], div, CmuOm;
   int i, j, k;
   double beta, Sij[ndir][ndir], Omij[ndir][ndir], xiw, fb;
   double SAa, SAb, SAc;

   for (icell=0;icell<ncell;icell++) {
      phi = setvar(icell);
      dep = setdep(icell);
      CelldblPnt = setdblcell(icell);

      memset(Src, 0.0, nvar * sizeof(double));

      // Momentum Sources ----------------------------------------------------
      if (eos != ICOUPLED) {
         for (idir=0;idir<ndir;idir++)
            Src[Uvar+idir] = - dep[PGrd+idir];
      }

      // Turbulence production --------------
      if (Kvar!=-1) {
         for (idir=0;idir<ndir;idir++)
            Src[Uvar+idir]-= 2.0 / 3.0 * dep[Rdep] * dep[KGrd+idir];
            //Src[Uvar+idir]-= 2.0 / 3.0 * dep[MTdep] * divergence( dep);

         stress( phi, dep, tau);
         Src[Kvar] = prodk( dep, tau); 
      }

      // turbulence Sources --------------------------------------------------
      switch (turb) {
      //    Spalart Allmaras ----------------
      case SA:
         // NASA #################################
         /*SAa = cb1 * ( 1.0 - ft2( phi, dep)) * Shat( phi, dep) * phi[MTvar];

         SAb = ( cw1 * fw( phi, dep) - cb1 / (kappa * kappa) * ft2( phi, dep)) * \
                    ( phi[MTvar] * phi[MTvar]) / ( dep[Ddep] * dep[Ddep]);

         SAc = cb2 / sigma * scalarp( ndir, &dep[MTGrd], &dep[MTGrd]);*/
         // Deck #################################
         SAa = cb1 * Shat( phi, dep) * phi[MTvar];

         SAb = cw1 * fw( phi, dep) * \
             ( phi[MTvar] * phi[MTvar]) / ( dep[Ddep] * dep[Ddep]);

         SAc = 0.0; //cb2 / sigma * scalarp( ndir, &dep[MTGrd], &dep[MTGrd]);

         Src[MTvar] = dep[Rdep] * (SAa - SAb + SAc);
         break;
      //    Kw Wilcox 2006 ------------------
      case KW06:
         strain( dep, Sij);
         rotation( dep, Omij);
         div = divergence( dep);

         for (i=0;i<ndir;i++)
            Sij[i][i] -= 0.5 * div;

         xiw = 0.0;
         for (i=0;i<ndir;i++) {
         for (j=0;j<ndir;j++) {
         for (k=0;k<ndir;k++)
            xiw += Omij[i][j] * Omij[j][k] * Sij[k][i]; }}

         CmuOm = Cmu * phi[OMvar];
         xiw = fabs( xiw / (CmuOm * CmuOm * CmuOm));

         fb = ( 1.0 + 85.0*xiw) / ( 1.0 + 100.0*xiw);
         beta = beta0 * fb;

         Src[Kvar] -= Cmu   * dep[Rdep]  * phi[Kvar]  * phi[OMvar];
         Src[OMvar] = alpha * phi[OMvar] / phi[Kvar]  * prodk( dep, tau);
         Src[OMvar]-= beta  * dep[Rdep]  * phi[OMvar] * phi[OMvar];
         if ( scalarp( ndir, &dep[KGrd], &dep[OMGrd]) > 0.0)
            Src[OMvar]+= 0.125 * dep[Rdep] / phi[OMvar] * \
                                 scalarp( ndir, &dep[KGrd], &dep[OMGrd]);
         break;
      //    Menter KwSST --------------------
      case KWSST:
         strain( dep, Sij);
         rotation( dep, Omij);

         alpha = blend( F1(phi,dep), gamma1, gamma2);
         beta  = blend( F1(phi,dep), beta1, beta2);

         Src[Kvar] -= Cmu   * dep[Rdep]  * phi[Kvar]  * phi[OMvar];
         Src[OMvar] = alpha * dep[Rdep]  / dep[MTdep] * prodk( dep, tau);
         Src[OMvar]-= beta  * dep[Rdep]  * phi[OMvar] * phi[OMvar];
         Src[OMvar]+= 2.0 * (1.0-F1(phi,dep)) * sigw2 * dep[Rdep] / phi[OMvar] * \
                                 scalarp( ndir, &dep[KGrd], &dep[OMGrd]);
         break;
      }

      getrhs( Rhs, icell, CelldblPnt[CellVol], Src);
   }
}
