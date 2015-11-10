#include "boundary.h"

static char* MonFile = "OUTPUT/monitors.dat";

void moninit( void)
{
   FILE *MonPnt;
   int imon;
   char buffer[3];

   if (nmon>0) {

      MonPnt = fopen( MonFile, "w");

      fprintf(MonPnt, "##############################\n");
      fprintf(MonPnt, "#  EriCa Wall Monitors file  #\n");
      fprintf(MonPnt, "##############################\n");

      for (imon=0;imon<nmon;imon++) {
         int isurf= monsurf[imon];
         int ivar = monphi[imon];
         int iave = monave[imon];
         int inorm= monorm[imon];

         fprintf(MonPnt," BND = %03d                   |", isurf);
      }
      fprintf(MonPnt,"\n");

      for (imon=0;imon<nmon;imon++) {
         int isurf= monsurf[imon];
         int ivar = monphi[imon];
         int iare = monare[imon];
         int iave = monave[imon];
         int inorm= monorm[imon];

         switch (iare) {
         case -1:
            strcpy( buffer, "AA");
            break;
         case 0:
            strcpy( buffer, "Ax");
            break;
         case 1:
            strcpy( buffer, "Ay");
            break;
         case 2:
            strcpy( buffer, "Az");
            break;
         default:
            FatalError("Wrong iare in monitor");
            break;
         }

         switch (iave) {
         case 0:
            fprintf(MonPnt," SUM( %s %s)        |", buffer, VarName( ivar));
            break;
         case 1:
            if (ivar>=0)
               fprintf(MonPnt," SUM( rho U %s %s)  |", buffer, VarName( ivar));
            else
               fprintf(MonPnt," SUM( rho U %s )    |", buffer);
            break;
         case 2:
            fprintf(MonPnt," SUM( mu %s D%s/Dn) |", buffer, VarName( ivar));
            break;
         default:
            FatalError("Wrong iave in monitor");
            break;
         }
      }
      fprintf(MonPnt,"\n");

      for (imon=0;imon<nmon;imon++) {
         int isurf= monsurf[imon];
         int ivar = monphi[imon];
         int iare = monare[imon];
         int iave = monave[imon];
         int inorm= monorm[imon];

         switch (inorm) {
         case 0:
            fprintf(MonPnt,"                             |");
            break;
         case 1:
            fprintf(MonPnt," SUM(  %s )                  |", buffer);
            break;
         case 2:
            fprintf(MonPnt," SUM( rho U %s)              |", buffer);
            break;
         default:
            FatalError("Wrong inorm in monitor");
            break;
         }
      }
      fprintf(MonPnt,"\n");

      fclose( MonPnt);
   }
}

void monitors( void)
{
   int imon, ighst, *GhstintPnt, icell, idir;
   double *GhstdblPnt, mu_i[nvar], areatot, AA;
   double alfa, *Ghstphi, *Ghstdep, *phi, *dep;
   double coefftg[ndir], Du_tauDn, coeff, phivar;
   FILE *MonPnt;

   if (nmon==0)
      return;

   MonPnt = fopen( MonFile, "a");

   for (imon=0;imon<nmon;imon++) {
      int ivar = monphi[imon];
      int iare = monare[imon];
      int iave = monave[imon];
      int inorm= monorm[imon];
      mon[imon] = 0.0;

      double norm = 0.0;

      for (ighst=0;ighst<nghst;ighst++) {
         GhstdblPnt = setdblghst(ighst);
         GhstintPnt = setintghst(ighst);
         Ghstphi    = setvar( ghst2gbl( ighst));
         Ghstdep    = setdep( ghst2gbl( ighst));
         phi        = setvar( GhstintPnt[GhstFidx]);
         dep        = setdep( GhstintPnt[GhstFidx]);

         if (GhstintPnt[GhstRef] == monsurf[imon]) {

            AA = sqrt( scalarp( ndir, &GhstdblPnt[GhstSrf], &GhstdblPnt[GhstSrf]));

            for (idir=0;idir<ndir;idir++)
               coefftg[idir] = 1.0 - GhstdblPnt[GhstSrf+idir] / AA; 

            setmuvar( mu_i, GhstdblPnt[GhstAlfa], dep, Ghstdep, phi, Ghstphi);

            switch (iare) {
               case -1:
                  AA = sqrt( scalarp( ndir, &GhstdblPnt[GhstSrf], &GhstdblPnt[GhstSrf]));
                  break;
               case 0:
               case 1:
               case 2:
                  AA = GhstdblPnt[GhstSrf+iare];
                  break;
            }

            switch (iave) {
            case 0:
               mon[imon] += AA * GhstdblPnt[GhstWallPhi+ivar];
               break;
            case 1:
               if (ivar>=0)
                  phivar = GhstdblPnt[GhstWallPhi+ivar];
               else
                  phivar = 1;

               mon[imon] += dep[Rdep] * scalarp( ndir, &GhstdblPnt[GhstSrf], &GhstdblPnt[GhstWallPhi+Uvar]) * phivar; 
               break;
            case 2:
               coeff = 1.0;

               if (ivar>=Uvar && ivar<Uvar+ndir)
                  coeff = coefftg[ivar-Uvar];

               mon[imon] += mu_i[ivar] * AA * GhstdblPnt[GhstWallDPDn+ivar] * coeff;
               break;
            }

            switch (inorm) {
            case 0:
               break;
            case 1:
               norm += AA;
               break;
            case 2:
               mon[imon] += dep[Rdep] * scalarp( ndir, &GhstdblPnt[GhstSrf], &GhstdblPnt[GhstWallPhi+Uvar]); 
               break;
            }
         }
      }

      if (inorm!=0)
         mon[imon] = mon[imon] / norm;

      fprintf( MonPnt, "       %21.8f |", mon[imon]);
   }

   fprintf( MonPnt, "\n");
   fclose( MonPnt);
}
