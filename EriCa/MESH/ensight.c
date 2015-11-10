#include "mesh.h"
#include "meshset.h"
#include "../EQNTNS/turb.h"

static double *timevalues;
double *nodexyz;
int    *c2n;
static FILE   *CasePnt;
int static iplot; //Initialized to zero authomatically
int static ElmTypCSR[6];

// Storing the number of elements divided by each type in CSR format
//   0 - tetra4
//   1 - pyra5
//   2 - penta6
//   3 - hexa8
//   4 - polyN
void ElmTypCSRstart()
{
   int icell;

   for (icell=0;icell<ncell;icell++) {
      if (nodadj[icell]<=4)
         ElmTypCSR[1]++;
      if (nodadj[icell]<=5)
         ElmTypCSR[2]++;
      if (nodadj[icell]<=6)
         ElmTypCSR[3]++;
      if (nodadj[icell]<=8)
         ElmTypCSR[4]++;
      ElmTypCSR[5]++;
   }
}

void printCEI( int itmst, int ntmst)
{
   if ( itmst % nplot != 0 && itmst != ntmst)
      return;

   fprintf(logfile, "\n --> EnSight plot\n\n");

   // Creation of CaseFile
   CasePnt  = fopen("ENSIGHT/eriCa.ENS.case","w");
   fprintf(CasePnt , "FORMAT\n");
   fprintf(CasePnt , "type: ensight gold\n\n");
   fprintf(CasePnt , "GEOMETRY\n");
   fprintf(CasePnt , "model: MESH.geo\n\n");
   fprintf(CasePnt , "VARIABLE\n");

   timevalues = realloc(timevalues, (iplot+1) * sizeof(double));

   if (itmst == 0)
      timevalues[iplot] = 0.0;
   else
      timevalues[iplot] = timst + timevalues[iplot-1];

   // Variable plots --------------------------------------------
   ensvar( iplot, "Pressure", 1, Pvar, 0, NULL);
   ensvar( iplot, "Velocity", ndir, Uvar, 0, NULL);
   ensvar( iplot, "Temperature", 1, Tvar, 0, NULL);
   ensvar( iplot, "TKE", 1, Kvar, 0, NULL);
   ensvar( iplot, "Omega", 1, OMvar, 0, NULL);

   char ScalName[9];
   int ivar;
   for (ivar=SCAvar;ivar<SCAvar+NSCAvar;ivar++) {
      sprintf(ScalName, "Scalar%2.2d", ivar - SCAvar);
      ensvar( iplot, ScalName, 1, ivar, 0, NULL);
   }

   // Dependants plots ------------------------------------------
   ensvar( iplot, "Density", 1, Rdep, 1, NULL);

   ensvar( iplot, "TimeStep", 1, TIMEdep , 1, NULL);
   ensvar( iplot, "Courant", 1, CFLdep , 1, NULL);

   ensvar( iplot, "Distance", 1, Ddep, 1, NULL);

   ensvar( iplot, "LamVisc", 1, MUdep, 1, NULL);
   ensvar( iplot, "TurbVisc", 1, MTdep, 1, NULL);

#if _DBG == 10
   // Debugging mode: plotting gradients ------------------------

   ensvar( iplot, "GradP", ndir, PGrd, 1, NULL);
   ensvar( iplot, "GradU", ndir, UGrd, 1, NULL);

   if (ndir>1)
      ensvar( iplot, "GradV", ndir, UGrd+1*ndir, 1, NULL);
   
   if (ndir>2)
      ensvar( iplot, "GradW", ndir, UGrd+2*ndir, 1, NULL);

   for (ivar=SCAvar;ivar<SCAvar+NSCAvar;ivar++) {
      sprintf(ScalName, "GrdScl%2.2d", ivar - SCAvar);
      ensvar( iplot, ScalName, ndir, SCAGrd, 1, NULL);
   }

   // PV-Coupling debugging -------------
   ensvar( iplot, "Pprime", 1, Ppdep, 1, NULL);
   ensvar( iplot, "Uprime", ndir, Updep, 1, NULL);
   ensvar( iplot, "Divergence", 1, Divdep, 1, NULL);
   ensvar( iplot, "Div-SIMPLE", 1, DivSdep, 1, NULL);
   ensvar( iplot, "DD-SIMPLE", ndir, DDdep, 1, NULL);

   // Turbulence debugging --------------
   ensvar( iplot, "Vorticity", 1, Vrtdep, 1, NULL);

   ensvar( iplot, "Prodk", 1, Pkdep, 1, NULL);
   ensvar( iplot, "GradTKE", ndir, KGrd, 1, NULL);
   ensvar( iplot, "GradOmega", ndir, OMGrd, 1, NULL);

   ensvar( iplot, "Nu-tilde", 1, MTvar, 0, NULL);
   ensvar( iplot, "GrdNu-tilde", ndir, MTGrd, 1, NULL);

   if (MTvar!=-1) {
      ensvar( iplot,"fw", 1, MTvar, 2, fw);
      ensvar( iplot,"Ft2", 1, MTvar, 2, ft2);
      ensvar( iplot,"Fv1", 1, MTvar, 2, fv1);
      ensvar( iplot,"Fv2", 1, MTvar, 2, fv2);
      ensvar( iplot,"Shat", 1, MTvar, 2, Shat);
   }

   ensvar( iplot, "PHIdist", 1, PHIDdep, 1, NULL);
#endif

   fprintf(CasePnt,"\nTIME\n");
   fprintf(CasePnt,"time set:              1\n");
   fprintf(CasePnt,"number of steps:       %d\n", iplot+1);
   fprintf(CasePnt,"filename start number: 0\n");
   fprintf(CasePnt,"filename increment:    1\n");
   fprintf(CasePnt,"time values:  ");

   int itmval;
   for (itmval=0;itmval<=iplot;itmval++)
      fprintf(CasePnt," %lf\n", timevalues[itmval]);

   fclose(CasePnt);
   iplot++;
}

//            Counter    String length + String           
//             dimens    variable   CharDep
void ensvar(int count, char *String,          \
            int ndim, int ivar, int vd,       \
            double (*funcp)(double *, double *))
{
   int icell, ipart, idim;
   double *phi, *dep;

   int StrLng = StringLength( String);

   char FName[23+StrLng+2], ENScase[48+StrLng+2], Suffix[5+2];

   if (ivar==-1)
      return;

   ipart = 0;

   FILE  *writeENS;
   strcpy(FName,"ENSIGHT/EriCa.ENS.");
   strcat(FName,String);
   sprintf(Suffix, ".%04d", count);
   strcat(FName,Suffix);

   writeENS = fopen(FName,"wb");
   char buffer[81];

   // File HEADER -----------------------------------------
   strcpy( buffer, "no description");
   fwrite( buffer, 80 * sizeof(char), 1, writeENS);

   strcpy( buffer, "part");
   fwrite( buffer, 80 * sizeof(char), 1, writeENS);

   ipart++;
   fwrite( &ipart, sizeof(int), 1, writeENS);

   int itype;
   for (itype = 0; itype < 5; itype++) {

      int istart, iend;
      istart = ElmTypCSR[itype];
      iend   = ElmTypCSR[itype+1];

      if (iend == istart)
         continue;

      switch (itype) {
      case 0:
         strcpy( buffer, "tetra4");
         break;
      case 1:
         strcpy( buffer, "pyramid5");
         break;
      case 2:
         strcpy( buffer, "penta6");
         break;
      case 3:
         strcpy( buffer, "hexa8");
         break;
      case 4:
         strcpy( buffer, "nfaced");
         break;
      }
      fwrite( buffer, 80 * sizeof(char), 1, writeENS);

      switch (vd) {
      case 0: // Variables
         for (idim=0;idim<ndim;idim++) {
            for (icell=istart;icell<iend;icell++) {
               phi = setvar( icell);
               float phivar = phi[ivar+idim];
               fwrite( &phivar, sizeof(float), 1, writeENS);
            }
         }
         break;
      case 1: // Dependents
         for (idim=0;idim<ndim;idim++) {
            for (icell=istart;icell<iend;icell++) {
               dep = setdep( icell);
               float phivar = dep[ivar+idim];
               fwrite( &phivar, sizeof(float), 1, writeENS);
            }
         }
         break;
      case 2: // VirtualFunction calculation
         if (ndim!=1 || funcp == NULL)
            FatalError("Virtual Function error in ENSVAR");

         for (icell=istart;icell<iend;icell++) {
            phi = setvar( icell);
            dep = setdep( icell);
            float phitmp = (*funcp)( phi, dep);
            float phivar = phitmp;
            fwrite( &phivar, sizeof(float), 1, writeENS);
         }
         break;
      }
      

      if (ndir == 2) {
         if (ndim > 1) {
            for (icell=0;icell<ncell;icell++) {
               float phivar = 0.0;
               fwrite( &phivar, sizeof(float), 1, writeENS);
            }
         }
      }
   }

   fclose(writeENS);

   // CaseFILE --------------------------------------------
   if (ndim == 1) 
      strcpy(ENScase,"scalar per element:      ");
   else
      strcpy(ENScase,"vector per element:      ");

   strcat(ENScase,String);
   strcat(ENScase,"   ");
   strcpy(FName,"EriCa.ENS.");
   strcat(FName,String);
   sprintf(Suffix, ".****");
   strcat(FName,Suffix);
   strcat(ENScase,FName);
   fprintf(CasePnt , "%s\n", ENScase);
}
