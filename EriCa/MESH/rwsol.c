#include "mesh.h"

void writeSOL( int itmst)
{
   if ((itmst+1) % nwrite != 0 && itmst != ntmst)
      return;

   int icell, ivar;
   double *CelldblPnt, *phi, *dep;
   FILE *solfile;

   fprintf(logfile, "\n --> Restart file\n\n");

   solfile = fopen("SYSTEM/EriCa.sol","wb");

   int nndir = ndir;

   fwrite( &nndir, sizeof(int), 1, solfile);
   fwrite( &ncell, sizeof(int), 1, solfile);
   fwrite( &nvar , sizeof(int), 1, solfile);

   // Writing cell centres coordinates ------
   for (icell=0;icell<ncell;icell++) {
      CelldblPnt = setdblcell( icell);
      fwrite( &CelldblPnt[CellXc], sizeof(double), ndir, solfile);
   }

   // Writing cell variables ----------------
   for (ivar=0;ivar<nvar;ivar++) {
      // Variable name ---------
      fwrite( VarName(ivar) , sizeof(char), 11, solfile);

      // Variable value --------
      for (icell=0;icell<ncell;icell++) {
         phi = setvar( icell);
         fwrite( &phi[ivar], sizeof(double), 1, solfile);
      }

      // Variable derivative ---
      for (icell=0;icell<ncell;icell++) {
         dep = setdep( icell);
         fwrite( &dep[PGrd+ivar*ndir], sizeof(double), ndir, solfile);
      }
   }

   fclose( solfile);
}

// Temporary read function (does not consider the gradients ..)
void readSOL( void)
{
   char flag = 0;
   int ierr;
   FILE *solfile;

   solfile = fopen("SYSTEM/EriCa.sol","rb");

   if (solfile==NULL)
      return;

   int nndir;
   ierr = fread( &nndir, sizeof(int), 1, solfile);
   if (nndir!=ndir)
      flag = 1;

   int nncell;
   ierr = fread( &nncell, sizeof(int), 1, solfile);
   if (nncell!=ncell)
      flag = 2;

   if (flag!=0) {
      printf("\nError READ file, code: %c\n", flag);
      FatalError(" ");
   }

   int nnvar;
   ierr = fread( &nnvar, sizeof(int), 1, solfile);

   // Coordinates reading ------------------- [to improve ...]
   double *coo;
   coo  = (double *) malloc( nncell * ndir * sizeof(double));
   ierr = fread( coo, sizeof(double), nncell*ndir, solfile);
   free(coo);

   // Variables reading ---------------------
   char   varnm[12];
   double *vvar, *grd, *phi;
   vvar = (double *) malloc( nncell * sizeof(double));
   grd  = (double *) malloc( nncell * ndir * sizeof(double));

   int ivvar, ivar, icell;

   for (ivvar=0;ivvar<nnvar;ivvar++) {
      // Reading info -----------
      ierr = fread( varnm, sizeof(char), 11, solfile);
      ierr = fread( vvar, sizeof(double), nncell, solfile);
      ierr = fread(  grd, sizeof(double), nncell*ndir, solfile);

      // Checking var presence --
      flag = 0;
      for (ivar=0;ivar<nvar;ivar++)
         if (strncmp( varnm, VarName(ivar), 11 )==0) {
            flag = 1;
            break;
         }

      if (flag==0)
         continue;

      // Overwriting var --------
      resvar[ivar] = 1;

      for (icell=0;icell<ncell;icell++) {
         phi = setvar( icell);
         phi[ivar] = vvar[icell];
      }
   }

   free(vvar);
   free( grd);

   fclose (solfile);
}
