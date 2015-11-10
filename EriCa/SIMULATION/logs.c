#include "simulation.h"
#include "../MESH/mesh.h"

FILE *logfile;
int ScreenSize = 39, LogSize = 48;

void start(void )
{
   logfile = fopen("EriCa.log","w");

   ScreenStars();
   printf("*        EriCa 0.0.alpha.%1dD           *\n",ndir);
   ScreenStars();
   if (_DBG != 0) 
      printf("DEBUGGING MODE = %1d\n",_DBG);
   if (__DBG != 0)
      printf("MariSa DEBUG   = %1d\n",__DBG);

   printf("\nReading inputs       ...........  ");
   fflush(stdout);

   LogStars();
   fprintf(logfile,"*             EriCa 0.0.alpha.%1dD               *\n",ndir);
   LogStars();
   if (_DBG != 0) 
      fprintf(logfile,"\nDEBUGGING MODE = %1d\n",_DBG);
   if (__DBG != 0)
      fprintf(logfile,"\nMariSa DEBUG   = %1d\n",__DBG);

   // Folder creation
   mkdir("ENSIGHT", S_IRWXU);
   mkdir("OUTPUT", S_IRWXU);
   mkdir("SYSTEM", S_IRWXU);
   mkdir("OUTPUT", S_IRWXU);

   if (_DBG != 0 || __DBG != 0)
      mkdir("DEBUG", S_IRWXU);
}

void end( void)
{
   printf("\n\n");
   ScreenStars();
   printf("*        EriCa 0.0.alpha.%1dD           *\n",ndir);
   ScreenStars();

   LogStars();
   fprintf(logfile,"*             EriCa 0.0.alpha.%1dD               *\n",ndir);
   LogStars();

   fclose(logfile);
   fflush(stdout);
}

void finish( int itmst, clock_t time)
{

   clock_t diff = clock() - time;
   int msec = diff * 1000 / CLOCKS_PER_SEC;

   int hrs =   (msec)/ 3600000;
   int min =  ((msec)% 3600000) / 60000;
   int sec = (((msec)% 3600000) % 60000) / 1000;

   fprintf(logfile,"\nTotal iterations time = %02dh.%02dm.%02ds\n",hrs,min,sec);

   printCEI( itmst, ntmst);

   writeSOL( ntmst);

   end();
}

void FatalError( char *message)
{
   printf("\n @@ ERROR:\n %s", message);
   fprintf(logfile, "\n @@ ERROR:\n %s\n", message);

   end();
   exit(1);
}

void summary(void )
{
// Log File
   LogSection( "SETTINGS");
   fprintf(logfile,"General Settings\n");
   fprintf(logfile," Nranks:            %d\n", nranks);
   fprintf(logfile," # time step:       %d\n", ntmst);
   if (cfl < 0.0) 
      fprintf(logfile," Transient Deltat:  %lf\n", timst);
   else
      fprintf(logfile," Steady Courant:    %lf\n", cfl);

   fprintf(logfile,"\nNumeric Schemes\n");
   fprintf(logfile," Implicit level:    %5.2lf%%\n", 100*theta);
   fprintf(logfile," Num scheme:        ");
   switch (ialgh) {
   case 0:
      fprintf(logfile,"SOR relaxation\n");
      break;
   case 1:
      fprintf(logfile,"ConjGrad\n");
      break;
   case 2:
      fprintf(logfile,"BiCGStab\n");
      break;
   default:
      FatalError("Wrong selection ALGH");
   }

   fprintf(logfile," Preconditioning:   ");
   switch (iprec) {
   case 1:
      fprintf(logfile,"Jacoby\n");
      break;
   default:
      fprintf(logfile,"NO preconditioning\n");
      break;
   }

   fprintf(logfile," # iterations:      %d\n", niter);
   fprintf(logfile," Tollerance:        %lf\n", toll);
   fprintf(logfile," Convective scheme: ");
   switch (iconv) {
   case 0:
      fprintf(logfile,"Linear UpWind\n");
      break;
   case 1:
      fprintf(logfile,"Central Difference\n");
      break;
   default:
      FatalError("Wrong selection CONV");
   }

   fprintf(logfile,"\nFluid Model\n");
   fprintf(logfile," Fluid type:        ");
   switch (eos) {
   case ISIMPLE:
      fprintf(logfile,"Incompr. SIMPLE\n");
      break;
   case ISIMPLEC:
      fprintf(logfile,"Incompr. SIMPLEC\n");
      break;
   case IPISO:
      fprintf(logfile,"Incompr. PISO\n");
      FatalError("PISO does not work yet");
      break;
   case ICOUPLED:
      fprintf(logfile,"Incompr. COUPLED\n");
      break;
   case IDBIG:
      fprintf(logfile,"Ideal gas DENS-BASED\n");
      break;
   case IFROZENP:
      fprintf(logfile,"PRES-BASED Frozen Pressure\n");
      break;
   default:
      FatalError("Wrong selection EOS");
   }

   fprintf(logfile," Turbulence model:  ");
   switch (turb) {
   case INVISCID:
      fprintf(logfile,"Inviscid\n");
      break;
   case LAMINAR:
      fprintf(logfile,"Laminar\n");
      break;
   case DISTANCE:
      fprintf(logfile,"Laminar + Distance\n");
      break;
   case SA:
      fprintf(logfile,"Spalart-Allmaras\n");
      break;
   case KW06:
      fprintf(logfile,"K-w 2006\n");
      break;
   case KWSST:
      fprintf(logfile,"K-w SST\n");
      break;
   default:
      FatalError("Wrong selection TURB");
   }

   fflush(logfile);
}

void meshscreen( void)
{
   printf("done!\n");
   printf("\nProcessing database  ...........  ");
   fflush(stdout);
}

// Visualization of a loading bar with ETA
// IT current iteration, NIT total number of iterations, TIME reference time
void loadingbar(int it, int nit, clock_t time)
{
   if (nit < 25) 
      return;

   if ( (it+1) % (nit/25) != 0 && it !=0 )
      return;

   if ( it == 0)
      printf("0%%          50%%        100%%         ETA\n");
 
   int jj;
   printf("\r[");

   for (jj=0;jj<(it+1)/(nit/25);jj++)
      printf("#");

   for (jj=(it+1)/(nit/25);jj<25;jj++)
      printf(" ");

   printf("]");

   int hrs, min, sec, msec;
   clock_t diff;

   diff = clock() - time;
   msec = diff * 1000 / CLOCKS_PER_SEC / (it+1);

   hrs =   ((nit-it-1) * msec)/ 3600000;
   min =  (((nit-it-1) * msec)% 3600000) / 60000;
   sec = ((((nit-it-1) * msec)% 3600000) % 60000) / 1000;

   printf(" %02dh.%02dm.%02ds", hrs, min, sec);

   fflush(stdout);
}

// Iterations output -----------------------------------
clock_t timei;
void iterating(int itmst, int niterr, double tolll, clock_t time)
{
   clock_t diff, timej;
   int msec, ivar;

   timej = clock();

   if (itmst==0) 
      diff = timej - time;
   else 
      diff = timej - timei;

   timei = timej;

   msec = diff * 1000 / CLOCKS_PER_SEC;

   fprintf(logfile, "    KIRiteration = %4d, Res = %f\n", niterr, tolll);

   // Variables Ranges ---------
   fprintf(logfile, "\n  ___________________________________________ \n");
   fprintf(logfile, " | Variable    |      Minimum |      Maximum |\n");
   fprintf(logfile, " |-------------|--------------|--------------|\n");

   for (ivar=0;ivar<nvar;ivar++)
      Range2table( VarName(ivar), ivar, 0);

   Range2table( "TimeStp", TIMEdep, 1);
   Range2table( "Courant", CFLdep, 1);

   fprintf(logfile, "  -------------------------------------------\n");

   fprintf(logfile, "\nTimeStep %d - CPU time %d.%03d s\n\n", itmst+1, msec/1000, msec%1000);
   LogStars();
   fprintf(logfile,"\n");
   fflush(logfile);

   loadingbar( itmst, ntmst, time);

   writeSOL( itmst);
}

void inisum(void )
{
   int ivar, ibnd, idep;

   LogSection("INITIALIZATION");

   fprintf(logfile,"Variables:\n");
   for (ivar=0;ivar<nvar;ivar++)
      if (resvar[ivar]==0)
         fprintf(logfile, " %s    Var = %15.8lf\n",VarName(ivar), varini[ivar]);
      else
         fprintf(logfile, " %s    From restart file\n",VarName(ivar));

   fprintf(logfile,"\nDependents:\n");
   for (idep=0;idep<ndepini;idep++)
      fprintf(logfile, " %s    Dep = %15.8lf\n",DepName(idep), depini[idep]);

   LogSection("BOUNDARIES");
   
   for (ibnd=0;ibnd<nbnd;ibnd++) {
      fprintf(logfile,"\nBnd # %d:\n", ibnd);
      for (ivar=0;ivar<nvar;ivar++)
         fprintf(logfile, " %s  Typ = %2d  Bnd = %11.2lf\n",VarName(ivar), \
                        bndtype[ibnd*nvar+ivar], bnd[ibnd*nvar+ivar]);
   }

   printf("done!\n");

   distance();

   LogSection("ITERATING");
   fflush(logfile);

   printf("\nIterations progress ... \n");
   fflush(stdout);

   free(resvar);
   free(varini);
   free(depini);
   free(bndtype);
   free(bnd);
}

void meminfo( struct CSRMat *Mat, struct RHS *b, struct CSRMat *MatGrd)
{
   LogSection("ALLOCATIONS");
   fprintf(logfile,"Mesh sizes\n");
   fprintf(logfile," Ncell:  %d\n", ncell);
   fprintf(logfile," Nghst:  %d\n", nghst);
   fprintf(logfile," NtotC:  %d\n", ntot );
   fprintf(logfile," Nface:  %d\n", nface);
   fprintf(logfile," Nnz:    %d\n", Matnnz(Mat));

   fprintf(logfile,"\nCardinalities\n");
   fprintf(logfile," Nvar:   %d\n", nvar);
   fprintf(logfile," Ndep:   %d\n", ndep);

   double pow2, tmp, totAlloc;
   totAlloc = 0.0;
   pow2 = pow(2.0,20.0);

   fprintf(logfile,"\nMain arrays\n");
   tmp = nvar * (ntot+ncell) * sizeof(double);
   totAlloc += tmp;
   fprintf(logfile," Variables:   %9.3f MB\n",tmp/pow2);
   tmp = ndep * ntot * sizeof(double);
   totAlloc += tmp;
   fprintf(logfile," Dependents:  %9.3f MB\n",tmp/pow2);

   fprintf(logfile,"\nGeometry arrays\n");
   tmp = ((ncell) * sizeof(int));
   totAlloc += tmp;
   fprintf(logfile," CSR mesh:    %9.3f MB\n",tmp/pow2);
   tmp = ncell * (ncelldbl * sizeof(double) + ncellint * sizeof(int));
   totAlloc += tmp;
   fprintf(logfile," Cells data:  %9.3f MB\n",tmp/pow2);
   tmp = nface * (nfacedbl * sizeof(double) + nfaceint * sizeof(int));
   totAlloc += tmp;
   fprintf(logfile," Faces data:  %9.3f MB\n",tmp/pow2);
   tmp = nghst * (nghstdbl * sizeof(double) + nghstint * sizeof(int));
   totAlloc += tmp;
   fprintf(logfile," Ghsts data:  %9.3f MB\n",tmp/pow2);

   fprintf(logfile,"\nMatrixes\n");
   tmp = MatMem( Mat);
   totAlloc += tmp;
   fprintf(logfile," Matrix:      %9.3f MB\n",tmp/pow2);
   tmp = MatMem( MatGrd);
   totAlloc += tmp;
   fprintf(logfile," MatrixGrd:   %9.3f MB\n",tmp/pow2);
   tmp = RhsMem( b);
   totAlloc += tmp;
   fprintf(logfile," Rhs:         %9.3f MB\n",tmp/pow2);

   fprintf(logfile,"___________________________\n");
   fprintf(logfile," TOTAL:       %9.3f MB\n",totAlloc/pow2);

   fflush(logfile);
}
