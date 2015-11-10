#include "simulation.h"

int  ncell, nghst, ntot, nface;

// FLUID VARIABLES -----------------------------------------------
int  turb=-1, eos=-1, restart=1, *resvar;
int  nvar, ndep, ndepini;
int  Pvar=-1, Uvar=-1, Tvar=-1, Kvar=-1, OMvar=-1, Dvar=-1, MTvar=-1;
int  SCAvar=-1, NSCAvar=0;
int  Rdep=-1, MUdep=-1, MTdep=-1, Ddep=-1;
int  PGrd=-1, UGrd=-1, TGrd=-1, KGrd=-1, OMGrd=-1, MTGrd=-1, SCAGrd=-1, CFLdep=-1;
int  SCHMIDTdep=-1, TIMEdep = -1;
double *varini = NULL, *depini = NULL, Cp = -1.0, kfld = -1.0, R;
#if _DBG == 10
int  Ppdep = -1, Updep = -1, Divdep = -1, Pkdep = -1, Vrtdep = -1;
int  DivSdep = -1, PHIDdep = -1, DDdep = -1;
#endif

void fldini(void )
{
   int ierr;
   file_c_manager fluid = { "./\0", "fluid.fld\0"};

   ierr = fluid_get_scalars( fluid, &NSCAvar);

   ierr = fluid_get_fluidtype( fluid, &eos, &turb);
   simstencil();

   ierr = fluid_get_fluidconst( fluid, &Cp, &kfld);

   varini = (double *) malloc( nvar * sizeof(double));
   ierr   = fluid_get_varini( fluid, &varini, nvar);

   depini = (double *) malloc( ndepini * sizeof(double));
   ierr = fluid_get_depini( fluid, &depini, ndepini);

   restart = fluid_get_restart( fluid);
   resvar  = (int *) calloc( nvar, sizeof(int));

   // CHECKING ----------------------------------------------
   if (depini == NULL || varini == NULL || turb ==-1 || eos == -1) {
      printf("\nERROR: fluid not properly defined!\n");
      printf("It is compulsory to define:\n");
      printf("turbulence\nfluidtype\nvarini\ndepini\n");
      FatalError(" ");
   }

   if ( eos == IDBIG) {
      if ( Cp < 0.0 || kfld < 0.0 )
         FatalError("fluidconst: Cp and k not defined");

      R = kfld / (kfld - 1.0) * Cp; 
   }
}


int  iconv = 0, ialgh = 0, iprec = 0, nwrite = 100;
double theta = 1.0, timst = -1.0, toll = 0.001, tollp = 0.000001, cfl = -1.0;
int  niter = 50, niterp = 1000, ntmst = 50, nplot = 10, nranks = 1;

void simini(void )
{
   int   ierr;
   file_c_manager simul = { "./\0", "simul.sim\0"};

   ierr = sim_get_convscheme( simul, &iconv);

   ierr = sim_get_steady( simul, &cfl, &ntmst);
   ierr = sim_get_transient( simul, &timst, &ntmst);

   ierr = sim_get_solver( simul, &ialgh, &iprec, &toll, &niter);

   ierr = sim_get_printCEI( simul, &nplot);

   ierr = sim_get_nranks( simul, &nranks);
#ifdef _OPENMP
   omp_set_num_threads(nranks);
#endif

   ierr = sim_get_timescheme( simul, &theta);

   ierr = sim_get_psimple( simul, &tollp, &niterp);

   ierr = sim_get_saveSOL( simul, &nwrite);

   // CHECKING ----------------------------------------------
   if (cfl < 0.0 && timst < 0.0)
      FatalError("transient/steady: choose one of them");

   if (eos > COMPR_b) {
      if (theta!=0.0)
         FatalError("timescheme: Density-Based works only if explicit");
   }
}

int  *bndtype, nbnd;
double *bnd;
int  *monsurf, *monphi, *monave, *monorm, *monare, nmon;
double *mon;

void bndini(void )
{
   int ierr, ibnd, nmon_bnd, *tmp_mon, imon, iimon;
   file_c_manager boundary = { "./\0", "boundary.bnd\0", \
                               "./SYSTEM/\0", "Bnd.dat\0"};

   // Number of BND ...
   ierr = bnd_get_nbnd( boundary, &nbnd);
   //   ... and allocation
   bndtype  = realloc(bndtype, nbnd * nvar * sizeof(int ));
   bnd      = realloc(bnd, nbnd * nvar * sizeof(double));

   nmon = 0;

   for (ibnd=0;ibnd<nbnd;ibnd++) {
      ierr = bnd_get_bndtype( boundary, ibnd, &bndtype[ibnd*nvar], nvar);

      if (ierr==1)
         FatalError("BNDTYPE missing");

      ierr = bnd_get_bndvar( boundary, ibnd, &bnd[ibnd*nvar], nvar);

      if (ierr==1)
         FatalError("BNDVAR missing");

      ierr = bnd_get_monitor( boundary, ibnd, &tmp_mon, &nmon_bnd);
 
      if (nmon_bnd>0) {
         imon    = nmon;
         iimon   = 0;

         nmon   += nmon_bnd; 
         monsurf = realloc( monsurf, nmon * sizeof(int));
         monphi  = realloc( monphi , nmon * sizeof(int));
         monave  = realloc( monave , nmon * sizeof(int));
         monare  = realloc( monare , nmon * sizeof(int));
         monorm  = realloc( monorm , nmon * sizeof(int));
         mon     = realloc( mon    , nmon * sizeof(double));

         do {
            monsurf[imon] = ibnd;
            monphi[imon]  = tmp_mon[iimon*4+0];
            monare[imon]  = tmp_mon[iimon*4+1];
            monave[imon]  = tmp_mon[iimon*4+2];
            monorm[imon]  = tmp_mon[iimon*4+3];

            imon++;
            iimon++;
         } while (imon < nmon);
      }
   }

   // CHECKING ---------------------------------------------
   for (imon=0;imon<nmon;imon++) {
      if (monphi[imon]==Pvar) {
         if (monare[imon] == -1)
            FatalError("monitor: Pressure must be assigned a direction");
      }

      if (monphi[imon]>=Uvar && monphi[imon]<Uvar+ndir)
         if (monare[imon] != -1)
            FatalError("monitor: Velocity must be assigned to Area = -1");

      if (monphi[imon]>=nvar)
         FatalError("monitor: Wrong IVAR");

      if (monphi[imon]==-1 && monave[imon]!=1)
         FatalError("monitor: IVAR = -1 only if AVE = 1");
   }
}
