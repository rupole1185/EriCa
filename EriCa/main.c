#include "EQNTNS/eqntns.h"

void main()
{
   int itmst, ierr, niterr;
   double tolll, *sol;
   struct CSRMat Mat;
   struct RHS    b;

   start();

// ReadMesh -----------------------------------------
   readmesh();

// Main SIMULATION settings -------------------------
   fldini();                   //fluid initialization
   simini();                   //simulations settings
   bndini();                   //boundary  conditions
   summary();

// MESH ---------------------------------------------
   meshscreen();
   createmesh();    //mesh creation, allocation, init
   gradinit();      //gradient matrix init (just GEO)
   moninit();

// MATs alloc ---------------------------------------
   int pttrn[nvar][nvar];
   pattern( pttrn);
   matini(ncell, nvar, &Mat);
   matptrn( &Mat, pttrn);
   if (theta==0.0)
      matdiag( &Mat);
   else
      matalc( &Mat, dnnz);
   rhsini(ncell, nvar, &b);

   meminfo( &Mat, &b, &MatGrd);

   sol = (double *) malloc(nvar*ncell* \
                    sizeof(double ));   // Sol  array
// ITERATIONS ---------------------------------------
   inisum();
   init_turb_const();

   clock_t time;
   time = clock();

   for (itmst=0;itmst<ntmst;itmst++) {
      matstart( 0.0, 1, &Mat);
      rhsstart( 0.0, &b);

      //Update BC, deps, Grads, etc -----------------
      boundary();

      gradients( nvar, setvar( 0), PGrd, ndep, setdep( 0));

      depupdt();

      //OutPuts -------------------------------------
      printCEI( itmst, ntmst);

      monitors();

      //Fluxes & sources ----------------------------
      convflx( &Mat, &b );
      viscflx( &Mat, &b );

      timestp( &Mat, &b );
      sources( &Mat, &b );

      //System definition ---------------------------
      pv_adjustment( &Mat, &b, sol);

      //Preconditioning
      switch (iprec) {
      case 1:
         pcjacb( &Mat);
         break;
      default:
         break;
      }

      tolll = toll;
      niterr= niter;

      //KirSolver Method
      switch (ialgh) {
      case 0:
         SORmethd( &Mat, &b, sol, 1.0, &tolll, &niterr);
         break;
      case 1:
         ConjGrad( &Mat, &b, sol, &tolll, &niterr);
         break;
      case 2:
         BiCGStab( &Mat, &b, sol, &tolll, &niterr);
         break;
      }

      //Checking Convergence ------------------------
      if (isnan(tolll))
         FatalError("Solution diverged");

      //Update solution -----------------------------
      update( &Mat, &b, sol);

      //Log file information ------------------------
      iterating( itmst, niterr, tolll, time);

      //User interaction ----------------------------
      user_order();
   }

   matdel( &Mat);
   matdel( &MatGrd);
   rhsdel( &b);

   finish( itmst, time);
}
