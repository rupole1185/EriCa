#include "eqntns.h"

// function that returns the cell distance from walls
//      (these are defined for velocities BNDType = -1)
double *phi_d, *grad_d;
int     dloop = 50, nitd = 100;
double  told = 0.0001, thetad = 1.0;

// Distance equation based on CFD-online website

void dist_teq( struct CSRMat *Mat_d, struct RHS *Rhs_d)
{
   int    iface, *FaceintPnt, acell, bcell, idir, icell;
   double *FacedblPnt, *AA, mu[1], alfa, Flux[1];
   double *CellAdblPnt, *CellBdblPnt, dx[ndir], Src[1], F_i[1][1];
   double Fluxa[1], Fluxb[1], *CelldblPnt, phi_a[1], phi_b[1];

   // Artificial Viscosity
   mu[0] = 1.0;

   for (iface=0;iface<nface;iface++) {
      //Setting Pointers
      FaceintPnt = setintface( iface);
      FacedblPnt = setdblface( iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      AA   = &FacedblPnt[FaceSrf];

      alfa = FacedblPnt[FaceAlfa];

      //Cell pointers 
      CellAdblPnt = setdblcell( acell);
      CellBdblPnt = setdblcell( bcell);

      // Phia & Phib
      phi_a[0]  = phi_d[acell];
      phi_b[0]  = phi_d[bcell];

      pntsdist( &CellBdblPnt[CellXc], &CellAdblPnt[CellXc], dx);

      //Initialization of F_i and mu_i
      memset( F_i, 0.0, 1 * 1 * sizeof(double));

      // DIFFUSION ---------------------------------------
      viscscheme( 1, mu, AA, dx, F_i);

      phi_a[0] *= -1.0;

      matvecdstd( 1, F_i, phi_a, Fluxa);
      matvecdstd( 1, F_i, phi_b, Fluxb);
      summvec( 1, Fluxa, Fluxb, Flux);

      //      Explicit
      getrhs( Rhs_d, acell, -1.0 * (1.0-thetad), Flux);
      getrhs( Rhs_d, bcell,  1.0 * (1.0-thetad), Flux);

      //      Implicit
      getnnz( Mat_d, acell, acell, -1.0 * thetad, F_i);

      if (bcell<ncell) 
         getnnz( Mat_d, acell, bcell,  1.0 * thetad, F_i);
      else 
         getrhs( Rhs_d, acell, -1.0 * thetad, Fluxb);

      if (acell<ncell) 
         getnnz( Mat_d, bcell, acell,  1.0 * thetad, F_i);
      else 
         getrhs( Rhs_d, bcell,   1.0 * thetad, Fluxa);

      getnnz( Mat_d, bcell, bcell, -1.0 * thetad, F_i);
   }

   for (icell=0;icell<ncell;icell++) {
      CelldblPnt = setdblcell(icell);

      // SOURCE -----------------------------------------
      Src[0] = + 1.0;

      getrhs( Rhs_d, icell, CelldblPnt[CellVol], Src);
   }
}

void dist_bc( void)
{
   int ighst, icell, *GhstintPnt;

   for (ighst=0;ighst<nghst;ighst++) {
      GhstintPnt = setintghst( ighst);

      icell = GhstintPnt[GhstFidx];

      phi_d[ ghst2gbl( ighst)] = phi_d[ icell];

      if ( isitwall( ighst)) 
         phi_d[ ghst2gbl( ighst)] = 0.0;
   }
}

void calc_dist()
{
   int icell, idir;
   double grad2, *dep;

   for (icell=0;icell<ncell;icell++) {
      dep = setdep( icell);

      grad2 = scalarp( ndir, &grad_d[ ndir*icell], &grad_d[ ndir*icell]);

      dep[Ddep] = sqrt( grad2 + 2.0 * phi_d[icell] ) - sqrt( grad2);

#if _DBG == 10
      dep[PHIDdep] = phi_d[icell];
#endif
   }
}

void distance( void)
{
   int iloop, icell;
   int    niterr;
   double tolll, *dep, cnst;
   struct CSRMat Mat_d;
   struct RHS    Rhs_d;

   if ( Ddep == -1)
      return;

   printf("\nWall distance ...\n");
   fflush(stdout);
   LogSection("WALL DISTANCE");

   phi_d = (double *) calloc( ntot, sizeof(double));
   grad_d= (double *) calloc( ndir*ntot, sizeof(double));

   for (icell=0;icell<ncell;icell++) 
      phi_d[icell]= 1.0;

   matini( ncell, 1, &Mat_d);
   matalc( &Mat_d, dnnz);
   rhsini( ncell, 1, &Rhs_d);

   clock_t time;
   time = clock();

   for (iloop=0;iloop<dloop;iloop++) {
      loadingbar( iloop, dloop, time);

      matstart( 0.0, 0, &Mat_d);
      rhsstart( 0.0, &Rhs_d);

      dist_bc();
      dist_teq( &Mat_d, &Rhs_d);

      niterr = nitd;
      tolll  = told;
      SORmethd( &Mat_d, &Rhs_d, phi_d, 1.0, &tolll, &niterr);

      if (isnan(tolll))
         FatalError(" Distance diverged");

      fprintf(logfile, "    Dist  iter   = %4d, Res = %f\n", niterr,  tolll);
      fflush(logfile);
   }

   gradients( 1, phi_d, 0, ndir, grad_d);
   calc_dist();

   matdel( &Mat_d);
   rhsdel( &Rhs_d);
   free( phi_d);
   free( grad_d);

   printf("\n");
}
