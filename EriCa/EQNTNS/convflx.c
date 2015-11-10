#include "eqntns.h"

//          artif coeff    Flux    Flux + diff variable      output
void convscheme(int iconv, int nvar, double multp, double Fi, double Flux[nvar][nvar])
{
   int ivar;
   double Ffi;

   switch (iconv) {
   // UpWind Difference Scheme - UDS
   case 0:
      Ffi = multp * fmax(Fi,0.0);

      for (ivar=0;ivar<nvar;ivar++)
         Flux[ivar][ivar] = Ffi;

      break;
   // Central difference scheme - CD
   case 1:
      Ffi = multp * 0.5 * Fi;

      for (ivar=0;ivar<nvar;ivar++)
         Flux[ivar][ivar] = Ffi;

      break;
   default :
      FatalError("Wrong ConvScheme selected!");
      break;
   }
}

// Convective fluxes
void convflx(struct CSRMat *Mat, struct RHS *Rhs )
{
   int iface, acell, bcell, *FaceintPnt, idir, ivar;
   double *AA, u_i[ndir], rho_a, rho_b, rho_i, *FacedblPnt, *varPnt, *dep_a;
   double phi_a[nvar], phi_b[nvar], F_a[nvar][nvar], F_b[nvar][nvar], *dep_b;
   double Flux_a[nvar], Flux_b[nvar], Flux[nvar], alfa, mssflw;

   for (iface=0;iface<nface;iface++)  {

      //Setting Pointers
      FaceintPnt = setintface( iface);
      FacedblPnt = setdblface( iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      //Initialization of F_a and F_b
      memset(F_a, 0.0, nvar*nvar* sizeof(double));
      memset(F_b, 0.0, nvar*nvar* sizeof(double));

      /*Variables to convect*/
      varPnt = setvar( acell);
      dep_a  = setdep( acell);
      memcpy( phi_a, varPnt, nvar * sizeof(double));
      rho_a = dep_a[Rdep];

      varPnt = setvar( bcell);
      dep_b  = setdep( bcell);
      memcpy( phi_b, varPnt, nvar * sizeof(double));
      rho_b  = dep_b[Rdep];

      //Convective fluxes calculation:
      // Coeff * Phi_a + Coeff * Phi_b
      // Coeff = max(F_i, ..., 0) depends on the convective scheme

      //Convective fluxes: Rho * u * A * Phi ---------------------------
      alfa  = FacedblPnt[FaceAlfa];

      rho_i   = facevar( alfa, rho_a, rho_b);
      AA      = &FacedblPnt[FaceSrf];

      for (idir=0;idir<ndir;idir++) 
         u_i[idir] = facevar( alfa, \
                     phi_a[Uvar+idir] + scalarp( ndir, &dep_a[UGrd+ndir*idir], &FacedblPnt[FaceDxA]), \
                     phi_b[Uvar+idir] + scalarp( ndir, &dep_b[UGrd+ndir*idir], &FacedblPnt[FaceDxB]));

      mssflw = rho_i * scalarp( ndir, u_i, AA);

      //conv returns the convective scheme
      convscheme( iconv, nvar, 1.0,       mssflw, F_a);
      convscheme( iconv, nvar,-1.0,-1.0 * mssflw, F_b);

      //////------------------------ PRESS-VEL COUPLING
      presflx(phi_a, F_a);
      presflx(phi_b, F_b);
      /////--------------------------------------------

      // Flux calculation
      matvecdstd(nvar, F_a, phi_a, Flux_a);
      matvecdstd(nvar, F_b, phi_b, Flux_b);
      summvec(nvar, Flux_a, Flux_b, Flux);

      // Function to add 2 RHS (it is the explicit part!!)
      getrhs( Rhs, acell, -1.0 * (1.0 - theta), Flux);
      getrhs( Rhs, bcell,  1.0 * (1.0 - theta), Flux);

      // A CELL fluxes ------------------------------------------------
      getnnz( Mat, acell, acell, 1.0 * theta, F_a);

      if (bcell<ncell) 
         getnnz( Mat, acell, bcell, 1.0 * theta, F_b);
      else 
         getrhs( Rhs, acell, -1.0 * theta, Flux_b);

      // B CELL fluxes ------------------------------------------------
      if (acell<ncell) 
         getnnz( Mat, bcell, acell, -1.0 * theta, F_a);
      else 
         getrhs( Rhs, bcell, 1.0 * theta, Flux_a);

      getnnz( Mat, bcell, bcell, -1.0 * theta, F_b);

      // Non-orthogonality --------------------------------------------
      // It is like an explicit flux where (Var * Dx) replaces phi
      for (ivar=0;ivar<nvar;ivar++) {
         phi_a[ivar] = scalarp( ndir, &dep_a[PGrd+(ivar*ndir)], &FacedblPnt[FaceDxA]);
         phi_b[ivar] = scalarp( ndir, &dep_b[PGrd+(ivar*ndir)], &FacedblPnt[FaceDxB]);
      }

      presflx( phi_a, F_a);
      presflx( phi_b, F_b);
 
      matvecdstd(nvar, F_a, phi_a, Flux_a);
      matvecdstd(nvar, F_b, phi_b, Flux_b);
      summvec(nvar, Flux_a, Flux_b, Flux);

      getrhs( Rhs, acell, -1.0, Flux);
      getrhs( Rhs, bcell,  1.0, Flux);
   }
}
