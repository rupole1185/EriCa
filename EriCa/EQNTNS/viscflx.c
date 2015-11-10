#include "eqntns.h"

void setmuvar( double mu_i[nvar], double alfa, double depa[ndep], double depb[ndep], \
                                               double phia[nvar], double phib[nvar])
{
   int ivar;

   for (ivar=0;ivar<nvar;ivar++) {

      if (ivar == Pvar) 
         mu_i[ivar] = 0.0;
      else if (ivar>=Uvar && ivar<Uvar+ndir) {
         mu_i[ivar] = facevar( alfa, depa[MUdep], depb[MUdep]);
         if (MTdep!=-1)
            mu_i[ivar]+= facevar( alfa, depa[MTdep], depb[MTdep]); }
      else if (ivar==Tvar)
         mu_i[ivar] = facevar( alfa, depa[MUdep], depb[MUdep]);
      else if (ivar==MTvar) {
         mu_i[ivar] = facevar( alfa, depa[MUdep] + depa[Rdep] * phia[MTvar], \
                                     depb[MUdep] + depb[Rdep] * phib[MTvar]);
         mu_i[ivar]/= sigma; }
      else if (ivar==Kvar) {
         mu_i[ivar] = facevar( alfa, depa[MUdep], depb[MUdep]);
         mu_i[ivar]+= sigk * facevar( alfa, depa[MTdep], depb[MTdep]); }
      else if (ivar==OMvar) {
         mu_i[ivar] = facevar( alfa, depa[MUdep], depb[MUdep]);
         mu_i[ivar]+= sigw * facevar( alfa, depa[MTdep], depb[MTdep]); }
      else if (ivar>=SCAvar && ivar < SCAvar+NSCAvar) {
         mu_i[ivar] = facevar( alfa, depa[MUdep], depb[MUdep])  /  \
                      facevar( alfa, depa[SCHMIDTdep+ivar-SCAvar], \
                                     depb[SCHMIDTdep+ivar-SCAvar]) ; }
      else 
         FatalError("Mu method not defined");
   }
}

void viscscheme( int nvar, double mu_i[nvar], double AA[ndir], \
                 double dx[ndir], double F_i[nvar][nvar])
{
   int ivar, idir;

   // only for FULLY ORTHOGONAL compon.=> face axis passes 
   //                                     though cells centres
   for (ivar=0;ivar<nvar;ivar++)
         F_i[ivar][ivar] -= mu_i[ivar] * scalarp( ndir, AA, AA) / \
                                         scalarp( ndir, AA, dx); 
}

// Viscous fluxes
void viscflx(struct CSRMat *Mat, struct RHS *Rhs )
{
   int iface, ivar, acell, bcell, *FaceintPnt, idir;
   double *AA, *FacedblPnt, *varPnt, *dep_a, *dep_b;
   double phi_a[nvar], phi_b[nvar], F_i[nvar][nvar];
   double Flux_a[nvar], Flux_b[nvar], Flux[nvar], alfa;
   double *CellAdblPnt, *CellBdblPnt, dx[ndir];
   double mu_i[nvar];

   if (MUdep == -1)
      return;

   for (iface=0;iface<nface;iface++)  {
      //Setting Pointers
      FaceintPnt = setintface( iface);
      FacedblPnt = setdblface( iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      //Cell pointers 
      CellAdblPnt = setdblcell( acell);
      CellBdblPnt = setdblcell( bcell);

      //Initialization of F_i and mu_i
      memset(F_i , 0.0, nvar * nvar * sizeof(double));

      /*Variables*/
      varPnt = setvar( acell);
      dep_a  = setdep( acell);
      memcpy( phi_a, varPnt, nvar * sizeof(double));

      varPnt = setvar( bcell);
      dep_b  = setdep( bcell);
      memcpy( phi_b, varPnt, nvar * sizeof(double));

      alfa  = FacedblPnt[FaceAlfa];

      // To better understand the ViscousFlux check the formula 11.22
      //  of the Versteeg book (page 332)

      //Cell centres distance
      pntsdist( &CellBdblPnt[CellXc], &CellAdblPnt[CellXc], dx);

      AA   = &FacedblPnt[FaceSrf];

      // Flux Matrix -----------
      setmuvar( mu_i, alfa, dep_a, dep_b, phi_a, phi_b);
      viscscheme( nvar, mu_i, AA, dx, F_i);

      // Note: in continuity equation there is no diffusion
      phi_a[Pvar] = 0.0; 
      phi_b[Pvar] = 0.0;

      for (ivar=0;ivar<nvar;ivar++) 
         phi_a[ivar] *= -1.0;

      // Flux calculation
      matvecdstd(nvar, F_i, phi_a, Flux_a);
      matvecdstd(nvar, F_i, phi_b, Flux_b);
      summvec(nvar, Flux_a, Flux_b, Flux);

      // Function to add 2 RHS (it is the explicit part!!)
      getrhs( Rhs, acell, -1.0 * (1.0 - theta), Flux);
      getrhs( Rhs, bcell,  1.0 * (1.0 - theta), Flux);

      // A CELL fluxes ------------------------------------------------
      getnnz( Mat, acell, acell, -1.0 * theta, F_i);

      if (bcell<ncell) 
         getnnz( Mat, acell, bcell,  1.0 * theta, F_i);
      else 
         getrhs( Rhs, acell, -1.0 * theta, Flux_b);

      // B CELL fluxes ------------------------------------------------
      if (acell<ncell) 
         getnnz( Mat, bcell, acell,  1.0 * theta, F_i);
      else 
         getrhs( Rhs, bcell,  1.0 * theta, Flux_a);
         // NB: this should be negative, but Flux_a has already been 
         //   multiplied by a negative value!!

      getnnz( Mat, bcell, bcell, -1.0 * theta, F_i);

      // Non-orthogonality --------------------------------------------
      // STANDARD MODE ################################################
      for (ivar=0;ivar<nvar;ivar++) {
         phi_a[ivar] = scalarp( ndir, &dep_a[PGrd+(ivar*ndir)], &FacedblPnt[FaceDxA]);
         phi_b[ivar] = scalarp( ndir, &dep_b[PGrd+(ivar*ndir)], &FacedblPnt[FaceDxB]);
      }

      for (ivar=0;ivar<nvar;ivar++) 
         phi_a[ivar] *= -1.0;

      phi_a[Pvar] = 0.0; 
      phi_b[Pvar] = 0.0;

      matvecdstd(nvar, F_i, phi_a, Flux_a);
      matvecdstd(nvar, F_i, phi_b, Flux_b);
      summvec(nvar, Flux_a, Flux_b, Flux);

      getrhs( Rhs, acell, -1.0, Flux);
      getrhs( Rhs, bcell,  1.0, Flux);

      // MATHUR & MURTHY ##############################################
      /*double grdavg[ndir];

      for (ivar=0;ivar<nvar;ivar++) {
         for (idir=0;idir<ndir;idir++)
            grdavg[idir] = facevar( alfa, dep_a[PGrd+(ivar*ndir)+idir], dep_b[PGrd+(ivar*ndir)+idir]);

         Flux[ivar] = scalarp( ndir, grdavg, AA);
         Flux[ivar]-= scalarp( ndir, grdavg, dx) * scalarp( ndir, AA, AA) / scalarp( ndir, AA, dx);
         Flux[ivar]*= mu_i[ivar];
      }

      getrhs( Rhs, acell,  1.0, Flux);
      getrhs( Rhs, bcell, -1.0, Flux);*/
   }
}
