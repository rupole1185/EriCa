#include "eqntns.h"

// Coupled alghorithm by Darwish

// Pressure treatment
void coupled_press( struct CSRMat *Mat, struct RHS *Rhs)
{
   int iface, acell, bcell, *FaceintPnt, idir;
   double *AA, *FacedblPnt, *dep_a, alfa, p_i, u_i[ndir], dx[ndir];
   double *phi_a, *phi_b, F_i[nvar][nvar], *dep_b, div_exp[nvar];
   double force[nvar], rho_i, F_a[nvar][nvar], F_b[nvar][nvar];
   double *dd_a, *dd_b, A_a, A_b, F_pa, F_pb, p_avg[ndir];
   double *CellAdblPnt, *CellBdblPnt;

   calcdd( Mat);

   for (iface=0;iface<nface;iface++)  {

      //Setting Pointers
      FaceintPnt = setintface( iface);
      FacedblPnt = setdblface( iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      //Variables
      phi_a  = setvar( acell);
      dep_a  = setdep( acell);

      phi_b  = setvar( bcell);
      dep_b  = setdep( bcell);

      alfa  = FacedblPnt[FaceAlfa];
      AA    = &FacedblPnt[FaceSrf];

      // Pressure effect on momentum ------------------
      memset(F_i, 0.0, nvar*nvar* sizeof(double));
      memset(force, 0.0, nvar* sizeof(double));

      p_i = facevar( alfa, phi_a[Pvar], phi_b[Pvar]);

      for (idir=0;idir<ndir;idir++) {
         // Explicit force
         force[Uvar+idir] = p_i * AA[idir];

         // Implicit force
         F_i[Uvar+idir][Pvar] = AA[idir];
      }

      getrhs( Rhs, acell, -1.0 * (1.0 - theta), force);
      getrhs( Rhs, bcell,  1.0 * (1.0 - theta), force);

      // A CELL fluxes -------------
      getnnz( Mat, acell, acell, 1.0 * theta * alfa, F_i);

      if (bcell<ncell) 
         getnnz( Mat, acell, bcell, 1.0 * theta * (1.0 - alfa), F_i);
      else 
         getrhs( Rhs, acell, -1.0 * theta, force);

      // B CELL fluxes -------------
      if (acell<ncell) 
         getnnz( Mat, bcell, acell, -1.0 * theta * alfa, F_i);
      else 
         getrhs( Rhs, bcell, 1.0 * theta, force);

      getnnz( Mat, bcell, bcell, -1.0 * theta * (1.0 - alfa), F_i);

      // Non-orthogonality ---------
      p_i = facevar( alfa, scalarp( ndir, &dep_a[PGrd], &FacedblPnt[FaceDxA]), \
                           scalarp( ndir, &dep_b[PGrd], &FacedblPnt[FaceDxB]));

      for (idir=0;idir<ndir;idir++)
         force[Uvar+idir] = p_i * AA[idir];

      getrhs( Rhs, acell, -1.0, force);
      getrhs( Rhs, bcell,  1.0, force);

      // Continuity equation reconstruction -----------
      // NB: momentum equation is 100% implicit!!!

      if (acell<ncell) {
         dd_a = setsol( acell, ndir, dd);
         A_a  = fabs( scalarp( ndir, dd_a, AA) / sqrt( scalarp( ndir, AA, AA)));
      }

      if (bcell<ncell) {
         dd_b = setsol( bcell, ndir, dd);
         A_b  = fabs( scalarp( ndir, dd_b, AA) / sqrt( scalarp( ndir, AA, AA)));
      }

      if (acell>=ncell)
         A_a = A_b;

      if (bcell>=ncell)
         A_b = A_a;

      CellAdblPnt = setdblcell( acell);
      CellBdblPnt = setdblcell( bcell);

      F_pa = CellAdblPnt[CellVol] / A_a;
      F_pb = CellBdblPnt[CellVol] / A_b;

      // Term 1: Velocity divergence
      memset(F_a, 0.0, nvar*nvar* sizeof(double));
      memset(F_b, 0.0, nvar*nvar* sizeof(double));
      memset(div_exp, 0.0, nvar* sizeof(double));

      rho_i = facevar( alfa, dep_a[Rdep], dep_b[Rdep]);

      for (idir=0;idir<ndir;idir++) {
         F_a[Pvar][Uvar+idir] +=        alfa  * rho_i * AA[idir]; 
         F_b[Pvar][Uvar+idir] += (1.0 - alfa) * rho_i * AA[idir]; 
      }

      // BC for divergence ----
      if (acell>=ncell)
         div_exp[Pvar] -= rho_i * alfa * scalarp( ndir, &phi_a[Uvar], AA);

      if (bcell>=ncell)
         div_exp[Pvar] -= rho_i * (1.0 - alfa) * scalarp( ndir, &phi_b[Uvar], AA);

      // Non-orthogonality ----
      for (idir=0;idir<ndir;idir++)
         u_i[idir] = facevar( alfa, \
              scalarp( ndir, &dep_a[UGrd+ndir*idir], &FacedblPnt[FaceDxA]), \
              scalarp( ndir, &dep_b[UGrd+ndir*idir], &FacedblPnt[FaceDxB]));

      div_exp[Pvar] -= rho_i * scalarp( ndir, u_i, AA); 

      // Term 2: Implicit pressure gradient
      pntsdist( &CellBdblPnt[CellXc], &CellAdblPnt[CellXc], dx);

      F_a[Pvar][Pvar] += rho_i * facevar( alfa, F_pa, F_pb) * scalarp( ndir, AA, AA) / \
                                                              scalarp( ndir, AA, dx);

      F_b[Pvar][Pvar] -= rho_i * facevar( alfa, F_pa, F_pb) * scalarp( ndir, AA, AA) / \
                                                              scalarp( ndir, AA, dx);

      // BC for pressure implicit gradient
      if (acell>=ncell)
         div_exp[Pvar] -= rho_i * facevar( alfa, F_pa, F_pb) * phi_a[Pvar] * scalarp( ndir, AA, AA) / \
                                                                             scalarp( ndir, AA, dx);

      if (bcell>=ncell)
         div_exp[Pvar] += rho_i * facevar( alfa, F_pa, F_pb) * phi_b[Pvar] * scalarp( ndir, AA, AA) / \
                                                                             scalarp( ndir, AA, dx);

      // Non-orthogonality ---
      p_i = scalarp( ndir, &dep_b[PGrd], &FacedblPnt[FaceDxB]) - \
            scalarp( ndir, &dep_a[PGrd], &FacedblPnt[FaceDxA]);
      
      div_exp[Pvar] += rho_i * facevar( alfa, F_pa, F_pb) * p_i * scalarp( ndir, AA, AA) / \
                                                                  scalarp( ndir, AA, dx); 

      // Term 3: Explicit pressure gradient
      for (idir=0;idir<ndir;idir++)
         p_avg[idir] = facevar( alfa, F_pa, F_pb) * facevar( 0.5,  dep_a[PGrd+idir], dep_b[PGrd+idir]);

      div_exp[Pvar] -= rho_i * scalarp( ndir, p_avg, AA);

      // ASSIGNING TO CELLS ----
      getrhs( Rhs, acell,  1.0, div_exp);
      getrhs( Rhs, bcell, -1.0, div_exp);

      // A CELL fluxes -------------
      getnnz( Mat, acell, acell,  1.0, F_a);

      getnnz( Mat, acell, bcell,  1.0, F_b);

      // B CELL fluxes -------------
      getnnz( Mat, bcell, acell, -1.0, F_a);

      getnnz( Mat, bcell, bcell, -1.0, F_b);
   }

   free(dd);
}

// Variable update
void coupled_upd(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   int icell, istart, ivar;
   double *phi, *phi_new;

   for (icell=0;icell<ncell;icell++) {
      phi     = setvar( icell);
      phi_new = setsol( icell, nvar, sol);

      for (ivar=0;ivar<nvar;ivar++)
         phi[ivar] = phi_new[ivar];
   }
}
