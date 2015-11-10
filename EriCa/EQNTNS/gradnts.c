#include "eqntns.h"

//GRADIENTs CALCULATION:
// Grad(phi) = (SUM phi_face * Area ) / Vol
//The formula above defines the gradient according to
// the divergence theorem.

// INPUTS : nvar    -> number of variables to calculate the gradient
//          varPNT  -> array with all the variables

// OUTPUTS: grdPNT  -> array whith all the computed gradients
//          grdSTART-> starting location of the gradient
//          nnvar   -> number of dependents in the gradient array

// Reference article:
// Discretization and parallel performance of an unstructured finite volume Navier-Stokes solver
// S.A. Karimian and A. G. Straatman
//
// IT HAS BEEN RECENTLY UPDATED IN ORDER TO TAKE BETTER ADVANTAGE OF THE DATA-STRUCTURE
// Vol * GRAD_Acell( PHI) = SUM_NEIGH( alfa * ( PHI_Acell + GRAD_Acell(PHI) * DX_Acell) + \
                               (1.0 - alfa) * ( PHI_Neigh + GRAD_Neigh(PHI) * DX_Neigh)) * Area

struct CSRMat MatGrd;
double c_a =  1.0;
double c_b = -1.0;

// Routine to fill the matrix (it is just geometry) ---------------
void gradinit( void)
{
   int icell, idir, jdir, iface, acell, bcell;
   int *FaceintPnt, *GhstintPnt;
   double *FacedblPnt, Rmat_a[ndir][ndir], Rmat_b[ndir][ndir];
   double Vmat[ndir][ndir], *CelldblPnt, alfa;

   matini( ncell, ndir, &MatGrd);
   matalc( &MatGrd, dnnz);

   memset(Vmat, 0.0, ndir*ndir* sizeof(double));

   // Matrix is filled in a separate loop since it is just geometry
   matstart( 0.0, 1, &MatGrd);

   //Diagonal part of the matrix
   for (icell=0;icell<ncell;icell++) {
      CelldblPnt = setdblcell(icell);

      for (idir=0;idir<ndir;idir++)
         Vmat[idir][idir] = CelldblPnt[CellVol];

      getnnz( &MatGrd, icell, icell, 1.0, Vmat);
   }

   //Loop among the faces
   for (iface=0;iface<nface;iface++) {
      FacedblPnt = setdblface(iface);
      FaceintPnt = setintface(iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      // Basic idea:
      // alfa ( GRAD_ACELL( phi) * DX_ACELL) + ( 1 - alfa) ( GRAD_BCELL( phi) * DX_BCELL)
      alfa  = FacedblPnt[FaceAlfa];

      for (idir=0;idir<ndir;idir++) { 
      for (jdir=0;jdir<ndir;jdir++) { 
         Rmat_a[idir][jdir] = FacedblPnt[FaceDxA+jdir];
         Rmat_b[idir][jdir] = FacedblPnt[FaceDxB+jdir];
      }
      }

      for (idir=0;idir<ndir;idir++) {
      for (jdir=0;jdir<ndir;jdir++) {
         Rmat_a[idir][jdir] *= FacedblPnt[FaceSrf+idir];
         Rmat_b[idir][jdir] *= FacedblPnt[FaceSrf+idir];
      }
      }

      getnnz( &MatGrd, acell, acell, -        alfa  * c_a, Rmat_a);

      getnnz( &MatGrd, acell, bcell, - (1.0 - alfa) * c_a, Rmat_b);

      getnnz( &MatGrd, bcell, acell, -        alfa  * c_b, Rmat_a);

      getnnz( &MatGrd, bcell, bcell, - (1.0 - alfa) * c_b, Rmat_b);
   }
}

// Iterative gradients calculation --------------------------------
void gradients( int nvar, double varPNT[ncell*nvar], \
                int grdSTART, int nnvar, double grdPNT[ncell*nnvar] )
{
   int icell, ivar, idir, iface, acell, bcell;
   int *FaceintPnt, ighst, *GhstintPnt;
   double *FacedblPnt, *phi_a, *phi_b, phi_i;
   double alfa, Rhs[ndir], *grads, *dep, *dep_g;
   struct RHS    b;

   // Mat Alloc for segregated gradients --------------------------
   rhsini( ncell, ndir, &b);

   // Gradient computation for all variables ----------------------
   for (ivar=0;ivar<nvar;ivar++) {
      rhsstart( 0.0, &b);

      grads = (double *) calloc( ndir*ncell, sizeof(double));

      //Loop among the faces
      for (iface=0;iface<nface;iface++) {
         FacedblPnt = setdblface(iface);
         FaceintPnt = setintface(iface);
   
         //A-B cells
         acell = FaceintPnt[FaceAcl];
         bcell = FaceintPnt[FaceBcl];
   
         phi_a = setsol(acell, nvar, varPNT);
         phi_b = setsol(bcell, nvar, varPNT);
   
         alfa  = FacedblPnt[FaceAlfa];

         //Variable at the face
         phi_i = facevar( alfa, phi_a[ivar], phi_b[ivar]);

         for (idir=0;idir<ndir;idir++)
            Rhs[idir] = phi_i * FacedblPnt[FaceSrf+idir];

         getrhs( &b, acell, c_a, Rhs);
         getrhs( &b, bcell, c_b, Rhs);
      }

      int    nitergrd = 1000;
      double tollgrd = 0.0001;
      SORmethd( &MatGrd, &b, grads, 1.0, &tollgrd, &nitergrd);
      //ConjGrad( &MatGrd, &b, grads, &tollgrd, &nitergrd);
      //BiCGStab( &MatGrd, &b, grads, &tollgrd, &nitergrd);

      fprintf(logfile, "    GradRecons %d/%d Iter = %4d, Res = %f\n", ivar+1, nvar, nitergrd, tollgrd);

      for (icell=0;icell<ncell;icell++) {
         dep = setsol(icell, nnvar, grdPNT);

         for (idir=0;idir<ndir;idir++)
            dep[grdSTART+(ivar*ndir)+idir] = grads[icell*ndir+idir];
      }

      free( grads); 

      //Defining Gradients for GHS cells -------------------------
      for (ighst=0;ighst<nghst;ighst++) {
         GhstintPnt = setintghst(ighst);
         icell      = GhstintPnt[GhstFidx];

         dep_g = setsol( ghst2gbl( ighst) , nnvar, grdPNT);
         dep   = setsol(            icell , nnvar, grdPNT);

         for (idir=0;idir<ndir;idir++)
            dep_g[grdSTART+(ivar*ndir)+idir] = \
                      dep[grdSTART+(ivar*ndir)+idir];
      }
   }

   rhsdel( &b);
}
