#include "eqntns.h"

double *dd, *p_grd, *p_sol;

// Calculation of the D coefficient for pressure equation
void calcdd( struct CSRMat *Mat)
{
   int icell, Ncol, *col, irow, icol, j, idir;
   double *elem, *dd_i;

   dd     = (double *) calloc(ncell*ndir, sizeof(double));

   for (icell=0;icell<ncell;icell++) {
      dd_i = setsol( icell, ndir, dd);

      for (idir=0;idir<ndir;idir++) {
         extrow( icell, Uvar+idir, &irow, &Ncol, &col, &elem, Mat);

         for (j=0;j<Ncol;j++) {
            icol = col[j];

            if (icol==irow)
               dd_i[idir] += elem[j];
            else if (eos==ISIMPLEC)
               dd_i[idir] -= elem[j];
         }
#if _DBG == 10
         if (dd_i[idir] == 0.0) {
            printf("Error dd!!! icell =  %d, dir = %d\n", icell, idir);
            FatalError(" ");
         }
#endif
      }
   }
}

// Continuity equation is used to find out the pressure corrections
void simple_peqn( void)
{
   double c_a, c_c;
   int    iface, acell, bcell, *FaceintPnt, idir, nniterp;
   double *FacedblPnt, *AA, *CellAdblPnt, *CellBdblPnt, A_a, A_b;
   double  F_i[1], Fp_i[1][1], u_f[ndir], u_a[ndir], u_b[ndir];
   double *phi_a, *dep_a, *phi_b, *dep_b, rho_i, *dd_a, *dd_b, p_sum[ndir];
   double F_pa, F_pb, dx[ndir], dxi, e_xi[ndir], alfa, u_fp, ttollp, ttmmpp;
   struct CSRMat Mat_p;
   struct RHS    Rhs_p;

   p_sol = (double *) calloc( ntot, sizeof(double));

   c_a = 1.0;
   c_c = 1.0;

   matini(ncell, 1, &Mat_p);
   matalc(&Mat_p, dnnz);
   rhsini(ncell, 1, &Rhs_p);

// Rileggersi: Versteeg pagg. 338 e seguenti ...
//    geometrical issues on Collocated grids
//   
// Face interpolation has three different contributions here reported
// u_f = f(U_a,U_b) + f(P_a, P_b) + f(grad(P_a), grad(P_b))
//          1             2                 3
// The first therm is straight forward to compute and is implicit
// The second and the third are explicit at a given time step  

   for (iface=0;iface<nface;iface++) {
      FaceintPnt = setintface( iface);
      FacedblPnt = setdblface( iface);

      //A-B cells
      acell = FaceintPnt[FaceAcl];
      bcell = FaceintPnt[FaceBcl];

      alfa = FacedblPnt[FaceAlfa];
      AA   = &FacedblPnt[FaceSrf];

      //Variables of the two cells
      phi_a = setvar( acell);
      dep_a = setdep( acell);
      phi_b = setvar( bcell);
      dep_b = setdep( bcell);

      //Cell dbl Pnt
      CellAdblPnt = setdblcell( acell);
      CellBdblPnt = setdblcell( bcell);

      rho_i = facevar( alfa, dep_a[Rdep], dep_b[Rdep]);

      // These are the coefficients averaged in all the directions
      //    For ghosts cells the correspondent coefficient
      //    of the fluid cell is considered
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

      // u* contribution to continuity equation -------------------
      //  1st term 
      //   we use u* to obtain the pressure correction
      for (idir=0;idir<ndir;idir++) {
         u_a[idir] = phi_a[Uvar+idir] + scalarp( ndir, &dep_a[UGrd+ndir*idir], &FacedblPnt[FaceDxA]);
         u_b[idir] = phi_b[Uvar+idir] + scalarp( ndir, &dep_b[UGrd+ndir*idir], &FacedblPnt[FaceDxB]);
      }

      for (idir=0;idir<ndir;idir++)
         u_f[idir] = facevar( alfa, u_a[idir], u_b[idir]);

      F_i[0] = rho_i * scalarp( ndir, AA, u_f);

      //  2nd and 3rd therms: face normal velocities from RhieChow
      F_pa = CellAdblPnt[CellVol] / A_a;
      F_pb = CellBdblPnt[CellVol] / A_b;

      pntsdist( &CellBdblPnt[CellXc], &CellAdblPnt[CellXc], dx);

      dxi = sqrt( scalarp( ndir, dx, dx) );

      for (idir=0;idir<ndir;idir++)
         e_xi[idir] = dx[idir] / dxi;

      double dxn = scalarp( ndir, dx, AA) / sqrt( scalarp( ndir, AA, AA));

      double p_a = phi_a[Pvar] + scalarp( ndir, &dep_a[PGrd], &FacedblPnt[FaceDxA]);
      double p_b = phi_b[Pvar] + scalarp( ndir, &dep_b[PGrd], &FacedblPnt[FaceDxB]);

      /////////////////////////////////////////// VERSTEEG
      u_fp  = facevar( alfa, F_pa, F_pb) * (p_b - p_a) / dxn;

      for (idir=0;idir<ndir;idir++)
         p_sum[idir] = facevar( alfa, F_pa * dep_a[PGrd+idir], F_pb * dep_b[PGrd+idir]);
         //p_sum[idir] = facevar( alfa, F_pa, F_pb) * facevar( 0.5, dep_a[PGrd+idir], dep_b[PGrd+idir]);

      u_fp -= scalarp( ndir, p_sum, e_xi);
      /////////////////////////////////////////// DARWISH
      /*for (idir=0;idir<ndir;idir++)
         p_sum[idir] = -1.0 * facevar( 0.5, dep_a[PGrd+idir], dep_b[PGrd+idir]);

      for (idir=0;idir<ndir;idir++)
         p_sum[idir]+= ( p_b - p_a) / dxn * e_xi[idir];

      u_fp = facevar( alfa, F_pa, F_pb) * scalarp( ndir, p_sum, AA) / sqrt( scalarp( ndir, AA, AA));*/
      ///////////////////////////////////////////

      // There is a sign inconsistency between:
      //    MATHUR AND MURTHY's "Pressure based method for unstructured mesh"
      //    VERSTEEG's "An Introduction to C.F.D."
      F_i[0] -= rho_i * u_fp * sqrt( scalarp( ndir, AA, AA));

      // A CELL fluxes ----------
      getrhs( &Rhs_p, acell,  c_a, F_i);

      // B CELL fluxes ----------
      getrhs( &Rhs_p, bcell, -c_a, F_i);

      // p' contribution to continuity equation -------------------
      // We need to know the sum of all the convective/diffusion
      // terms for the velocity
      double Cnst = rho_i * facevar( alfa, F_pa, F_pb);

      Fp_i[0][0]  = Cnst * scalarp( ndir, AA, AA) / scalarp( ndir, dx, AA);

      // A CELL fluxes ----------          //LE BC qui sono DIFFERENTI!!!
      getnnz( &Mat_p, acell, acell, -c_c, Fp_i);

      if (bcell<ncell) 
         getnnz( &Mat_p, acell, bcell, c_c, Fp_i);
      else 
         getnnz( &Mat_p, acell, acell, c_c * ghst2fld( gbl2ghst( bcell), Pvar), Fp_i);

      // B CELL fluxes ----------
      if (acell<ncell)
         getnnz( &Mat_p, bcell, acell, c_c, Fp_i);
      else 
         getnnz( &Mat_p, bcell, bcell, c_c * ghst2fld( gbl2ghst( acell), Pvar), Fp_i);

      getnnz( &Mat_p, bcell, bcell, -c_c, Fp_i);
   }

#if _DBG == 10
   int icell;
   double *dep, *dd_i;

   for (icell=0;icell<ncell;icell++) {
      dep = setdep( icell);

      dep[DivSdep] = extrhs( icell, 0, &Rhs_p);

      dd_i = setsol( icell, ndir, dd);

      for (idir=0;idir<ndir;idir++)
         dep[DDdep+idir] = dd_i[idir];
   }
#endif

   // Solution -------------
   nniterp = niterp;
   ttollp  = tollp;

//   BiCGS seems to be the best method among the ones implemented ...
//   Cases in future may confirm this assumption ...
//   SORmethd( &Mat_p, &Rhs_p, p_sol, 1.0, &ttollp, &nniterp);
//   pcjacb( &Mat_p);
//   ConjGrad( &Mat_p, &Rhs_p, p_sol, &ttollp, &nniterp);
   BiCGStab( &Mat_p, &Rhs_p, p_sol, &ttollp, &nniterp);

   fprintf(logfile, "    Press iter   = %4d, Res = %f\n", nniterp,  ttollp);

   matdel( &Mat_p);
   rhsdel( &Rhs_p);
}

void simple_bc( void)
{
   int ighst, *GhsintPnt, icell;

   for (ighst=0;ighst<nghst;ighst++) {
      GhsintPnt = setintghst(ighst);
      icell     = GhsintPnt[GhstFidx];
      
      p_sol[ghst2gbl(ighst)] = p_sol[icell] * ghst2fld( ighst, Pvar);
   }

   // p' gradients
   p_grd = (double *) calloc( ntot*ndir, sizeof(double));

   gradients( 1, p_sol, 0, ndir, p_grd);
}

void simple_upd(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell*nvar])
{
   int icell, ivar, idir;
   double *phi, *dep, *pgrd, *CelldblPnt, *dd_i;

   // Solution update -------------------
   // (all variables except pressure)
   for(icell=0;icell<ncell;icell++) {
      phi = setvar( icell);

      memcpy( &phi[Uvar], &sol[icell*nvar+Uvar], (nvar-1)*sizeof(double));
   }

   if ( eos == IFROZENP)
      return;

   // Not sure if these two instructions are relevant ...
   boundary();
   gradients( nvar, setvar( 0), PGrd, ndep, setdep( 0));

   // Velocity coefficients -------------
   calcdd( Mat);

   // Pressure equation resolution ------
   simple_peqn();

   // SIMPLE BC -------------------------
   simple_bc();

   // PresVel corrections ---------------
   for(icell=0;icell<ncell;icell++)  {
      phi        = setvar(icell);
      CelldblPnt = setdblcell(icell);
  
      //Pressure update P = P + p'
      phi[Pvar] += 0.5 * p_sol[icell];

      //Velocity update U = u* + u'
      dd_i = setsol( icell, ndir, dd);
      pgrd = setsol( icell, ndir, p_grd);

      for (idir=0;idir<ndir;idir++)
         phi[Uvar+idir] -= CelldblPnt[CellVol] * pgrd[idir] / dd_i[idir];

#if _DBG == 10
      // Adding P' and U' to CEI plot
      dep        = setdep(icell);
      dep[Ppdep] = p_sol[icell];
      for (idir=0;idir<ndir;idir++)
         dep[Updep+idir] = -CelldblPnt[CellVol] * pgrd[idir] / dd_i[idir];
#endif
   }

   // Memory release --------------------
   free(p_sol);
   free(p_grd);
   free(dd);
}

// Pressure Equation is solved separately from the main
// system. This is why all the variables are set to zero
// except the TIME matrix [for computational reason]
void simple_conf(double phi[nvar], double F_i[nvar][nvar])
{
   phi[Pvar]       = 0.0;
   F_i[Pvar][Pvar] = 0.0; 
}

// Time step advancing matrix is identity for SIMPLE method
void simple_tmmt( double phi[nvar], double dep[ndep], double F_i[nvar][nvar])
{
   int ivar;
   
   for (ivar=0;ivar<nvar;ivar++) 
      F_i[ivar][ivar] = dep[Rdep];
}
