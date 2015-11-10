#include "mesh.h"
#include "meshset.h"

int *xadj, *adjncy, *vlmadj, *nodadj;
double *dirncy, *voladj, *xcadj, *xcfncy, gspc[ndir];

// INPUT FILE LAYOUT -----------------------------------------
// NDIR
// Ncell
// XADJ + ADJNCY
// XCADJ
// DIRNCY
// XCFNCY
// VOLADJ
// VLMADJ
// NODADJ

//CREATE graph -----------------------------------------------
// Descriptions of basic values:
// CELL DATA ---------------------
// vlmadj[ncells] = cell volume group
// voladj[ncells] = cell volumes
// xcadj[ndir*ncell] = cell centres
// xadj[ncells+1] = CSR structure
// CONNECTIVITY DATA -------------
// nconn              = xadj[ncell+1] - xadj[0]
// adjncy[nconn]      = connected cells
// dirncy[ndir*nconn] = face normals
// xcfncy[ndir*nconn] = face centres

void readmesh(void )
{
   int  ierr;
   FILE *GeoFile;

   GeoFile = fopen("SYSTEM/Mesh.erc","rb");

   if (GeoFile == NULL)
      FatalError("No mesh file found!");

#if _DBG == 100
   FILE *GeoOut;
   GeoOut = fopen("DEBUG/Mesh.EriCa.erc","w");
#endif

   // NDIR -------------------------------------------
   int nndir;
   ierr = fread( &nndir, sizeof(int), 1, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "NDIR\n%d\n", nndir);
#endif

   // NCELLS -----------------------------------------
   ierr = fread( &ncell, sizeof(int), 1, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "NCELLS\n%d\n", ncell);
#endif

   // XADJ, ADJNCY -----------------------------------
   xadj = (int *) malloc((ncell+1) * sizeof(int));
   ierr = fread( xadj, sizeof(int), ncell+1, GeoFile);

   int nconn;
   nconn  = xadj[ncell] - xadj[0];
   adjncy = (int *) malloc(nconn * sizeof(int));
   ierr   = fread( adjncy, sizeof(int), nconn, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "XADJ\n");
   int icell, iconn, idir;
   for (icell=0;icell<ncell+1;icell++)
      fprintf( GeoOut, "%d ", xadj[icell]);

   fprintf( GeoOut, "\nADJNCY\n");
   for (iconn=xadj[0];iconn<xadj[ncell];iconn++)
      fprintf( GeoOut, "%d ", adjncy[iconn]);
#endif

   // XCADJ ------------------------------------------
   xcadj = (double *) malloc(nndir*ncell * sizeof(double));
   ierr  = fread( xcadj, sizeof(double), nndir*ncell, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "\nXCADJ\n");

   for (icell=0;icell<ncell;icell++) {
      for (idir=0;idir<nndir;idir++)
         fprintf( GeoOut, "%f ", xcadj[icell*nndir+idir]);
      fprintf( GeoOut, "\n");
   }
#endif

   // DIRNCY -----------------------------------------
   dirncy = (double *) malloc(nndir*nconn * sizeof(double));
   ierr   = fread( dirncy, sizeof(double), nndir*nconn, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "DIRNCY\n");

   for (iconn=xadj[0];iconn<xadj[ncell];iconn++) {
      for (idir=0;idir<nndir;idir++)
         fprintf( GeoOut, "%f ", dirncy[iconn*nndir+idir]);
      fprintf( GeoOut, "\n");
   }
#endif

   // XCFNCY -----------------------------------------
   xcfncy = (double *) malloc(nndir*nconn * sizeof(double));
   ierr   = fread( xcfncy, sizeof(double), nndir*nconn, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "XCFNCY\n");

   for (iconn=xadj[0];iconn<xadj[ncell];iconn++) {
      for (idir=0;idir<nndir;idir++)
         fprintf( GeoOut, "%f ", xcfncy[iconn*nndir+idir]);
      fprintf( GeoOut, "\n");
   }
#endif

   // VOLADJ -----------------------------------------
   voladj = (double *) malloc(ncell * sizeof(double));
   ierr   = fread( voladj, sizeof(double), ncell, GeoFile);

#if _DBG == 100
   fprintf( GeoOut, "VOLADJ\n");

   for (icell=0;icell<ncell;icell++)
      fprintf( GeoOut, "%f\n", voladj[icell]);
#endif

   // VLMADJ -----------------------------------------
   vlmadj = (int *) malloc(ncell * sizeof(int));
   ierr   = fread( vlmadj, sizeof(int), ncell, GeoFile);
   
#if _DBG == 100
   fprintf( GeoOut, "VLMADJ\n");

   for (icell=0;icell<ncell;icell++)
      fprintf( GeoOut, "%d\n", vlmadj[icell]);
#endif

   if (feof( GeoFile))
      FatalError(" Old Mesh version!");

   // NODADJ -----------------------------------------
   nodadj = (int *) malloc(ncell * sizeof(int));
   ierr   = fread( nodadj, sizeof(int), ncell, GeoFile);
   
#if _DBG == 100
   fprintf( GeoOut, "NODADJ\n");

   for (icell=0;icell<ncell;icell++)
      fprintf( GeoOut, "%d\n", nodadj[icell]);
#endif

   fclose( GeoFile);
#if _DBG == 100
   fclose( GeoOut);
#endif
}

//Create Face ----------------------------------------------
void createface(int iface, int acell, int bcell,      \
                double xcf[ndir], double norm[ndir]) 
{
   int *FaceintPnt;
   double *FacedblPnt;

   FaceintPnt = setintface(iface);
   FacedblPnt = setdblface(iface);

   // Face 2 Cell connectivity
   FaceintPnt[FaceAcl] = acell;
   FaceintPnt[FaceBcl] = bcell;

   // Face normal/Area
   memcpy(&FacedblPnt[FaceSrf], norm, ndir*sizeof(double));

   // This part calculate the face-projections of
   // the two cell centres -> non orthogonal meshes

   // Norm is the Area and not its normal but
   // this does not affect the results as
   // t is obtained as function of the Area 
   // instead of Norm
   double ta, tb, *CellAdblPnt, *CellBdblPnt;

   CellAdblPnt = setdblcell(acell);
   CellBdblPnt = setdblcell(bcell);

   // P = Xcf + t * norm  (P -> point where we ought to interpolate)
   ta = givet( &FacedblPnt[FaceSrf], &CellAdblPnt[CellXc], xcf);
   // N.B. tb is always positive because the Norm goes outside cellA
   tb = givet( &FacedblPnt[FaceSrf], &CellBdblPnt[CellXc], xcf);

   // Remember "regola della leva inversa"
   FacedblPnt[FaceAlfa] = tb/(tb-ta);

   // Dx = P - Xc = Xcf + t * norm - Xc
   setintP( &FacedblPnt[FaceDxA], &CellAdblPnt[CellXc], xcf, ta, &FacedblPnt[FaceSrf]);
   setintP( &FacedblPnt[FaceDxB], &CellBdblPnt[CellXc], xcf, tb, &FacedblPnt[FaceSrf]);

   // Gradient calculation variables ------------
   FacedblPnt[FaceAlfaGrd] = 0.5;

   int idir;
   for (idir=0;idir<ndir;idir++)
      FacedblPnt[FaceDxGrd+idir] = xcf[idir] - 0.5 * (CellAdblPnt[CellXc+idir] + CellBdblPnt[CellXc+idir]);

   if (bcell>=ncell) {
      FacedblPnt[FaceAlfaGrd] = 0.0;

      for (idir=0;idir<ndir;idir++)
         FacedblPnt[FaceDxGrd+idir] = 0.0;
   }
}

//Create Cell ----------------------------------------------
void createcell(int icell, double vol, double xc[ndir],    \
                           int volmod, double displ[ndir] )
{
   int    *CellintPnt;
   double *CelldblPnt;

   CellintPnt = setintcell(icell);
   CelldblPnt = setdblcell(icell);

   // VolModel
   CellintPnt[CellVlm] = volmod;

   // Cell centre
   memcpy(&CelldblPnt[CellXc], xc, ndir*sizeof(double));

   // Cell volume
   CelldblPnt[CellVol] = vol;

   if (icell<ncell) return;

   // GhostCells have a displacement to apply
   memcpy(&CelldblPnt[CellXc], displ, ndir*sizeof(double));
}

//Create GhostCell -----------------------------------------
void createghst(int ighst, int *conn, int icell,   \
                double xcf[ndir], double norm[ndir])
{
   int *GhstintPnt, surf, *CellintPnt;
   double *GhstdblPnt, *CelldblPnt;

   surf = -*conn - 100;

   //Adding the ghost cell to the stencil
   *conn = ighst + ncell;

   GhstintPnt = setintghst(ighst);
   GhstdblPnt = setdblghst(ighst);

   GhstintPnt[GhstFidx] = icell;
   GhstintPnt[GhstRef]  = surf-1;

   //Copying BND types and values in the Ghst arrays
   memcpy(&GhstintPnt[GhstType],&bndtype[(surf-1)*nvar],nvar*sizeof(int));
   memcpy(&GhstdblPnt[GhstBnd] ,&bnd[(surf-1)*nvar]    ,nvar*sizeof(double));

   //Geometric quantities: the ghost cell is orthogonal to the face
   //                      and with the same volume of the GhstFcell
   CellintPnt = setintcell( icell);
   CelldblPnt = setdblcell( icell);

   // Ghost cell is orthogonal to BND face and
   // its volume is the same as the GhstFcell
   double xcghst[ndir], t;
   int idir;

   // XcGhst Opposed to Fluid Centre Projection
   // Old definition of t ...
   //t = - 1.0 * givet( norm, &CelldblPnt[CellXc], xcf);
   t = 0.0;

   for (idir=0;idir<ndir;idir++)
      xcghst[idir] = t * norm[idir] + xcf[idir];

   double volghst = CelldblPnt[CellVol];

   createcell( *conn, volghst, &CelldblPnt[CellXc], 0, xcghst);

   // Calculation of Alfa => same as for faces
   // P = Xcf + t * norm  (P -> point where we ought to interpolate)
   double ta = givet( norm, &CelldblPnt[CellXc], xcf);
   // N.B. tb is always positive because the Norm goes outside cellA
   double tb = givet( norm, xcghst, xcf);

   // N.B. by FaceConvention and considering Ghst cells are always
   //      the BCELL of the face stencil, alfa is considered from
   //      the fluid cell
   GhstdblPnt[GhstAlfa] = tb/(tb-ta);

   setintP( &GhstdblPnt[GhstDx], &CelldblPnt[CellXc], xcf, ta, norm);

   // NB: issue with sign conventions. This number is negative because
   //     usually the flux is considered positive when going inside the system
   //     But this convention may not be general ...
   GhstdblPnt[GhstDxAlfa] = ta * sqrt( scalarp( ndir, norm, norm));

   // STORING inward normal (from fld to wall)
   for (idir=0;idir<ndir;idir++)
      GhstdblPnt[GhstSrf+idir] = norm[idir];

   // Distance between the two interpolation points
   GhstdblPnt[GhstDist] = (tb - ta) * sqrt( scalarp( ndir, norm, norm));

#if _DBG == 10
   //Temporary: ghost cells should know all their BC!!!
   if (surf > nbnd) {
      printf("\nError: not enough BND defined %d %d\n", surf, nbnd);
      printf("Fluid Cell # = %d\n", icell);
      FatalError("Problem in BC generation");
   }
#endif
}
