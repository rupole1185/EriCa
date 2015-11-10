#include "mesh.h"
#include "meshset.h"

//Geometric calculation for faces --------------------------
double givet( double norm[ndir], double xc[ndir], double xcf[ndir])
{
   double xx = 0.0;
   int idir;

   for (idir=0;idir<ndir;idir++)
      xx += (xc[idir] - xcf[idir]) * norm[idir];

   return xx / scalarp( ndir, norm, norm);
}

void setintP( double Dx[ndir], double xc[ndir], double xcf[ndir], \
              double t, double norm[ndir])
{
   int idir;

   for (idir=0;idir<ndir;idir++)
      Dx[idir] = xcf[idir] + t * norm[idir] - xc[idir];
}

//Counting    cardinalities --------------------------------
void createstencil(void )
{
   //Initialization of ghost cell number and faces
   nghst  = 0;
   nface  = 0;

   //Cell properties -----------------
   //double
   ncelldbl = 0;
   CellXc   = ncelldbl;
   ncelldbl+= ndir;
   CellVol  = ncelldbl++;

   //integer
   ncellint = 0;
   CellVlm  = ncellint++;

   //Ghst properties ----------------
   //double
   nghstdbl = 0;
   GhstAlfa = nghstdbl++;
   GhstBnd  = nghstdbl;
   nghstdbl+= nvar;
   GhstDx   = nghstdbl;
   nghstdbl+= ndir;
   GhstDxAlfa = nghstdbl++;
   GhstWallPhi = nghstdbl;
   nghstdbl+= nvar;
   GhstWallDPDn = nghstdbl;
   nghstdbl+= nvar;
   GhstSrf  = nghstdbl;
   nghstdbl+= ndir;
   GhstDist = nghstdbl++;

   //integer
   nghstint = 0;
   GhstFidx = nghstint++;
   GhstRef  = nghstint++;
   GhstType = nghstint;
   nghstint+= nvar;

   //Face properties -----------------
   //double
   nfacedbl = 0;
   FaceSrf  = nfacedbl;
   nfacedbl+= ndir;
   FaceAlfa = nfacedbl++;
   FaceDxA  = nfacedbl;
   nfacedbl+= ndir;
   FaceDxB  = nfacedbl;
   nfacedbl+= ndir;
   FaceDxGrd= nfacedbl;
   nfacedbl+= ndir;
   FaceAlfaGrd = nfacedbl++;

   //integer
   nfaceint = 0;
   FaceAcl  = nfaceint++;
   FaceBcl  = nfaceint++;
}

void writemeshdata(void )
{
   int icell, j, istart, iend, *CelintPnt, *FaceintPnt, iface, idir;
   double *CeldblPnt, *locvar, *FacedblPnt;
   FILE *MeshPnt;

   // -------------------------------------------------------
   FILE *filePnt;
   filePnt = fopen("DEBUG/Fxc.dat", "w");
 
   for (icell=0;icell<ncell;icell++) {
      CeldblPnt = setdblcell(icell);

      for (idir=0;idir<ndir;idir++)
         fprintf(filePnt," %lf", CeldblPnt[CellXc+idir]);

      fprintf(filePnt," \n");
   }
   fclose(filePnt);

   filePnt = fopen("DEBUG/Gxc.dat", "w");

   for (icell=ncell;icell<ntot;icell++) {
      CeldblPnt = setdblcell(icell);

      for (idir=0;idir<ndir;idir++)
         fprintf(filePnt," %lf", CeldblPnt[CellXc+idir]);

      fprintf(filePnt," \n");
   }
   fclose(filePnt);

   MeshPnt = fopen("DEBUG/MeshInfo.dat","w");

   fprintf(MeshPnt,"MESH data for current simulation:\n");
   fprintf(MeshPnt,"Current mesh FLD = %d\n", ncell);
   fprintf(MeshPnt,"Current mesh GHS = %d\n", nghst);
   fprintf(MeshPnt,"Current mesh FAC = %d\n", nface);
   fprintf(MeshPnt,"\n\n");

   fprintf(MeshPnt,"CELL data for current simulation:\n\n");

   // -------------------------------------------------------
   for (icell=0;icell<ncell;icell++) {
      CelintPnt = setintcell(icell);
      CeldblPnt = setdblcell(icell);

      fprintf(MeshPnt,"---------------------------------------------\n");
      fprintf(MeshPnt,"             Data set - Cell[%d]\n", icell);
      for (idir=0;idir<ndir;idir++)
         fprintf(MeshPnt,"Xc[%d] = %lf\n", idir+1, CeldblPnt[CellXc+idir]); 

      fprintf(MeshPnt,"Vol= %lf\n", CeldblPnt[CellVol]); 

      istart = xadj[icell];
      iend   = xadj[icell+1];

      fprintf(MeshPnt,"Nearby cell # = %d\n List: \n", iend-istart);
      for (j=istart;j<iend;j++) {
         fprintf(MeshPnt,"\n Cell= %d \n Norm = ", adjncy[j]);
         for (idir=0;idir<ndir;idir++)
            fprintf(MeshPnt," %e", dirncy[j*ndir+idir]);
      }
      fprintf(MeshPnt,"\n\n\n");
   }

   fprintf(MeshPnt,"FACE data for current simulation:\n\n");

   for (iface=0;iface<nface;iface++) {
      FaceintPnt = setintface(iface);
      FacedblPnt = setdblface(iface);
      
      fprintf(MeshPnt,"---------------------------------------------\n");
      fprintf(MeshPnt,"            Iface = %d\n", iface);
      fprintf(MeshPnt,"Acell = %d, Bcell = %d\n", FaceintPnt[FaceAcl],  \
                            FaceintPnt[FaceBcl]);

      for (idir=0;idir<ndir;idir++)
         fprintf(MeshPnt,"AA[%d] = %f\n", idir+1, FacedblPnt[FaceSrf+idir]);
   }

   fclose(MeshPnt);
}
