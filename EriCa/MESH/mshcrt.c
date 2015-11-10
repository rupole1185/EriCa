#include "mesh.h"
#include "meshset.h"

//Cell pointers
//integer
int CellVlm;
//double
int CellVol, CellXc;

//Ghost pointers
//integer
int GhstFidx, GhstRef, GhstType;
//double
int GhstDxAlfa, GhstDx, GhstAlfa, GhstBnd, \
    GhstWallPhi, GhstWallDPDn, GhstSrf, GhstDist;

//Face pointers
//integer
int FaceAcl, FaceBcl;
//double
int FaceSrf, FaceAlfa, FaceDxA, FaceDxB, FaceDxGrd, FaceAlfaGrd;

//MESH info: local to this module
double *celldbldata, *ghstdbldata, *facedbldata;
int    *cellintdata, *ghstintdata, *faceintdata;
int    ncelldbl, ncellint, nfacedbl, nfaceint, nghstdbl, nghstint;

void createmesh(void )
{
   int icell, istart, iend, j;
   double tmp[ndir];

   //MESH creation ----------------------------------------------
   //Creation of the stencil ------------------------------------
   createstencil();
   
   //SPEEDUPS AND CELL DEFINITION -------------------------------
   // Cell allocations
   cellintdata = (int *)    malloc(ncell * ncellint * sizeof(int));
   celldbldata = (double *) malloc(ncell * ncelldbl * sizeof(double));

   //To avoid slow code its is recomended to count and then
   //allocate! This saves a lot of memory movements!
   for (icell=0;icell<ncell;icell++) {
      istart = xadj[icell];
      iend   = xadj[icell+1];

      createcell(icell, voladj[icell], &xcadj[icell*ndir], vlmadj[icell], tmp);

      for (j=istart;j<iend;j++) {
         if (adjncy[j] < 0)
	    nghst++;

         if (adjncy[j] > icell || adjncy[j] < 0)
	    nface++;
      }
   }

   ntot = ncell + nghst;

   // Expanding Cell arrays
   cellintdata = realloc(cellintdata, ntot * ncellint * sizeof(int));
   celldbldata = realloc(celldbldata, ntot * ncelldbl * sizeof(double));

   // Ghst allocations
   ghstintdata = (int *)    malloc(nghst * nghstint * sizeof(int));
   ghstdbldata = (double *) malloc(nghst * nghstdbl * sizeof(double));

   // Face allocations 
   faceintdata = (int *)    malloc(nface * nfaceint * sizeof(int));
   facedbldata = (double *) malloc(nface * nfacedbl * sizeof(double));

   nghst = 0;
   nface = 0;

   //GHOST AND FACE DEFINITION ----------------------------------
   for (icell=0;icell<ncell;icell++){
      istart = xadj[icell];
      iend   = xadj[icell+1];

      for (j=istart;j<iend;j++) {

         //Ghost cells
         if (adjncy[j] < 0) 
           createghst(nghst++, &adjncy[j], icell, &xcfncy[j*ndir], &dirncy[j*ndir]);

         //Faces
         if (adjncy[j] > icell)
       	    createface(nface++, icell, adjncy[j], &xcfncy[j*ndir], &dirncy[j*ndir]);
      }
   }

#if _DBG == 100
   writemeshdata();
#endif

   //Final Operations: release memory ---------------------------
   free( dirncy);
   free( voladj);
   free(  xcadj);
   free( xcfncy);
   free( vlmadj);

   dnnzcount();
   ElmTypCSRstart();

   free( xadj);
   free( adjncy);
   free( nodadj);

   //VARIABLES ALLOCATION ---------------------------------------
   meshalloc();

   if (restart==0)
      readSOL();
}

//Allocating mesh variables --------------------------------
double *var, *dep;

void meshalloc(void )
{
   int icell, ivar;
   double *phi_i, *dep_i;

   //Var and Dep allocation
   var  = (double *) malloc(nvar*ntot * sizeof(double));
   dep  = (double *) malloc(ndep*ntot * sizeof(double));

   //Initialization of Var to USER defined VARINI
   for (icell=0;icell<ntot;icell++) {
      phi_i = setvar( icell);
      dep_i = setdep( icell);

      memcpy( phi_i, varini, nvar*sizeof(double));
      memcpy( dep_i, depini, ndepini*sizeof(double));
   }
}

//Allocating Matrix ---------------------------------------
int *dnnz;

void dnnzcount()
{
   int icell, istart, iend, j;
   dnnz   = (int *) calloc(ncell, sizeof(int));

   for (icell=0;icell<ncell;icell++) {
      dnnz[icell] = 1;

      istart = xadj[icell];
      iend   = xadj[icell+1];

      for (j=istart;j<iend;j++) {
         if (adjncy[j]<ncell) 
            dnnz[icell] ++;
      }
   }
}
