#include "../SIMULATION/simulation.h"
#include "marisa.h"

//MSHCRT ---------------------------------------
extern int *dnnz;

//Cell pointers
extern int ncellint, ncelldbl;
//Ghst pointers
extern int nghstint, nghstdbl;
//Face pointers
extern int nfaceint, nfacedbl;

//Cell pointers
//integer  VolumeModel
extern int CellVlm;
//double   CellVolume, CellCentre (x,y,z)
extern int CellVol   , CellXc;

//Ghost pointers
//integer  FldCellIDX, Reference, BCType (0 Dirichlet, 1 Newman)
extern int GhstFidx  ,   GhstRef, GhstType;
//double               FldCellCentre to Projection to FaceNormal 
extern int GhstDxAlfa, GhstDx, \
//             Alfa,   BCVar,    Wall Var,   Wall Flux, Surf. Area
           GhstAlfa, GhstBnd, GhstWallPhi, GhstWallDPDn, GhstSrf, GhstDist;

//Face pointers
//integer  Acell  , Bcell
extern int FaceAcl, FaceBcl;
//double   Area (xyz), Alfa,  Dxof Acell, Dxof Bcell
extern int FaceSrf, FaceAlfa, FaceDxA, FaceDxB, FaceDxGrd, FaceAlfaGrd;

//MSH generation -------------------------------
void readmesh( void);
void meshcreate(void );

//MSHPNT ---------------------------------------
//Set Face Pointers
int* setintface(int iface);
double* setdblface(int iface);

//Set Cell Pointers
int* setintcell(int icell);
double* setdblcell(int icell);

//Set ghost Pointers
int* setintghst(int ighost);
double* setdblghst(int ighost);

//Set variables Pointers
double* setvar(int icell);
double* setdep(int icell);
double* setsol(int icell, int nnvar, double *pntr);

//ENSIGHT print ---------------------
void printCEI( int itmst, int ntmst);
void writeSOL( int itmst);
