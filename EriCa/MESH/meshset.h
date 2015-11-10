#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../SIMULATION/simulation.h"

//MSHCRT ---------------------------------------
//MESH GRAPH
extern int *xadj, *adjncy, *vlmadj, *nodadj;
extern double *dirncy, *voladj, *xcadj, *xcfncy;
extern double gspc[ndir];

//MESH info: local to this module
extern double *celldbldata, *facedbldata, *ghstdbldata;
extern int    *cellintdata, *faceintdata, *ghstintdata;

//MESH variables
extern double *var, *dep;

//MSHGRPH -------------------------------------
void meshalloc(void );
void dnnzcount(void );

void createcell(int icell, double vol, double xc[ndir], int volmod, double displ[ndir]);
void createghst(int ighst, int *conn, int icell, double xcf[ndir], double norm[ndir]);
void createface(int iface, int acell, int bcell, double xcf[ndir], double norm[ndir]);  

void createstencil(void );
double givet( double norm[ndir], double xc[ndir], double xcf[ndir]);
void setintP( double Dx[ndir], double xc[ndir], double xcf[ndir], \
              double t, double norm[ndir]);
//ENSIGHT GEOMETRY ---------------------------
void ensvar(int count, char *String,            \
            int ndim, int ivar, int vd,         \
            double (*funcp)( double *, double *));
