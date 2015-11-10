#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../BOUNDARY/boundary.h"
#include "turb.h"

//FLUXES ----------------------
void convflx(struct CSRMat *Mat, struct RHS *Rhs );
void viscflx(struct CSRMat *Mat, struct RHS *Rhs );
void sources(struct CSRMat *Mat, struct RHS *Rhs);

//pres-vel coupling subroutines
extern double *dd, *p_grd, *p_sol;

void presflx(double phi[nvar], double F_i[nvar][nvar]);
void timematrix(double phi[nvar], double dep[ndep], double F_i[nvar][nvar]);
void copyvar(double sol[ncell*nvar]);
void update(struct CSRMat *Mat, struct RHS *Rhs, double sol[ncell]);
void pattern(int pattern[nvar][nvar]);

//GRADIENTS --------------------
extern struct CSRMat MatGrd;

void gradinit( void);

void gradients( int nvar, double varPNT[ncell*nvar], \
                int grdSTART, int nnvar, double grdPNT[ncell*nnvar] );

//THERMODYNAMIC ----------------
void calcrho(void );
void calccfl(void );
void calcdis(void );
void calcmt(void );
double divergence( double dep[ndep]);
double drdp( double phi[nvar]);
double drdt( double phi[nvar]);
double enthalpy( double phi[nvar]);
double dhdp( double phi[nvar]);
double dhdt( double phi[nvar]);
#if _DBG == 10
void calcdbgdep(void );
#endif

//MATHEMATICS ------------------
void pntsdist( double xxa[ndir], double xxb[ndir], double dx[ndir]);
double facevar(double alfa, double var_a, double var_b);
