#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../BOUNDARY/boundary.h"

// CONSTANTS ------------------
static double kappa= 0.41, Cmu = 0.09;
//    KW06
static double sigk = 0.6, sigw = 0.5, Clim = 7.0 / 8.0;
static double alpha= 13.0 / 25.0, beta0 = 0.0708;
//    KWSST
static double sigk1 = 0.85, sigw1 = 0.5,   beta1 = 0.075;
static double sigk2 = 1.0,  sigw2 = 0.856, beta2 = 0.0828;
static double a1 = 0.31, gamma1, gamma2;
//    SA
static double cb1 = 0.1355, cb2 = 0.622, sigma = 2.0 / 3.0;
static double cw1, cw2 = 0.3, cw3 = 2.0, cv1 = 7.1; 
static double ct3 = 1.2, ct4 = 0.5;

// FUNCTIONS ---------------
void init_turb_const( void);
double blend(double alfa, double vara, double varb);
void stress( double phi[nvar], double dep[ndep], double tij[ndir][ndir]);
void strain( double dep[ndep], double Sij[ndir][ndir]);
double prodk( double dep[ndep], double tij[ndir][ndir]);
double matsum( double mat[ndir][ndir]);

// Kw routines
double CDkw( double phi[nvar], double dep[ndep]);
double F1( double phi[nvar], double dep[ndep]);
double F2( double phi[nvar], double dep[ndep]);

// SA routines
double fw( double phi[nvar], double dep[ndep]);
double Shat( double phi[nvar], double dep[ndep]);
double fv1( double phi[nvar], double dep[ndep]);
double fv2( double phi[nvar], double dep[ndep]);
double ft2( double phi[nvar], double dep[ndep]);
