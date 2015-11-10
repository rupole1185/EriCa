#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <maya.h>
#include <sys/types.h>
#include <sys/stat.h>

// COMPILATION PARAMETERS -----------------

// 1D - 2D - 3D
#define ndir 3

// DEBUGGING COMPILATION
// _DBG =  0 standard  mode
// _DBG = 10 standard debugging
// _DBG =100 debugging matrixes
#define _DBG 10

// SIMINI --------------------------------
extern int    ncell, nghst, ntot, nface;

// Values of turb
#define INVISCID     0
#define LAMINAR_b   49
#define LAMINAR     50
#define DISTANCE_b  89
#define DISTANCE    95
#define TURBULENT_b 99
#define SA         110
#define KW06       201
#define KWSST      210

// Values of eos
#define ISIMPLE      0
#define ISIMPLEC     5
#define IPISO       10
#define ICOUPLED    50
#define IFROZENP    95
#define COMPR_b     99
#define IDBIG      100

extern int  turb, eos, restart, *resvar;
extern int  nvar, ndep, ndepini;
extern int  Pvar, Uvar, Tvar, Kvar, OMvar, MTvar, SCAvar, NSCAvar;
extern int  Rdep, MUdep, MTdep, Ddep, SCHMIDTdep;
extern int  PGrd, UGrd, TGrd, KGrd, OMGrd, MTGrd, SCAGrd, CFLdep, TIMEdep;
#if _DBG == 10
extern int  Ppdep, Updep, Divdep, Pkdep, Vrtdep, DivSdep, PHIDdep, DDdep;
#endif
extern double *varini, *depini, Cp, kfld, R;

void fldini(void );

extern int  iconv, ialgh, iprec;
extern double theta, timst, toll, tollp, cfl;
extern int  niter, niterp, ntmst, nplot, nranks, nwrite;

void simini(void );

extern int    *bndtype, nbnd;
extern double *bnd;
extern int    *monsurf, *monphi, *monave, *monorm, *monare, nmon;
extern double *mon;

void bndini(void );

// LOG WRITE ----------------------------
extern FILE *logfile;
extern clock_t timei;
int ScreenSize, LogSize;

void FatalError( char *message);

//TOOLS ---------------------------------
const char* foptions( int ichoice);
const char* soptions( int ichoice);
const char* boptions( int ichoice);
const char* VarName( int ivar);
const char* DepName( int idep);
