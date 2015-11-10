#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// DEBUGGING MODE for SOLVER ---------------------
// 0 - Default mode
// 1 - Standard debugging mode
//-1 - Timing mode
#define _DBG 0

// SOLVER TYPE-----
// 0 -> LOCAL
// 1 -> PETSC
#define _SOLVER 0

#if _SOLVER == 1
#include <petsc.h>
#endif

// OPENMP -----------------------------------------
#ifdef _OPENMP
#include <omp.h>
#endif

/*! \addtogroup ALGEBRA 
 *  @{
 */

#ifndef _MAT_ALG
#define _MAT_ALG

//!Sparse matrix structure
typedef struct CSRMat {
     int  *AI, *AJ, *ptrn, *AIivar, *AJivar;
     double *A ; 
     int  dim, nvar, nnz, nptrn;

     struct CSRMat *PrcPtr;
#if _SOLVER == 1
     Mat PETScMat;
#endif
} SparseMat;

//!RHS structure
typedef struct RHS  {
     double *b ;
     int  dim, nvar;

#if _SOLVER == 1
     Vec PETScVec;
#endif
} RHSVector;
#endif

// MATRIX ALGEBRA ******************************************************
//   DENSE MATRIX / VECTORS -----------------------
//! Function to calculate the scalar product among two different arrays
double scalarp(int  dim, double a[dim], double b[dim]);

//! Subroutine to calculate the Matrix-Vector product
void matvecdstd(int  dim, double mat[dim][dim], double vec_in[dim], \
                                                double vec_out[dim]);

//! Summ of two vectors
void summvec(int  dim, double Vec1[dim], double Vec2[dim], double Sum[dim]);

//! Calculate the inverse of a given matrix (Pseudo inverse if rectangular)
void invmat(int m, int  n, double matrix[m][n]);

//! LU factorization
void LU( int n, double mat[n][n], double L[n][n] , double U[n][n]);

//   SPARSE MATRIX --------------------------------
//! Matrix vector product form
void matvecd(struct CSRMat *Mat, double Vec[Mat->dim*Mat->nvar], \
                                 double prod[Mat->dim*Mat->nvar]);

//! Subroutine to calculate the residual
void updtres(struct CSRMat *Mat, double Vec[Mat->dim*Mat->nvar], \
             double x[Mat->dim*Mat->nvar], double res[Mat->dim*Mat->nvar]);

// MAT MANAGER **********************************************************
//   DENSE MATRIX / VECTORS -----------------------
//! Returns the sum of the elements of a Vector
int  sveciel(int  dim, int  Vec[dim]);

//! Returns the sum of the elements of a Vector
double sveciel_dbl(int  dim, double Vec[dim]);

//   SPARSE MATRIX --------------------------------
//! CSR Matrix cardinalities
void matini(int  dim, int  nvar, struct CSRMat *Mat);

//! CSR Matrix pattern
void matptrn( struct CSRMat *Mat, int ptrn[Mat->nvar][Mat->nvar]);

/*CSR Matrix allocation*/
//! Diagonal
void matdiag( struct CSRMat *Mat);
//! User dnnz
void matalc( struct CSRMat *Mat, int dnnz[Mat->dim]);

//! CSR Matrix deallocation
void matdel(struct CSRMat *Mat);

//! CSR right_hand_side allocation
void rhsini(int  dim, int  nvar, struct RHS *B);

//! CSR right_hand_side deallocation
void rhsdel(struct RHS *B);

//! CSR initialization to zero
void matstart(double init, int opt, struct CSRMat *Mat); 

//! RHS initialization to zero
void rhsstart(double init, struct RHS *B);

//! CSR Matrix new element
void getnnz(struct CSRMat *Mat, int  row, int  col, \
            double coeff, double add[Mat->nvar][Mat->nvar]);

//! RHS vector element inserting
void getrhs(struct RHS *B, int  row, double coeff, double add[B->nvar]);

//!  Extracting row of given MATRIX 
void extrow( int  idim, int  ivar, int  *i, int  *ncol, \
             int  **col, double **elem, struct CSRMat *Mat);

//!  Extract element from the rhs 
double extrhs( int idim, int ivar, struct RHS *Rhs);

// QUERIES **************************************************************
//   DENSE MATRIX ---------------------------------
//   SPARSE MATRIX --------------------------------
extern int nprt;
static char FName[18+10];
const char * fname(char Name[10]);

//! Print MATLAB file for sparse matrix
void prtmat(struct CSRMat *Mat, struct RHS *B);

//! Subroutine to check the MATRIX
void checkmat(struct CSRMat *Mat);

//! Returns the number of bytes stored in Matrix
double MatMem(struct CSRMat *Mat);

//! Returns the number of bytes stored in RHS   
double RhsMem(struct RHS *B);

//! Returns the nnz of a Matrix
int Matnnz(struct CSRMat *Mat);

//! Returns the nconn of a Matrix
int Matnconn(struct CSRMat *Mat);

// PRECONDITIONING ******************************************************
//! Jacoby preconditioner initialization
void pcjacb(struct CSRMat *Mat);

//! Generation of the K1 term: Precon K = K1 * K2
void createK1( struct CSRMat *Prc);

// KIRSOLVER METHODS and SOR METHOD *************************************
//!  CG 
void ConjGrad(struct CSRMat *Mat, struct RHS *B, \
              double x[Mat->dim*Mat->nvar], double *tollr, int  *iterr);

//!  BCGSTAB 
void BiCGStab(struct CSRMat *Mat, struct RHS *B, \
              double x[Mat->dim*Mat->nvar], double *tollr, int  *iterr);

//!  SOR 
void SORmethd(struct CSRMat *Mat, struct RHS *B, \
              double x[Mat->dim*Mat->nvar], double omega, \
              double *tollr, int  *iterr);

/*! @}
 */

