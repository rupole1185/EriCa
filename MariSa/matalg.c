#include "marisa.h"

//Subroutine to calculate the Matrix-Vector product-------------------
/*Matrix is in SPARSE form*/
void matvecd( struct CSRMat *Mat, double Vec[Mat->dim*Mat->nvar], \
                                  double prod[Mat->dim*Mat->nvar])
{ 
   int     i,ja,dim,nvar,j;
   double  prodi;

   nvar = Mat->nvar;
   dim  = Mat->dim;

#pragma omp parallel for default(none) private( i, ja, prodi) \
   shared(dim, nvar, Mat, prod, Vec)
   for (i=0;i<dim*nvar;i++)  {
      prodi = 0.0;

      for (ja=Mat->AI[i];ja<Mat->AI[i+1];ja++) /*Scalar product*/
         prodi += Mat->A[ja] * Vec[Mat->AJ[ja]]; 

      prod[i] = prodi;
   }
}

//Subroutine to calculate the Vec - Matrix * x ---------
void updtres( struct CSRMat *Mat, double Vec[Mat->dim*Mat->nvar], \
              double x[Mat->dim*Mat->nvar], double res[Mat->dim*Mat->nvar])
{ 
   int     i, dim;

   dim = Mat->dim * Mat->nvar;

   matvecd( Mat, x, res);

   for (i=0;i<dim;i++)
      res[i] = Vec[i] - res[i];
}
