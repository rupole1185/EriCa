#include "marisa.h"

/*General initialization of Preconditioners*/
// AI and AJ taken from the Reference Matrix
// AA initialized to zero
void pcstrt(struct CSRMat *Mat)
{
   Mat->PrcPtr = malloc( sizeof( struct CSRMat));

   matini( Mat->dim * Mat->nvar, 1, Mat->PrcPtr);

   matdiag( Mat->PrcPtr);
}

void pcdel( struct CSRMat *Mat)
{
   if (Mat->PrcPtr == NULL)
      return;

   matdel( Mat->PrcPtr);
   Mat->PrcPtr = NULL;
}

/*Jacoby preconditioner initialization*/
/*Diagonal elements of Mat are copied into the Prc structure*/
void pcjacb(struct CSRMat *Mat)
{
   int     i, j, nvar, dim;
   double tmp[1][1];

   if (Mat->PrcPtr == NULL)
      pcstrt( Mat);

   nvar = Mat->nvar;
   dim  = Mat->dim;

   matstart( 0.0, 1, Mat->PrcPtr);

   for (i=0;i<dim*nvar;i++) {
      for (j=Mat->AI[i];j<Mat->AI[i+1];j++) {
         if (i==Mat->AJ[j])  {

            /* the inverse of the matrix is stored */
            tmp[0][0] = 1.0 / Mat->A[j];
            getnnz(Mat->PrcPtr,i,i,1.0,tmp);  }}}
}

/*Generation of the K1 term: Precon K = K1 * K2*/
void createK1( struct CSRMat *Prc)
{
   int     i, j, nvar, dim;
   double tmp[1][1];

   pcstrt( Prc);

   nvar = Prc->nvar;
   dim  = Prc->dim;

   for (i=0;i<dim*nvar;i++) {
      for (j=Prc->AI[i];j<Prc->AI[i+1];j++) {
         if (i==Prc->AJ[j])  {

            /* the sqrt inverse of the matrix is stored */
            tmp[0][0] = sqrt(Prc->A[j]);
            getnnz(Prc->PrcPtr,i,i,1.0,tmp);  }}}
}
