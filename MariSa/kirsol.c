#include "marisa.h"
#include <fenv.h>

double *res, *p, *Ap, *z, *res_old, *z_old;

/* Conjugate Gradient */
/* INPUT: Mat, RHS, x0, Toll_max, Niter
  OUTPUT: x, Toll_fin, Niter_fin         */
void ConjGrad(struct CSRMat *Mat, struct RHS *B, double x[Mat->dim*Mat->nvar], double *tollr, int  *iterr)
{
   int  i, iter = 0, dim;
   double toll, alpha, beta;

#if _DBG == 1
   checkmat( Mat);
   prtmat( Mat, B);

   FILE *DbgPnt;

   DbgPnt = fopen(fname("SolCon"), "w");
   fprintf(DbgPnt, "Iterations Residual\n");
#endif

   /*This is the real dimension of the problem*/
   dim = Mat->dim * Mat->nvar;

   res = (double *) malloc(dim * sizeof(double)); 
   p   = (double *) malloc(dim * sizeof(double)); 
   Ap  = (double *) malloc(dim * sizeof(double)); 
   z   = (double *) malloc(dim * sizeof(double)); 
   res_old = (double *) malloc(dim * sizeof(double)); 
   z_old   = (double *) malloc(dim * sizeof(double)); 

   if (B->dim!=Mat->dim || B->nvar!=Mat->nvar) {
      printf("Error Bdim != Matdim!!\n");
      printf(" Matdim %d Matnvar %d\n", Mat->dim, Mat->nvar);
      printf(" B  dim %d B  nvar %d\n", B->dim, B->nvar);
      return; }

   /*initial residual*/
   updtres( Mat, B->b, x, res);

   if ( scalarp( dim, res, res) != 0.0) {

      if (Mat->PrcPtr != NULL) 
         matvecd( Mat->PrcPtr, res, z);
      else
         memcpy(z, res, dim* sizeof(double)); 

      memcpy(p, z, dim* sizeof(double)); 

       /*iterative loop*/
      do {
         matvecd( Mat, p, Ap);

         alpha = scalarp(dim, res, z) / scalarp(dim, p, Ap);
 
         memcpy(res_old, res, dim* sizeof(double)); 

         for (i=0;i<dim;i++) {
            x[i]   = x[i]   + alpha * p[i];
            res[i] = res[i] - alpha * Ap[i];  }

         toll = scalarp(dim, res, res);

         memcpy(z_old, z, dim* sizeof(double));
 
         if (Mat->PrcPtr != NULL) 
            matvecd( Mat->PrcPtr, res, z);
         else 
            memcpy(z, res, dim* sizeof(double)); 

         beta = scalarp(dim, z, res) / scalarp(dim,z_old,res_old);

         for (i=0;i<dim;i++)
            p[i]  = z[i] + beta * p[i];

#if _DBG == 1
         fprintf(DbgPnt, "%d  %f\n", iter+1, sqrt(toll));
#endif
      } while (++iter != *iterr && sqrt(toll) > *tollr);
   }

   *iterr = iter;
   *tollr = sqrt(toll);

#if _DBG == 1
   fclose(DbgPnt);

   DbgPnt = fopen(fname("Solution"), "w");
   for (i=0;i<dim;i++)
      fprintf(DbgPnt,"%f\n", x[i]);
   fclose(DbgPnt);
#endif

   free(res); 
   free(p); 
   free(Ap); 
   free(z); 
   free(res_old); 
   free(z_old); 
}

double *y, *resh, *v, *s, *t, *aux1, *aux2;

/* Biconjucate gradient stabilized */
/* INPUT: Mat, RHS, x0, Toll_max, Niter
  OUTPUT: x, Toll_fin, Niter_fin         */
void BiCGStab(struct CSRMat *Mat, struct RHS *B, double x[Mat->dim*Mat->nvar], double *tollr, int  *iterr)
{
   int     i, iter = 0, dim;
   double toll, rho[2], omega[2], beta, alpha;
   struct CSRMat *K1;

#if _DBG == 1
   checkmat( Mat);
   prtmat( Mat, B);

   FILE *DbgPnt;

   DbgPnt = fopen(fname("SolCon"), "w");
   fprintf(DbgPnt, "Iterations Residual\n");
#endif

   /*This is the real dimension of the problem*/
   dim = Mat->dim * Mat->nvar;

   res = (double *) malloc(dim * sizeof(double)); 
   resh= (double *) malloc(dim * sizeof(double)); 
   p   = (double *) malloc(dim * sizeof(double)); 
   Ap  = (double *) malloc(dim * sizeof(double)); 
   z   = (double *) malloc(dim * sizeof(double)); 
   y   = (double *) malloc(dim * sizeof(double)); 
   v   = (double *) malloc(dim * sizeof(double)); 
   s   = (double *) malloc(dim * sizeof(double)); 
   t   = (double *) malloc(dim * sizeof(double)); 
   aux1= (double *) malloc(dim * sizeof(double)); 
   aux2= (double *) malloc(dim * sizeof(double)); 
   res_old = (double *) malloc(dim * sizeof(double)); 
   z_old   = (double *) malloc(dim * sizeof(double)); 

   if (B->dim!=Mat->dim || B->nvar!=Mat->nvar) {
      printf("Error Bdim != Matdim!!\n");
      printf(" Matdim %d Matnvar %d\n", Mat->dim, Mat->nvar);
      printf(" B  dim %d B  nvar %d\n", B->dim, B->nvar);
      return; }


   if (Mat->PrcPtr!= NULL) 
      createK1( Mat->PrcPtr);

   /*initial residual*/
   updtres( Mat, B->b, x, res);

   if (scalarp( dim, res, res) != 0.0) {

      memcpy(resh, res, dim* sizeof(double)); 

      rho[1]   = 1.0;
      alpha    = 1.0;
      omega[1] = 1.0;

      for (i=0;i<dim;i++)  {
         p[i]   = 0.0;
         v[i]   = 0.0; }

      do {
         rho[0] = scalarp(dim, resh, res);

         beta = (rho[0]/rho[1]) * (alpha/omega[1]);

         for (i=0;i<dim;i++) 
            p[i] = res[i] + beta*(p[i] - omega[1] * v[i]);

         if (Mat->PrcPtr != NULL) 
            matvecd( Mat->PrcPtr, p, y);
         else
            memcpy(y, p, dim* sizeof(double)); 

         matvecd( Mat, y, v);

         alpha = rho[0] / scalarp(dim, resh, v);
      
         for (i=0;i<dim;i++)
            s[i] = res[i] - alpha * v[i];

         if (Mat->PrcPtr != NULL) 
            matvecd( Mat->PrcPtr, s, z);
         else
            memcpy(z, s, dim* sizeof(double)); 

         matvecd( Mat, z, t);

         /*In case of precondition ...*/
         if (Mat->PrcPtr != NULL) {
            matvecd( Mat->PrcPtr->PrcPtr, t, aux1);
            matvecd( Mat->PrcPtr->PrcPtr, s, aux2);
 
            omega[0] = scalarp(dim, aux1, aux2) / scalarp(dim, aux1, aux1); 
            }
         else 
            omega[0] = scalarp(dim, t, s) / scalarp(dim, t, t);

         for (i=0;i<dim;i++) {
            x[i] = x[i] + alpha * y[i] + omega[0] * z[i];
            res[i] = s[i] - omega[0] * t[i];   }
      
         toll = scalarp(dim, res, res);

         rho[1]   = rho[0];
         omega[1] = omega[0];

#if _DBG == 1
         fprintf(DbgPnt, "%d  %f\n", iter+1, sqrt(toll));
#endif
      } while (++iter != *iterr && sqrt(toll) > *tollr);
   }

   *iterr = iter;
   *tollr = sqrt(toll);

#if _DBG == 1
   fclose(DbgPnt);

   DbgPnt = fopen(fname("Solution"), "w");
   for (i=0;i<dim;i++)
      fprintf(DbgPnt,"%f\n", x[i]);
   fclose(DbgPnt);
#endif

   if (Mat->PrcPtr!= NULL) 
      pcdel( Mat->PrcPtr->PrcPtr);

   free(res); 
   free(resh); 
   free(p); 
   free(Ap); 
   free(z); 
   free(y); 
   free(v); 
   free(s); 
   free(t); 
   free(aux1); 
   free(aux2); 
   free(res_old); 
   free(z_old); 
}

/* Successive Over Relaxation */
/* INPUT: Mat, RHS, x0, Toll_max, Niter
  OUTPUT: x, Toll_fin, Niter_fin         */
void SORmethd(struct CSRMat *Mat, struct RHS *B, double x[Mat->dim*Mat->nvar], double omega, double *tollr, int  *iterr)
{
   int     i, j,iter = 0, dim, istart, iend, imin, imax, incr;
   double toll, aii;
   FILE *DbgPnt;

#if _DBG == 1
   checkmat( Mat);
   prtmat( Mat, B);

   DbgPnt = fopen(fname("SolCon"), "w");
   fprintf(DbgPnt, "Iterations Residual\n");
#endif

   /*This is the real dimension of the problem*/
   dim = Mat->dim * Mat->nvar;

   res = (double *) malloc(dim * sizeof(double)); 

   if (B->dim!=Mat->dim || B->nvar!=Mat->nvar) {
      printf("Error Bdim != Matdim!!\n");
      printf(" Matdim %d Matnvar %d\n", Mat->dim, Mat->nvar);
      printf(" B  dim %d B  nvar %d\n", B->dim, B->nvar);
      return; }

#pragma omp parallel default(none) private( j, i, aii) \
shared( DbgPnt, iter, toll, dim, res, Mat, x, omega, B, iterr, tollr)
{
   do {
#pragma omp for
      for (i=0;i<dim;i++) {
         res[i] = 0.0;
         aii    = 0.0;

         for (j=Mat->AI[i];j<Mat->AI[i+1];j++) {
            res[i] += Mat->A[j] * x[Mat->AJ[j]];

            if (Mat->AJ[j] == i) 
               aii = Mat->A[j];
         }

         res[i] = (B->b[i] - res[i]) / aii;
         x[i]   = x[i] + omega * res[i];
      }
      // N.B.: technically speaking RES is not the residual.
      // This shold be obtained by the following:
      // updtres( Mat, B->b, x, res);
      // which increases the execution time a lot!!

#pragma omp master
{
      toll = sqrt(scalarp(dim, res, res));
#if _DBG == 1
      fprintf(DbgPnt, "%d  %f\n", iter+1, toll);
#endif
      ++iter;
}
#pragma omp barrier
   } while (iter != *iterr && toll > *tollr);
}

   *iterr = iter;
   *tollr = toll;

#if _DBG == 1
   fclose(DbgPnt);

   DbgPnt = fopen(fname("Solution"), "w");
   for (i=0;i<dim;i++)
      fprintf(DbgPnt,"%f\n", x[i]);
   fclose(DbgPnt);
#endif

   free(res); 
}
