#include "marisa.h"

//Function to calculate the scalar product among two different arrays*/------------------------
double scalarp(int  dim, double a[dim], double b[dim])
{
   int  i;
   double prod;

   prod = 0.0;

   for (i=0;i<dim;i++)
      prod += a[i]*b[i];

   return prod;
}

//Standard Matrix vector product
void matvecdstd(int  dim, double mat[dim][dim], double vec_in[dim], \
                                                double vec_out[dim])
{ 
   int     i;

   for (i=0;i<dim;i++) 
      vec_out[i] = scalarp(dim,&mat[i][0],vec_in); 
}

//Returns the sum of the elements of a integer Vector*/---------
int  sveciel(int  dim, int  Vec[dim])
{ 
   int  sum, i;

   sum = 0;

   for (i=0;i<dim;i++)
      sum += Vec[i];

   return sum;
}

//Returns the sum of the elements of a double Vector*/---------
double sveciel_dbl(int  dim, double Vec[dim])
{ 
   int   i;
   double sum;

   sum = 0.0;

   for (i=0;i<dim;i++)
      sum += Vec[i];

   return sum;
}

//Summ of two vectors-----------------------------------
void summvec(int  dim, double Vec1[dim], double Vec2[dim], double Sum[dim])
{
   int  i;

   for (i=0;i<dim;i++)
      Sum[i] = Vec1[i] + Vec2[i];
}

//Inverse of a matrix--------------------------------------------
void invmat(int m, int  n, double mat[m][n])
{
    if (m==n) {
       double ratio,a, matrix[n][2*n];
       int  i, j, k;

       for(i = 0; i < n; i++){
           for(j = 0; j <   n; j++)
               matrix[i][j] = mat[i][j];
   
           for(j = n; j < 2*n; j++){
               if(i==(j-n))
                   matrix[i][j] = 1.0;
               else
                   matrix[i][j] = 0.0;
           }
       }
       for(i = 0; i < n; i++){
           for(j = 0; j < n; j++){
               if(i!=j){
                   ratio = matrix[j][i]/matrix[i][i];
                   for(k = 0; k < 2*n; k++){
                       matrix[j][k] -= ratio * matrix[i][k];
                   }
               }
           }
       }
       for(i = 0; i < n; i++){
           a = matrix[i][i];
           for(j = 0; j < 2*n; j++){
               matrix[i][j] /= a;
           }
       }

       for(i = 0; i < n; i++){
           for(j = 0; j < n; j++)
               mat[i][j] = matrix[i][n+j];
       }
    }
}

// LU factorization ---------------------------------------------
void LU( int n, double mat[n][n], double L[n][n] , double U[n][n])
{
// No row exchanges
// Square matrix
   int k , i , j ;

   for ( k = 0; k<n; k++) {
      L[k][k] = 1;
      for ( i = k+1; i<n; i++) {
         L[i][k]= mat[i][k] / mat[k][k];
         // mat[i][k] = L[i][k];
         for ( j = k+1; j<n; j++) 
            mat[i][j] = mat[i][j] - L[i][k] * mat[k][j];
      }
      for ( j = k; j<n ; j++) 
         U[k][j] = mat[k][j];
   }
}
