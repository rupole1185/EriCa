#include "marisa.h"

//CSR MATRIX preallocation --------------------------------------------
void matini(int  dim, int  nvar, struct CSRMat *Mat)
{
   Mat->nvar   = nvar;
   Mat->dim    = dim;
   Mat->nptrn  = nvar*nvar;
   Mat->PrcPtr = NULL;

   Mat->AI = NULL;
   Mat->AJ = NULL;

   Mat->ptrn   = (int  *) calloc( nvar*nvar , sizeof(int ));
   memset( Mat->ptrn, 0, nvar*nvar* sizeof(int));
}

// Matrix Pattern: In case a similar patter is present
//   in all blocks it is possible to pass it in order to achieve
//   even better memory and efficiency performancies.
//   1 - Element not to be allocated
//   0 - Element to allocate
void matptrn( struct CSRMat *Mat, int ptrn[Mat->nvar][Mat->nvar])
{
   int ivar, iivar;

   for ( ivar=0; ivar < Mat->nvar; ivar++) {
      for ( iivar=0; iivar < Mat->nvar; iivar++) {
         if (ptrn[ivar][iivar] != 0) 
            Mat->ptrn[ivar*Mat->nvar+iivar] = 1; 
      }
   }

   Mat->nptrn -= sveciel( Mat->nvar * Mat->nvar, Mat->ptrn);
}

// Matrix allocation: case of diagonal matrix
void matdiag( struct CSRMat *Mat)
{
   int i, *dnnz, dim;

   dim = Mat->dim;
   
   dnnz   = (int *) calloc(dim, sizeof(int));

   for (i=0;i<dim;i++)
      dnnz[i] = 1;

   matalc( Mat, dnnz);

   free( dnnz);
}

// Matrix allocation: dnnz represents the non zero in each row
// NB: pattern are used to minimize the allocation space
void matalc( struct CSRMat *Mat, int dnnz[Mat->dim])
{
   int     i, ivar, iivar, dnnzt, nnzp, dim, nvar, nnzptrn[Mat->nvar];

   nvar = Mat->nvar;
   dim  = Mat->dim;

   /*Total number of non zeros*/
   dnnzt = sveciel(dim, dnnz);

   Mat->nnz = dnnzt;

   if (Mat->AI==NULL)
      Mat->AI = (int  *) calloc(((dim * nvar) + 1), sizeof(int ));
   if (Mat->AJ==NULL)
      Mat->AJ = (int  *) calloc(dnnzt * Mat->nptrn, sizeof(int ));

   Mat->A  = (double *)  calloc(dnnzt * Mat->nptrn, sizeof(double));

   Mat->AIivar = (int *) calloc(nvar+1 , sizeof(int));
   Mat->AJivar = (int *) calloc(Mat->nptrn , sizeof(int));

   //Ivaridx definition
   for (ivar=0;ivar<nvar;ivar++) {
      Mat->AIivar[ivar+1] = Mat->AIivar[ivar];

   for (iivar=0;iivar<nvar;iivar++)
      if (Mat->ptrn[ivar*nvar+iivar] == 0) {
         Mat->AJivar[Mat->AIivar[ivar+1]] = iivar;
         Mat->AIivar[ivar+1]++;   }
   }

   if (Mat->AIivar[nvar]!=Mat->nptrn) {
      printf("Error Matrix allocations\n");
      exit(1);
   }

   //Nnz contained in a single block
   for (ivar=0;ivar<nvar;ivar++)
      nnzptrn[ivar] = nvar - sveciel( nvar, &Mat->ptrn[ivar*nvar]);

   //AI is allocated according to the pattern given
   nnzp = 0; 
   for (i=0;i<dim;i++)  {
      for (ivar=0;ivar<nvar;ivar++) {
         Mat->AI[i*nvar+ivar] = nnzp;
         nnzp += (nnzptrn[ivar]*dnnz[i]);
      }
   }
   Mat->AI[dim*nvar] = nnzp;

#if _DBG == 0
   free(Mat->ptrn);
#endif

   matstart(0.0, 0, Mat);

#if _DBG == 1
   printf("\nMatrix info:");
   printf("DIM  = %d\n", Mat->dim);
   printf("NVAR = %d\n", Mat->nvar);
   printf("NNZ  = %d\n", Mat->nnz);

   printf("PatternScheme\n");
   for ( ivar=0; ivar < Mat->nvar; ivar++) {
      for ( iivar=0; iivar < Mat->nvar; iivar++) {
            printf(" %d", Mat->ptrn[ivar*Mat->nvar+iivar]);
      }
      printf("\n");
   }
#endif
}

void matdel( struct CSRMat *Mat)
{
   free(Mat->AI);
   free(Mat->AJ);
   free(Mat->A);
}

//CSR RHS preallocation ----------------------------------------------
void rhsini(int  dim, int  nvar, struct RHS *B)
{
   int     i, ivar;

   B->nvar     = nvar;
   B->dim      = dim;

   B->b    = (double *) calloc(dim*nvar , sizeof(double));

   rhsstart( 0.0, B); 
}

void rhsdel( struct RHS *B)
{
   free(B->b);
}

//CSR MATRIX initialization to given value ---------------------------
void matstart(double init, int opt, struct CSRMat *Mat)
{
   int  i, nnz, dim, nvar;

   dim  = Mat->dim;
   nvar = Mat->nvar;

   nnz = Mat->AI[dim*nvar]; 

   //Resetting A Matrix
   if (opt == 0) {
      for (i=0;i<nnz;i++)
         Mat->AJ[i] = -1;
   }

   for (i=0;i<nnz;i++)
      Mat->A[i] = init;
}

void rhsstart(double init, struct RHS *B)
{
   int  dim, nvar, i;

   dim  = B->dim;
   nvar = B->nvar;

   //Resetting RHS
   for (i=0;i<dim*nvar;i++)
      B->b[i] = init;
}

//CSR MATRIX element inserting ----------------------------------------
void getnnz(struct CSRMat *Mat, int  row, int  col, double coeff, double add[Mat->nvar][Mat->nvar])
{
   int iptrn, i, j, istart, iend, ivar, iivar, nvar, flag=0, dim;
   int ivars, ivare, jvar;

   if (coeff==0.0) return;

   dim = Mat->dim;

   if (row<0 || row>=dim) return;
   if (col<0 || col>=dim) return;

   nvar = Mat->nvar;

   for (ivar=0;ivar<nvar;ivar++) {
      istart = Mat->AI[row*nvar+ivar];
      iend   = Mat->AI[row*nvar+ivar+1];

      ivars = Mat->AIivar[ivar];
      ivare = Mat->AIivar[ivar+1];

      for (jvar=ivars;jvar<ivare;jvar++) {
         iivar = Mat->AJivar[jvar];
   
         for (i=istart;i<iend;i++)  {
            j = Mat->AJ[i];

            if (j<0 || j == col*nvar+iivar) {
               Mat->AJ[i] = col*nvar+iivar;
               Mat->A[i] +=(coeff * add[ivar][iivar]);
#if _DBG == 1 
               flag++;
#endif
               istart = i+1;
               break;
            }
         } 
      }
   }

#if _DBG == 1
   /*Passare a variabile ierr!!!*/
   if (flag!=Mat->nptrn) {
      printf("ERROR: wrong number of variables! %d %d %d\n", flag, Mat->nptrn, Mat->nvar);
      printf("Row = %d, Col = %d\n", row, col);
   }

   // Checking conformity with pattern
   flag = 0;

   for (ivar=0;ivar<nvar;ivar++) {
   for (jvar=0;jvar<nvar;jvar++) {
      if ((Mat->ptrn[ivar*nvar+jvar] == 1) && \
           add[ivar][jvar] != 0.0)
         flag = 1;
   }
   }

   if (flag == 1) {
      printf("\n Matrix differs from pattern!\n");

      printf("\n ADD MAT\n");
      for (ivar=0;ivar<nvar;ivar++) {
         for (jvar=0;jvar<nvar;jvar++) {
            printf(" %f", add[ivar][jvar]);
         }
         printf("\n");
      }

      printf("\n PATTERN\n");
      for (ivar=0;ivar<nvar;ivar++) {
         for (jvar=0;jvar<nvar;jvar++) {
            printf(" %d", Mat->ptrn[ivar*nvar+jvar]);
         }
         printf("\n");
      }
      exit(1);
   }

#endif
}


/*RHS vector element inserting*/
void getrhs(struct RHS *B, int  row, double coeff, double add[B->nvar])
{
   int  nvar,ivar;

   if (row<0 || row>=B->dim) return;

   nvar = B->nvar;

   for (ivar=0;ivar<nvar;ivar++) 
      B->b[row*nvar+ivar] += coeff * add[ivar];
}

/* Extract row from the matrix */
void extrow( int  idim, int  ivar, int  *i, int  *ncol, int  **col, double **elem, struct CSRMat *Mat)
{
   *i = idim * Mat->nvar + ivar;

   *ncol = Mat->AI[*i+1] - Mat->AI[*i];

   *col  = &(Mat->AJ[Mat->AI[*i]]);
   *elem = &(Mat->A[Mat->AI[*i]]);
}

/* Extract element from the rhs */
double extrhs( int idim, int ivar, struct RHS *Rhs)
{
   return Rhs->b[idim * Rhs->nvar+ivar];
}
