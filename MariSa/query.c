#include "marisa.h"

int nprt = 0;

// Subroutine that returns the name of the FOLDER/FILEXXX.dat
const char * fname(char Name[10])
{
   char Suffix[5+2];

   strcpy(FName,"DEBUG/MrS");
   strcat(FName,Name);
   sprintf(Suffix, ".%04d", nprt);
   strcat(FName,Suffix);
   sprintf(Suffix, ".dat");
   strcat(FName,Suffix);

   return FName;
}

//MATRIX CHECKS -------------------------------------------------------
//Print  MATLAB file
void prtmat(struct CSRMat *Mat, struct RHS *B)
{
   int  i, j;
   FILE     *writeOUT;

   nprt++;

   writeOUT = fopen(fname("Matr"),"w"); 
   fprintf(writeOUT, "Row, Col, AIJ\n");

   fprintf(writeOUT, "Nvar = %d\n", Mat->nvar);

   for (i=0; i<Mat->dim*Mat->nvar; i++)  {
      for (j=Mat->AI[i]; j<Mat->AI[i+1]; j++) 
         fprintf(writeOUT, "%d %d %15.10E\n", i+1,Mat->AJ[j]+1,Mat->A[j]);
   }
   
   fclose(writeOUT);

   writeOUT = fopen(fname("Rhs"),"w"); 
   fprintf(writeOUT, " BI\n");

   fprintf(writeOUT, "Nvar = %d\n", B->nvar);

   for (i=0; i<B->dim*B->nvar; i++) 
      fprintf(writeOUT, "%15.10E\n", B->b[i]);
   
   fclose(writeOUT);
}

/*Subroutine to check the MATRIX*/
void checkmat(struct CSRMat *Mat)
{
   int     idim, ivar, nvar, dim;
   int     irow, ncol, *col, iel;
   double *elem, aii;

   /*This is the real dimension of the problem*/
   dim = Mat->dim;
   nvar= Mat->nvar;

   for (idim=0;idim<dim;idim++) {
      for (ivar=0;ivar<nvar;ivar++) {
         extrow( idim, ivar, &irow, &ncol, &col, &elem, Mat);

         for (iel=0;iel<ncol;iel++)
            if (irow==col[iel])
               aii = fabs(elem[iel]);

         for (iel=0;iel<ncol;iel++) {
            if (fabs(elem[iel])>aii) {
               printf("\nRow %d, var %d, is not DD\n", idim, ivar); 
               printf("Col %d, var %d, may have a problem\n", col[iel]/nvar, col[iel]%nvar);
               printf("VarAii = %f, VarAij = %f\n", aii, elem[iel] );
               exit(1);
            }
         }
      }
   }
}

// Memory enquiries -----------------
double MatMem( struct CSRMat *Mat)
{
   int nnz, dim;
   double size;

   dim = Mat->dim * Mat->nvar;
   nnz = Mat->AI[dim];

   size = (dim+1) * sizeof(int) + nnz * (sizeof(double)+sizeof(int)) + \
                                  (Mat->nvar+1+Mat->nptrn)*sizeof(int);

   if (Mat->PrcPtr != NULL)
      size += MatMem( Mat->PrcPtr);

   return size;
}

double RhsMem( struct RHS *B)
{
   int dim;
   double size;

   dim = B->dim * B->nvar;

   size = dim * sizeof(double);

   return size;
}

int Matnnz( struct CSRMat *Mat)
{
   return Mat->nnz;
}

int Matnconn( struct CSRMat *Mat)
{
   return Mat->AI[Mat->dim*Mat->nvar];
}
