#include "eqntns.h"

void timestp(struct CSRMat *Mat, struct RHS *Rhs )
{
   int icell, *cellintPnt;
   double Mm[nvar][nvar], Flux[nvar], cnst, *phi, *dep, *celldblPnt;

   for (icell=0;icell<ncell;icell++) {
      cellintPnt = setintcell( icell);
      celldblPnt = setdblcell( icell);

      //Initialization of Mm
      memset(Mm, 0.0, nvar * nvar * sizeof(double));
    
      /*Variables to use*/
      phi = setvar( icell);
      dep = setdep( icell);

      cnst = celldblPnt[CellVol] / dep[TIMEdep];

      //Time step matrix to assemble according to the alghoritm
      timematrix( phi, dep, Mm);
 
      //Adding to RHS
      matvecdstd(nvar, Mm, phi, Flux);
      getrhs( Rhs, icell, cnst, Flux);

      //Adding to Mat
      getnnz( Mat, icell, icell, cnst, Mm);
   }
}
