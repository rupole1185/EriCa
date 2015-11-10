#include "simulation.h"

// ERICA DATA STRUCTURE
// Ncells -> Number of Cells
// Nvar   -> Number of Variables
//
// All variables are stored in single arrays with the following structure:
// VARIABLEARRAY[      1] - VARIABLEARRAY[   NVAR] -> 1st cell variables
// VARIABLEARRAY[ NVAR+1] - VARIABLEARRAY[ 2*NVAR] -> 2nd cell variables
// and so on and so for ...

void simstencil()
{
   //Setting total number of independent variables: -------------------
   nvar = 0;

   Pvar = nvar++;
   Uvar = nvar;
   nvar+= ndir;

   if (eos > COMPR_b)
      Tvar = nvar++;

   if (turb==SA)
      MTvar = nvar++;

   if (turb>TURBULENT_b && turb/100==2)
      Kvar = nvar++;

   if (turb/100==2 && turb<250)
      OMvar = nvar++;

   if (NSCAvar!=0) {
      SCAvar = nvar;
      nvar  += NSCAvar;
   }

   // Defining dependent variables --------------------------------------
   ndep    = 0;
   ndepini = 0;

   Rdep = ndep++;
   ndepini++;

   if (turb > LAMINAR_b) {
      MUdep = ndep++;
      ndepini++;

      if (NSCAvar!=0) {
         SCHMIDTdep = ndep;
         ndep      += NSCAvar;
         ndepini   += NSCAvar; }
   }

   if (turb > DISTANCE_b && turb != KW06) 
      Ddep = ndep++;

   if (turb > TURBULENT_b) 
      MTdep = ndep++;

   //Defining dependent variables [GRADIENTS]
   PGrd = ndep;
   ndep+= ndir;

   UGrd = ndep;
   ndep+= ndir*ndir;

   if (eos > COMPR_b) {
      TGrd = ndep;
      ndep+= ndir;}

   if (MTvar!=-1) {
      MTGrd = ndep;
      ndep += ndir;}

   if (Kvar!=-1) {
      KGrd  = ndep;
      ndep += ndir;}

   if (OMvar!=-1) {
      OMGrd = ndep;
      ndep += ndir;}

   if (SCAvar!=-1) {
      SCAGrd = ndep;
      ndep += ndir*NSCAvar;
   }

   // Checks gradient allocations
   if (ndep != PGrd+(nvar*ndir))
      FatalError("Error Gradient Allocation");

   CFLdep = ndep++;
   TIMEdep= ndep++;

#if _DBG == 10
   switch (eos) {
   case ISIMPLE:
   case ISIMPLEC:
   case IPISO:
      Ppdep  = ndep++;
      Updep  = ndep;
      ndep  += ndir;
      DivSdep= ndep++;
      DDdep  = ndep;
      ndep  += ndir;
      break;
   case IFROZENP:
      break;
   case IDBIG:
      break;
   }

   Vrtdep = ndep++;

   Divdep = ndep++;

   if (Kvar!=-1)
      Pkdep  = ndep++;

   if (Ddep!=-1)
      PHIDdep = ndep++;
#endif
}

const char* boptions( int ichoice)
{
   static char optname[12];

   switch (ichoice) {
   case 0:
      strncpy(optname, "bndtype    ", 11);
      break;
   case 1:
      strncpy(optname, "bndvar     ", 11);
      break;
   case 2:
      strncpy(optname, "monitor    ", 11);
      break;
   }

   return optname;
}
