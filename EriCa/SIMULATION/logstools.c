#include "simulation.h"
#include "../MESH/mesh.h"

void ScreenStars(void )
{
   int i;
   for (i=0;i<ScreenSize;i++)
      printf("*");
   printf("\n");
}

void LogStars(void )
{
   int i;
   for (i=0;i<LogSize;i++)
      fprintf(logfile,"*");
   fprintf(logfile,"\n");
}

const char* VarName( int ivar)
{
   static char varname[12];

   if (ivar == Pvar)
      strncpy(varname, "Pressure   ", 11);
   else if (ivar==Uvar)
      strncpy(varname, "VelocityX  ", 11);
   else if (ivar==Uvar+1 && ndir > 1)
      strncpy(varname, "VelocityY  ", 11);
   else if (ivar==Uvar+2 && ndir > 2)
      strncpy(varname, "VelocityZ  ", 11);
   else if (ivar== Tvar)
      strncpy(varname, "Temperature", 11);
   else if (ivar== Kvar)
      strncpy(varname, "TKE        ", 11);
   else if (ivar== OMvar)
      strncpy(varname, "Omega      ", 11);
   else if (ivar== MTvar)
      strncpy(varname, "Nu-tilde   ", 11);
   else if (ivar>= SCAvar && ivar < SCAvar+NSCAvar)
      sprintf(varname, "Scalar%2.2d   ", ivar - SCAvar);
   else
      strncpy(varname, "ERROR!!    ", 11);

   return varname;
}

const char* DepName( int idep)
{
   static char depname[12];

   if (idep == Rdep)
      strncpy(depname, "Density    ", 11);
   else if (idep == MUdep)
      strncpy(depname, "Viscosity  ", 11);
   else if (idep>= SCHMIDTdep && idep < SCHMIDTdep+NSCAvar)
      sprintf(depname, "Schmidt%2.2d  ", idep - SCHMIDTdep);
   else
      strncpy(depname, "ERROR!!    ", 11);

   return depname;
}

void LogSection( char *SecName)
{
   fprintf(logfile,"\n\n");
   fprintf(logfile,"%s", SecName);
   fprintf(logfile," ");
   int i = 0, init = StringLength( SecName);
   for(i=init+1;i<LogSize;i++)
      fprintf(logfile,"-");
   fprintf(logfile,"\n");
}

int StringLength( char *string)
{
   int i=0;
   do {
      i++;
   } while (string[i] != '\0');

   return i;
}

// Printing Max and Min into a table
void Range2table( char *String, int idep, int vd) 
{
   if (idep < 0)
      return;

   double * (*setphidep)( int);

   switch (vd) {
   case 0:
      setphidep = &setvar;
      break;
   case 1:
      setphidep = &setdep;
      break;
   default:
      break;
   }

   double min, max, *phidep;
   int icell;

   min = 1.E10;
   max =-1.E10;

   for (icell=0;icell<ncell;icell++) {
      phidep = setphidep( icell);

      if (min > phidep[idep])
         min = phidep[idep];

      if (max < phidep[idep])
         max = phidep[idep];
   }

   fprintf(logfile, " | %-11s | %12.3f | %12.3f |\n", String, min, max);
}
