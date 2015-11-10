#include "simulation.h"

int file_exists(const char * filename)
{
   FILE *file;

   if (file = fopen(filename, "r"))
   {
      fclose(file);
      return 1;
   }

   return 0;
}

void user_order(void )
{
   // Modifying SIMULATION settings ----------------------------
   if (file_exists("MODIFY.sim")) {
      fprintf(logfile,"\n @@@ USER ORDER: simulation changes\n");
      fprintf(logfile,"\n@@@@@@@@@@@@@@@@@ NEW SETTINGS @@@@@@@@@@@@@@@@@");
      printf("\r");
      simini();
      summary();
      fprintf(logfile,"\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
      fprintf(logfile,"\n\n @@@ simulation continues ...\n\n");
      remove("MODIFY.sim");
   }

   // Stopping simulation --------------------------------------
   if (file_exists("STOP.sim")) {
      fprintf(logfile,"\n @@@ USER ORDER: stopping simulation\n\n");
      ntmst = 0;
      remove("STOP.sim");
   }

   fflush(stdout);
   fflush(logfile);
}
