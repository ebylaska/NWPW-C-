#ifdef ATL_NCPU

#include "atlas_misc.h"
#include "assert.h"

void PrintUsage(char *nam)
{
   fprintf(stderr, "\nUSAGE: %s [-o <outfile>]\n", nam);
   exit(-1);
}

void GetFlags(int nargs, char **args, FILE **fpout)
{
   int i;

   *fpout=stdout;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-') PrintUsage(args[0]);
      switch(args[i][1])
      {
      case 'o':
         *fpout = fopen(args[++i], "w");
         assert(*fpout);
         break;
      default:
         PrintUsage(args[0]);
      }
   }
}

void getLaunchOrder(int P, int *lo)
{
   int i, j, k, stop, dest;

   for (i=0; (1<<i) < P; i++)
   lo[0] = 0;
   k = 1;
   for (i--; i >= 0; i--)
   {
      stop = k;
      for (j=0; j < stop; j++)
      {
         dest = lo[j] + (1<<i);
         if (dest < P)
            lo[k++] = dest;
         if (k == P)
            return;
      }
   }
}
main(int nargs, char **args)
{
   FILE *fpout;
   int i, j, k;
   #if ATL_NCPU > 0
      int lo[ATL_NCPU];
   #endif
   GetFlags(nargs, args, &fpout);

   fprintf(fpout, "#ifndef ATLAS_NTHREADS_H\n   #define ATLAS_NTHREADS_H\n\n");
/*
 * I presently build Antoine's pthread implementation even on windows for
 * comparison purposes.  Need to get rid of 00 when this is no longer the
 * case.
 */
   fprintf(fpout, "/* Get rid of 00 if you don't want to build pthreads */\n");
   fprintf(fpout,
      "   #ifndef ATL_WINTHREADS00\n      #include \"pthread.h\"\n   #endif\n");
   #if ATL_NCPU != 0
      fprintf(fpout, "   #define ATL_NTHREADS %d\n", ATL_NCPU);
      for (i=0; (1<<i) < ATL_NCPU; i++);
      fprintf(fpout, "   #define ATL_NTHRPOW2 %d\n", i);
      getLaunchOrder(ATL_NCPU, lo);
      fprintf(fpout, "   #ifdef ATL_LAUNCHORDER\n");
      fprintf(fpout, "       static int ATL_launchorder[%d] = {0", ATL_NCPU);
      for (i=1; i < ATL_NCPU; i++)
         fprintf(fpout, ",%d", lo[i]);
      fprintf(fpout, "};\n   #endif\n");
   #else
      fprintf(fpout, "   #define ATL_NTHREADS 4\n");
      fprintf(fpout, "   #define ATL_NTHRPOW2 2\n")
      fprintf(fpout, "   #ifdef ATL_LAUNCHORDER\n");
      fprintf(fpout, "       static int ATL_launchorder[4] = {0,2,1,3}\n");
      fprintf(fpout, "   #endif\n");
   #endif
   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
   exit(0);
}
#endif
