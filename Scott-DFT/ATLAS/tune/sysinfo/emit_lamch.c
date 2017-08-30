#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int sComputeRound(void)
/*
 * Blind translation of netlib LAPACK LAMCH's rounding computation
 * RETURNS: 1 if numbers are correctly rounded, 0 if they are truncated
 */
{
   volatile float a, b, c, f;
   int rnd=0;
   b = a = 1.0;
   do
   {
      a *= 2.0;
      c = a + b;
      c = c - a;
   }
   while (c == 1.0);
   b = 0.5*FLT_RADIX;
   c = .01*(-FLT_RADIX);
   f = b + c;
   c = f + a;
   rnd = (c == a);

   c = 0.01*(FLT_RADIX);
   f = b + c;
   c = f + a;
   if (rnd && c == a)
      rnd = 0;
   return(rnd);
}

float sComputeSafmin(void)
/*
 * BFI translation of netlib LAPACK LAMCH's safmin calc, adapted to use float.h
 * RETURNS: LAMCH's safmin
 */
{
   volatile float small;
   small = 1.0/(pow(FLT_RADIX, FLT_MAX_EXP-2)*(4.0-2.0*FLT_EPSILON));
   if (small >= FLT_MIN)
      return(small*(1.0+0.5*FLT_EPSILON));
   return(FLT_MIN);
}

void emit_slamch(char *path)
{
   FILE *fpout;
   char *name;
   volatile float f;
   int len = 16;

   if (path)
   {
      len += strlen(path);
      name = malloc(len);
      assert(name);
      sprintf(name, "%s/atlas_slamch.h", path);
      fpout = fopen(name, "w");
      assert(fpout);
      free(name);
   }
   else
      fpout = stdout;
   fprintf(fpout, "#ifndef ATLAS_SLAMCH_H\n");
   fprintf(fpout, "   #define ATLAS_SLAMCH_H\n\n");

   fprintf(fpout, "#define ATL_slaMANTDIG     %d\n", FLT_MANT_DIG);
   fprintf(fpout, "#define ATL_slaMINEXP      %d\n", FLT_MIN_EXP);
   fprintf(fpout, "#define ATL_slaMAXEXP      %d\n", FLT_MAX_EXP);
   fprintf(fpout, "#define ATL_slaBASE        %d\n", FLT_RADIX);
   f = 0.5;
   f *= FLT_EPSILON;
   fprintf(fpout, "#define ATL_slaEPSILON     %60.53e\n", f);
   f = 0.5 * FLT_RADIX;
   f *= FLT_EPSILON;
   fprintf(fpout, "#define ATL_slaPRECISION   %60.53e\n", f);
   fprintf(fpout, "#define ATL_slaUNDERTHRESH %60.53e\n", FLT_MIN);
   f = pow(FLT_RADIX, FLT_MAX_EXP-2)*(4.0-2.0*FLT_EPSILON);
   fprintf(fpout, "#define ATL_slaOVERTHRESH  %60.53e\n", f);
   fprintf(fpout, "#define ATL_slaSAFMIN      %60.53e\n",
           sComputeSafmin());
   fprintf(fpout, "#define ATL_slaROUND       %d\n", sComputeRound());

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

int dComputeRound(void)
/*
 * Blind translation of netlib LAPACK LAMCH's rounding computation
 * RETURNS: 1 if numbers are correctly rounded, 0 if they are truncated
 */
{
   volatile double a, b, c, f;
   int rnd=0;
   b = a = 1.0;
   do
   {
      a *= 2.0;
      c = a + b;
      c = c - a;
   }
   while (c == 1.0);
   b = 0.5*FLT_RADIX;
   c = .01*(-FLT_RADIX);
   f = b + c;
   c = f + a;
   rnd = (c == a);

   c = 0.01*(FLT_RADIX);
   f = b + c;
   c = f + a;
   if (rnd && c == a)
      rnd = 0;
   return(rnd);
}

double dComputeSafmin(void)
/*
 * BFI translation of netlib LAPACK LAMCH's safmin calc, adapted to use float.h
 * RETURNS: LAMCH's safmin
 */
{
   volatile double small;
   small = 1.0/(pow(FLT_RADIX, DBL_MAX_EXP-2)*(4.0-2.0*DBL_EPSILON));
   if (small >= DBL_MIN)
      return(small*(1.0+0.5*DBL_EPSILON));
   return(DBL_MIN);
}

void emit_dlamch(char *path)
{
   FILE *fpout;
   char *name;
   volatile double f;
   int len = 16;

   if (path)
   {
      len += strlen(path);
      name = malloc(len);
      assert(name);
      sprintf(name, "%s/atlas_dlamch.h", path);
      fpout = fopen(name, "w");
      assert(fpout);
      free(name);
   }
   else
      fpout = stdout;
   fprintf(fpout, "#ifndef ATLAS_DLAMCH_H\n");
   fprintf(fpout, "   #define ATLAS_DLAMCH_H\n\n");

   fprintf(fpout, "#define ATL_dlaMANTDIG     %d\n", DBL_MANT_DIG);
   fprintf(fpout, "#define ATL_dlaMINEXP      %d\n", DBL_MIN_EXP);
   fprintf(fpout, "#define ATL_dlaMAXEXP      %d\n", DBL_MAX_EXP);
   fprintf(fpout, "#define ATL_dlaBASE        %d\n", FLT_RADIX);
   f = 0.5;
   f *= DBL_EPSILON;
   fprintf(fpout, "#define ATL_dlaEPSILON     %60.53e\n", f);
   f = 0.5 * FLT_RADIX;
   f *= DBL_EPSILON;
   fprintf(fpout, "#define ATL_dlaPRECISION   %60.53e\n", f);
   fprintf(fpout, "#define ATL_dlaUNDERTHRESH %60.53e\n", DBL_MIN);
   f = pow(FLT_RADIX, DBL_MAX_EXP-2)*(4.0-2.0*DBL_EPSILON);
   fprintf(fpout, "#define ATL_dlaOVERTHRESH  %60.53e\n", f);
   fprintf(fpout, "#define ATL_dlaSAFMIN      %60.53e\n",
           dComputeSafmin());
   fprintf(fpout, "#define ATL_dlaROUND       %d\n", dComputeRound());

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

int main (int nargs, char **args)
{
   char *path = "res/";
   if (nargs > 1)
      path = args[1];
   emit_dlamch(path);
   emit_slamch(path);
   return(0);
}
