/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 *                    (C) Copyright 2009 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),sysinfo.h))

static char *resfile=NULL;
static FILE *fpres=NULL;
   #ifdef TREAL
      return(((1.0e-6 * M)*(2.0*N+1.0))/time);
   #else
      return((((6.0*M)*(N+1.0) + (2.0*M)*N)*1.0e-6)/time);
   #endif
}

void Times2Flops(ATL_INT M, ATL_INT N, ATL_INT ntim, double *mf)
/*
 * Converts time to MFLOP
 */
{
   int i;

   for (i=0; i < ntim; i++)
      mf[i] = Time2Flop(M, N, mf[i]);
}

static double mysum(ATL_CINT N, double *d)
{
   int i;
   double sum;

   sum = d[0];
   for (i=1; i < N; i++)
      sum += d[i];
   return(sum);
}

#define ATL_UGEMV Mjoin(Mjoin(Mjoin(PATL,gemvT_a1_x1),_b0),_y1)
void ATL_UGEMV(ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A,
               ATL_CINT lda, const TYPE *X, ATL_CINT incX,
               const SCALAR beta, TYPE *Y, ATL_CINT incY);
static void mvsimN(
   size_t celts,    /* # of elts in cache size we are blocking for */
   size_t pgelts,   /* # of elts on a virtual mem page (best guess) */
   ATL_CINT xu,     /* unrolling on X's loop by this kernel */
   ATL_CINT yu,     /* unrolling on Y's loop by this kernel */
   enum ATLAS_TRANS TA,
   ATL_CINT M,
   ATL_CINT N,
   const TYPE *A,
   ATL_CINT lda,
   const TYPE *X,
   const TYPE *beta,
   TYPE *Y)
/*
 * This routine assumes the Notranspose case, where we read X in the outer
 * loop, and do a axpy with each column of A into Y in the inner loop,
 * and we therefore cut M in order encourage cache reuse of Y
 */
{
   ATL_INT Mp, m, i;
   #ifdef TREAL
      #define BETA *beta
      const TYPE alpha = ATL_rone;
   #else
      const TYPE alpha[2] = {ATL_rone, ATL_rzero};
      #define BETA beta
   #endif
/*
 * Compute where to cut M in order to get to reuse the vector in
 * the celts-length cache
 */
   Mp = (celts - 2*xu) / (2*xu+1);
   if (Mp)
   {
/*
 *    Keep the M partition a multiple of the page size if that doesn't drastically
 *    change it.  This will tend to minimize the cost of moving row-wise by often
 *    avoiding an extra page load that isn't used much.
 */
      if (celts > pgelts+pgelts)
         celts = (celts/pgelts)*pgelts;
   }
   else Mp = M;
   for (i=0; i < M; i += Mp)
   {
      m = M - i;
      m = (m > Mp) ? Mp : m;
      ATL_UGEMV(m, N, alpha, A+(i SHIFT), lda, X, 1, BETA, Y+(i SHIFT), 1);
   }
}

static void mvsimT(
   size_t celts,    /* # of elts in cache size we are blocking for */
   size_t pgelts,   /* # of elts on a virtual mem page (best guess) */
   ATL_CINT xu,     /* unrolling on X's loop by this kernel */
   ATL_CINT yu,     /* unrolling on Y's loop by this kernel */
   enum ATLAS_TRANS TA,
   ATL_CINT M,
   ATL_CINT N,
   const TYPE *A,
   ATL_CINT lda,
   const TYPE *X,
   const TYPE *beta,
   TYPE *Y)
/*
 * This routine assumes the transpose case, where we write Y in the outer
 * loop, and apply X to each column of A in the inner loop, and we therefore
 * cut M in order encourage cache reuse of X
 */
{
   ATL_INT Mp, m, i;
   #ifdef TREAL
      #define BETA *beta
      const TYPE alpha = ATL_rone;
   #else
      const TYPE alpha[2] = {ATL_rone, ATL_rzero};
      #define BETA beta
   #endif
/*
 * Compute where to cut M in order to get to reuse the vector in
 * the celts-length cache
 */
   Mp = (celts - 2*yu) / (2*yu+1);
   if (Mp)
   {
/*
 *    Keep the M partition a multiple of the page size if that doesn't drastically
 *    change it.  This will tend to minimize the cost of moving row-wise by often
 *    avoiding an extra page load that isn't used much.
 */
      if (celts > pgelts+pgelts)
         celts = (celts/pgelts)*pgelts;
   }
   else Mp = M;
   for (i=0; i < M; i += Mp)
   {
      m = M - i;
      m = (m > Mp) ? Mp : m;
//      ATL_UGEMV(m, N, alpha, A+(i SHIFT), lda, X+(i SHIFT), 1, BETA, Y, 1);
      ATL_UGEMV(N, m, alpha, A+(i SHIFT), lda, X+(i SHIFT), 1, BETA, Y, 1);
   }
}

double mvtime_OC(
   int nreps,           /* number of reps to do for one timing sample */
   ATL_INT flushelts,   /* size of area to flush to avoid cache reuse */
   ATL_INT celts,       /* # of elts for cache blocking by MV driver */
   ATL_INT pgelts,      /* guess for virtual mem page size -- used by MV driver */
   enum ATLAS_TRANS TA,
   ATL_INT M,           /* # of rows of array A */
   ATL_INT N,           /* # of cols of array A */
   ATL_INT lda,         /* leading dim */
   TYPE *beta,
   int xu,              /* unrolling on X */
   int yu,              /* unrolling on Y */
   int FAa,             /* if (FA. = 0) enforce no alignment */
   int MAa,             /* else force op to be aligned to at least FA bytes */
   int FAx,             /* if MA. != 0, disallow op to be aligned to MA. bytes */
   int MAx,
   int FAy,
   int MAy)
/*
 * Times the kernel for out-of-cache (where flushelts sets the cache that it
 * is not allowed to be in) use.
 * RETURNS: elapsed time in seconds to average repitition of indicated problem.
 */
{
   #ifdef TREAL
      TYPE NONE = -1.0;
   #else
      TYPE NONE[2] = {-1.0, 0.0};
   #endif
   double t0, t1;
   TYPE *A, *X, *Y, *a, *x, *y;
   void *vmem;
   void (*mvsim)(size_t celts, size_t pgelts, ATL_CINT xu, ATL_CINT yu,
                 enum ATLAS_TRANS TA, ATL_CINT M, ATL_CINT N, const TYPE *A,
                 ATL_CINT lda, const TYPE *X, const TYPE *beta, TYPE *Y);
   ATL_INT Aelts, Xelts, Yelts, setspan, ygap, xgap, agap, pregap, setsz, nsets;
   ATL_INT i, j;
   size_t ptr_st;
   int maxalign;

   mvsim = (TA == AtlasNoTrans || TA == AtlasConj) ? mvsimN : mvsimT;
   if (MAx)
      assert(MAx != FAx);
   if (MAy)
      assert(MAy != FAy);
   if (MAa)
      assert(MAa != FAa);
/*
 * Find basic length of each operand in elements
 */
   Aelts = lda * N;
   Xelts = (TA == AtlasNoTrans || TA == AtlasConj) ? N : M;
   Yelts = (TA == AtlasNoTrans || TA == AtlasConj) ? M : N;
/*
 * Map memory so that we can enforce all required alignments while moving
 * through memory; mem starts with maxalign-aligned memory, so that we can
 * guarantee all further alignments
 */
   maxalign = (FAx >= FAa) ? FAx : FAa;
   maxalign = (maxalign >= FAy) ? maxalign : FAy;
   if (MAx | MAy | MAa)
   {
      maxalign = (MAx >= MAa) ? MAx : MAa;
      maxalign = (maxalign >= MAy) ? maxalign : MAy;
   }
   if (MAx)
   {
      j = (FAx) ? FAx : ATL_sizeof;
      for (i=0; (i % j != 0 || i%MAx == 0); i += ATL_sizeof);
      pregap = i;
   }
   else pregap = 0;
   xgap = ATL_MulBySize(Xelts);
   if (FAy || MAy)
   {
      j = (FAy) ? FAy : ATL_sizeof;
      if (MAy)
         for (i=pregap+xgap; (i%j != 0 || i%MAy == 0); i += ATL_sizeof);
      else
         for (i=pregap+xgap; (i%j != 0); i += ATL_sizeof);
      xgap = i - pregap;
   }
   ygap = ATL_MulBySize(Yelts);
   if (FAa || MAa)
   {
      j = (FAa) ? FAa : ATL_sizeof;
      if (MAa)
         for (i=pregap+xgap+ygap; (i%j != 0 || i%MAa == 0); i += ATL_sizeof);
      else
         for (i=pregap+xgap+ygap; (i%j != 0); i += ATL_sizeof);
      ygap = i - pregap - xgap;
   }
   agap = ATL_MulBySize(Aelts);

   if (maxalign)
   {
      j = pregap;
      for (i=pregap+xgap+ygap+agap; i%maxalign != 0; i++);
      agap = i-xgap-ygap;
   }
   setspan = xgap + ygap + agap;
   assert(setspan%ATL_sizeof == 0);
   setsz = ATL_MulBySize(M+N+M*N);
   nsets = (ATL_MulBySize(flushelts)+setsz-1)/setsz;
   if (!nsets)
      nsets = 1;
   vmem = malloc(maxalign + nsets*setspan);
   assert(vmem);
   if (maxalign)   /* start maxaligned to guarantee all alignments */
      for (ptr_st = (size_t)vmem; ptr_st%maxalign; ptr_st++);
   else ptr_st = (size_t) vmem;
   X = (TYPE*) (ptr_st + pregap);
   Y = (TYPE*) (ptr_st + pregap + xgap);
   A = (TYPE*) (ptr_st + pregap + xgap + ygap);
/*
 * Set ptrs to last set in memory
 */
   setspan /= ATL_sizeof;
   a = A += (nsets-1) * setspan;
   x = X += (nsets-1) * setspan;
   y = Y += (nsets-1) * setspan;
   for (i=nsets; i; i--)
   {
      #define DEBUG_FA
      #ifdef DEBUG_FA
         if (FAa)
            assert(((size_t)a)%FAa == 0);
         if (FAx)
            assert(((size_t)x)%FAx == 0);
         if (FAy)
            assert(((size_t)y)%FAy == 0);
         if (MAa)
            assert(((size_t)a)%MAa != 0);
         if (MAx)
            assert(((size_t)x)%MAx != 0);
         if (MAy)
            assert(((size_t)y)%MAy != 0);
      #endif
      Mjoin(PATL,gegen)(Yelts, 1, y, Yelts, M);
      Mjoin(PATL,gegen)(Xelts, 1, x, Xelts, N+127*50+77);
      if (i&1)
         Mjoin(PATL,scal)(Xelts, NONE, x, 1);
      Mjoin(PATL,gegen)(M, N, A, lda, N*M+513*7+90);
      a -= setspan; x -= setspan; y -= setspan;
   }
   a = A; x = X; y = Y;

   j=0;
   t0 = time00();
   for (i=nreps; i; i--)
   {
      mvsim(celts, pgelts, xu, yu, TA, M, N, a, lda, x, beta, y);
      if (++j < nsets) { a -= setspan; x -= setspan; y -= setspan; }
      else  { a = A; x = X; y = Y; j=0; }
   }
   t1 = time00();
   free(vmem);
   t1 = (t1-t0) / (1.0*nreps);
   if (verb)
      fprintf(stdout, "   M=%d, N=%d, lda=%d, nreps=%d, time=%e, mflop=%.2f\n",
              M, N, lda, nreps, t1, Time2Flop(M, N, t1));
   return(t1);
}

void DoTimes(int verb, ATL_INT flshelts, ATL_INT celts, ATL_INT pgelts,
             ATL_INT ntim, ATL_INT nrep, ATL_INT xu, ATL_INT yu,
             enum ATLAS_TRANS TA, ATL_INT M, ATL_INT N, ATL_INT lda, TYPE *beta,
             int FAa, int MAa, int FAx, int MAx, int FAy, int MAy)
{
   double *times;
   int i;

   times = malloc(ntim * sizeof(double));
   assert(times);

#ifdef TREAL
   fprintf(stdout,
           "GEMV: M=%d, N=%d, lda=%d, AF=[%d,%d,%d], AM=[%d,%d,%d], beta=%e:\n",
           M, N, lda, FAa, FAx, FAy, MAa, MAx, MAy, *beta);
#else
   fprintf(stdout,
      "GEMV: M=%d, N=%d, lda=%d, AF=[%d,%d,%d], AM=[%d,%d,%d], beta=[%e,%e]:\n",
           M, N, lda, FAa, FAx, FAy, MAa, MAx, MAy, *beta, beta[1]);
#endif
   for (i=0; i < ntim; i++)
      times[i] = mvtime_OC(nrep, flshelts, celts, pgelts, TA, M, N, lda, beta,
                           xu, yu, FAa, MAa, FAx, MAx, FAy, MAy);
   SortDoubles(ntim, times);
   Times2Flops(M, N, ntim, times);
   if (fpres)
   {
      #if defined(PentiumCPS) || defined(WALL)
         fprintf(fpres, "%d 1\n", ntim);
      #else
         fprintf(fpres, "%d 0\n", ntim);
      #endif
      for (i=0; i < ntim; i++)
         fprintf(fpres, "%le\n", times[i]);
      fclose(fpres);
   }
   fprintf(stdout, "NSAMPLES=%d, MAX=%.2f, MIN=%.2f, AVG=%.2f, MED=%.2f\n",
           ntim, times[0], times[ntim-1], mysum(ntim, times)/ntim,
           times[ntim>>1]);
   free(times);
}

void PrintUsage(char *name, char *arg, int i)
{
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' on %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s [flags], where flags are:\n", name);
   fprintf(stderr, "   -v <#> : set verbosity level\n");
   fprintf(stderr, "   -C <#> : set flushsz = # (kbytes)\n");
   fprintf(stderr, "   -p <#> : set pagesz = # (kbytes)\n");
   fprintf(stderr,
       "   -G <#> : set GEMV cache size (for blocking) to # (kbytes)\n");
   fprintf(stderr, "   -x <#> : unrolling for X in kernel is #\n");
   fprintf(stderr, "   -y <#> : unrolling for Y in kernel is #\n");
   fprintf(stderr, "   -m <#> : set # of rows of matrix to #\n");
   fprintf(stderr, "   -n <#> : set # of cols of matrix to #\n");
   fprintf(stderr, "   -l <#> : set leading dimension of array A to #\n");
   fprintf(stderr, "   -F <#> : do at least # MFLOPS for each timing interval\n");
   fprintf(stderr, "   -f <file> : output timing summary in <file>; if file exists read & report\n");
   fprintf(stderr, "   -A n/t/c/z : set transpose (z = Conj, NoTrans)\n");
   fprintf(stderr,
           "   -r <#> : do # repetitions of the call for each timing interval\n");
   fprintf(stderr,
      "   -# <#> : report # timings (each interval may have multiple calls)\n");
   fprintf(stderr,
"   -F[x,y,a] <#> : if(# > 0) -> force op to be aligned to at least # bytes\n");
   fprintf(stderr,
"                   if(# < 0) -> force op to be aligned to < # bytes.\n");
   fprintf(stderr, "   -b <beta> : 2 floats for complex, one for real.\n");
   exit(i ? i : -1);
}

void GetFlags(int nargs, char **args, int *verb,
              ATL_INT *flushelts, ATL_INT *celts, ATL_INT *pgelts,
              ATL_INT *xu, ATL_INT *yu, ATL_INT *ntim, ATL_INT *nrep,
              enum ATLAS_TRANS *TA, ATL_INT *m, ATL_INT *n, ATL_INT *lda,
              TYPE *beta,
              int *FAa, int *MAa, int *FAx, int *MAx, int *FAy, int *MAy)
{
   double mfF=ATL_nkflop/1000.0, flops;
   ATL_INT j, h;
   int i;
   char ch;

   *verb = 1;
   #ifdef ATL_PAGESZ
      *pgelts = ATL_DivBySize(ATL_PAGESZ);
   #else
      *pgelts = 4*ATL_DivBySize(1024);
   #endif
   *celts = 0.75*ATL_L1elts;
   *flushelts = 8*1024*ATL_DivBySize(1024);
   *xu = *yu = 1;
   *m = 800;
   *n = 200;
   *nrep = *lda = 0;
   *ntim = 3;
   *FAa = *MAa = *FAx = *MAx = *FAy = *MAy = 0;
   *beta = 1.0;
   #ifdef TCPLX
      beta[1] = 0.0;
   #endif

   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "No '-' preceeding flag!", i);
      switch(args[i][1])
      {
      case 'f' :  /* set resfile output */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -f ", i-1);
         resfile = args[i];
         break;
      case 'v' :  /* set verbosity */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -v ", i-1);
         *verb = atoi(args[i]);
         break;
      case 'G' :  /* set GEMV blocking cache size in KB */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -G ", i-1);
         j = atoi(args[i]);
         *celts = j*ATL_DivBySize(1024);
         break;
      case 'A' :  /* set transpose */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -A ", i-1);
         ch = args[i][0];
         if (ch == 't' || ch == 'T')
            *TA = AtlasTrans;
         else if (ch == 'c' || ch == 'C')
            *TA = AtlasConjTrans;
         else if (ch == 'z' || ch == 'Z')
            *TA = AtlasConj;
         else
            *TA = AtlasNoTrans;
         break;
      case 'C' :  /* set flushsz in KB */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -C ", i-1);
         j = atoi(args[i]);
         *flushelts = j*ATL_DivBySize(1024);
         break;
      case 'p' :  /* set pagesz in KB */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -p ", i-1);
         j = atoi(args[i]);
         *pgelts = j*ATL_DivBySize(1024);
         break;
      case 'x' :  /* set xu */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -x ", i-1);
         *xu = atoi(args[i]);
         break;
      case 'y' :  /* set yu */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -y ", i-1);
         *yu = atoi(args[i]);
         break;
      case 'm' :  /* set M */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -m ", i-1);
         *m = atoi(args[i]);
         break;
      case 'n' :  /* set N */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -n ", i-1);
         *n = atoi(args[i]);
         break;
      case 'l' :  /* set lda */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -l ", i-1);
         *lda = atoi(args[i]);
         break;
      case 'a' : /* alias for setting alpha in r1ktime */
      case 'b' : /* set beta */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -b ", i-1);
         *beta = atof(args[i]);
         #ifdef TCPLX
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -b ", i-1);
            beta[1] = atof(args[i]);
         #endif
         break;
      case 'F' :  /* set nrep by specifying MFLOPS, or force alignment */
         ch = args[i][2];
         if (!ch)   /* specifying MFLOPS */
         {
            if (++i >= nargs)
               PrintUsage(args[0], "out of flags in -F ", i-1);
            j = atoi(args[i]);
            mfF = j;
         }
         else
         {
            if (ch != 'a' && ch != 'y' && ch != 'x')
               PrintUsage(args[0], args[i], i);
            if (++i >= nargs)
               PrintUsage(args[0], args[i-1], i-1);
            j = atoi(args[i]);
            if (j < 0)
            {
               if (ch == 'a')
                  *MAa = -j;
               else if (ch == 'y')
                  *MAy = -j;
               else if (ch == 'x')
                  *MAx = -j;
            }
            else
            {
               if (ch == 'a')
                  *FAa = j;
               else if (ch == 'y')
                  *FAy = j;
               else if (ch == 'x')
                  *FAx = j;
            }
         }
         break;
      case 'r' :  /* set nrep directly as integer */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -r ", i-1);
         *nrep = atoi(args[i]);
         break;
      case '#' :  /* set number of timings to report */
         if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -# ", i-1);
         *ntim = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], args[i], i);
      }
   }
   if (!(*nrep))
   {
      flops = Time2Flop(*m, *n, 1.0) * 1000.0;  /* Get kiloFLOPS in GEMV */
      *nrep = (mfF+flops-1)/flops;
      if (*nrep < 1) *nrep = 1;
   }
   if (!(*lda))
      *lda = *m + 8;
}
int main(int nargs, char **args)
{
   ATL_INT flushelts, celts, pgelts, xu, yu, ntim, nrep, m, n, lda;
   int FAa, MAa, FAx, MAx, FAy, MAy;    /* Force & Max align for ops */
   int verb;
   enum ATLAS_TRANS TA;
   double *dres;
   #ifdef TREAL
      TYPE beta;
   #else
      TYPE beta[2];
   #endif

   GetFlags(nargs, args, &verb, &flushelts, &celts, &pgelts, &xu, &yu, &ntim,
            &nrep, &TA, &m, &n, &lda, SADD beta,
            &FAa, &MAa, &FAx, &MAx, &FAy, &MAy);
   if (resfile)
   {
      dres = ReadResultsFile(1, resfile);
      if (dres)
      {
         fprintf(stdout, "TIMINGS READ IN FROM '%s':\n", resfile);
         PrintResultsFromFile(stdout, dres);
         free(dres);
         exit(0);
      }
      fpres = fopen(resfile, "w");
      assert(fpres);
   }
   DoTimes(verb, flushelts, celts, pgelts, ntim, nrep, xu, yu, TA, m, n, lda,
           SADD beta, FAa, MAa, FAx, MAx, FAy, MAy);
   exit(0);
}
