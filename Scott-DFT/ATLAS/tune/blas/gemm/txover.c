/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 * Copyright (C) 2010, 2009, 2010, 2009 R. Clint Whaley
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
#include <math.h>
#include <limits.h>
#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_f77.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),sysinfo.h))

#define dumb_seed(iseed_) srand(iseed_)
#ifndef RAND_MAX  /* rather dangerous non-ansi workaround */
   #define RAND_MAX ((unsigned long)(1<<30))
#endif
#define dumb_rand() ( 0.5 - ((double)rand())/((double)RAND_MAX) )


___main(){}
__main(){}
MAIN__(){}
_MAIN_(){}

double time00();
int Mjoin(PATL,tgemm_p)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,tgemm)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
int Mjoin(PATL,gemm_serial)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   Mjoin(PATL,gemm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   return(1);
}

typedef int (*GEMMPTR)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
#define NSHAPES 7
static char *shapesuff[NSHAPES] =
   {"SQmnk", "SmnLk", "SmkLn", "SnkLm", "SmLnk", "SnLmk", "SkLmn"};
enum ATLAS_MATSHAPE
{SquareMNK=0,  /* square, including recursive where all dim cut */
 ShortMN_LK=1, ShortMK_LN=2, ShortNK_LM=3,/* recursive shapes, one dim uncut */
 ShortM_LNK=4, ShortN_LMK=5, ShortK_LMN=6};/* 1-D staticly blocked shapes */
int ATL_FINDP = ATL_NCPU;
int ATL_OLDFINDP = 1;
static int *ATL_Ps;
static int ATL_nP;
static double ATL_CONMFLOP, ATL_SHPMFLOP=0;

double mmcase
(
   enum ATLAS_TRANS TA,                 /* Transpose setting for A */
   enum ATLAS_TRANS TB,                 /* Transpose setting for B */
   size_t flushsz,                      /* size of flush area */
   GEMMPTR gemm,                        /* ptr to matmul to time */
   double mflopF,                       /* how many mflops to force */
   ATL_CINT M, ATL_CINT N, ATL_CINT K,  /* matrix dimensions */
   int *np  /* number of procs returned by gemm call */
)
{
   #ifdef TCPLX
      const TYPE alpha[2] = {ATL_rone, ATL_rzero},
                 beta[2] = {ATL_rnone, ATL_rzero};
   #else
      const TYPE alpha = ATL_rone, beta = ATL_rnone;
   #endif
   ATL_INT lda, ldb, ldc, asize, bsize;
   size_t setsz, nsets, incsz, nrep, i, j;
   double t0, t1;
   TYPE *A, *B, *C, *mp;
   char *sp;
   static int cnt=0;

   lda = (TA == AtlasNoTrans) ? M : K;
   ldb = (TB == AtlasNoTrans) ? K : N;
   ldc = M;
   lda += 300;
   ldb += 300;
   ldc += 300;
   setsz = M*N + M*K + K*N;
   nsets = (ATL_DivBySize(flushsz)+setsz-1)/setsz;
   if (nsets < 1)
      nsets = 1;
   asize = (TA == AtlasNoTrans) ? lda*K : lda*M;
   bsize = (TB == AtlasNoTrans) ? ldb*N : ldb*K;
   incsz = ldc*N + asize + bsize;
   do
   {
      mp = malloc(nsets*ATL_MulBySize(incsz));
      if (!mp)
         nsets--;
   }
   while (!mp && nsets > 1);
   assert(mp);
   t0 = (mflopF*1000000.0)/(2.0*M*N*K);
   nrep = t0;
   if (t0-nrep > 0.45)
      nrep++;
   if (nrep < 1) nrep = 1;
   Mjoin(PATL,gegen)(nsets*incsz, 1, mp, nsets*incsz, M*N+lda);

   t0 = time00();
   for (i=0; i < nrep; i++)
   {
      A = ((i/nsets)*nsets == i) ? mp : A + incsz;
      B = A + asize;
      C = B + bsize;
      *np = gemm(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   t1 = time00() - t0;
   free(mp);
   if (gemm == Mjoin(PATL,tgemm_p))
   {
      sp = "parall";
      t0 = (*np > 1) ? (2.0*M*N*K*nrep) / (1000000.0*t1) : 0.0;
   }
   else
   {
      sp = "serial";
      t0 = (2.0*M*N*K*nrep) / (1000000.0*t1);
   }
   printf("%6d %6.6s %6d %6d %6d %6d %9.2f\n", ++cnt, sp, *np, M, N, K, t0);
   return(t0);
}

int FindBestNP
(
   enum ATLAS_TRANS TA,                 /* Transpose setting for A */
   enum ATLAS_TRANS TB,                 /* Transpose setting for B */
   size_t flushsz,                      /* size of flush area */
   double mflopF,                       /* how many mflops to force */
   ATL_CINT M, ATL_CINT N, ATL_CINT K,  /* matrix dimensions */
   int P                                /* max # of procs to try */
)
/*
 * This routine tries using a varying number of processors to decompose
 * the problem.  It tries all powers of 2 between 1 and NP.
 * RETURNS: optimal number of processors to use for this problem
 * NOTE: give serial code 3% advantage, since parallel code should be
 *       noticably better, and is prone to large variance
 */
{
   int j, p, npB;
   double mfB, mf, mfS;

   assert(P > 1);
   npB = ATL_FINDP = P;
   mfB = mmcase(TA, TB, flushsz, Mjoin(PATL,tgemm_p), mflopF, M, N, K, &p);
   mfS = (ATL_SHPMFLOP > 0.0) ? ATL_SHPMFLOP : ATL_CONMFLOP;
   if (p < P || mfB < mfS)
   {
      mf = mmcase(TA, TB, flushsz, Mjoin(PATL,gemm_serial), mflopF, M, N, K,
                  &p) * 1.03;   /* give serial code 3% advantage */
      ATL_SHPMFLOP = mf;
      if (mf > mfB) { npB = 1; mfB = mf; }
   }
   for (j=2; j < P; j *= 2);
   for (j >>= 1; j >= 2; j >>= 1)
   {
      ATL_FINDP = j;
      mf = mmcase(TA, TB, flushsz, Mjoin(PATL,tgemm_p), mflopF, M, N, K, &p);
      if (mf > mfB) { npB = j; mfB = mf; }
      if (mf > (j>>1)*mfS)
         break;
   }
   return(npB);
}

void GetDims
(
   enum ATLAS_MATSHAPE shape,  /* matrix shape */
   int restS,   /* size for the remaining dimenions to be restricted to */
   int growS,   /* size for the growing dimension(s) to be set to */
   int *M, int *N, int *K
)
{
   switch(shape)
   {
    case ShortMN_LK:
       *M = *N = restS;
       *K = growS;
       return;
    case ShortMK_LN:
       *M = *K = restS;
       *N = growS;
       return;
    case ShortNK_LM:
       *N = *K = restS;
       *M = growS;
       return;
    case ShortM_LNK:
       *M = restS;
       *N = *K = growS;
       return;
    case ShortN_LMK:
       *N = restS;
       *M = *K = growS;
       return;
    case ShortK_LMN:
       *K = restS;
       *M = *N = growS;
       return;
    default:   /* square */
      *M = *N = *K = growS;
   }
}

int GetXover
(
   size_t flushsz, double mflopF,
   enum ATLAS_TRANS TA,                 /* Transpose setting for A */
   enum ATLAS_TRANS TB,                 /* Transpose setting for B */
   enum ATLAS_MATSHAPE shape,
   ATL_CINT minN,   /* smallest size for growing dimension(s) */
   ATL_CINT maxN,   /* largest size for growing dimensions(s) */
   ATL_CINT NN,     /* size for non-growing dimension(s) */
   ATL_CINT P       /* max number of processors to use */
)
/*
 * This routine finds point between minN and maxN where the optimal number
 * of processors goes down.  It uses recursive halving, and assumes that a
 * given high value of P will win for all N > i, if it wins for i.
 * RETURNS: size where higher P first beats all other # procs between
 *          [minN,maxN]; 0 if P processors never wins.
 */
{
   ATL_INT M, N, K, j, mymax, mymin, gap, k;
   int np, itmp;

/*
 * If we can't ever use all processors, then no crossover for this P
 */
   GetDims(shape, NN, maxN, &M, &N, &K);
   np = FindBestNP(TA, TB, flushsz, mflopF, M, N, K, P);
   if (np < P)
      return(0);
/*
 * Make sure parallel gemm isn't already faster on smallest case
 */
   GetDims(shape, NN, minN, &M, &N, &K);
   np = FindBestNP(TA, TB, flushsz, mflopF, M, N, K, P);
   if (np == P)
      return(minN);
/*
 * Begin recursive halving search
 */
   mymax = maxN;
   mymin = minN;
   do
   {
      j = ((mymax-mymin)>>1) + mymin;
      GetDims(shape, NN, j, &M, &N, &K);
      np = FindBestNP(TA, TB, flushsz, mflopF, M, N, K, P);
      if (np == P)  /* full P still winning at this size */
         mymax = j;
      else          /* P lost, try larger problems*/
         mymin = j;
      if (j < 100 || shape == SquareMNK)
         gap = 8;
      else if (j < 1000)
         gap = 16;
      else if (j < 2000)
         gap = 32;
      else
         gap = 128;
   }
   while (mymax-mymin > gap);
   return((np == P) ? j : mymax);
}

int *GetXoversForAllP
(
   size_t flushsz, double mflopF,
   enum ATLAS_TRANS TA,                 /* Transpose setting for A */
   enum ATLAS_TRANS TB,                 /* Transpose setting for B */
   enum ATLAS_MATSHAPE shape,
   ATL_CINT minN,   /* smallest size for growing dimension(s) */
   ATL_CINT maxN0,  /* largest size for growing dimensions(s) */
   ATL_CINT NN      /* size for non-growing dimension(s) */
)
/*
 * RETURNS: a log2(NCPU) length array; each ith entry is the crossover point
 *          where using p=2^i provided better performance than any lesser
 *          amount of parallelism.  A 0 means that this number of processors
 *          was not useful.  The last entry is always NCPU, even if this
 *          is not a power of 2
 */
{
   static int *Ps, *xos;
   static int log2p=0;
   int i, j, maxN=maxN0;

   if (log2p == 0)
   {
      for (log2p=1; (1<<log2p) < ATL_NCPU; log2p++);
      Ps = malloc(sizeof(int)*log2p*2);
      xos = Ps + log2p;
      Ps[log2p-1] = ATL_NCPU;
      for (j=log2p-2; j >= 0; j--)
         Ps[j] = 2<<j;
   }
   for (j=log2p-1; j >= 0; j--)
   {
      if (maxN > minN)
      {
         xos[j] = GetXover(flushsz, mflopF, TA, TB, shape, minN, maxN,NN,Ps[j]);
         if (xos[j])
            maxN = xos[j];
      }
      else
        xos[j] = 0;
   }
   return(xos);
}

int GetShapesMaxN(enum ATLAS_TRANS TA, enum ATLAS_TRANS TB,
                  enum ATLAS_MATSHAPE shape)

{
   switch(shape)
   {
    case ShortMN_LK:
       if (TA == AtlasNoTrans || TB == AtlasTrans)
       {
          #ifdef TREAL
             return(8000);
          #elif defined(SCPLX)
             return(8000);
          #else
             return(4000);
          #endif
       }
    case ShortMK_LN:
    case ShortNK_LM:
       #ifdef TREAL
          return(10000);
       #elif defined(SCPLX)
          return(8000);
       #else
          return(6000);
       #endif
    case ShortM_LNK:
    case ShortN_LMK:
    case ShortK_LMN:
       #ifdef SREAL
          return(6000);
       #elif defined(SCPLX)
          return(4000);
       #elif defined(DCPLX)
          return(4000);
       #else
          return(6000);
       #endif
    default:   /* square */
       return(1000);
   }
   return(2000);
}

ATL_INT *GetXoversByShape
(
   size_t flushsz,              /* size of cache to flush */
   double mflopF,               /* mflops to force for timing resolution */
   enum ATLAS_TRANS TA,         /* Transpose setting for A */
   enum ATLAS_TRANS TB,         /* Transpose setting for B */
   enum ATLAS_MATSHAPE shape    /* shape of the matrix */
)
/*
 * RETURNS: log2(p)*9 array of crossover points for different processor counts.
 *          The contiguous dim is log2(P) crossover points, and tracks the
 *          crossover points for p=2, 4, 8, P (doubling the number of cores
 *          at each step, until the maximum is reached).  The noncontiguous
 *          dimensions gives the size of the restricted dimensions, where
 *          are set to 1<<i, 0 <= i < 9;
 */
{
   ATL_INT *xos, i, k, j, *ix, log2p, minN, maxN, *ip;

   printf("\n\nFINDING XOVERS FOR %s\n", shapesuff[shape]);
   printf(" COUNT  WHICH      P      M      N      K     MFLOP\n");
   printf("====== ====== ====== ====== ====== ====== =========\n");

   minN = Mjoin(PATL,GetNB)();
   maxN = GetShapesMaxN(TA, TB, shape);

   for (log2p=1; (1<<log2p) < ATL_NCPU; log2p++);
   ix = calloc(9*log2p, sizeof(ATL_INT));
   assert(ix);
   for (i=8; i >= 0; i--)
   {
      ATL_SHPMFLOP = 0.0;
/*
 *    Get crossover points for all power-of-2 P
 */
      xos = GetXoversForAllP(flushsz, mflopF, TA, TB, shape, minN, maxN, 1<<i);
      ip = ix + log2p*i;
      j = 0;
      for (k=0; k < log2p; k++)
      {
         if (xos[k])
         {
            if (j)
               j = Mmin(j, xos[k]);
            else
               j = xos[k];
         }
         ip[k] = xos[k];
      }
/*
 *    New minN is the crossover of this time, since we go from biggest to
 *    smallest on other dims
 */
      if (j)
         minN = j;
      if (!j)       /* if we can't thread this size, can't thread smaller */
         break;
   }
   return(ix);
}

void FindXovers
(
   size_t flushsz,              /* size of cache to flush */
   double mflopF,               /* mflops to force for timing resolution */
   FILE *fpout                  /* file to print output to */
)
{
   int i, k, t, u, l2p;
   ATL_INT *xo;
   enum ATLAS_TRANS tr[2] = {AtlasNoTrans, AtlasTrans}, TA, TB;
   char ctr[2] = {'N', 'T'}, cta, ctb;

   for (l2p=1; (1<<l2p) < ATL_NCPU; l2p++);
   fprintf(fpout, "/* This file generated by %s\n */\n", __FILE__);
   fprintf(fpout, "#ifndef ATL_TXOVER_H\n   #define ATL_TXOVER_H 1\n\n");
   fprintf(fpout, "   #define ATL_PDIM %d\n", l2p);
   for (i=0; i < NSHAPES; i++)
   {
      for (t=0; t < 2; t++)
      {
         TA = tr[t];
         cta = ctr[t];
         for (u=0; u < 2; u++)
         {
            TB = tr[u];
            ctb = ctr[u];
            if (i == ShortMK_LN)
            {
               if (TA == AtlasTrans)
               {
                  fprintf(fpout,
                          "#define ATL_tmmT%c_%s_XO ATL_tmmN%c_%s_XO\n",
                          ctb, shapesuff[i], ctb, shapesuff[i]);
                  continue;
               }
            }
            else if (i == ShortMK_LN)  /* TA not important */
            {
               if (TA == AtlasTrans)
               {
                  fprintf(fpout,
                          "#define ATL_tmmT%c_%s_XO ATL_tmmN%c_%s_XO\n",
                          ctb, shapesuff[i], ctb, shapesuff[i]);
                  continue;
               }
            }
            else if (i == ShortNK_LM)  /* TB not important */
            {
               if (TB == AtlasTrans)
               {
                  fprintf(fpout,
                          "#define ATL_tmm%cT_%s_XO ATL_tmm%cN_%s_XO\n",
                          cta, shapesuff[i], cta, shapesuff[i]);
                  continue;
               }
            }
/*
 *          If M == N, then NN case is roughly equivalent to TT case
 */
            else if (i == SquareMNK || i == ShortMN_LK || i == ShortK_LMN)
            {
               if (TA == AtlasTrans && TB == AtlasTrans)
               {
                  fprintf(fpout,
                          "#define ATL_tmmTT_%s_XO ATL_tmmNN_%s_XO\n",
                          shapesuff[i], shapesuff[i]);
                  continue;
               }
            }
            if (i == SquareMNK)
            {
               xo = GetXoversForAllP(flushsz, mflopF, TA, TB, i,
                                      Mjoin(PATL,GetNB)(),
                                      GetShapesMaxN(TA, TB, i), 0);
               fprintf(fpout, "static const int ATL_tmm%c%c_%s_XO[%d] = \n   {",
                       cta, ctb, shapesuff[i], l2p);
               for (k=0; k < l2p-1; k++)
                  fprintf(fpout, "%ld,", xo[k]);
               fprintf(fpout, "%d};\n", xo[l2p-1]);
            }
            else
            {
               xo = GetXoversByShape(flushsz, mflopF, TA, TB, i);
               fprintf(fpout, "static const int ATL_tmm%c%c_%s_XO[%d] = \n   {",
                       cta, ctb, shapesuff[i], l2p*9);
               for (k=0; k < 9*l2p-1; k++)
                  fprintf(fpout, "%ld,", xo[k]);
               fprintf(fpout, "%d};\n", xo[9*l2p-1]);
               free(xo);
            }
         }
      }
      fprintf(fpout, "static const int *ATL_tmm_%s_XO[4] =\n", shapesuff[i]);
      fprintf(fpout, "{ATL_tmmNN_%s_XO, ATL_tmmNT_%s_XO,\n",
              shapesuff[i], shapesuff[i]);
      fprintf(fpout, " ATL_tmmTN_%s_XO, ATL_tmmTT_%s_XO};\n",
              shapesuff[i], shapesuff[i]);
   }
   fprintf(fpout, "\n#endif /* end ifndef ATL_TXOVER_H */\n");
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);

   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -F <mflop> : force <mflops> of timed computation\n");
   fprintf(stderr, "   -f <flushKB> : flush at least this mem in LRU timers\n");
   fprintf(stderr, "   -o <file> : write header file to <file>\n");

   exit(ierr ? ierr : -1);
}

FILE *GetFlags(int nargs, char **args, int *flushKB, int *mflopF)
{
   int i;
   char *fout;
   FILE *fpout;
   #ifdef L2SIZE
      *flushKB = L2SIZE;
   #else
      *flushKB = 4*1024;
   #endif
   *mflopF = 1;
   fout = Mstr(Mjoin(res/atlas_,Mjoin(Mjoin(Mjoin(PRE,tXover_),ATL_NCPU),p.h)));
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'o':   /* -o <file> *f (++i >= nargs) */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         fout = args[i];
         break;
      case 'f':                         /* -f <flushKB> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *flushKB = atoi(args[i]);
         break;
      case 'F':                         /* -F <mflop> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *mflopF = atoi(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   fpout = fopen(fout, "w");
   assert(fpout);
   return(fpout);
}

main(int nargs, char **args)
{
   int flushKB = 4*1024, mflopF=1;
   FILE *fpout;
   int i;

#if ATL_NCPU > 1
   fpout = GetFlags(nargs, args, &flushKB, &mflopF);
   for (i=1; (1<<i) < ATL_NCPU; i++);
   ATL_nP = i;
   ATL_Ps = malloc(ATL_nP*sizeof(int));
   assert(ATL_Ps);
   ATL_Ps[ATL_nP-1] = ATL_NCPU;
   for (i=ATL_nP-2; i >= 0; i--)
      ATL_Ps[i] = (2<<i);
   ATL_CONMFLOP = mmcase(AtlasNoTrans, AtlasNoTrans, flushKB,
                         Mjoin(PATL,gemm_serial), mflopF, 1000, 1000, 1000, &i);
   FindXovers(flushKB*1024, mflopF, fpout);
#endif
   exit(0);
}
