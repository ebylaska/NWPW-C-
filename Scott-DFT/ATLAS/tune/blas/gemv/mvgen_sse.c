/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 *                    (C) Copyright 2000 R. Clint Whaley
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

/*
 * NOTE: do prefetch by having small program that reads output files, and
 *       simply schedules the prefetch inst within the code, by searching
 *       for particular inst in standard output.  Can print markers in
 *       generated code ID regions of interest
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "atlas_misc.h"
#include "atlas_r1parse.h"

#define NUMFLAGS        8
char *F_exp[NUMFLAGS] =
{
NULL,
NULL,
"X&Y are copied into all legal alignments",
"X&Y type-aligned, LDA*sizeof multiple of 16",
"A&X 16-byte aligned",
"preload nmu*nu rows of A ahead",
"preload nmu*nu elts of X ahead",
"Schedule X loads before compute block",
"Schedule 1st row of A loads before compute block",
"Generate single precision, not double",
"Use A ptrs rather than lda to index columns of A"
};

#define F_PRELOADA      2    /* preload mu*nmu elts of A ahead */
#define F_PRELOADX      3    /* preload mu*nmu elts of X ahead */
#define F_LDXSEP        4    /* load X before starting computation */
#define F_LDASEP        5    /* load 1st row of A before starting computation */

#define SET_FLAG(bits_, flg_, val_) \
{\
   if (val_) (bits_) |= (1<<(flg_)); \
   else (bits_) &= ~(1<<(flg_)); \
}

#define TEST_FLAG(bits_, flg_) ( (bits_) & (1<<(flg_)) )


int GetPower2(int n)
{
   int pwr2, i;

   if (n == 1) return(0);
   for (pwr2=0, i=1; i < n; i <<= 1, pwr2++);
   if (i != n) pwr2 = 0;
   return(pwr2);
}

#define ShiftThresh 2
char *GetDiv(int N, char *inc)
/*
 * Given a runtime variable whose name is in the string inc that you want
 * to divide by the compile-time constant N, produces the appropriate shift
 * (N is power of 2) or division (N not a power of 2).
 */
{
   static char ln[256];
   int pwr2 = GetPower2(N);
   if (N == 1) sprintf(ln, "%s", inc);
   else if (pwr2) sprintf(ln, "((%s) >> %d)", inc, pwr2);
   else sprintf(ln, "((%s) / %d)", inc, N);
   return(ln);
}

char *GetMul(int N, char *inc)
/*
 * let inc be a runtime variable that you wish to multiply by the
 * compile-time constant N.  This routine attempts to use at most
 * ShiftThresh shifts and adds instead of using a multiply.  If more
 * than ShiftThresh adds are required, just uses multiply as normal.
 */
{
   static char ln0[256];
   char ln[256];
   char *p=ln;
   int i, n=N, iPLUS=0;

   if (n == 0)
   {
      ln[0] = '0';
      ln[1] = '\0';
   }
   while(n > 1)
   {
      for (i=0; n >= (1<<i); i++);
      if ( (1 << i) > n) i--;
      if (iPLUS++) *p++ = '+';
      sprintf(p, "((%s) << %d)", inc, i);
      p += strlen(p);
      n -= (1 << i);
   }
   if (n == 1)
   {
      if (iPLUS++) *p++ = '+';
      sprintf(p, "%s", inc);
   }
   if (iPLUS > ShiftThresh) sprintf(ln0, "(%d*(%s))", N, inc);
   else if (iPLUS) sprintf(ln0, "(%s)", ln);
   else sprintf(ln0, "%s", ln);
   return(ln0);
}

char *GetAadd(
   int flag,
   int I,       /* constant value to go to unrolled row */
   int J,       /* constant value to go to unrolled col */
   int UseI)    /* Should we add "+i" to ptr ref? */
{
   static char ln[256];
   char *sp = (UseI) ? "+i" : "";

   if (TEST_FLAG(flag, R1F_APTRS))  /* use ptrs to index columns */
   {
      if (I)
         sprintf(ln, "A%d+%d%s", J, I, sp);
      else if (UseI)
         sprintf(ln, "A%d+i", J);
      else
         sprintf(ln, "A%d", J);
   }
   else  /* use lda to index columns */
   {
      if (J)
      {
         if (I)
            sprintf(ln, "A+lda%d+%d%s", J, I, sp);
         else
            sprintf(ln, "A+lda%d%s", J, sp);
      }
      else
      {
         if (I)
            sprintf(ln, "A+%d%s", I, sp);
         else if (UseI)
            sprintf(ln, "A%s", sp);
         else
            sprintf(ln, "A");
      }
   }
   return(ln);
}

char **GetLDAs(int nu)
{
   static char **ldas=NULL;
   static int nu0=0;
   int j;

   if (nu0 == nu)
      return(ldas);
   if (nu0)
   {
      for (j=0; j < nu0; j++)
         free(ldas[j]);
      free(ldas);
      ldas = NULL;
   }
   nu0 = nu;
   if (nu)
   {
      ldas = malloc(nu*sizeof(char*));
      assert(ldas);
      ldas[0] = malloc(sizeof(char));
      assert(ldas[0]);
      ldas[0][0] = '\0';
      for (j=1; j < nu; j++)
      {
         ldas[j] = malloc(8*sizeof(char));
         assert(ldas[j]);
         ldas[j][0] = '+'; ldas[j][1] = 'l'; ldas[j][2] = 'd'; ldas[j][3] = 'a';
         sprintf(ldas[j]+4, "%d", j);
      }
   }
   return(ldas);
}

void GetLdStPre(int flag, char **Xld, char **Ald, char **Ast, char *pre)
{
   if (TEST_FLAG(flag, R1F_ALGLDA) | TEST_FLAG(flag, R1F_ALLALIGNXY))
   {
      *Xld = "_mm_load_p";
      *Ald = "_mm_load_p";
      *Ast = "_mm_store_p";
   }
   else
   {
      *Xld = "_mm_loadu_p";
      *Ald = "_mm_loadu_p";
      *Ast = "_mm_storeu_p";
   }
   *pre = (TEST_FLAG(flag, R1F_SINGLE)) ? 's' : 'd';
}

void gen_sMUxNU0(FILE *fpout, char *spc, int flag, int mu, int nu, int I0)
{
   int incM;
   int i, j;
   char pre;

   pre = TEST_FLAG(flag, R1F_SINGLE) ? 's' : 'd';
   incM = (pre == 'd') ? 2 : 4;
   for (i=0; i < mu; i++)
   {
      for (j=0; j < nu; j++)
      {
         if (!j)
            fprintf(fpout, "%sx%d = _mm_load_s%c(X+i+%d);\n", spc, i, pre, i+I0);
         if (!i)
         {
            fprintf(fpout, "%sy%d = _mm_load_s%c(%s);\n",
                    spc, j, pre, GetAadd(flag, I0, j, 1));
            fprintf(fpout, "%sy%d = _mm_mul_s%c(y%d, x0);\n", spc, j, pre, j);
         }
         else
         {
            fprintf(fpout, "%sa%d_%d = _mm_load_s%c(%s);\n",
                    spc, i*incM, j, pre, GetAadd(flag, i+I0, j, 1));
            fprintf(fpout, "%sa%d_%d = _mm_mul_s%c(a%d_%d, x%d);\n",
                    spc, i*incM, j, pre, i*incM, j, i*incM);
            fprintf(fpout, "%sy%d = _mm_add_s%c(y%d, a%d_%d);\n\n",
                    spc, j, pre, j, i, j);
         }
      }
   }
}

void gen_sMUxNU(FILE *fpout, char *spc, int flag, int mu, int nu, int I0)
{
   int i, j, incM;
   char pre;

   pre = TEST_FLAG(flag, R1F_SINGLE) ? 's' : 'd';
   incM = (pre == 'd') ? 2 : 4;
   for (i=0; i < mu; i++)
   {
      for (j=0; j < nu; j++)
      {
         if (!j)
            fprintf(fpout, "%sx%d = _mm_load_s%c(X+i+%d);\n", spc, i, pre, i+I0);
         fprintf(fpout, "%sa%d_%d = _mm_load_s%c(%s);\n",
                 spc, i*incM, j, pre, GetAadd(flag, i+I0, j, 1));
         fprintf(fpout, "%sa%d_%d = _mm_mul_s%c(a%d_%d, x%d);\n",
                 spc, i*incM, j, pre, i*incM, j, i*incM);
         fprintf(fpout, "%sy%d = _mm_add_s%c(y%d, a%d_%d);\n\n",
                 spc, j, pre, j, i, j);
      }
   }
}

void gen_MUxNU0(FILE *fpout, char *spc, int flag, int mu, int nu, int I0)
{
   char *Xld, *Ald, *Ast;
   char pre;
   int i, j, incM;

   GetLdStPre(flag, &Xld, &Ald, &Ast, &pre);
   incM = (pre == 'd') ? 2 : 4;
   assert(mu%incM == 0);

   for (i=0; i < mu; i += incM)
   {
      for (j=0; j < nu; j++)
      {
         if (!j)
            fprintf(fpout, "%sx%d = %s%c(X+i+%d);\n",
                    spc, i, Xld, pre, i+I0);
         if (!i)
         {
            fprintf(fpout, "%sy%d = %s%c(%s);\n",
                    spc, j, Ald, pre, GetAadd(flag, I0, j, 1));
            fprintf(fpout, "%sy%d = _mm_mul_p%c(y%d, x0);\n", spc, j, pre, j);
         }
         else
         {
            fprintf(fpout, "%sa%d_%d = %s%c(%s);\n",
                    spc, i, j, Ald, pre, GetAadd(flag, i+I0, j, 1));
            fprintf(fpout, "%sa%d_%d = _mm_mul_p%c(a%d_%d, x%d);\n",
                    spc, i, j, pre, i, j, i);
            fprintf(fpout, "%sy%d = _mm_add_p%c(y%d, a%d_%d);\n\n",
                    spc, j, pre, j, i, j);
         }
      }
   }
}

void gen_MUxNU(FILE *fpout, char *spc, int flag, int mu, int nu, int I0)
{
   char *Xld, *Ald, *Ast;
   char pre;
   int i, j, incM;

   GetLdStPre(flag, &Xld, &Ald, &Ast, &pre);
   incM = (pre == 'd') ? 2 : 4;
   assert(mu%incM == 0);

   for (i=0; i < mu; i += incM)
   {
      fprintf(fpout, "%sx%d = %s%c(X+i+%d);\n", spc, i, Xld, pre, i+I0);
      for (j=0; j < nu; j++)
      {
         fprintf(fpout, "%sa%d_%d = %s%c(%s);\n",
                 spc, i, j, Ald, pre, GetAadd(flag, i+I0, j, 1));
         fprintf(fpout, "%sa%d_%d = _mm_mul_p%c(a%d_%d, x%d);\n",
                 spc, i, j, pre, i, j, i);
         fprintf(fpout, "%sy%d = _mm_add_p%c(y%d, a%d_%d);\n\n",
                 spc, j, pre, j, i, j);
      }
   }
}

void genMloop(
   FILE *fpout, /* file to print to */
   char *spc,   /* string with indentation spaces */
   int flag,
   char *mstart,/* usually "0", but will be "M8" for UR=8 cleanup */
   char *mend,  /* ending clause for i */
   int nmu,     /* # of reps of mu register block in M unroll */
   int mu,      /* unrolling on M dimension */
   int nu)      /* unrolling on N dimension */
{
   int i, j;
   int MU = nmu * mu;

   fprintf(fpout, "\n%sfor (%s; i < %s; i += %d)\n", spc, mstart, mend, MU, MU);
   fprintf(fpout, "%s{/* ----- BEGIN M-LOOP BODY ----- */\n", spc);
   for (i=0; i < nmu; i++)
   {
      fprintf(fpout, "%s   /* --- BEGIN MUxNU UNROLL %d --- */\n", spc, i);
      gen_MUxNU(fpout, spc-3, flag, mu, nu, i*mu);
      fprintf(fpout, "%s   /* --- END MUxNU UNROLL %d --- */\n", spc, i);
   }
   fprintf(fpout, "%s}/* ----- END M-LOOP BODY ----- */\n", spc);
}

void genNloopBody(
   FILE *fpout, /* file to print to */
   char *spc,   /* string with indentation spaces */
   int flag,
   int nmu,     /* # of reps of mu register block in M unroll */
   int mu,      /* unrolling on M dimension */
   int nu)      /* unrolling on N dimension */
{
   int i, j, MU = nmu*mu;
   char mstart[8], mend[8];
   char pre;

   pre = TEST_FLAG(flag, R1F_SINGLE) ? 's' : 'd';

   sprintf(mstart, "i=MAp");
   sprintf(mend, "M%d", MU);
   fprintf(fpout, "%sif (MAp != 1)\n", spc);
   fprintf(fpout, "%s{/* peel to zero Y */\n", spc);
   fprintf(fpout, "%s   i=0;\n", spc);
   gen_MUxNU0(fpout, spc-3, flag, 2, nu, 0);
   fprintf(fpout, "%s} /* end zero Y peel */\n", spc);
   fprintf(fpout, "%selse /* if (MAp == 1)*/\n", spc);
   fprintf(fpout, "%s{/* peel to force X/A alignment, zero Y */\n", spc);
   fprintf(fpout, "%s   i=0;\n", spc);
   gen_sMUxNU0(fpout, spc-3, flag, 1, nu, 0);
   fprintf(fpout, "%s} /* end force-align/zeroY peel */\n", spc);
   genMloop(fpout, spc, flag, mstart, mend, nmu, mu, nu);
   if (MU > 2)
   {
      sprintf(mstart, "i=M%d", MU);
      fprintf(fpout, "%sif (M != M%d)\n", spc, MU);
      fprintf(fpout,"%s{/* ----- BEGIN VECTOR UNROLL M CLEANUP ----- */\n",spc);
      spc -= 3;
//      if (mu > 6)
         genMloop(fpout, spc, flag, mstart, "M2", 1, 2, nu);
//      else  /* put case statement in here eventually */
//      {
//         gen_MUxNU(fpout, spc-3, flag, mu, nu, i*mu);
//      }
   }
   fprintf(fpout, "%sif (M != M2)\n", spc);
   fprintf(fpout, "%s{/* ----- BEGIN SCALAR M CLEANUP ----- */\n", spc);
   gen_sMUxNU(fpout, spc-3, flag, 1, nu, 0);
   fprintf(fpout, "%s}/* ----- END SCALAR M CLEANUP ----- */\n", spc);

   if (MU > 2)
   {
      spc += 3;
      fprintf(fpout, "%s}/* ----- END VECTOR UNROLL M CLEANUP ----- */\n", spc);
   }
/*
 * Handle post-loop reduction & store to Y
 */
   for (j=0; j < nu/2; j++)
   {
      fprintf(fpout, "%s_my_hadd_pd(y%d, y%d);\n", spc, j*2, j*2+1);
      fprintf(fpout, "%s#if defined(BETAX) || defined(BETA1)\n", spc);
      fprintf(fpout, "%s   a0_%d = _mm_load_pd(Y+j+%d);\n", spc, j, j*2);
      fprintf(fpout,
              "%s   #ifdef BETAX\n%s     a0_%d = _mm_mul_pd(a0_%d, vbeta);\n",
              spc, spc, j, j);
      fprintf(fpout, "%s   #endif\n", spc);
      fprintf(fpout, "%s   y%d = _mm_add_pd(y%d, a0_%d);\n", spc, j*2, j*2, j);
      fprintf(fpout, "%s#endif\n", spc);
      fprintf(fpout, "%s_mm_store_pd(Y+j+%d, y%d);\n", spc, j*2, j*2);
   }
   if (nu%2)
   {
      fprintf(fpout, "%s_my_hadd_pd(y%d, y%d);\n", spc, j*2, j*2);
      fprintf(fpout, "%s#if defined(BETAX) || defined(BETA1)\n", spc);
      fprintf(fpout, "%s   a0_%d = _mm_load_sd(Y+j+%d);\n", spc, j, j*2);
      fprintf(fpout,
              "%s   #ifdef BETAX\n%s      a0_%d = _mm_mul_sd(a0_%d, vbeta);\n",
              spc, spc, j, j);
      fprintf(fpout, "%s   #endif\n", spc);
      fprintf(fpout, "%s   y%d = _mm_add_sd(y%d, a0_%d);\n", spc, j*2, j*2, j);
      fprintf(fpout, "%s#endif\n", spc);
      fprintf(fpout, "%s_mm_store_sd(Y+j+%d, y%d);\n", spc, j*2, j*2);
   }
}

void genNloop(
   FILE *fpout, /* file to print to */
   char *spc,   /* string with indentation spaces */
   int flag,
   char *nstart,/* usually "0", but will be "M8" for UR=8 cleanup */
   char *nend,  /* ending clause for j */
   int nmu,     /* # of reps of mu register block in M unroll */
   int mu,      /* unrolling on M dimension */
   int nu)      /* unrolling on N dimension */
{
   int j;
   if (TEST_FLAG(flag, R1F_APTRS))
   {
      fprintf(fpout, "\n%sfor (j=%s; j < %s; j += %d",
              spc, nstart, nend, nu);
      for (j=0; j < nu; j++)
         fprintf(fpout, ", A%d += lda%d", j, nu);
      fprintf(fpout, ")\n");
   }
   else
      fprintf(fpout, "\n%sfor (j=%s; j < %s; j += %d, A += lda%d)\n",
              spc, nstart, nend, nu, nu);
   fprintf(fpout, "%s{/* BEGIN N-LOOP UR=%d */\n", spc, nu);
   genNloopBody(fpout, spc-3, flag, nmu, mu, nu);
   fprintf(fpout, "%s}/* END N-LOOP UR=%d */\n", spc, nu);
}

#define NUNDEF 3
static char *udefs[NUNDEF] = {"MA", "MAp", "A0"};
void genMV(
   FILE *fpout, /* file to print to */
   char *spc,   /* string with indentation spaces */
   char *rout,  /* routine name */
   int flag,    /* bit field with flags controlling generation options */
   int nmu,     /* # of reps of mu register block in M unroll */
   int mu,      /* unrolling on M dimension */
   int nu)      /* unrolling on N dimension */
{
   int i, j;
   int MU = nmu*mu, incM;
   char nbnd[8];
   incM = TEST_FLAG(flag, R1F_SINGLE) ? 4 : 2;

   fprintf(fpout, "%s#include <xmmintrin.h>\n", spc);
   fprintf(fpout, "%s#include \"atlas_misc.h\"\n", spc);
   fprintf(fpout, "%s#define _my_hadd_pd(dst, src) \\\n", spc);
   fprintf(fpout, "%s   __asm__ __volatile__ (\"haddpd %%2, %%0\" : \"=x\"(dst) : \"0\" (dst), \"x\"(src))\n", spc);
   fprintf(fpout, "\n%svoid %s\n", spc, rout);
   fprintf(fpout, "%s   (ATL_CINT N, ATL_CINT M, const SCALAR alpha,\n", spc);
   fprintf(fpout, "%s    TYPE *A, ATL_CINT lda1, const TYPE *X, ATL_CINT incX,\n",
           spc);
   fprintf(fpout, "%s    const SCALAR beta, TYPE *Y, ATL_CINT incY)\n", spc);
   fprintf(fpout, "%s{/* BEGIN GEMV: nMU=%d, MU=%d, NU=%d */\n",spc, nmu, mu, nu);
   spc -= 3;
   fprintf(fpout, "%sATL_INT i, j;\n", spc);
   if (!TEST_FLAG(flag, R1F_SINGLE))
      fprintf(fpout,"%sATL_CINT MAp = ((((size_t)A)&0xF) || M==1) ? 1 : 2;\n", spc);
   else
      fprintf(fpout, "%sATL_CINT MAp = ((size_t)A)&0xF ? (( (((((size_t)A)+15)>>4)<<4) - ((size_t)A) )/sizeof(TYPE)) : 4;\n", spc);
   fprintf(fpout, "%sATL_CINT MA = M - MAp;\n", spc);
   if (TEST_FLAG(flag, R1F_APTRS))
   {
      fprintf(fpout, "%s#define A0 A\n", spc);
      if (nu)
      {
         fprintf(fpout, "%sconst TYPE *A1=A0+lda1", spc);
         for (j=2; j < nu; j++)
            fprintf(fpout, ", *A%d=A%d+lda1", j, j-1);
         fprintf(fpout, ";\n");
      }
   }
   fprintf(fpout, "%sATL_CINT", spc);
   fprintf(fpout, " M%d=(%s)+MAp", MU, GetMul(MU, GetDiv(MU, "MA")));
   if (MU != 2)
      fprintf(fpout, ", M2=((MA>>1)<<1)+MAp");
   if (nu > 1)
   {
      if ((nu-1) & nu)
         fprintf(fpout, ", N%d=(N/%d)*%d", nu, nu, nu);
      else
         fprintf(fpout, ", N%d=%s", nu, GetMul(nu, GetDiv(nu, "N")));
   }
   for (j=2; j <= nu; j++)
   {
      if (!TEST_FLAG(flag, R1F_APTRS))
         fprintf(fpout, ", lda%d=lda%d+lda1", j, j-1);
   }
   if (TEST_FLAG(flag, R1F_APTRS))
      fprintf(fpout, ", lda%d=%s", nu, GetMul(nu, "lda1"));
   fprintf(fpout, ";\n");
   fprintf(fpout, "%s#ifdef BETAX\n%s   __m128d vbeta;\n%s#endif\n", spc,spc,spc);
   fprintf(fpout, "%s__m128d x0", spc);
   for (j=1; j < mu; j++)
      fprintf(fpout, ", x%d", j);
   for (j=0; j < nu; j++)
      fprintf(fpout, ", y%d", j);
   for (j=0; j < nu; j++)
      for (i=0; i < mu; i += incM)
         fprintf(fpout, ", a%d_%d", i, j);
   fprintf(fpout, ";\n");

   fprintf(fpout, "%sif (!M || !N) return;\n", spc);
   fprintf(fpout, "%s#ifdef BETAX\n%s   vbeta = _mm_load1_pd(&beta);\n%s#endif\n",
           spc, spc, spc);
   if (nu > 1)
      sprintf(nbnd, "N%d", nu);
   else
      sprintf(nbnd, "N");
   genNloop(fpout, spc, flag, "0", nbnd, nmu, mu, nu);
   if (nu > 1)  /* later on, use case statement rather than NU=1 cleanup */
   {
      genNloop(fpout, spc, flag, nbnd, "N", nmu, mu, 1);
   }
   spc += 3;
   fprintf(fpout, "%s}/* END GEMV: nMU=%d, MU=%d, NU=%d */\n",
           spc, nmu, mu, nu);
   for (i=0; i < NUNDEF; i++)
      fprintf(fpout, "%s#ifdef %s\n%s   #undef %s\n%s#endif\n",
              spc, udefs[i], spc, udefs[i], spc);
}

void PrintUsage(char *name, char *arg, int i)
{
   int j;
   if (i > 0)
      fprintf(stderr, "BAD ARG '%s' on %dth FLAG\n", arg, i);
   fprintf(stderr, "USAGE: %s [flags], where flags are:\n", name);
   fprintf(stderr, "   -M <#> : repeat mu unroll # times in loop\n");
   fprintf(stderr, "   -m <#> : unroll (wt reg blking) M loop by #\n");
   fprintf(stderr, "   -n <#> : unroll&jam N loop by #\n");
   fprintf(stderr, "   -S <flag#> [0/1] : unset/set flag.  Flags include:\n");
   for (j=0; j < NUMFLAGS; j++)
      fprintf(stderr, "      %2d : %s\n", j, R1F_exp[j]);
   fprintf(stderr, "   -F <#> : set flag bitfield to #\n");
   exit(i ? i : -1);
}

int GetFlags(int nargs, char **args, int *NMU, int *MU, int *NU)
{
   int i, flag=0, j, set;

   *NMU = *NU = 1;
   *MU = 2;

   flag = 0;
   SET_FLAG(flag, R1F_ALGLDA, 1);
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], "No '-' preceeding flag!", i);
      switch(args[i][1])
      {
      case 'S':
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -S ", i-1);
         set = atoi(args[i]);
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -S ", i-1);
         j = atoi(args[i]);
         SET_FLAG(flag, set, j);
        break;
      case 'F':
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -S ", i-1);
         flag = atoi(args[i]);
         break;
      case 'M':
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -M ", i-1);
         *NMU = atoi(args[i]);
        break;
      case 'm':
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -m ", i-1);
         *MU = atoi(args[i]);
        break;
      case 'n':
        if (++i >= nargs)
            PrintUsage(args[0], "out of flags in -n ", i-1);
         *NU = atoi(args[i]);
        break;
      default:
         PrintUsage(args[0], args[i], i);
      }
   }
   return(flag);
}

#define NSPCS 128
int main(int nargs, char **args)
{
   int i, nmu, mu, nu, flag;
   char spc[NSPCS];
   char pre;
   char name[1024];

   for (i=0; i < NSPCS; i++)
      spc[i] = ' ';
   spc[NSPCS-1] = '\0';
   flag = GetFlags(nargs, args, &nmu, &mu, &nu);
   pre = TEST_FLAG(flag, R1F_SINGLE) ? 's' : 'd';
   sprintf(name, "Mjoin(Mjoin(ATL_%cgemvT_a1_x1,BNM),_y1)", pre);
   genMV(stdout, spc+NSPCS-1, name, flag, nmu, mu, nu);
   exit(0);
}
