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
#include <string.h>
#include <assert.h>
#include "atlas_misc.h"
#include "atlas_r1parse.h"

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

static int Mylcm(const int M, const int N)
/*
 * Returns least common multiple (LCM) of two positive integers M & N by
 * computing greatest common divisor (GCD) and using the property that
 * M*N = GCD*LCM.
 */
{
   register int tmp, max, min, gcd=0;

   if (M != N)
   {
      if (M > N) { max = M; min = N; }
      else { max = N; min = M; }
      if (min > 0)  /* undefined for negative numbers */
      {
         do  /* while (min) */
         {
            if ( !(min & 1) ) /* min is even */
            {
               if ( !(max & 1) ) /* max is also even */
               {
                  do
                  {
                     min >>= 1;
                     max >>= 1;
                     gcd++;
                     if (min & 1) goto MinIsOdd;
                  }
                  while ( !(max & 1) );
               }
               do min >>=1 ; while ( !(min & 1) );
            }
/*
 *          Once min is odd, halve max until it too is odd.  Then, use
 *          property that gcd(max, min) = gcd(max, (max-min)/2)
 *          for odd max & min
 */
MinIsOdd:
            if (min != 1)
            {
               do  /* while (max >= min */
               {
                  max -= (max & 1) ? min : 0;
                  max >>= 1;
               }
               while (max >= min);
            }
            else return( (M*N) / (1<<gcd) );
            tmp = max;
            max = min;
            min = tmp;
         }
         while(tmp);
      }
      return( (M*N) / (max<<gcd) );
   }
   else return(M);
}

/*
 * For SYR and SYR2, generate a macro which does * a small NUxNU
 * triangular matrix so that GER kernel can be called
 * on rest of NU-wide panel.
 */
void UnrollSYR2
(
   FILE *fpout,         /* stream to print to */
   char *name,          /* name for macro */
   char pre,            /* precisition/type prefix */
   enum ATLAS_UPLO Uplo,
   int  nu              /* unroll factor */
)
/*
 * Real precision unroll of SYR2
 */
{
   int i, j;

   fprintf(fpout, "#define %s(A_, lda_, x_, y_) \\\n{ \\\n", name);
   fprintf(fpout, "   TYPE *aa=(A_); \\\n");
   fprintf(fpout, "   ATL_CINT lda0_ = 0");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", lda%d_ = lda%d_+(lda_)", i, i-1);
   fprintf(fpout, "; \\\n   const TYPE x0_=*(x_)");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", x%d_=(x_)[%d]", i, i);
   fprintf(fpout, "; \\\n   const TYPE y0_=*(y_)");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", y%d_=(y_)[%d]", i, i);
   fprintf(fpout, "; \\\n");
   if (Uplo == AtlasUpper)
   {
      for (j=0; j < nu; j++)
         for (i=0; i <= j; i++)
            fprintf(fpout, "   aa[lda%d_+%d] += x%d_*y%d_ + y%d_*x%d_; \\\n",
                    j, i, i, j, i, j);
   }
   else
   {
      for (j=0; j < nu; j++)
         for (i=j; i < nu; i++)
            fprintf(fpout, "   aa[lda%d_+%d] += x%d_*y%d_ + y%d_*x%d_; \\\n",
                    j, i, i, j, i, j);
   }
   fprintf(fpout, "}\n");
}

void UnrollHER2
(
   FILE *fpout,         /* stream to print to */
   char *name,          /* name for macro */
   char pre,            /* precisition/type prefix */
   enum ATLAS_UPLO Uplo,
   int  nu              /* unroll factor */
)
/*
 * Complex type unroll of HER2
 */
{
   int i, j;

   fprintf(fpout, "#define %s(A_, lda_, x_, y_, xt_, yt_) \\\n{ \\\n", name);
   fprintf(fpout, "   TYPE *aa=(A_); \\\n");
   fprintf(fpout, "   ATL_CINT lda0_ = 0");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", lda%d_ = lda%d_+(lda_)+(lda_)", i, i-1);
   fprintf(fpout, "; \\\n   const TYPE x0r=*(x_), x0i=(x_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", x%dr=(x_)[%d], x%di=(x_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", xt%dr=(xt_)[%d], xt%di=(xt_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n   const TYPE y0r=*(y_), y0i=(y_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", y%dr=(y_)[%d], y%di=(y_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n   const TYPE yt0r=*(yt_), yt0i=(yt_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", yt%dr=(yt_)[%d], yt%di=(yt_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n");
   if (Uplo == AtlasUpper)
   {
      for (j=0; j < nu; j++)
      {
         for (i=0; i < j; i++)
         {
            fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%dr-x%di*yt%di + y%dr*xt%dr-y%di*xt%di; \\\n",
                    j, 2*i, i, j, i, j, i, j, i, j);
            fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%di+x%di*yt%dr + y%dr*xt%di+y%di*xt%dr; \\\n",
                    j, 2*i+1, i, j, i, j, i, j, i, j);
         }
         fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%dr-x%di*yt%di + y%dr*xt%dr-y%di*xt%di; \\\n",
                 j, 2*j, j, j, j, j, j, j, j, j);
         fprintf(fpout, "   aa[lda%d_+%d] = 0.0; \\\n", j, 2*j+1);
      }
   }
   else
   {
      for (j=0; j < nu; j++)
      {
         fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%dr-x%di*yt%di + y%dr*xt%dr-y%di*xt%di; \\\n",
                 j, 2*j, j, j, j, j, j, j, j, j);
         fprintf(fpout, "   aa[lda%d_+%d] = 0.0; \\\n", j, 2*j+1);
         for (i=j+1; i < nu; i++)
         {
            fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%dr-x%di*yt%di + y%dr*xt%dr-y%di*xt%di; \\\n",
                    j, 2*i, i, j, i, j, i, j, i, j);
            fprintf(fpout,
      "   aa[lda%d_+%d] += x%dr*yt%di+x%di*yt%dr + y%dr*xt%di+y%di*xt%dr; \\\n",
                    j, 2*i+1, i, j, i, j, i, j, i, j);
         }
      }
   }
   fprintf(fpout, "}\n");
}
void UnrollSYR1
(
   FILE *fpout,         /* stream to print to */
   char *name,          /* name for macro */
   char pre,            /* precisition/type prefix */
   enum ATLAS_UPLO Uplo,
   int  nu              /* unroll factor */
)
/*
 * Real precision unroll of SYR1
 */
{
   int i, j;

   fprintf(fpout, "#define %s(A_, lda_, x_, y_) \\\n{ \\\n", name);
   fprintf(fpout, "   TYPE *aa=(A_); \\\n");
   fprintf(fpout, "   ATL_CINT lda0_ = 0");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", lda%d_ = lda%d_+(lda_)", i, i-1);
   fprintf(fpout, "; \\\n   const TYPE x0_=*(x_)");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", x%d_=(x_)[%d]", i, i);
   fprintf(fpout, "; \\\n   const TYPE y0_=*(y_)");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", y%d_=(y_)[%d]", i, i);
   fprintf(fpout, "; \\\n");
   if (Uplo == AtlasUpper)
   {
      for (j=0; j < nu; j++)
         for (i=0; i <= j; i++)
            fprintf(fpout, "   aa[lda%d_+%d] += x%d_*y%d_; \\\n",
                    j, i, i, j);
   }
   else
   {
      for (j=0; j < nu; j++)
         for (i=j; i < nu; i++)
            fprintf(fpout, "   aa[lda%d_+%d] += x%d_*y%d_; \\\n",
                    j, i, i, j);
   }
   fprintf(fpout, "}\n");
}

void UnrollHER1
(
   FILE *fpout,         /* stream to print to */
   char *name,          /* name for macro */
   char pre,            /* precisition/type prefix */
   enum ATLAS_UPLO Uplo,
   int  nu              /* unroll factor */
)
/*
 * Complex type unroll of HER1
 */
{
   int i, j;

   fprintf(fpout, "#define %s(A_, lda_, x_, xt_) \\\n{ \\\n", name);
   fprintf(fpout, "   TYPE *aa=(A_); \\\n");
   fprintf(fpout, "   ATL_CINT lda0_ = 0");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", lda%d_ = lda%d_+(lda_)+(lda_)", i, i-1);
   fprintf(fpout, "; \\\n   const TYPE x0r=*(x_), x0i=(x_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", x%dr=(x_)[%d], x%di=(x_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]");
   for (i=1; i < nu; i++)
      fprintf(fpout, ", xt%dr=(xt_)[%d], xt%di=(xt_)[%d]", i, 2*i, i, 2*i+1);
   fprintf(fpout, "; \\\n");
   if (Uplo == AtlasUpper)
   {
      for (j=0; j < nu; j++)
      {
         for (i=0; i < j; i++)
         {
            fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%dr-x%di*xt%di; \\\n",
                    j, 2*i, i, j, i, j);
            fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%di+x%di*xt%dr; \\\n",
                    j, 2*i+1, i, j, i, j, i, j, i, j);
         }
         fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%dr-x%di*xt%di; \\\n",
                 j, 2*j, j, j, j, j);

         fprintf(fpout, "   aa[lda%d_+%d] = 0.0; \\\n", j, 2*j+1);
      }
   }
   else
   {
      for (j=0; j < nu; j++)
      {
         fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%dr-x%di*xt%di; \\\n",
                 j, 2*j, j, j, j, j);
         fprintf(fpout, "   aa[lda%d_+%d] = 0.0; \\\n", j, 2*j+1);
         for (i=j+1; i < nu; i++)
         {
            fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%dr-x%di*xt%di; \\\n",
                    j, 2*i, i, j, i, j);
            fprintf(fpout, "   aa[lda%d_+%d] += x%dr*xt%di+x%di*xt%dr; \\\n",
                    j, 2*i+1, i, j, i, j);
         }
      }
   }
   fprintf(fpout, "}\n");
}

int FixMB(char pre, int mu, int mb)
/*
 * Makes sure MB is not a power of two, and that it won't mess up allignment
 * and that it is a multiple of mu
 */
{
   int MU;

/*
 * Find MU necessary to keep 16-byte alignment & take LCM wt kernel's mu
 */
   if (pre == 's')
      MU = 4;
   else if (pre == 'c' || pre == 'd')
      MU = 2;
   else
      MU = 1;
   MU = Mylcm(mu, MU);
/*
 * Force mb to multiple of mu
 */
   mb = (2*MU < mb) ? (mb/MU)*MU : MU;
/*
 * Don't use power of 2 (avoid worst-case column cache set conflict)
 */
   if (!(mb & (mb-1)))
      mb = (mb > 2*MU) ? mb-MU : mb;
   return(mb);
}

void s1hgen
(
   ATL_r1node_t *r1B,   /* standard 8-entry R1SUMM kernel list */
   int LVL,             /* 0:out-of-cache, 1: in-L1, 2: in-L2 */
   int L1Elts,          /* number of elements in L1 cache */
   char pre,
   char *path           /* path to generate header files in */
)
{
   ATL_r1node_t *r1r, *r1g;     /* restricted and general kernels to use */
   FILE *fpout;
   int NU, imf, i;
   char ln[1024];
   char PRE = toupper(pre);

/*
 * Select kernels based on LVL
 */
   if (!LVL)
   {
      imf = 0;
      r1r = r1B;
   }
   else if (LVL == 2)
   {
      imf = 3;
      r1r = r1B->next->next;
   }
   else if (LVL == 1)
   {
      imf = 4;
      r1r = r1B->next->next->next->next;
   }
   r1g = r1r->next;
   if (r1r->ID == r1g->ID)
      r1r = NULL;
/*
 * Don't support restricted kernel unless it is much faster than general,
 * if it requires changing the blocking
 */
   else if (r1r->XU != r1g->XU && r1r->mflop[imf] < 1.05*r1g->mflop[imf])
     r1r = NULL;

   NU = r1r ? Mylcm(r1r->YU, r1g->YU) : r1g->YU;

   assert(strlen(path) < 1010);         /* the shame, the shame */
   if (LVL)
      sprintf(ln, "%s/atlas_%csyr_L%d.h", path, pre, LVL);
   else
      sprintf(ln, "%s/atlas_%csyr.h", path, pre);
   fpout = fopen(ln, "w");
   fprintf(fpout, "/*\n * This file generated on line %d of %s\n */\n",
           __LINE__, __FILE__);
   if (LVL)
   {
      fprintf(fpout,
              "#ifndef ATLAS_%cSYR_L%d_H\n   #define ATLAS_%cSYR_L%d_H\n\n",
              PRE, LVL, PRE, LVL);
      fprintf(fpout, "#include \"atlas_%cr1_L%d.h\"\n", pre, LVL);
   }
   else
   {
      fprintf(fpout,
              "#ifndef ATLAS_%cSYR_H\n   #define ATLAS_%cSYR_H\n\n",
              PRE, PRE);
      fprintf(fpout, "#include \"atlas_%cr1.h\"\n", pre);
   }

//   NU = (NU < 4) ? NU+NU : NU;
   fprintf(fpout, "\n#define ATL_s1NU %d\n", NU);
/*
 * Only out-of-cache needs blocking stuff; others always unblocked
 */
   if (!LVL)
   {
      fprintf(fpout, "\n");
      if (!r1g->CacheElts)
         fprintf(fpout, "#define ATL_NOBLOCK_S1 1\n");
      if (NU == r1g->YU)
         fprintf(fpout, "#define ATL_GetPartS1 ATL_GetPartR1\n");
      else
      {
         i = r1g->CacheElts;
         i = (i) ? (i-2*NU)/(2*NU+1) : (1<<29);
         i = FixMB(pre, r1g->XU, i);
         fprintf(fpout, "#define ATL_GetPartS1(A_, lda_, mb_, nb_) { (mb_) = %d; (nb_) = %d; }\n", i, NU);
      }
      if (r1r)
      {
         fprintf(fpout, "ATL_s1USERRESTRICTK 1\n");
         if (!r1r->CacheElts)
            fprintf(fpout, "#define ATL_NOBLOCK_S1r 1\n");
         if (NU == r1r->YU)
            fprintf(fpout, "#define ATL_GetPartS1r ATL_GetPartR1r\n");
         else
         {
            i = r1r->CacheElts;
            i = (i) ? (i-2*NU)/(2*NU+1) : (1<<29);
            i = FixMB(pre, r1g->XU, i);
            fprintf(fpout, "#define ATL_GetPartS1r(A_, lda_, mb_, nb_) { (mb_) = %d; (nb_) = %d; }\n", i, NU);
         }
      }
   }

   fprintf(fpout, "\n");
   if (pre == 's' || pre == 'd')
   {
      UnrollSYR1(fpout, "ATL_SYR1U_nu", pre, AtlasUpper, NU);
      UnrollSYR1(fpout, "ATL_SYR1L_nu", pre, AtlasLower, NU);
   }
   else
   {
      UnrollHER1(fpout, "ATL_HER1U_nu", pre, AtlasUpper, NU);
      UnrollHER1(fpout, "ATL_HER1L_nu", pre, AtlasLower, NU);
   }

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

void s2hgen
(
   ATL_r1node_t *r1B,   /* standard 8-entry R1SUMM kernel list */
   int L1Elts,          /* number of elements in L1 cache */
   char pre,
   char *path           /* path to generate header files in */
)
/*
 * Generates atlas_<pre>syr2.h file, given r1B
 */
{
   char ln[1024];
   FILE *fpout;
   ATL_r1node_t *r1p, *r1IC, *r1OC, *r1ICr, *r1OCr;
   int CacheElts, MU, NU, i;
   char PRE = toupper(pre);
/*
 * Fill in standard names for these files
 */
   PutKernNameInStr(r1B);
/*
 * Default to blocking for L1
 */
   r1OCr = r1B->next->next->next->next->next->next;
   r1OC = r1OCr->next;
   r1ICr = r1B->next->next->next->next;
   r1IC = r1ICr->next;
   CacheElts = r1OC->CacheElts;
/*
 * If best out-of-cache (OC) kernels use no blocking or L2 blocking, we
 * need to check if L2-blocking will be faster for SYR2
 */
   r1p = r1B->next;
   if (r1p->CacheElts == 0 || r1p->CacheElts > L1Elts)
   {
      double mfL1, mfL2;
      mfL1 = (r1OCr->mflop[0] + r1IC->mflop[4])/2.0;
      mfL2 = (r1p->CacheElts) ? r1p->mflop[1] : r1p->mflop[2];
      mfL2 = (mfL2+r1B->next->next->next->mflop[3])/2.0;
      if (mfL2 >= 1.02*mfL1)
      {
         r1OCr = r1B;
         r1OC  = r1OCr->next;
         r1ICr = r1B->next->next;
         r1IC  = r1ICr->next;
         CacheElts = (r1OC->CacheElts) ? r1OC->CacheElts : r1IC->CacheElts;
      }
   }
   if (r1ICr->ID == r1IC->ID)
      r1ICr = NULL;
   if (r1OCr->ID == r1OC->ID)
      r1OCr = NULL;
   MU = Mylcm(r1OC->XU, r1IC->XU);
   if (r1ICr)
      MU = Mylcm(MU, r1ICr->XU);
   if (r1OCr)
      MU = Mylcm(MU, r1OCr->XU);
   NU = Mylcm(r1OC->YU, r1IC->YU);
   if (r1ICr)
      NU = Mylcm(NU, r1ICr->YU);
   if (r1OCr)
      NU = Mylcm(NU, r1OCr->YU);


   sprintf(ln, "%s/atlas_%csyr2.h", path, pre);
   fpout = fopen(ln, "w");
   assert(fpout);
   fprintf(fpout, "/*\n * This file generated on line %d of %s\n */\n",
           __LINE__, __FILE__);
   fprintf(fpout,
           "#ifndef ATLAS_%cSYR2_H\n   #define ATLAS_%cSYR2_H\n\n",
           PRE, PRE);

   fprintf(fpout, "#include \"atlas_%cr1kernels.h\"\n", pre);

   fprintf(fpout, "#define ATL_s2CacheElts %d\n", CacheElts);
   fprintf(fpout, "#define ATL_s2MU %d\n", MU);
   fprintf(fpout, "#define ATL_s2NU %d\n", NU);
   fprintf(fpout, "#define ATL_R1OC %s\n", r1OC->str);
   fprintf(fpout, "#define ATL_R1OCr %s\n", r1OCr ? r1OCr->str : r1OC->str);
   fprintf(fpout, "#define ATL_R1IC %s\n", r1IC->str);
   fprintf(fpout, "#define ATL_R1ICr %s\n", r1ICr ? r1ICr->str : r1IC->str);
   i = FLAG_IS_SET(r1OC->flag, R1F_ALIGNX2A);
   i |= FLAG_IS_SET(r1IC->flag, R1F_ALIGNX2A);
   if (r1OCr)
      i |= FLAG_IS_SET(r1OCr->flag, R1F_ALIGNX2A);
   if (r1ICr)
      i |= FLAG_IS_SET(r1ICr->flag, R1F_ALIGNX2A);
   if (i)
      fprintf(fpout, "#define ATL_s2ALIGNX2A 1\n");

   if (r1OCr)
   {
      fprintf(fpout, "#define ATL_s2USERESTRICTK_OC 1\n");
      fprintf(fpout, "#define ATL_s2UseRestrictK_OC(M_, N_, A_, lda_) \\\n");
      assert(r1OCr->ldamul > 1);  /* only allowed restriction right now! */
      if (r1OCr->ldamul > 1)
         fprintf(fpout, "   (%s == ATL_sizeof*(lda_))",
                 GetMul(r1OCr->ldamul, GetDiv(r1OCr->ldamul,
                                              "ATL_sizeof*(lda_)")));
   }
   if (r1ICr)
   {
      fprintf(fpout, "#define ATL_s2USERESTRICTK_IC 1\n");
      fprintf(fpout, "#define ATL_s2UseRestrictK_IC(M_, N_, A_, lda_) \\\n");
      assert(r1ICr->ldamul > 1);  /* only allowed restriction right now! */
      if (r1ICr->ldamul > 1)
         fprintf(fpout, "   (%s == ATL_sizeof*(lda_))",
                 GetMul(r1ICr->ldamul, GetDiv(r1ICr->ldamul,
                                              "ATL_sizeof*(lda_)")));
   }
   i = (CacheElts-4*NU)/(2*NU+2);
   i = FixMB(pre, MU, i);
   fprintf(fpout, "\n#define ATL_GetPartS2(A_, lda_, mb_, nb_) { (mb_) = %d; (nb_) = ATL_s2NU; }\n\n", i);
   if (pre == 's' || pre == 'd')
   {
      UnrollSYR2(fpout, "ATL_SYR2U_nu", pre, AtlasUpper, NU);
      UnrollSYR2(fpout, "ATL_SYR2L_nu", pre, AtlasLower, NU);
   }
   else
   {
      UnrollHER2(fpout, "ATL_HER2U_nu", pre, AtlasUpper, NU);
      UnrollHER2(fpout, "ATL_HER2L_nu", pre, AtlasLower, NU);
   }

   fprintf(fpout, "\n#endif\n");
   fclose(fpout);
}

void PrintPrototype(FILE *fpout, char pre, char *rout, char *type, char *styp)
{
   fprintf(fpout, "void %s(ATL_CINT, ATL_CINT, const %s, const %s*, ATL_CINT, const %s*, ATL_CINT, %s*, ATL_CINT);\n",
           rout, styp, type, type, type);
}

void GenKernFiles(char pre, char *path, ATL_r1node_t *r1b)
/*
 * r1b is a list of rank-1 kernels that must be compiled (including those
 * needed to form SYR and SYR2).  This list should be unique (same kernel
 * not compiled twice).  r1b->str will have the routine name to give the
 * kernel during compilation.
 */
{
   ATL_r1node_t *r1p;
   char ln[2048];
   for (r1p = r1b; r1p; r1p = r1p->next)
   {
      if (r1p->genstr)   /* generate kernel if necessary */
      {
         assert(!system(r1p->genstr));
         sprintf(ln, "cp %s %s/%s\n", r1p->rout, path, r1p->str);
      }
      else
         sprintf(ln, "cp CASES/%s %s/%s.c\n", r1p->rout, path, r1p->str);
      if (system(ln))
      {
         fprintf(stderr, "FAILED: %s\n", ln);
         exit(-1);
      }
   }
}

void EmitMakefile(char pre, char *path, ATL_r1node_t *r1b)
/*
 * r1b is a list of rank-1 kernels that must be compiled (including those
 * needed to form SYR and SYR2).  This list should be unique (same kernel
 * not compiled twice).  r1b->str will have the routine name to give the
 * kernel during compilation.
 */
{
   ATL_r1node_t *r1p, *r1k;
   char *kern, *outf, *typD;
   FILE *fpout;
   int i, ialias=0;
   const char UPRE = (pre == 'z' || pre == 'd') ? 'D' : 'S';
   static char *aliased[16];

   assert(path);
   kern = (pre == 'z' || pre == 'c') ? "geru" : "ger";
   if (pre == 'd')
      typD = "DREAL";
   else if (pre == 's')
      typD = "SREAL";
   else if (pre == 'c')
      typD = "SCPLX";
   else if (pre == 'z')
      typD = "DCPLX";
   else
      assert(0);

   i = strlen(path);
   outf = malloc((i+10)*sizeof(char));
   assert(outf);
   strcpy(outf, path);
   strcpy(outf+i, "/Make_");   /* Make_<pre>r1 */
   outf[i+6] = pre;
   outf[i+7] = 'r';
   outf[i+8] = '1';
   outf[i+9] = '\0';
   fpout = fopen(outf, "w");
   assert(fpout);

   fprintf(fpout, "#\n#  This file generated at line %d of %s\n#\n",
           __LINE__, __FILE__);
   fprintf(fpout, "include Make.inc\n\n");
   fprintf(fpout, "R1CC = $(%cKC)\nR1FLAGS = $(CDEFS) $(%cKCFLAGS)",
            UPRE, UPRE);
   fprintf(fpout, " -D%s\n\n", typD);
   fprintf(fpout, "obj =");
   for (r1p=r1b; r1p; r1p = r1p->next)
      fprintf(fpout, " %s.o", r1p->str);
   fprintf(fpout, "\n");

   fprintf(fpout, "lib : %clib\n%clib : %cr1k.grd\n", pre, pre, pre);
   fprintf(fpout, "%cr1k.grd : $(obj)\n", pre);
   fprintf(fpout, "\t$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(obj)\n");
   fprintf(fpout, "\t$(RANLIB) $(ATLASlib)\n");
   fprintf(fpout, "\ttouch %cr1k.grd\n", pre);

   fprintf(fpout, "%cclean : clean\n", pre);
   fprintf(fpout, "clean :\n\trm -f $(obj) %cr1k.grd\n\n", pre);

   fprintf(fpout, "%ckilllib : killlib\n", pre);
   fprintf(fpout, "killlib : \n");
   fprintf(fpout, "\t$(ARCHIVER) d $(ATLASlib) $(obj)\n");
   fprintf(fpout, "\t$(RANLIB) $(ATLASlib)\n");
   fprintf(fpout, "killall : killlib clean\n");
   fprintf(fpout, "\t rm -f");
   for (r1p=r1b; r1p; r1p = r1p->next)
      fprintf(fpout, " %s.c", r1p->str);
   fprintf(fpout, "\n\n");
/*
 * Spit out build command for all surviving kernels
 */
   for (r1p=r1b; r1p; r1p = r1p->next)
   {
      fprintf(fpout, "%s.o : %s.c\n", r1p->str, r1p->str);
      if (r1p->comp)
         fprintf(fpout, "\t %s", r1p->comp);
      else
         fprintf(fpout, "\t $(R1CC)");
      fprintf(fpout, " -o %s.o -c -DATL_UGERK=%s", r1p->str, r1p->str);
      if (r1p->cflags)
         fprintf(fpout, " %s -D%s", r1p->cflags, typD);
      else
         fprintf(fpout, " $(R1FLAGS)");
      fprintf(fpout, " %s.c\n", r1p->str);
   }
   free(outf);
}

void r1khgen(char pre, char *path, ATL_r1node_t *r1b, char **aliases)
{
   char *ln;
   int i;
   FILE *fpout;
   ATL_r1node_t *r1p;
   char *styp, *type = (pre == 'd' || pre == 'z') ? "double" : "float";
   char PRE;

   PRE = toupper(pre);
   if (pre == 'd')
      styp = "double";
   else if (pre == 's')
      styp = "float";
   else if (pre == 'c')
      styp = "float*";
   else
      styp = "double*";

   i = strlen(path);
   ln = malloc(i+32*sizeof(char));
   sprintf(ln, "%s/atlas_%cr1kernels.h", path, pre);

   fpout = fopen(ln, "w");
   fprintf(fpout, "/*\n * This file generated on line %d of %s\n */\n",
           __LINE__, __FILE__);
   fprintf(fpout,
           "#ifndef ATLAS_%cR1KERNELS_H\n   #define ATLAS_%cR1KERNELS_H\n\n",
           PRE, PRE);

   for (r1p=r1b; r1p; r1p = r1p->next)
      PrintPrototype(fpout, pre, r1p->str, type, styp);
   if (aliases)
   {
      fprintf(fpout, "\n");
      for (i=0; aliases[i]; i += 2)
         fprintf(fpout, "#define %-24s %s\n", aliases[i], aliases[i+1]);

   }
   fprintf(fpout, "\n#endif /* end guard around atlas_%cr1kernels.h */\n", pre);

   fclose(fpout);
   free(ln);
}

void r1hgen(char pre, char *path, int LVL, ATL_r1node_t *kp, ATL_r1node_t *kpR)
{
   FILE *fpout;
   char PRE = toupper(pre);
   char *type = (pre == 'd' || pre == 'z') ? "double" : "float";
   char *styp, *sp;
   char gerk[32];
   int mb, DOPROTO=0;
   int irest=12, iconj=8;

   if (kpR && kpR->ID == kp->ID)   /* nuke restricted if same as normal */
      kpR = NULL;
   if (LVL < 0)
   {
      DOPROTO = 1;
      LVL = 0;
   }
/*
 * Name the kernel according to cache block level and data type
 */
   assert(LVL >= 0 && LVL <= 9);
   sprintf(gerk, "ATL_%cgerk_L%d_restrict", pre, LVL);
   gerk[irest] = '\0';
   sp = malloc(strlen(path) + 16);
   assert(sp);
   if (!LVL)
      sprintf(sp, "%s/atlas_%cr1.h", path, pre);
   else
      sprintf(sp, "%s/atlas_%cr1_L%d.h", path, pre, LVL);
   fpout = fopen(sp, "w");
   free(sp);
   if (pre == 'd')
      styp = "double";
   else if (pre == 's')
      styp = "float";
   else if (pre == 'c')
      styp = "float*";
   else
      styp = "double*";
   fprintf(fpout, "#ifndef ATLAS_%cR1_L%d_H\n#define ATLAS_%cR1_L%d_H\n\n",
           PRE, LVL, PRE, LVL);
   fprintf(fpout, "#include \"atlas_type.h\"\n\n");
   fprintf(fpout, "#define ATL_r1CacheElts %d\n", kp->CacheElts);
   fprintf(fpout, "#define ATL_r1MU %d\n", kp->XU);
   fprintf(fpout, "#define ATL_r1NU %d\n", kp->YU);
   if (FLAG_IS_SET(kp->flag, R1F_ALIGNX2A))
      fprintf(fpout, "#define ATL_r1ALIGNX2A 1\n");
   if (FLAG_IS_SET(kp->flag, R1F_FYU))
      fprintf(fpout, "#define ATL_r1NMUL %d\n", kp->YU);
   if (kp->CacheElts == 0)
      fprintf(fpout, "#define ATL_r1NOBLOCK\n");
   if (kp->CacheElts > 0)
   {
      mb = (kp->CacheElts - 2*kp->YU) / (2*kp->YU+1);   /* elts in cache */
      mb = FixMB(pre, kp->XU, mb);
   }
   else
      mb = 0;
   if (DOPROTO)
   {
      fprintf(fpout, "#define %s ATL_UGERK\n", gerk);
      PrintPrototype(fpout, pre, gerk, type, styp);
   }
   else
      fprintf(fpout, "#include \"atlas_%cr1kernels.h\"\n", pre);
/*
 * Print handles for kernel names for use by fixed code
 */
   fprintf(fpout, "#define ATL_GERK ATL_%cgerk_L%d\n", pre, LVL);
   fprintf(fpout, "#define ATL_GERKr ATL_%cgerk_L%d_restrict\n\n", pre, LVL);

   fprintf(fpout, "#define ATL_GetPartR1(A_, lda_, mb_, nb_) { (mb_) = %d; (nb_) = ATL_r1NU; }\n", mb);
   if (kpR)
   {
      fprintf(fpout, "#define ATL_r1USERESTRICTK 1\n");
      fprintf(fpout, "#define ATL_r1CacheEltsr %d\n", kpR->CacheElts);
      fprintf(fpout, "#define ATL_r1MUr %d\n", kpR->XU);
      fprintf(fpout, "#define ATL_r1NUr %d\n", kpR->YU);
      if (FLAG_IS_SET(kpR->flag, R1F_ALIGNX2A))
         fprintf(fpout, "#define ATL_r1ALIGNX2Ar 1\n");
      if (FLAG_IS_SET(kpR->flag, R1F_FYU))
         fprintf(fpout, "#define ATL_r1NMULr %d\n", kpR->YU);
      if (kpR->CacheElts == 0)
         fprintf(fpout, "#define ATL_r1NOBLOCKr\n");
      gerk[irest] = '_';
      if (DOPROTO)
         PrintPrototype(fpout, pre, gerk, type, styp);
      mb = (kpR->CacheElts - 2*kpR->YU) / (2*kpR->YU+1);
      mb = (mb >= kpR->XU) ? (mb/kpR->XU)*kpR->XU : kpR->XU;
      fprintf(fpout, "#define ATL_GetPartR1r(A_, lda_, mb_, nb_) { (mb_) = %d; (nb_) = ATL_r1NUr; }\n");
      fprintf(fpout, "#define ATL_r1UseRestrictK(M_, N_, A_, lda_) \\\n");
      assert(kpR->ldamul > 1);  /* only allowed restriction right now! */
      if (kpR->ldamul > 1)
         fprintf(fpout, "   (%s == ATL_sizeof*(lda_))",
                 GetMul(kpR->ldamul, GetDiv(kpR->ldamul, "ATL_sizeof*(lda_)")));
   }
   fprintf(fpout,
           "\n#endif  /* end protection around header file contents */\n");
   fclose(fpout);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set type/precision prefix \n");
   fprintf(stderr, "   -d <dir> : output files using path <dir>\n");
   fprintf(stderr, "   -F <file> : read kernel file & gen headers\n");
   fprintf(stderr, "    The following flags can be used if -F is not:\n");
   fprintf(stderr, "      -l <l1mul> : use l1mul*L1CacheSize for blocking\n");
   fprintf(stderr, "      -x <xu> : xu elements of X are accessed at once\n");
   fprintf(stderr, "      -y <yu> : yu elements of Y are accessed at once\n");
   fprintf(stderr, "      -f <iflag> : set the flag bitfield to iflag\n");
   exit(ierr ? ierr : -1);
}

void GetFlags(int nargs, char **args, char *PRE, char **FNAM, char **DIR,
              int *XU, int *YU, int *L1MUL, int *IFLAG)
{
   int i, k;
   char pre='d';

   *DIR = "./";
   *IFLAG = *XU = *YU = *L1MUL = 0;
   *FNAM = NULL;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'x':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *XU = atoi(args[i]);
         break;
      case 'y':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *YU = atoi(args[i]);
         break;
      case 'l':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *L1MUL = atoi(args[i]);
         break;
      case 'f':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *IFLAG = atoi(args[i]);
         break;
      case 'd':
         if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
         *DIR = args[i];
         break;
      case 'p':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        pre = tolower(args[i][0]);
        assert(pre == 's' || pre == 'd' || pre == 'z' || pre == 'c');
        break;
      case 'F':
        if (++i >= nargs)
            PrintUsage(args[0], i-1, NULL);
        *FNAM = args[i];
        break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
   *PRE = pre;
   if (*FNAM == NULL && (*XU == 0 || *YU == 0))
   {
      *FNAM = malloc(16*sizeof(char));
      sprintf(*FNAM, "res/%cR1SUMM", pre);
   }
}

main(int nargs, char **args)
{
   char *fnam, *path, *aliases[16];
   ATL_r1node_t *r1b, *r1p, *r1B;
   int i, xu, yu, l1mul, iflag, L1Elts;
   char pre;

   GetFlags(nargs, args, &pre, &fnam, &path, &xu, &yu, &l1mul, &iflag);

/*
 * If we just want simple tuning header, no need to read file for details
 */
   if (!fnam)
   {
      r1b = GetR1Node();
      r1b->next = NULL;
      r1b->XU = xu;
      r1b->YU = yu;
      r1b->flag = iflag;
      r1b->CacheElts = l1mul * GetL1CacheElts(pre);
      r1hgen(pre, path, -1, r1b, NULL);
      exit(0);
   }
/*
 * Otherwise, we should be doing a full-blown install; read in summary file
 */
   r1b = ReadR1File(fnam);
   SetAllR1TypeFlags(pre, r1b);
/*
 * Find out which are geniune kernels, and which are aliased
 */
   r1B = GetSortedUniqueR1Kerns(pre, r1b, aliases);
/*
 * Generate prototype file for all routines
 */
   r1khgen(pre, path, r1B, aliases);
/*
 * For each cache level, generate a header file which provides macros
 * describing the best GER kernels
 */
   r1p = r1b;
   for (i=0; i < 3; i++)
   {
      r1hgen(pre, path, i, r1p->next, r1p);
      r1p = r1p->next->next;
   }
/*
 * Generate Makefiles to compile all GER kernels (including those used by
 * SYR and SYR2).  These Makefiles & kernels will wind up in
 *   BLDdir/src/blas/ger/
 */
   EmitMakefile(pre, path, r1B);
/*
 * Get required .c kernel files
 */
   GenKernFiles(pre, path, r1B);
/*
 * Generate headers for SYR & SYR2
 */
   L1Elts = GetL1CacheElts(pre);
   r1p = r1b->next->next->next->next->next->next;  /* L1-blocked kernels */
   s2hgen(r1b, L1Elts, pre, path);
   for(i=0; i < 3; i++)
      s1hgen(r1b, i, L1Elts, pre, path);
   KillAllR1Nodes(r1b);
   KillAllR1Nodes(r1B);
   return(0);
}
