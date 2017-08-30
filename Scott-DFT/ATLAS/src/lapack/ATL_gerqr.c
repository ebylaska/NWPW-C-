/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 * Copyright (C) 2009 Siju Samuel
 *
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
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

#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_lapack.h"
#include "atlas_lvl3.h"
#if defined(ATL_NCPU) && ATL_NCPU > 1
   #include "atlas_tcacheedge.h"
#else
   #include "atlas_cacheedge.h"
#endif

#ifdef CacheEdge
   #if CacheEdge > 4194304 || CacheEdge == 0
      #define LA_CE 262144
   #else
      #define LA_CE 262144
   #endif
#else
   #define LA_CE 262144
#endif

int ATL_gerqr(int M, int N, TYPE *A, int LDA, TYPE  *TAU,
               TYPE *ws_LQ2, TYPE *ws_T, int LDT,
               TYPE *WORKM, int buildT)
/*
 * This is a recursive implementation of ATL_dgeqrf.c; it performs a QR
 * factorization of a panel (M > N) with a bottom level of ATL_geqr2. The
 * recursion is on columns only; it divides by 2 until it reaches a
 * stopping point; at which time it calls ATL_geqr2 to complete a sub-panel,
 * ATL_larft and ATL_larfb to propagate the results, etc.
 *
 * ATL_gerqr.c :
 * int ATL_gerqr(int M, int N, TYPE *A, int LDA, TYPE  *TAU,
 *               TYPE *ws_QR2, TYPE *ws_T, int LDT,
 *               TYPE *WORKM, int buildT)
 *      NOTE :   ATL_geqr2.c will get compiled to four precisions
 *               single precision real,      double precision real
 *               single precision complex,   double precision complex
 *  Purpose
 *  =======
 *
 *  ATL_gerqr computes a QR factorization of a real M-by-N matrix A:
 *  A = R * Q.
 *
 *  Arguments
 *  =========
 *
 *  M       (input) INTEGER
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) INTEGER
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (input/output)  array, dimension (LDA,N)
 *          On entry, the M-by-N matrix A.
 *          On exit, the elements on and above the diagonal of the array
 *          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
 *          upper triangular if m >= n); the elements below the diagonal,
 *          with the array TAU, represent the orthogonal matrix Q as a
 *          product of min(m,n) elementary reflectors (see Further
 *          Details).
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *          If LDA == 0, returns work sizes in *WORKV and *WORKM.
 *
 *  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
 *          The scalar factors of the elementary reflectors (see Further
 *          Details).
 *
 *  ws_RQ2  (workspace) workspace for rq2 factorization. To be allocated
 *          with space of max(M,N)
 *
 *  ws_T    (input/output).  Is the size of T matrix. To be allocated
 *          with a space of min(M,N) X min(M,N). If buildT flag is true,
 *          T is computed and populated as output. If buildT is false,
 *          T must not be used as output
 *
 *  LDT     (input) INTEGER
 *          The leading dimension of the array T.  LDT >= max(1,min(M,N)).
 *
 *
 *  WORKM   Work space matrix, double precision matrix at least M rows by
 *          N columns; the amount used by larft and larfb.
 *
 *  buildT  If non-zero, dgerqr will build in ws_T the complete T necessary
 *          for the original panel it is passed;
 *          such that Q= I - transpose(Y) * T * Y.
 *          If zero, ws_T will contain only those elements of T necessary to
 *          complete the panel.
 */
{

   int I, top, topMN,  bottom, INFO, IINFO, lbuilt, rbuilt;
   int LDA2 = LDA SHIFT;                        /* for complex LDA *2         */
   int LDT2 = LDT SHIFT;                        /* for complex LDT *2         */
   ATL_CINT minMN = Mmin(M, N);
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone, ATL_rzero};
   #else
      #define ONE ATL_rone
   #endif

   if (M < 1 || N < 1) return(0);               /* Nothing to do.             */

//   topMN = (minMN>>2)<<1;
//   bottom = minMN - topMN;
//   top = M - bottom;

//------------------------------------------------------TODO Check Beg
/*
 * Choose a smart recursive column partitioning based on M:
 */
   if (minMN >= NB+NB) /* big prob, put remainder on right */
   {
      topMN = ATL_MulByNB(ATL_DivByNB(minMN>>1));
      bottom = minMN - topMN;
      top  = M - bottom;
   }
   else  /* small prob, keep M mult of MU (MU more critical than NU) */
   {
      bottom = ((minMN>>1)/ATL_mmMU)*ATL_mmMU;
      topMN = minMN - bottom;
      top = M - bottom;
   }
//--------------------------------------------------------TODO Check End


/*
 * Stop recursion if problem cache contained, or no more benefit from L3BLAS
 */
   if (ATL_MulBySize(N)*minMN <= LA_CE || minMN <= 8 || !topMN || !bottom)
   {
      #if 0
      if (minMN >= 4)
      {
         Mjoin(PATL,gemoveT)(N, minMN, ONE, A, LDA, WORKM, N);
         ATL_geql2(N, minMN, WORKM, N, TAU, ws_LQ2);
         Mjoin(PATL,gemoveT)(minMN, N, ONE, WORKM, N, A, LDA);
      }
      else
      #endif
         ATL_gerq2(minMN, N, A+((M -minMN) SHIFT) , LDA, TAU, ws_LQ2);

      if (buildT  || M > minMN)
      {
         ATL_larft(LABackward, LARowStore, N, minMN,
                    A+((M -minMN) SHIFT),   LDA,
                    TAU, ws_T, LDT);            /* Build the T matrix.        */
      }

     if (M > minMN)
     {
        ATL_larfb(CblasRight, CblasNoTrans, /* From right,Not  Transposed, */
              LABackward, LARowStore,       /* Backward, Rowwise.          */
              M-minMN,                      /* Update area is top rows     */
              N,                            /* 'N' columns (remainder)     */
              minMN,                        /* V is 'bottom' columns.      */
              A+((M -minMN) SHIFT)  , LDA,  /* V array (inside A)          */
             (ws_T), LDT,                   /* Build T at T[top,top].      */
              A, LDA,                       /* Array to update.            */
              WORKM, M);                    /* Workspace, leading dim.     */

     }

      return(0);                            /* All okay.                   */
   }



   /*-------------------------------------------------------------*/
   /* On the bottom half, we use the same workspaces.             */
   /* Because we know we have a top hand side we must always      */
   /* build T, so we can multiply by Q before doing the top side. */
   /*-------------------------------------------------------------*/
   ATL_gerqr(bottom,N, (A+(top SHIFT)), LDA, (TAU+(topMN SHIFT)), ws_LQ2,
             ( ws_T+(topMN SHIFT)+topMN*LDT2), LDT,  /* Build T at T[top,top].*/
              WORKM, 1);

   /*---------------------------------------------------------------*/
   /* Now we must adjust the top hand side according to our T.    */
   /* We must apply H'                                             */
   /*---------------------------------------------------------------*/


   ATL_larfb(CblasRight, CblasNoTrans,       /* From right,Not  Transposed,  */
              LABackward, LARowStore,        /* Backward, Rowwise.           */
              top,                           /* Update area is top rows ,    */
              N,                             /* 'N' columns (remainder)      */
              bottom,                        /* V is 'bottom' columns.       */
              (A +(top SHIFT))  , LDA,       /* V array (inside A)           */
             ( ws_T+(topMN SHIFT)+topMN*LDT2), LDT, /* Build T at T[top,top].*/
              A, LDA,                        /* Array to update.             */
              WORKM, M);                     /* Workspace, leading dim.      */

   /*---------------------------------------------------------------*/
   /* On the top  half,                                             */
   /*---------------------------------------------------------------*/
   ATL_gerqr(top, N-bottom,               /* Shorter columns.       */
              (A), LDA,                   /* A at A[0,0].           */
              (TAU),                      /* TAU at TAU[top].       */
              ws_LQ2,                     /* No returned values.    */
             ( ws_T), LDT,                /* Build T at T[top,top]. */
              WORKM, buildT);             /* No returned values.    */

   /*---------------------------------------------------------------*/
   /* If we build T, the bottom side must be completely built, and  */
   /* the top side should be partially built. We need to fill in    */
   /* the lower  left  hand block, 'bottom' rows by 'top' columns.  */
   /* The formula is -T2 * (Y2 * Y1^T) * T1.                        */
   /* The routine is in ATL_larft.c.                                */
   /*---------------------------------------------------------------*/

   if (buildT )
   {
      ATL_larft_block(LABackward, LARowStore,N, minMN, minMN-bottom, bottom,
         A+((M -minMN) SHIFT), LDA, ws_T, LDT);      /* compute upper corner. */
   }

   return(0);
} /* END ATL_gerqr */

