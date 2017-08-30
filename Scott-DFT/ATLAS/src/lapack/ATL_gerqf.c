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

#ifdef  SREAL
     #define MYOPT LASreal
#endif
#ifdef  DREAL
    #define MYOPT  LADreal
#endif
#ifdef  SCPLX
    #define MYOPT  LAScplx
#endif
#ifdef  DCPLX
    #define MYOPT  LADcplx
#endif


int ATL_gerqf(ATL_CINT M, ATL_CINT N, TYPE  *A, ATL_CINT lda, TYPE  *TAU,
               TYPE *WORK, ATL_CINT LWORK)
/*-----------------------------------------------------------------------------
 * This is the C translation of the standard LAPACK Fortran routine:
 *      SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
 *
 * ATL_gerqf.c :
 * int ATL_gerqf(int M, int N, TYPE  *A, int LDA, TYPE  *TAU,
 *              TYPE *WORK, int LWORK)
 *
 *  Purpose
 *  =======
 *
 *  ATL_gerqf  computes a RQ factorization of a real/complex M-by-N matrix A:
 *  A = R * Q.
 *
 *  Compared to LAPACK, here, a recursive pannel factorization is implemented.
 *  Refer ATL_gerqr.c andd ATL_larft.c for details.
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
 *          On entry, the m by n matrix A.
 *          On exit, if m <= n, the upper triangle of the subarray
 *          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
 *          if m >= n, the elements on and above the (m-n)-th subdiagonal
 *          contain the m by n upper trapezoidal matrix R; the remaining
 *          elements, with the array TAU, represent the orthogonal matrix
 *          (unitary matrix incase of complex precision )  as a
 *          as a product of elementary reflectors (see Further
 *          Details).
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,M).
 *
 *  TAU     (output) array, dimension (min(M,N))
 *          The scalar factors of the elementary reflectors (see Further
 *          Details).
 *
 *  WORK    (workspace/output) array, dimension (MAX(1,LWORK))
 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *
 *  LWORK   (input) INTEGER
 *          The dimension of the array WORK.  LWORK >= max(1,N).
 *          For optimum performance LWORK >= N*NB, where NB is
 *          the optimal blocksize.
 *
 *          If LWORK = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the WORK array, returns
 *          this value as the first entry of the WORK array, and no error
 *          message related to LWORK is issued .
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *
 *  Further Details
 *  ===============
 *
 *  The matrix Q is represented as a product of elementary reflectors
 *
 *     Q = H(1) H(2) . . . H(k), where k = min(m,n).    (for Real precison)
 *     Q = H(1)' H(2)' . . . H(k)', where k = min(m,n). (for Complex Precison)
 *         (Note : Conjugate Transpose of H is taken above)
 *
 *  Each H(i) has the form
 *
 *     H(i) = I - tau * v * v'                 (for Real Precision)
 *     H(i) = I - tau * v * conjugate(v)'      (for Complex  Precision)
 *
 *  where tau is a real scalar, and v is a real vector with
 *  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in
 *  A(m-k+i,1:n-k+i-1), and tau in TAU(i).
 *----------------------------------------------------------------------------*/
{
   ATL_CINT minMN = Mmin(M, N), maxMN = Mmax(M, N);
   ATL_INT n, nb, j;
   TYPE  *ws_RQ2,  *ws_T, *ws_larfb;   /* Workspace for RQ2,T, larfb          */
   void *vp=NULL;

   nb = clapack_ilaenv(LAIS_OPT_NB, LAgeqrf, MYOPT+LALeft+LAUpper, M, N,-1,-1);

/*
 * If it is a workspace query, return the size of work required.
 *    wrksz = wrksz of ATL_larfb + ATL_larft + ATL_gerq2
 * RCW Q: Why can't LARFB & GERQ2 workspaces be overlapped?
 */
   if (LWORK < 0)
   {
      *WORK = ( M*nb + nb*nb + maxMN )  ;
      return(0);
   }
   else if (M < 1 || N < 1)  /* quick return if no work to do */
      return(0);
/*
 * If the user gives us too little space, see if we can allocate it ourselves
 */
   else if (LWORK < (M*nb + nb*nb + maxMN))
   {
      vp = malloc(ATL_MulBySize(M*nb + nb*nb + maxMN) + ATL_Cachelen);
      if (!vp)
         return(-7);
       WORK = ATL_AlignPtr(vp);
   }
/*
 * QL is transpose of RQ; Use this fact to go from row-major RQ to
 * col-major QL (typically more than 10% faster even wt transpose costs)
 */
   if (M == N && N >= 128)
   {
      Mjoin(PATL,sqtrans)(N, A, lda);
      n = ATL_geqlf(M, N, A, lda, TAU, WORK, LWORK);
      Mjoin(PATL,sqtrans)(N, A, lda);
      return(n);
   }

/*
 * Assign workspace areas for ATL_larft, ATL_gerq2, ATL_larfb
 */
   ws_T = WORK;                         /* T at begining of work */
   ws_RQ2 = WORK +(nb SHIFT)*nb;        /* After T Work space             */
   ws_larfb = ws_RQ2 + (maxMN SHIFT);   /* After workspace for T and RQ2  */

/*
 * Leave one iteration to be done outside loop, so we don't build T
 * Any loop iterations are therefore known to be of size nb (no partial blocks)
 */
   n = (minMN / nb) * nb;
   if (n == minMN)
      n -= Mmin(nb, minMN);  // I am not sure why miMN is taking here SIJU

   for (j=0; j < n; j += nb)
   {

/*
 *    Starts Factorization from the bottom row
 */
      ATL_assert(!ATL_gerqr(nb, N-j, A+((M -j -nb) SHIFT), lda,
                            TAU+((minMN -(j+nb)) SHIFT),
                            ws_RQ2, ws_T, nb, ws_larfb, 1));
      if (j+nb < M)  /* if there are more cols left on top to update them */
      {
/*
 *       ======================================================================
 *       Form the triangular factor of the block reflector
 *          H = H(i+ib-1) . . . H(i+1) H(i)
 *       After gerqr, ws_T contains 'T', the nb x nb triangular factor 'T'
 *       of the block reflector. It is an output used in the next call, dlarfb.
 *
 *       The ws_T array used above is an input to ATL_larfb; it is 'T' in
 *       that routine, and LDT x K (translates here to NB  x NB).
 *       WORK is an LDWORK x NB workspace in ATL_larfb(maxMN X NB).
 *       ======================================================================
 */

         ATL_larfb(CblasRight, CblasNoTrans, LABackward, LARowStore,
                   M-j-nb, N-j, nb,
                   A+((M - j -nb ) SHIFT),    /* Pointet to A in larfb        */
                   lda, ws_T, nb,
                   A,                        /* Pointer to C in larfb         */
                   lda, ws_larfb, M);
      }
   }

/*
 *  Build Last panel.  buildT  is passed as 0, as there is no need to build T
 */
   nb = minMN - n;       /* remianing nb                                      */

   ATL_assert(!ATL_gerqr(M-n, N-n, A, lda, TAU,
                         ws_RQ2, ws_T, nb, ws_larfb, 0));

   if (vp)
      free(vp);
   return(0);
} /* END ATL_gerqf */

