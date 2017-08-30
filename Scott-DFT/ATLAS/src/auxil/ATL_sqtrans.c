/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 * Copyright (C) 2009 R. Clint Whaley
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
#define NB 32
#define ATL_MulByNB(n_) ((n_)<<5)
#define ATL_DivByNB(n_) ((n_)>>5)

static void Mjoin(PATL,sqtrans0)(ATL_CINT N, TYPE *C, ATL_CINT ldc)
/*
 * Does an in-place transpose of a square matrix.
 * NOTE: this should only be used on small matrices, as it is not optimized
 *       for the TLB.
 */
{
   ATL_INT j;
/*
 * We will work by reflecting swapping columns & rows across diagonal,
 * starting from the last column, so that early cols are retained in cache
 */
   for (j=N-1; j; j--)
      Mjoin(PATL,swap)(j, C+ldc*(j SHIFT), 1, C+(j SHIFT), ldc);
}

void Mjoin(PATL,sqtrans)(ATL_CINT N, TYPE *C, ATL_CINT ldc)
/*
 * Does an in-place transpose of a square matrix.  This routine is blocked
 * to help with TLB
 */
{
   ATL_CINT Nnb = ATL_MulByNB(ATL_DivByNB(N)), Nr = N - Nnb;
   ATL_INT i, j;

   if (N < NB+NB)
   {
      Mjoin(PATL,sqtrans0)(N, C, ldc);
      return;
   }
/*
 * Loop in reverse order, so first part of matrix retained in cache
 */
   if (Nr)
   {
      for (i=0; i < Nnb; i += NB)
         Mjoin(PATL,geswapT)(NB, Nr, C+((Nnb*ldc+i)SHIFT), ldc,
                             C+((Nnb+i*ldc)SHIFT), ldc);
      Mjoin(PATL,sqtrans0)(Nr, C+((Nnb*(ldc+1))SHIFT), ldc);
   }
   for (j=Nnb-NB; j >= 0; j -= NB)
   {

      for (i=0; i < j; i += NB)
         Mjoin(PATL,geswapT)(NB, NB, C+((j*ldc+i)SHIFT), ldc,
                             C+((j+i*ldc)SHIFT), ldc);
      Mjoin(PATL,sqtrans0)(NB, C+((j*(ldc+1))SHIFT), ldc);
   }
}

/*
 * temporary F77 interfaces for testing
 */
Mjoin(PRE,sqtrans_)(int *N, TYPE *C, int *ldc)
{
   Mjoin(PATL,sqtrans)(*N, C, *ldc);
}
Mjoin(PRE,gemovet_)(int *M, int *N, TYPE *A, int *lda, TYPE *B, int *ldb)
{
   #ifdef TREAL
      #define ONE 1.0
   #else
      TYPE ONE[2] = {1.0, 0.0};
   #endif
   Mjoin(PATL,gemoveT)(*M, *N, ONE, A, *lda, B, *ldb);
}
