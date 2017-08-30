#define ATL_TUNING
/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 *                    (C) Copyright 1999 R. Clint Whaley
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
#include "atlas_lvl2.h"
#include "atlas_lvl3.h"
#if defined(ATL_INL1)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),r1_L1.h))
#elif defined(ATL_INL2)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),r1_L2.h))
#else
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),r1.h))
#endif


#ifdef TREAL
   #ifdef ATL_INL1
      #define ATL_ger Mjoin(PATL,ger_L1)
   #elif defined(ATL_INL2)
      #define ATL_ger Mjoin(PATL,ger_L2)
   #else
      #define ATL_ger Mjoin(PATL,ger)
   #endif
#else
   #ifdef Conj_
      #ifdef ATL_INL1
         #define ATL_ger Mjoin(PATL,gerc_L1)
      #elif defined(ATL_INL2)
         #define ATL_ger Mjoin(PATL,gerc_L2)
      #else
         #define ATL_ger Mjoin(PATL,gerc)
      #endif
   #else
      #ifdef ATL_INL1
         #define ATL_ger Mjoin(PATL,geru_L1)
      #elif defined(ATL_INL2)
         #define ATL_ger Mjoin(PATL,geru_L2)
      #else
         #define ATL_ger Mjoin(PATL,geru)
      #endif
   #endif
#endif
void ATL_ger(const int M, const int N, const SCALAR alpha0,
             const TYPE *X, const int incX, const TYPE *Y, const int incY,
             TYPE *A, const int lda)
/*
 * This dumb version just to get things working.  Later on, only copy
 * Y if we are blocking along X, or the kernel requires it, handle
 * restricted usages, see if we want to enforce alignment, etc.
 */
{
   void (*getX)(const int N, const SCALAR alpha, const TYPE *X,
                const int incX, TYPE *Y, const int incY);
   void *vp=NULL;
   size_t t1, t2;
   TYPE *x = (TYPE*)X, *y = (TYPE*)Y;
   const int ALPHA_IS_ONE = SCALAR_IS_ONE(alpha0);
   int APPLYALPHAX=0, COPYX, COPYY, ALIGNX2A=0;
   int incy=1, imb, mb, k, m, Nm=N, nr=0;
   #ifdef TREAL
      #define one ATL_rone
      TYPE alpha = alpha0;
   #else
      TYPE one[2] = {ATL_rone, ATL_rzero}, *alpha=(TYPE*)alpha0;
   #endif
   void (*gerk)(ATL_CINT, ATL_CINT, const SCALAR, const TYPE*, ATL_CINT,
                const TYPE*, ATL_CINT, TYPE*, ATL_CINT);

   if (M < 1 || N < 1 || SCALAR_IS_ZERO(alpha))
      return;
/*
 * ATLAS's GER kernels loop over M in inner loop, which is bad news if M is
 * very small.  Call code that requires no copy of A & B for these degenerate
 * cases
 */
   if (M < 16)
   {
      #ifdef Conj_
         Mjoin(PATL,gerck_Mlt16)(M, N, alpha, X, incX, Y, incY, A, lda);
      #else
         Mjoin(PATL,gerk_Mlt16)(M, N, alpha, X, incX, Y, incY, A, lda);
      #endif
      return;
   }
/*
 * For very small N, we can't afford the data copy, so call AXPY-based routine
 */
   if (N < 4)
   {
      #ifdef Conj_
         Mjoin(PATL,gerck_axpy)(M, N, alpha, X, incX, Y, incY, A, lda);
      #else
         Mjoin(PATL,gerk_axpy)(M, N, alpha, X, incX, Y, incY, A, lda);
      #endif
      return;
   }
   #ifndef ATL_R1NOBLOCK
      imb = mb = M;
   #else
      ATL_GetPartR1(A, lda, mb, i);
      imb = mb;
   #endif
   #ifdef ATL_r1ALIGNX2A
      ALIGNX2A = 1;
   #endif
   #ifdef ATL_r1NMUL
      Nm = (N/ATL_r1NU)*ATL_r1NU;
      nr = N - Nm;
   #endif
   gerk = ATL_GERK;
   #ifdef ATL_r1USERESTRICTK
      if (ATL_r1UseRestrictK(M, N, A, lda))
      {
         #ifndef ATL_R1NOBLOCKr
            imb = mb = M;
         #else
            ATL_GetPartR1r(A, lda, mb, i);
            imb = mb;
         #endif
         #ifdef ATL_r1ALIGNX2Ar
            ALIGNX2A = 1;
         #else
            ALIGNX2A = 0;
         #endif
         #ifdef ATL_r1NMULr
            Nm = (N/ATL_r1NUr)*ATL_r1NUr;
            nr = N - Nm;
         #elif defined(ATL_r1NMUL)
            Nm = N;
            nr = 0;
         #endif
         gerk = ATL_GERKr;
      }
   #endif
   #ifdef Conj_
      COPYY = 1;
   #else
      COPYY = (incY != 1);
   #endif
   COPYX = (incX != 1);
   if (ALIGNX2A && !COPYX)
   {
      t1 = (size_t) A;
      t2 = (size_t) X;
      COPYX = (t1-(t1/ATL_Cachelen)*ATL_Cachelen) !=
              (t2-(t2/ATL_Cachelen)*ATL_Cachelen);
   }
/*
 * If both vectors are copied or neither is, apply alpha to the shortest
 * vector (apply to Y on tie in order to match lapack).
 * If exactly one vector is copied, then apply alpha to the one we are already
 * copying, regardless of size (copy with scale is typically not noticably
 * slower than copying, despite the extra flops).
 */
   if (COPYX != COPYY)
      APPLYALPHAX = COPYX;
   else if (!COPYY && !COPYX)
   {
      APPLYALPHAX = M < N;
      if (!ALPHA_IS_ONE)
      {
         COPYX = APPLYALPHAX;
         COPYY = !APPLYALPHAX;
      }
   }
   else /* if (COPYY && COPYX) */
      APPLYALPHAX = M < N;

   if (COPYX | COPYY)
   {
      k = Mmax(imb, mb);
      k = Mmin(k, M);
      vp = malloc(ATL_MulBySize(COPYX*k+COPYY*N) + 2*ATL_Cachelen);
/*
 *    If we cannot allocate enough space to copy the vectors, give up and
 *    call the axpy-based implementation
 */
      if (!vp)
      {
         Mjoin(PATL,gerk_axpy)(M, N, alpha, X, incX, Y, incY, A, lda);
         return;
      }
      if (COPYY)
      {
         y = ATL_AlignPtr(vp);
         x = y + (N SHIFT);
         x = (ALIGNX2A) ? ATL_Align2Ptr(x, A) : ATL_AlignPtr(x);
         if (!APPLYALPHAX && !ALPHA_IS_ONE)  /* need to apply alpha to Y */
         {
            #ifdef Conj_
               Mjoin(PATL,moveConj)(N, alpha, Y, incY, y, 1);
            #else
               Mjoin(PATL,cpsc)(N, alpha, Y, incY, y, 1);
            #endif
            alpha = one;
         }
         else  /* do not apply alpha */
         #ifdef Conj_
            Mjoin(PATL,copyConj)(N, Y, incY, y, 1);
         #else
            Mjoin(PATL,copy)(N, Y, incY, y, 1);
         #endif
      }
      else if (ALIGNX2A)
         x = ATL_Align2Ptr(vp, A);
      else
         x = ATL_AlignPtr(vp);
   }
   getX = (COPYX) ? Mjoin(PATL,cpsc) : NULL;
   m = M;
   imb = (imb) ? imb : mb;
   imb = Mmin(imb, m);
   do
   {
/*
 *    Copy X if necessary
 */
      if (getX)
         getX(imb, alpha, X, incX, x, 1);
      else
         x = (TYPE*) X;
/*
 *    Call optimized kernel (can be restricted or general)
 */
      gerk(imb, Nm, one, x, 1, y, incy, A, lda);
/*
 *    Some kernels require N%NU=0; if so nr is remainder, do cleanup with axpy
 */
      if (nr)
         Mjoin(PATL,gerk_axpy)(imb, nr, one, x, 1, y+(Nm SHIFT), 1,
                               A+lda*(Nm SHIFT), lda);
      A += imb SHIFT;
      X += (imb*incX)SHIFT;
      m -= imb;
      imb = Mmin(m,mb);
   }
   while(m);
   if (vp)
      free(vp);
}
