#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_reflvl2.h"
#include "atlas_reflevel2.h"
#if defined(ATL_INL1)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr_L1.h))
   #define ATL_syr Mjoin(PATL,syr_L1)
   #define ATL_gerK  Mjoin(PATL,gerk_L1)
   #define ATL_gerKr Mjoin(PATL,gerk_L1_restrict)
#elif defined(ATL_INL2)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr_L2.h))
   #define ATL_syr Mjoin(PATL,syr_L2)
   #define ATL_gerK  Mjoin(PATL,gerk_L2)
   #define ATL_gerKr Mjoin(PATL,gerk_L2_restrict)
#else
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr.h))
   #define ATL_syr Mjoin(PATL,syr)
   #define ATL_gerK  Mjoin(PATL,gerk_L0)
   #define ATL_gerKr Mjoin(PATL,gerk_L0_restrict)
#endif

typedef void (*ATL_gerk_t)(ATL_CINT, ATL_CINT, const SCALAR, const TYPE*,
                           ATL_CINT, const TYPE*, ATL_CINT, TYPE*, ATL_CINT);

#ifndef ATL_S1NX
   #define ATL_S1NX 32
#endif
Mjoin(PATL,syr_kU)
(
   ATL_gerk_t gerk,             /* func ptr to selected GER kernel */
   ATL_CINT N,                  /* size of prob to solve */
   const TYPE alpha,            /* alpha */
   const TYPE *x,               /* vector X -- may have alpha applied */
   const TYPE *xt,              /* X^T */
   TYPE *A,                     /* symmetric matrix, A = A + x*xt */
   ATL_CINT lda                 /* row stride of A */
)
{
   ATL_INT nx=Mmin(ATL_S1NX,N), j;

   j = N - nx;
   j = (j/ATL_s1NU)*ATL_s1NU;
   if (j != N-nx)
      nx += N-nx-j;
   Mjoin(PATL,refsyr)(AtlasUpper, nx, alpha, xt, 1, A, lda);
   for (j=nx; j < N; j += ATL_s1NU)
   {
      gerk(j, ATL_s1NU, ATL_rone, x, 1, xt+j, 1, A+j*lda, lda);
      ATL_SYR1U_nu(A+j*(lda+1), lda, x+j, xt+j);
   }
}

Mjoin(PATL,syr_kL)
(
   ATL_gerk_t gerk,             /* func ptr to selected GER kernel */
   ATL_CINT N,                  /* size of prob to solve */
   const TYPE alpha,            /* alpha */
   const TYPE *x,               /* vector X -- may have alpha applied */
   const TYPE *xt,              /* X^T */
   TYPE *A,                     /* symmetric matrix, A = A + x*xt */
   ATL_CINT lda                 /* row stride of A */
)
{
   ATL_INT nx=Mmin(ATL_S1NX,N), i, NN;

   i = N - nx;
   i = (i/ATL_s1NU)*ATL_s1NU;
   if (i != N-nx)
      nx += N-nx-i;
   NN = N - nx;
   for (i=0; i < NN; i += ATL_s1NU)
   {
      ATL_SYR1L_nu(A, lda, x, xt);
      gerk(N-i-ATL_s1NU, ATL_s1NU, ATL_rone, x+ATL_s1NU, 1, xt, 1,
           A+ATL_s1NU, lda);
      A += ATL_s1NU*(lda+1);
      xt += ATL_s1NU;
      x += ATL_s1NU;
   }
   Mjoin(PATL,refsyr)(AtlasLower, nx, alpha, xt, 1, A, lda);
}

void ATL_syr(const enum ATLAS_UPLO Uplo, ATL_CINT N, const TYPE alpha,
               const TYPE *X, ATL_CINT incX, TYPE *A, ATL_CINT lda)
{
   void *vp=NULL;
   TYPE *x, *xt, *xx=(TYPE*)X;
   #ifdef ATL_s1USERRESTRICTK
      ATL_gerk_t gerk;
   #else
      #define gerk ATL_gerK
   #endif
   ATL_INT MB, NB, mb, nb, Nmb, i, incx=incX;
   int COPYX=0, ALIGNX2A=0;
   const int ALPHA_IS_ONE=(alpha == ATL_rone);

   if (N < 1 || (alpha == ATL_rzero))
      return;
/*
 * For very small problems, avoid overhead of func calls & data copy
 */
   if (N < 50)
   {
      Mjoin(PATL,refsyr)(Uplo, N, alpha, X, incX, A, lda);
      return;
   }
/*
 * Determine which kernel to use, and whether we need to copy vectors
 */
   COPYX = (incX != 1);
   #ifdef ATL_s1USERESTRICTK
      gerk = ATL_s1UseRestrictK(M, N, A, lda) ? ATL_gerKr : ATL_gerK;
      if (gerk == ATL_gerKr)
      {
         #ifdef R1F_ALIGNX2Ar
            COPYX = COPYX | (ATL_Align2Ptr(X, A) != X);
            ALIGNX2A = 1;
         #endif
         #ifdef ATL_NOBLOCK_S1r
            MB = NB = N;
         #else
            ATL_GetPartS1r(rA, lda, MB, NB);
         #endif
      }
      else
      {
         #if R1F_ALIGNX2A
            COPYX = COPYX | (ATL_Align2Ptr(X, A) != X);
            ALIGNX2A = 1;
         #endif
         #ifdef ATL_NOBLOCK_S1
            MB = NB = N;
         #else
            ATL_GetPartS1(rA, lda, MB, NB);
         #endif
      }
   #else
      #ifdef R1F_ALIGNX2A
         ALIGNX2A = 1;
         COPYX = COPYX | (ATL_Align2Ptr(X, A) != X);
      #endif
      #ifdef ATL_NOBLOCK_S1
         MB = NB = N;
      #else
         ATL_GetPartS1(rA, lda, MB, NB);
      #endif
   #endif
   if (MB > N || MB < 1)
      MB = N;
   if (nb > N || nb < 1)
      nb = N;
   if (ALPHA_IS_ONE)  /* X need not be distinct from Xt */
   {
      if (COPYX)
      {
         vp = malloc(ATL_MulBySize(N)+ATL_Cachelen);
         if (!vp)
         {
            Mjoin(PATL,refsyr)(Uplo, N, alpha, X, incX, A, lda);
            return;
         }
         x = xt = ALIGNX2A ? ATL_Align2Ptr(vp, A) : ATL_AlignPtr(vp);
         Mjoin(PATL,copy)(N, X, incX, x, 1);
         COPYX = 0;
      }
      else
         x = xt = (TYPE*) X;
   }
   else if (incX == 1)          /* apply alpha to X, orig vec Xt */
   {
      COPYX = 1;
      xt = (TYPE*) X;
      vp = malloc(ATL_MulBySize(MB)+ATL_Cachelen);
      if (!vp)
      {
         Mjoin(PATL,refsyr)(Uplo, N, alpha, X, incX, A, lda);
         return;
      }
      x = ALIGNX2A ? ATL_Align2Ptr(vp, A) : ATL_AlignPtr(vp);
   }
   else                         /* must copy both X & Xt, apply alpha to x */
   {
      COPYX = 1;
      vp = malloc(ATL_MulBySize(MB+N)+2*ATL_Cachelen);
      if (!vp)
      {
         Mjoin(PATL,refsyr)(Uplo, N, alpha, X, incX, A, lda);
         return;
      }
      x = ALIGNX2A ? ATL_Align2Ptr(vp, A) : ATL_AlignPtr(vp);
      xt = x + MB;
      xt = ATL_AlignPtr(xt);
      Mjoin(PATL,copy)(N, X, incX, xt, 1);
/*
 *    Set it up so that we copy from contiguous vector, not original X
 */
      xx = xt;
      incx = 1;
   }
   Nmb = ((N-1)/MB)*MB;
   if (Uplo == AtlasUpper)
   {
      for (i=0; i < Nmb; i += MB)
      {
         if (COPYX)
            Mjoin(PATL,cpsc)(MB, alpha, xx+i*incx, incx, x, 1);
         Mjoin(PATL,syr_kU)(gerk, MB, alpha, x, xt+i, A+i*(lda+1), lda);
         gerk(MB, N-i-MB, ATL_rone, x, 1, xt+i+MB, 1, A+(MB+i)*lda+i, lda);
         if (!COPYX)
            x += MB;
      }
      mb = N - Nmb;
      if (COPYX)
         Mjoin(PATL,cpsc)(mb, alpha, xx+Nmb*incx, incx, x, 1);
      Mjoin(PATL,syr_kU)(gerk, mb, alpha, x, xt+Nmb, A+Nmb*(lda+1), lda);
   }
   else         /* Uplo == AtlasLower */
   {
      mb = N - Nmb;
      if (COPYX)
         Mjoin(PATL,cpsc)(mb, alpha, xx, incx, x, 1);
      Mjoin(PATL,syr_kL)(gerk, mb, alpha, x, xt, A, lda);
      for (i=mb; i < N; i += MB)
      {
         if (COPYX)
            Mjoin(PATL,cpsc)(MB, alpha, xx+i*incx, incx, x, 1);
         else
            x += mb;
         gerk(MB, i, ATL_rone, x, 1, xt, 1, A+i, lda);
         Mjoin(PATL,syr_kL)(gerk, MB, alpha, x, xt+i, A+i*(lda+1), lda);
         mb = MB;
      }
   }

   if (vp)
     free(vp);
}
