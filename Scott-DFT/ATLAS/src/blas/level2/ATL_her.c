#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_reflvl2.h"
#include "atlas_reflevel2.h"
#if defined(ATL_INL1)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr_L1.h))
   #define ATL_her Mjoin(PATL,her_L1)
   #define ATL_gerK  Mjoin(PATL,gerk_L1)
   #define ATL_gerKr Mjoin(PATL,gerk_L1_restrict)
#elif defined(ATL_INL2)
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr_L2.h))
   #define ATL_her Mjoin(PATL,her_L2)
   #define ATL_gerK  Mjoin(PATL,gerk_L2)
   #define ATL_gerKr Mjoin(PATL,gerk_L2_restrict)
#else
   #include Mstr(Mjoin(Mjoin(atlas_,PRE),syr.h))
   #define ATL_her Mjoin(PATL,her)
   #define ATL_gerK  Mjoin(PATL,gerk_L0)
   #define ATL_gerKr Mjoin(PATL,gerk_L0_restrict)
#endif

typedef void (*ATL_gerk_t)(ATL_CINT, ATL_CINT, const SCALAR, const TYPE*,
                           ATL_CINT, const TYPE*, ATL_CINT, TYPE*, ATL_CINT);

#ifndef ATL_S1NX
   #define ATL_S1NX 32
#endif
Mjoin(PATL,her_kU)
(
   ATL_gerk_t gerk,             /* func ptr to selected GER kernel */
   ATL_CINT N,                  /* size of prob to solve */
   const TYPE alpha,            /* alpha */
   const TYPE *x,               /* input vector X */
   const TYPE *xh,              /* alpha*X^H */
   TYPE *A,                     /* hermitian matrix, A = A + x*xh */
   ATL_CINT lda                 /* row stride of A */
)
{
   ATL_INT nx=Mmin(ATL_S1NX,N), j;
   TYPE one[2] = {ATL_rone, ATL_rzero};
   ATL_CINT lda2 = lda+lda;

   j = N - nx;
   j = (j/ATL_s1NU)*ATL_s1NU;
   if (j != N-nx)
      nx += N-nx-j;
   Mjoin(PATL,refher)(AtlasUpper, nx, alpha, x, 1, A, lda);
   for (j=nx; j < N; j += ATL_s1NU)
   {
      gerk(j, ATL_s1NU, one, x, 1, xh+j+j, 1, A+j*lda2, lda);
      ATL_HER1U_nu(A+j*(lda2+2), lda, x+j+j, xh+j+j);
   }
}

Mjoin(PATL,her_kL)
(
   ATL_gerk_t gerk,             /* func ptr to selected GER kernel */
   ATL_CINT N,                  /* size of prob to solve */
   const TYPE alpha,            /* alpha */
   const TYPE *x,               /* input vector X */
   const TYPE *xh,              /* alpha*X^H */
   TYPE *A,                     /* hermitian matrix, A = A + x*xh */
   ATL_CINT lda                 /* row stride of A */
)
{
   ATL_INT nx=Mmin(ATL_S1NX,N), i, NN;
   ATL_CINT lda2 = lda+lda;
   const TYPE one[2] = {ATL_rone, ATL_rzero};

   i = N - nx;
   i = (i/ATL_s1NU)*ATL_s1NU;
   if (i != N-nx)
      nx += N-nx-i;
   NN = N - nx;
   for (i=0; i < NN; i += ATL_s1NU)
   {
      ATL_HER1L_nu(A, lda, x, xh);
      gerk(N-i-ATL_s1NU, ATL_s1NU, one, x+ATL_s1NU+ATL_s1NU, 1, xh, 1,
           A+ATL_s1NU+ATL_s1NU, lda);
      A += ATL_s1NU*(lda2+2);
      xh += ATL_s1NU+ATL_s1NU;
      x += ATL_s1NU+ATL_s1NU;
   }
   Mjoin(PATL,refher)(AtlasLower, nx, alpha, x, 1, A, lda);
}

void ATL_her(const enum ATLAS_UPLO Uplo, ATL_CINT N, const TYPE alpha,
               const TYPE *X, ATL_CINT incX, TYPE *A, ATL_CINT lda)
{
   const TYPE one[2] = {ATL_rone, ATL_rzero}, calpha[2] = {alpha, ATL_rzero};
   ATL_CINT lda2 = lda+lda, incx = incX+incX;
   void *vp=NULL;
   TYPE *x, *xh;
   #ifdef ATL_s1USERRESTRICTK
      ATL_gerk_t gerk;
   #else
      #define gerk ATL_gerK
   #endif
   ATL_INT MB, NB, mb, nb, Nmb, i;
   int COPYX=0, ALIGNX2A=0;
   const int ALPHA_IS_ONE=(alpha == ATL_rone);

   if (N < 1 || (alpha == ATL_rzero))
      return;
/*
 * For very small problems, avoid overhead of func calls & data copy
 */
   if (N < 50)
   {
      Mjoin(PATL,refher)(Uplo, N, alpha, X, incX, A, lda);
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
   i = N + COPYX*MB;
   vp = malloc(ATL_MulBySize(i)+2*ATL_Cachelen);
   if (!vp)
   {
      Mjoin(PATL,refher)(Uplo, N, alpha, X, incX, A, lda);
      return;
   }
   xh = ATL_AlignPtr(vp);
   if (COPYX)
   {
      x = xh + N+N;
      x = ALIGNX2A ? ATL_Align2Ptr(x, A) : ATL_AlignPtr(x);
   }
   else
      x = (TYPE*) X;
   if (ALPHA_IS_ONE)
      Mjoin(PATL,copyConj)(N, X, incX, xh, 1);
   else
      Mjoin(PATL,moveConj)(N, calpha, X, incX, xh, 1);
   Nmb = ((N-1)/MB)*MB;
   if (Uplo == AtlasUpper)
   {
      for (i=0; i < Nmb; i += MB)
      {
         if (COPYX)
            Mjoin(PATL,copy)(MB, X+i*incx, incX, x, 1);
         Mjoin(PATL,her_kU)(gerk, MB, alpha, x, xh+i+i, A+i*(lda2+2), lda);
         gerk(MB, N-i-MB, one, x, 1, xh+((i+MB)<<1), 1, A+(MB+i)*lda2+i+i, lda);
         if (!COPYX)
            x += MB+MB;
      }
      mb = N - Nmb;
      if (COPYX)
         Mjoin(PATL,copy)(mb, X+Nmb*incx, incX, x, 1);
      Mjoin(PATL,her_kU)(gerk, mb, alpha, x, xh+Nmb+Nmb, A+Nmb*(lda2+2), lda);
   }
   else         /* Uplo == AtlasLower */
   {
      mb = N - Nmb;
      if (COPYX)
         Mjoin(PATL,copy)(mb, X, incX, x, 1);
      Mjoin(PATL,her_kL)(gerk, mb, alpha, x, xh, A, lda);
      for (i=mb; i < N; i += MB)
      {
         if (COPYX)
            Mjoin(PATL,copy)(MB, X+i*incx, incX, x, 1);
         else
            x += mb+mb;
         gerk(MB, i, one, x, 1, xh, 1, A+i+i, lda);
         Mjoin(PATL,her_kL)(gerk, MB, alpha, x, xh+i+i, A+i*(lda2+2), lda);
         mb = MB;
      }
   }

   if (vp)
     free(vp);
}
