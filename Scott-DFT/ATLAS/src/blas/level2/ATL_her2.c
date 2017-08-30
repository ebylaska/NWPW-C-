#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_reflvl2.h"
#include "atlas_reflevel2.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),syr2.h))

#ifndef ATL_S2NX
  #define ATL_S2NX 28
#endif
void Mjoin(PATL,her2U_k)
(
   ATL_CINT N,
   const TYPE *X,
   const TYPE *Y,
   const SCALAR alpha,  /* Need orig alpha to pass to ref HER2 */
   const TYPE *Xt,
   const TYPE *Yt,
   TYPE *A,
   ATL_CINT lda
)
{
   ATL_INT i, m, nx=(ATL_S2NX<=N) ? ATL_S2NX : N, nr = N-nx;
   #ifdef TREAL
      #define one ATL_rone
   #else
      TYPE one[2] = {ATL_rone, ATL_rzero};
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      void (*gerk_IC)(ATL_CINT, ATL_CINT, const double, const double*,
                      ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
   #else
      #define gerk_IC ATL_R1IC
   #endif
   #ifdef ATL_s2USERESTRICTK_OC
      void (*gerk_OC)(ATL_CINT, ATL_CINT, const double, const double*,
                      ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
   #else
      #define gerk_OC ATL_R1OC
   #endif

   i = (nr/ATL_s2NU)*ATL_s2NU;
   if (nr != i)
   {
      nx += nr - i;
      nr = i;
   }
   nx = (nx <= N) ? nx : N;
   Mjoin(PATL,refher2)(AtlasUpper, nx, alpha, X, 1, Y, 1, A, lda);
   if (nx == N)
      return;
   A += lda*nx SHIFT;
   Yt += nx SHIFT;
   Xt += nx SHIFT;
   #ifdef ATL_s2USERESTRICTK_OC
      gerk_OC = ATL_s2UseRestrictK_OC(M, N, A, lda) ? ATL_R1OCr : ATL_R1OC;
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      gerk_IC = ATL_s2UseRestrictK_IC(M, N, A, lda) ? ATL_R1ICr : ATL_R1IC;
   #endif

   m = nx;
   do
   {
      gerk_OC(m, ATL_s2NU, one, X, 1, Yt, 1, A, lda);
      gerk_IC(m, ATL_s2NU, one, Y, 1, Xt, 1, A, lda);
      ATL_HER2U_nu(A+(m SHIFT), lda, X+(m SHIFT), Y+(m SHIFT), Xt, Yt);
      Xt += ATL_s2NU SHIFT;
      Yt += ATL_s2NU SHIFT;
      A += (ATL_s2NU*lda)SHIFT;
      m += ATL_s2NU;
      nr -= ATL_s2NU;
   }
   while (nr);
}

void Mjoin(PATL,her2L_k)
(
   ATL_CINT N,          /* length of X, Y, and both dim of A */
   const TYPE *X,       /* N-length vector X */
   const TYPE *Y,       /* N-length vector Y */
   const SCALAR alpha,  /* Need orig alpha to pass to ref HER2 */
   const TYPE *Xt,      /* hermitian transpose of X with conj(alpha) applied */
   const TYPE *Yt,      /* hermitian transpose of Y with alpha applied */
   TYPE *A,             /* ldaxN matrix */
   ATL_CINT lda         /* row stride */
)
{
   ATL_INT i, j, m, nn, nx=(ATL_S2NX<=N) ? ATL_S2NX : N, nr = N-nx;
   #ifdef TREAL
      #define one ATL_rone
   #else
      TYPE one[2] = {ATL_rone, ATL_rzero};
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      void (*gerk_IC)(ATL_CINT, ATL_CINT, const double, const double*,
                      ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
   #endif
   #ifdef ATL_s2USERESTRICTK_OC
      void (*gerk_OC)(ATL_CINT, ATL_CINT, const double, const double*,
                      ATL_CINT, const double*, ATL_CINT, double*, ATL_CINT);
   #endif

   #ifdef ATL_s2USERESTRICTK_OC
      gerk_OC = ATL_s2UseRestrictK_OC(M, N, A, lda) ? ATL_R1OCr : ATL_R1OC;
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      gerk_IC = ATL_s2UseRestrictK_IC(M, N, A, lda) ? ATL_R1ICr : ATL_R1IC;
   #endif
   i = (nr/ATL_s2NU)*ATL_s2NU;
   if (nr != i)
   {
      nx += nr - i;
      nr = i;
   }
   if (nx < N)
   {
      nn = N - nx;
      m = N;
      for (j=0; j < nn; j += ATL_s2NU)
      {
         ATL_HER2L_nu(A, lda, X, Y, Xt, Yt);
         X += ATL_s2NU+ATL_s2NU;
         Y += ATL_s2NU+ATL_s2NU;
         A += ATL_s2NU+ATL_s2NU;
         m -= ATL_s2NU;
         gerk_OC(m, ATL_s2NU, one, X, 1, Yt, 1, A, lda);
         gerk_IC(m, ATL_s2NU, one, Y, 1, Xt, 1, A, lda);
         A += lda*(ATL_s2NU+ATL_s2NU);
         Xt += ATL_s2NU+ATL_s2NU;
         Yt += ATL_s2NU+ATL_s2NU;
      }
   }
   else
      nx = N;
   Mjoin(PATL,refher2)(AtlasLower, nx, alpha, X, 1, Y, 1, A, lda);
}

void Mjoin(PATL,her2)(const enum ATLAS_UPLO Uplo, ATL_CINT N,
                      const SCALAR alpha, const TYPE *X, ATL_CINT incX,
                      const TYPE *Y, ATL_CINT incY, TYPE *A, ATL_CINT lda)
{
   TYPE *x, *y;
   TYPE *xh, *yh;        /* hermitian transposes */
   void *vp=NULL;
   #ifdef TCPLX
      const TYPE one[2] = {ATL_rone, ATL_rzero}, alconj[2]={*alpha, -alpha[1]};
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      void (*gerk_IC)(ATL_CINT, ATL_CINT, const TYPE, const TYPE*,
                      ATL_CINT, const TYPE*, ATL_CINT, TYPE*, ATL_CINT);
   #endif
   #ifdef ATL_s2USERESTRICTK_OC
      void (*gerk_OC)(ATL_CINT, ATL_CINT, const TYPE, const TYPE*,
                      ATL_CINT, const TYPE*, ATL_CINT, TYPE*, ATL_CINT);
   #endif
   ATL_INT MB, mb, nb, i, j;
   int COPYX, COPYY;
   const int ALPHA_IS_ONE = SCALAR_IS_ONE(alpha);

   if (N < 1 || SCALAR_IS_ZERO(alpha))
      return;
   ATL_GetPartS2(A, lda, MB, i);
   MB = (MB/ATL_s2MU)*ATL_s2MU;
   #ifdef ATL_s2ALIGNX2A
      COPYY = COPYX = 1;
   #else
      COPYX = (incX != 1);
      COPYY = (incY != 1);
   #endif
   i = (COPYX+COPYY)*N + N+N;
   vp = malloc(ATL_MulBySize(i) + 4*ATL_Cachelen);
   if (!vp)
   {
      Mjoin(PATL,refher2)(Uplo, N, alpha, X, incX, Y, incY, A, lda);
      return;
   }
   xh = ATL_AlignPtr(vp);
   yh = xh + N+N;
   yh = ATL_AlignPtr(yh);
   if (COPYY)
   {
      y = yh + N+N;
      #ifdef ATL_s2ALIGNX2A
         y = ATL_Align2Ptr(y, A);
      #else
         y = ATL_AlignPtr(y);
      #endif
      Mjoin(PATL,copy)(N, Y, incY, y, 1);
   }
   else
      y = (TYPE*) Y;
   if (COPYX)
   {
      x = COPYY ? y+N+N : yh+N+N;
      #ifdef ATL_s2ALIGNX2A
         x = ATL_Align2Ptr(x, A);
      #else
         x = ATL_AlignPtr(x);
      #endif
      Mjoin(PATL,copy)(N, X, incX, x, 1);
   }
   else
      x = (TYPE*) X;
   Mjoin(PATL,moveConj)(N, alpha, y, 1, yh, 1);   /* yh = alpha*conj(y) */
   Mjoin(PATL,moveConj)(N, alconj, x, 1, xh, 1);  /* xh = conj(alpha*x) */
   #ifdef ATL_s2USERESTRICTK_OC
      gerk_OC = ATL_s2UseRestrictK_OC(M, N, A, lda) ? ATL_R1OCr : ATL_R1OC;
   #endif
   #ifdef ATL_s2USERESTRICTK_IC
      gerk_IC = ATL_s2UseRestrictK_IC(M, N, A, lda) ? ATL_R1ICr : ATL_R1IC;
   #endif
   if (Uplo == AtlasUpper)
   {
      mb = mb <= N ? mb : N;
      for (i=0; i < N; i += MB)
      {
         mb = N-i;
         mb = Mmin(mb, MB);
         Mjoin(PATL,her2U_k)(mb, x+i+i, y+i+i, alpha, xh+i+i, yh+i+i,
                             A+((i+i*lda)<<1), lda);
         for (j=i+mb; j < N; j += ATL_s2NU)
         {
            nb = N-j;
            nb = Mmin(nb, ATL_s2NU);
            gerk_OC(mb, nb, one, x+i+i, 1, yh+j+j, 1,
                    A+((i+j*lda)<<1), lda);
            gerk_IC(mb, nb, one, y+i+i, 1, xh+j+j, 1,
                    A+((i+j*lda)<<1), lda);
         }
      }
   }
   else
   {
      for (i=0; i < N; i += MB)
      {
         mb = N-i;
         mb = Mmin(mb, MB);
         for (j=0; j < i; j += ATL_s2NU)
         {
            nb = N-j;
            nb = Mmin(nb, ATL_s2NU);
            gerk_OC(mb, nb, one, x+i+i, 1, yh+j+j, 1,
                    A+((i+j*lda)<<1), lda);
            gerk_IC(mb, nb, one, y+i+i, 1, xh+j+j, 1,
                    A+((i+j*lda)<<1), lda);
         }
         Mjoin(PATL,her2L_k)(mb, x+i+i, y+i+i, alpha, xh+i+i, yh+i+i,
                             A+i*((lda+1)SHIFT), lda);
      }
   }
   if (vp)
      free(vp);
}
