#include "atlas_misc.h"
#include "atlas_level1.h"
#include "atlas_reflvl2.h"
#include "atlas_reflevel2.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),syr2.h))

#ifndef ATL_S2NX
  #define ATL_S2NX 28
#endif
void Mjoin(PATL,syr2U_k)
(
   ATL_CINT N,
   const TYPE *X,
   const TYPE *Y,
   TYPE *A,
   ATL_CINT lda
)
{
   const TYPE *Xt=X, *Yt=Y;
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
   Mjoin(PATL,refsyr2)(AtlasUpper, nx, ATL_rone, X, 1, Y, 1, A, lda);
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
      ATL_SYR2U_nu(A+m, lda, X+m, Yt);
      Xt += ATL_s2NU SHIFT;
      Yt += ATL_s2NU SHIFT;
      A += (ATL_s2NU*lda)SHIFT;
      m += ATL_s2NU;
      nr -= ATL_s2NU;
   }
   while (nr);
}

void Mjoin(PATL,syr2L_k)
(
   ATL_CINT N,          /* length of X, Y, and both dim of A */
   const TYPE *X,       /* N-length vector X */
   const TYPE *Y,       /* N-length vector Y */
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
         ATL_SYR2L_nu(A, lda, X, Y);
         m -= ATL_s2NU;
         A += ATL_s2NU;
         gerk_OC(m, ATL_s2NU, one, X+ATL_s2NU, 1, Y, 1, A, lda);
         gerk_IC(m, ATL_s2NU, one, Y+ATL_s2NU, 1, X, 1, A, lda);
         A += lda*ATL_s2NU;
         X += ATL_s2NU;
         Y += ATL_s2NU;
      }
   }
   else
      nx = N;
   Mjoin(PATL,refsyr2)(AtlasLower, nx, ATL_rone, X, 1, Y, 1, A, lda);
}

void Mjoin(PATL,syr2)(const enum ATLAS_UPLO Uplo, ATL_CINT N,
                      const SCALAR alpha, const TYPE *X, ATL_CINT incX,
                      const TYPE *Y, ATL_CINT incY, TYPE *A, ATL_CINT lda)
{
   TYPE *x, *y;
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
   int APPLYALPHATOX=1, COPYX, COPYY;
   const int ALPHA_IS_ONE = SCALAR_IS_ONE(alpha);

   if (N < 1 || SCALAR_IS_ZERO(alpha))
      return;
   ATL_GetPartS2(A, lda, MB, i);
   MB = (MB/ATL_s2MU)*ATL_s2MU;
   #ifdef ATL_s2ALIGNX2A
      COPYX = 1;
   #else
      COPYX = (incX != 1);
   #endif
   COPYY = (incY != 1);
   if (COPYX == 0 && COPYY == 0)        /* Neither presently copied */
      APPLYALPHATOX = COPYX = 1;
   else if (COPYX == COPYY)             /* both are being copied */
      APPLYALPHATOX = 0;
   else
      APPLYALPHATOX = COPYX;
   i = (COPYX+COPYY)*N;
   if (i)
   {
      vp = malloc(ATL_MulBySize(i) + 2*ATL_Cachelen);
      if (!vp)
      {
         Mjoin(PATL,refsyr2)(Uplo, N, alpha, X, incX, Y, incY, A, lda);
         return;
      }
      if (COPYY)
      {
         y = COPYX ? (void*)(((char*)vp)+ATL_MulBySize(N)+ATL_Cachelen) : vp;
         y = ATL_AlignPtr(y);
         if (!APPLYALPHATOX && !ALPHA_IS_ONE)
            Mjoin(PATL,cpsc)(N, alpha, Y, incY, y, 1);
         else
            Mjoin(PATL,copy)(N, Y, incY, y, 1);
      }
      else
         y = (TYPE*)Y;
      if (COPYX)
      {
         #ifdef ATL_s2ALIGNX2A
            x = ATL_Align2Ptr(vp, A);
         #else
            x = ATL_AlignPtr(vp);
         #endif
         if (APPLYALPHATOX && !ALPHA_IS_ONE)
            Mjoin(PATL,cpsc)(N, alpha, X, incX, x, 1);
         else
            Mjoin(PATL,copy)(N, X, incX, x, 1);

      }
      else
         x = (TYPE*) X;
   }
   else
   {
      x = (TYPE*) X;
      y = (TYPE*) Y;
   }
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
         Mjoin(PATL,syr2U_k)(mb, x+i, y+i, A+i+i*lda, lda);
         for (j=i+mb; j < N; j += ATL_s2NU)
         {
            nb = N-j;
            nb = Mmin(nb, ATL_s2NU);
            gerk_OC(mb, nb, one, x+i, 1, y+j, 1, A+i+j*lda, lda);
            gerk_IC(mb, nb, one, y+i, 1, x+j, 1, A+i+j*lda, lda);
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
            gerk_OC(mb, nb, one, x+i, 1, y+j, 1, A+i+j*lda, lda);
            gerk_IC(mb, nb, one, y+i, 1, x+j, 1, A+i+j*lda, lda);
         }
         Mjoin(PATL,syr2L_k)(mb, x+i, y+i, A+i*(lda+1), lda);
      }
   }
   if (vp)
      free(vp);
}
