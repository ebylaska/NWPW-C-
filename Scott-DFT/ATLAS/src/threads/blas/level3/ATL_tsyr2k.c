#include "atlas_misc.h"
#include "atlas_threads.h"
#include "atlas_tlvl3.h"

void Mjoin(PATL,tvsyApAt)(const enum ATLAS_UPLO Uplo, ATL_CINT N, const void *A,
                          ATL_CINT lda, const void *beta, void *C, ATL_CINT ldc)
{
   Mjoin(PATL,syApAt)(Uplo, N, A, lda, SVVAL((const TYPE*)beta), C, ldc);
}

void Mjoin(PATL,tsyr2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   ATL_SYR2K_t sy;
   #ifdef TREAL
      const TYPE ONE = ATL_rone, ZERO = ATL_rzero;
   #else
      const TYPE ONE[2]={ATL_rone, ATL_rzero}, ZERO[2]={ATL_rzero,ATL_rzero};
   #endif

   if (N < 1)
      return;
   if (SCALAR_IS_ZERO(alpha) || K < 1)
   {
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,trscal)(Uplo, N, N, beta, C, ldc);
      return;
   }
   sy.alpha2 = sy.alpha = SADD alpha;
   sy.beta  = SADD beta;
   sy.one = SADD ONE;
   sy.zero = SADD ZERO;
   sy.tvgemm = Mjoin(PATL,tvgemm);
   sy.tvApAt = Mjoin(PATL,tvsyApAt);
   sy.K = K;
   sy.lda = lda;
   sy.ldb = ldb;
   sy.ldc = ldc;
   sy.eltsh = Mjoin(PATL,shift);
   sy.Uplo = Uplo;
   sy.trans = Trans;
   if (Trans == AtlasNoTrans)
   {
      sy.TA = AtlasNoTrans;
      sy.TB = AtlasTrans;
      sy.TA2 = AtlasTrans;
      sy.TB2 = AtlasNoTrans;
   }
   else
   {
      sy.TA = AtlasTrans;
      sy.TB = AtlasNoTrans;
      sy.TA2 = AtlasNoTrans;
      sy.TB2 = AtlasTrans;
   }
   sy.nb = Mjoin(PATL,GetNB)();
   ATL_tvsyr2k_rec(&sy, N/sy.nb, N%sy.nb, A, B, C);
}
