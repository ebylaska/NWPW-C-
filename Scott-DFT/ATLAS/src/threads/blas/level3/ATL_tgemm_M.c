#include "atlas_misc.h"
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
/*
 * prototype the typeless tGEMM helper routines
 */
void ATL_DoWorkMM(ATL_LAUNCHSTRUCT_t *lp, void *vp);
int ATL_StructIsInitMM(void *vp);
int ATL_thrdecompMM_K
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC);
int ATL_thrdecompMM_N
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC);
int ATL_thrdecompMM_M
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC);
void Mjoin(PATL,InitTMMNodes)
   (const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB, const TYPE *alpha,
    const TYPE *beta, const TYPE *one, const TYPE *zero,
    ATL_thread_t *btp, ATL_TMMNODE_t *ptmms);
int Mjoin(PATL,tgemm_M)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
                       const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   ATL_thread_t tp[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   ATL_TMMNODE_t mms[ATL_NTHREADS];
   int i, np, DividedK=0;
   #ifdef TREAL
      TYPE ONE=ATL_rone, ZERO=ATL_rzero;
   #else
      TYPE ONE[2] = {ATL_rone, ATL_rzero}, ZERO[2] = {ATL_rzero, ATL_rzero};
   #endif

   if (M < 1 || N < 1)
      return;
   if (K < 1 || SCALAR_IS_ZERO(alpha))
   {
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return;
   }
   Mjoin(PATL,InitTMMNodes)(TA, TB, SADD alpha, SADD beta, SADD ONE,
                            SADD ZERO, tp, mms);
   np = ATL_thrdecompMM_M(mms, TA, TB, M/MB, M%MB, N/NB, N%NB, K/KB, K%KB,
                          A, lda, B, ldb, C, ldc, ATL_NTHREADS, 0, 0);
#ifdef DEBUG
fprintf(stderr, "np=%d\n\n", np);
#endif
   if (np < 2)
   {
      Mjoin(PATL,gemm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return(1);
   }
   ls.opstruct = (char*) mms;
   ls.opstructstride = (int) ( ((char*)(mms+1)) - (char*)mms );
   ls.OpStructIsInit = ATL_StructIsInitMM;
   ls.CombineOpStructs = NULL;
   ls.DoWork = ATL_DoWorkMM;
   ls.rank2thr = tp;
   for (i=0; i < ATL_NTHREADS; i++)
   {
      tp[i].vp = &ls;
      tp[i].rank = i;
   }
   ATL_thread_start(tp, 0, ATL_log2tlaunch, tp);
   ATL_thread_join(tp);
   return(np);
}
