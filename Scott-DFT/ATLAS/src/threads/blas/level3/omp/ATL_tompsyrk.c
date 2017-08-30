#include "atlas_misc.h"
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
#include "atlas_omplvl3.h"

/*
 * Prototype functions in ATL_Xtsyrk
 */
int ATL_IsInitSYRK_M(void *vp);
void ATL_DoWorkSYRK_M(ATL_LAUNCHSTRUCT_t *lp, void *vp);
int ATL_tsyrkdecomp_M
   (ATL_TSYRK_M_t *syp, const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA,
    ATL_CINT N, ATL_CINT K, const void *alpha, const void *A, ATL_CINT lda,
    const void *beta, void *C, ATL_CINT ldc, ATL_CINT nb, const int mu,
    const int eltsh, const enum ATLAS_TRANS TB, double minmf,
    void (*gemmK)(ATL_CINT, ATL_CINT, ATL_CINT, const void*, const void *,
                  ATL_CINT,const void*, ATL_CINT, const void*, void*, ATL_CINT),
    void (*tvsyrk)(const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CINT,
                   ATL_CINT, const void*, const void*, ATL_CINT, const void*,
                   void*, ATL_CINT));

int ATL_tsyrkdecomp_K
   (ATL_TSYRK_K_t *psyrk,
    void (*syrkK)(const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CINT,
                  ATL_CINT, const void*, const void*, ATL_CINT, const void*,
                  void*, ATL_CINT),
    const int eltsh, const int nb, const void *zero, const void *one,
    const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT Kblks, const int kr,
    const void *alpha, const void *A, ATL_CINT lda,
    const void *beta, void *C, ATL_CINT ldc);
void ATL_DoWorkSYRK_K(ATL_LAUNCHSTRUCT_t *lp, void *vp);
int ATL_IsInitSYRK_K(void *vp);

/*
 * Prototype functions in ATL_Xtompsyrk
 */
void ATL_tompsyrk_K_rec(ATL_TSYRK_K_t *syp, int np,  ATL_CINT Nblks,
                        ATL_CINT nr, ATL_CINT K, const void *A, void *C);

/*
 *prototype function for ATL_tsyrk.c
 */
void Mjoin(PATL,CombineStructsSYRK)(void *vme, void *vhim);
void Mjoin(PATL,tvsyrk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const void *alpha, const void *A, ATL_CINT lda,
    const void *beta, void *C, ATL_CINT ldc);

void Mjoin(PATL,tompsyrk_K_rec)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc, ATL_CINT nb)
/*
 * This typed wrapper routine sets up type-specific data structures, and
 * calls the appropriate typeless recursive routine in order to recursively
 * cut N until workspace can be allocated, and then the K-dimension will be
 * threaded.  During the recursion, parallel performance is acheived by
 * calling the threaded GEMM.
 */
{
   ATL_CINT Nblks = N/nb, nr = N - nb*Nblks;
   ATL_TSYRK_K_t syp[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   ATL_thread_t tp[ATL_NTHREADS];
   #ifdef TCPLX
      TYPE ZERO[2] = {ATL_rzero, ATL_rzero}, ONE[2] = {ATL_rone, ATL_rzero};
   #else
      TYPE ZERO=ATL_rzero, ONE=ATL_rone;
   #endif
   int i;

   ls.opstructstride = (int) ( ((char*)(syp+1)) - (char*)syp );
   ls.OpStructIsInit = ATL_IsInitSYRK_K;
   ls.DoWork = ATL_DoWorkSYRK_K;
   ls.CombineOpStructs = Mjoin(PATL,CombineStructsSYRK);
   ls.rank2thr = tp;
   for (i=0; i < ATL_NTHREADS; i++)
   {
      tp[i].vp = &ls;
      tp[i].rank = i;
   }
   syp[0].lp = &ls;
   syp[0].Uplo = Uplo;
   syp[0].Trans = Trans;
   syp[0].TB = (Trans == AtlasNoTrans) ? AtlasTrans : AtlasNoTrans;
   syp[0].K = K;
   syp[0].alpha = SADD alpha;
   syp[0].beta = SADD beta;
   syp[0].zero = SADD ZERO;
   syp[0].one  = SADD ONE;
   syp[0].lda = lda;
   syp[0].ldc = ldc;
   syp[0].gemmT = Mjoin(PATL,tvompgemm);
   syp[0].tvsyrk = Mjoin(PATL,tvsyrk);
   syp[0].eltsh = Mjoin(PATL,shift);
   syp[0].nb = nb;
   ls.opstruct = (char*) syp;
   ATL_tompsyrk_K_rec(syp,Mjoin(PATL,tNumGemmThreads)(N, N>>1, K), Nblks, nr,
                      K, A, C);
}

static int ATL_tompsyrk_M
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS TA, ATL_CINT N,
    ATL_CINT K, const void *alpha, const TYPE *A, ATL_CINT lda,
    const void *beta, TYPE *C, ATL_CINT ldc)
{
   ATL_thread_t tp[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   ATL_TSYRK_M_t syp[ATL_NTHREADS];
   int i, p;
   int iam;
   ATL_LAUNCHSTRUCT_t *lp;
   p = ATL_tsyrkdecomp_M(syp, Uplo, TA, N, K, alpha, A, lda, beta, C, ldc,
                         MB, ATL_mmMU, Mjoin(PATL,shift),
                         (TA == AtlasNoTrans) ? AtlasTrans : AtlasNoTrans,
                         ATL_TGEMM_PERTHR_MF, (TA == AtlasNoTrans) ?
                         Mjoin(PATL,tsvgemmNT):Mjoin(PATL,tsvgemmTN),
                         Mjoin(PATL,tvsyrk));
   if (p < 2)
      return(0);
   ls.opstruct = (char*) syp;
   ls.opstructstride = (int) ( ((char*)(syp+1)) - (char*)syp );
   ls.OpStructIsInit = ATL_IsInitSYRK_M;
   ls.DoWork = ATL_DoWorkSYRK_M;
   ls.CombineOpStructs = NULL;
   ls.rank2thr = tp;
   #pragma omp parallel private(lp,iam)
   {
      iam = omp_get_thread_num();
      tp[iam].rank = iam;
      tp[iam].vp = &ls;
      lp = tp[iam].vp;
      if(lp->OpStructIsInit(lp->opstruct+lp->opstructstride*iam))
         lp->DoWork(lp,lp->opstruct+lp->opstructstride*iam);
   }
   return(p);
}

void Mjoin(PATL,tompsyrk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   #ifdef TREAL
      const TYPE ONE = ATL_rone, ZERO = ATL_rzero;
   #else
      const TYPE ONE[2]={ATL_rone, ATL_rzero}, ZERO[2]={ATL_rzero, ATL_rzero};
   #endif
   size_t nblksN;
   int i, np, nb;
   void Mjoin(PATL,ptsyrk)
      (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
       ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
       const SCALAR beta, TYPE *C, ATL_CINT ldc);

   if (!Mjoin(PATL,tNumGemmThreads)(N, N>>1, K))
      goto DOSERIAL;
   if (N < 1)
      return;
   if (SCALAR_IS_ZERO(alpha) || K < 1)
   {
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,trscal)(Uplo, N, N, beta, C, ldc);
      return;
   }

   nb = MB;
   if (K > (N<<ATL_NTHRPOW2) && (((size_t)N)*N*sizeof(TYPE) > ATL_PTMAXMALLOC))
   {
      Mjoin(PATL,tompsyrk_K_rec)(Uplo, Trans, N, K, alpha, A, lda, beta,
                               C, ldc, nb);
      return;
   }
   np = ATL_tompsyrk_M(Uplo, Trans, N, K, SADD alpha, A, lda,
                     SADD beta, C, ldc);
   if (np < 2)
   {
DOSERIAL:
      Mjoin(PATL,syrk)(Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
      return;
   }
}

