#include "atlas_misc.h"
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
#include <omp.h>
/*
 * prototype the typeless tGEMM helper routines
 */
void ATL_DoWorkMM(ATL_LAUNCHSTRUCT_t *lp, void *vp);
int ATL_StructIsInitMM(void *vp);
int ATL_thrdecompMM
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC);
void Mjoin(PATL,CombineStructsMM)(void *vme, void *vhim);
static void InitTMMNodes(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                         const TYPE *alpha, const TYPE *beta,
                         const TYPE *one, const TYPE *zero,
                         ATL_thread_t *btp, ATL_TMMNODE_t *ptmms)
{
   int i;
   void (*gemmK)(ATL_CINT, ATL_CINT, ATL_CINT, const void*, const void *,
                 ATL_CINT, const void*, ATL_CINT, const void*, void*, ATL_CINT);

   if (TA == AtlasNoTrans)
   {
#ifdef TCPLX
      if (TB == AtlasConjTrans)
         gemmK = Mjoin(PATL,tsvgemmNC);
      else
#endif
      gemmK = (TB == AtlasNoTrans)?Mjoin(PATL,tsvgemmNN):Mjoin(PATL,tsvgemmNT);
   }
#ifdef TCPLX
   else if (TA == AtlasConjTrans)
   {
      if (TB == AtlasNoTrans)
         gemmK = Mjoin(PATL,tsvgemmCN);
      else if (TB == AtlasConjTrans)
         gemmK = Mjoin(PATL,tsvgemmCC);
      else
         gemmK = Mjoin(PATL,tsvgemmCT);
   }
#endif
   else
   {
#ifdef TCPLX
      if (TB == AtlasConjTrans)
         gemmK = Mjoin(PATL,tsvgemmTC);
      else
#endif
      gemmK = (TB == AtlasNoTrans)?Mjoin(PATL,tsvgemmTN):Mjoin(PATL,tsvgemmTT);
   }
   for (i=0; i < ATL_NTHREADS; i++)
   {
      ptmms[i].mb = MB;
      ptmms[i].nb = NB;
      ptmms[i].kb = KB;
      ptmms[i].gemmK = gemmK;
      ptmms[i].eltsz = ATL_sizeof;
      ptmms[i].eltsh = Mjoin(PATL,shift);
      ptmms[i].K = 0;
      ptmms[i].nCw = 0;
      ptmms[i].rank = i;
      ptmms[i].alpha = (void*) alpha;
      ptmms[i].beta  = (void*) beta;
      ptmms[i].one = (void*) one;
      ptmms[i].zero  = (void*) zero;
      ptmms[i].Cinfp[0] = ptmms+i;
   }
}

void Mjoin(PATL,tompgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
                       const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   ATL_thread_t tp[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   ATL_TMMNODE_t mms[ATL_NTHREADS];
   int i, np;
   /*openmp needed variables*/
   int iam, abit, mask, src, dest;
   ATL_LAUNCHSTRUCT_t *lp;
   omp_lock_t omp_locks[ATL_NTHREADS]; /*used to lock threads for combine*/

/*
 *Orginal code from tgemm() unmodified
 */
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

/*
 * Don't thread if the compute/mem ratio is essentially like a Level 1 blas
 */
   i = Mmin(M,N);
   i = Mmin(i, K);
   if (i < 8)
      np = 1;
   else
   {
      InitTMMNodes(TA, TB, SADD alpha, SADD beta, SADD ONE, SADD ZERO, tp, mms);
      np = ATL_thrdecompMM(mms, TA, TB, M/MB, M%MB, N/NB, N%NB, K/KB, K%KB,
                           A, lda, B, ldb, C, ldc, ATL_NTHREADS, 0, 0);
   }
#ifdef DEBUG
fprintf(stderr, "np=%d\n\n", np);
#endif
   if (np < 2)
   {
      Mjoin(PATL,gemm)(TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      return;
   }
   ls.opstruct = (char*) mms;
   ls.opstructstride = (int) ( ((char*)(mms+1)) - (char*)mms );
   ls.OpStructIsInit = ATL_StructIsInitMM;
   ls.CombineOpStructs = Mjoin(PATL,CombineStructsMM);
   ls.DoWork = ATL_DoWorkMM;
   ls.rank2thr = tp;

/*
 *end unmodified code
 */

   /*
    *initialize openmp locks
    */
   for(i = 0; i < ATL_NTHREADS; i++)
   {
      omp_init_lock(&omp_locks[i]);
   }

   #pragma omp parallel private(lp,iam,mask,abit,i,src) shared(omp_locks)
   {
      /*
       *each thread intializes private variables
       */
      iam = omp_get_thread_num();
      omp_set_lock(&omp_locks[iam]);
      tp[iam].rank = iam;
      tp[iam].vp = &ls;
      lp = tp[iam].vp;
      /*
       *call the work routine
       */
      lp->DoWork(lp,lp->opstruct+lp->opstructstride*iam);

      /*
       *combines the answers in log2 steps majority of code taken from
       *log2launch
       */
      mask = 0;
      for( i = 0; i < ATL_NTHRPOW2; i++)
      {
         if(!(iam & mask))
         {
            abit = (1<<i);
            if(!(iam & abit))
            {
               src = iam ^ abit;
               if ( src < ATL_NTHREADS && (!lp->OpStructIsInit ||
                    lp->OpStructIsInit(lp->opstruct+lp->opstructstride*src)) )
               {
                  /*
                   *waits for the previous level to finish before combine
                   */
                  omp_set_lock(&omp_locks[src]);
                  lp->CombineOpStructs(lp->opstruct+lp->opstructstride*iam,
                                       lp->opstruct+lp->opstructstride*src);
                  omp_unset_lock(&omp_locks[src]);
               }
            }
            else
               /*
                *once finshed releases lock for current level
                */
               omp_unset_lock(&omp_locks[iam]);
         }
         mask |= abit;
      }
      omp_unset_lock(&omp_locks[iam]);  /*cleans up locks*/
   }
}
void Mjoin(PATL,tvompgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const void *alpha,
                       const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
                       const void *beta, void *C, ATL_CINT ldc)

/*
 *This void wrapper for tompgemm is used in some typeless structures
 */
{
   Mjoin(PATL,tompgemm)(TA, TB, M, N, K, SVVAL((const TYPE*)alpha), A, lda,
                        B, ldb, SVVAL((const TYPE*)beta), C, ldc);
}

