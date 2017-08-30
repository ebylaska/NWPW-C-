#include "atlas_misc.h"
#define ATL_LAUNCHORDER
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
#include "math.h"
#include "omp.h"

/*
 * Recursive decompisiton on trapazoidal-shaped matrix ($C$ after splitting)
 */
#ifndef ATL_MINL3THRFLOPS
   #ifdef ATL_TGEMM_ADDP
      #define ATL_MINL3THRFLOPS \
         (((2.0*ATL_TGEMM_ADDP)*ATL_TGEMM_ADDP)*ATL_TGEMM_ADDP)
   #else
      #define ATL_MINL3THRFLOPS (((2.0*MB)*NB)*KB)
   #endif
#endif

/*
 * Prototype functons in ATL_Xtsyrk
 */

int ATL_IsInitSYRK_K(void *vp);
int ATL_tsyrkdecomp_K
   (ATL_TSYRK_K_t *psyrk,
    void (*syrkK)(const enum ATLAS_UPLO, const enum ATLAS_TRANS, ATL_CINT,
                  ATL_CINT, const void*, const void*, ATL_CINT, const void*,
                  void*, ATL_CINT),
    int np, const int eltsh, const int nb, const void *zero, const void *one,
    const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT Kblks, const int kr,
    const void *alpha, const void *A, ATL_CINT lda,
    const void *beta, void *C, ATL_CINT ldc);

void ATL_tompsyrk_K(ATL_TSYRK_K_t *syp, int np, ATL_CINT N, ATL_CINT K,
                 const void *A, void *C)
{
   const int nb = syp->nb;
   int iam,i,abit,mask,src,dest;
   ATL_LAUNCHSTRUCT_t *lp;
   omp_lock_t omp_locks[ATL_NTHREADS];

   np = (Mmin(N,K) < 8) ? 1 :
        ATL_tsyrkdecomp_K(syp, syp->tvsyrk,np, syp->eltsh, nb, syp->zero,
                          syp->one, syp->Uplo, syp->Trans, N, K/nb, K%nb,
                          syp->alpha, A, syp->lda, syp->beta, C, syp->ldc);
   if (np < 2)
   {
      syp->tvsyrk(syp->Uplo, syp->Trans, N, K, syp->alpha, A, syp->lda,
                  syp->beta, C, syp->ldc);
      return;
   }
   for( i = 0; i < ATL_NTHREADS; i++ )
   {
      omp_init_lock(&omp_locks[i]);
   }
   #pragma omp parallel private(iam,src,lp,mask,abit,i)
   {
      iam=omp_get_thread_num();
       omp_set_lock(&omp_locks[iam]);
       lp = syp->lp;
       if(lp->OpStructIsInit(lp->opstruct+lp->opstructstride*iam))
          lp->DoWork(lp,lp->opstruct+lp->opstructstride*iam);
       mask = 0;
       for( i = 0; i < ATL_NTHRPOW2; i++)
       {
          if(!(iam & mask))
          {
             abit = (1<<i);
             if(!(iam & abit))
             {
                src = iam ^ abit;
                if( src < ATL_NTHREADS && (!lp->OpStructIsInit ||
                    lp->OpStructIsInit(lp->opstruct+lp->opstructstride*src)) )
                {
                   omp_set_lock(&omp_locks[src]);
                   lp->CombineOpStructs(lp->opstruct+lp->opstructstride*iam,
                                        lp->opstruct+lp->opstructstride*src);
                    omp_unset_lock(&omp_locks[src]);
                 }
              }
              else
                 omp_unset_lock(&omp_locks[iam]);
           }
           mask |= abit;
        }
       // omp_unset_lock(&omp_locks[iam]);
   }
}

void ATL_tompsyrk_K_rec(ATL_TSYRK_K_t *syp,int np, ATL_CINT Nblks, ATL_CINT nr,
                     ATL_CINT K, const void *A0, void *C00)
/*
 * This routine recurs on N until we can allocate the full NxN workspace,
 * at which point it stops the recursion and distributes K for parallel
 * operation
 */
{
   const enum ATLAS_TRANS TA = syp->Trans;
   ATL_CINT lda = syp->lda, ldc = syp->ldc, eltsh = syp->eltsh;
   ATL_CINT nb = syp->nb, N = Nblks*nb+nr;
   ATL_INT sz, nblksL, nblksR, nrL, nrR, nL, nR;
   const void *A1;
   void *C10, *C01, *C11;
/*
 * Stop recursion & call threaded SYRK if we can allocate workspace for all of C
 */
   sz = (N * N) << eltsh;
/*
 * Quit recurring if we can allocate space for C workspace and we can
 * no longer usefully split Nblks, or we can usefully split K
 */
   if (sz <= ATL_PTMAXMALLOC && (nb*ATL_NTHREADS < K || Nblks < ATL_NTHREADS))
   {
      ATL_tompsyrk_K(syp, np, Nblks*nb+nr, K, A0, C00);
      return;
   }
   nblksL = (Nblks+1)>>1;
   nblksR = Nblks - nblksL;
   if (nblksL >= nblksR)
   {
      nrL = nr;
      nrR = 0;
   }
   else
   {
      nrL = 0;
      nrR = nr;
   }

   nL = nblksL * nb + nrL;
   nR = nblksR * nb + nrR;
   if (syp->Uplo == AtlasUpper)
   {
      sz = nL<<eltsh;
      C01 = MindxT(C00,sz*ldc);
      A1 = (TA == AtlasNoTrans) ? MindxT(A0,sz) : MindxT(A0,sz*lda);
      C11 = MindxT(C01,sz);
      ATL_tompsyrk_K_rec(syp, np, nblksL, nrL, K, A0, C00);
      syp->gemmT(syp->Trans, syp->TB, nL, nR, K, syp->alpha, A0, lda, A1, lda,
                 syp->beta, C01, ldc);
      ATL_tompsyrk_K_rec(syp, np , nblksR, nrR, K, A1, C11);
   }
   else /* Lower triangular matrix */
   {
      sz = nL<<eltsh;
      C10 = MindxT(C00,sz);
      A1 = (TA == AtlasNoTrans) ? MindxT(A0,sz) : MindxT(A0,sz*lda);
      sz += (ldc*nL)<<eltsh;
      C11 = MindxT(C00,sz);
      ATL_tompsyrk_K_rec(syp, np, nblksL, nrL, K, A0, C00);
      syp->gemmT(syp->Trans, syp->TB, nR, nL, K, syp->alpha, A1, lda, A0, lda,
                 syp->beta, C10, ldc);
      ATL_tompsyrk_K_rec(syp, np, nblksR, nrR, K, A1, C11);
   }
}
