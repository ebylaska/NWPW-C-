#include "atlas_misc.h"
#define ATL_LAUNCHORDER         /* we want static ATL_launchorder array */
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
#include <omp.h>

/*
 *protypes for functions found in ATL_symm.c and ATL_hymm.c
 */
 void Mjoin(PATL,DoWorkHEMM)(ATL_LAUNCHSTRUCT_t *lp, void *vp);
 int Mjoin(PATL,StructIsInitHEMM)(void *vp);

static void ATL_omphemmL_rec
   (ATL_TSYMM_t *syp, ATL_CINT Mblks, ATL_CINT mr, ATL_CINT Nblks, ATL_CINT nr,
    const TYPE *A, const TYPE *B, TYPE *C)
{
   const TYPE *A10, *A01, *B10;
   TYPE *C10;
   const int nb = syp->nb;
   ATL_INT nbR, nbL, rR, rL, nL, nR;
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone,ATL_rzero};
   #else
      TYPE ONE = ATL_rone;
   #endif

   nbR = Mblks>>1;
   nbL = Mblks - nbR;
/*
 * Stop recursion once we are no longer getting parallelism
 */
   if (nbR*Nblks < ATL_TGEMM_XOVER)
   {
      Mjoin(PATL,hemm)(syp->side, syp->uplo, Mblks*nb+mr, syp->N,
                       SVVAL((TYPE*)syp->alpha), A, syp->lda, B, syp->ldb,
                       SVVAL((TYPE*)syp->beta), C, syp->ldc);
      return;
   }
   rL = (nbR == nbL) ? mr : 0;
   rR = mr - rL;
   nL = nbL*nb + rL;
   nR = nbR*nb + rR;
   B10 = B + (nL SHIFT);
   C10 = C + (nL SHIFT);
   ATL_omphemmL_rec(syp, nbL, rL, Nblks, nr, A, B, C);
   ATL_omphemmL_rec(syp, nbR, rR, Nblks, nr, A+(syp->lda+1)*(nL SHIFT), B10, C10);
   if (syp->uplo == AtlasLower)
   {
      A10 = A + (nL SHIFT);
      Mjoin(PATL,tompgemm)(AtlasConjTrans, AtlasNoTrans, nL, syp->N, nR,
                        SVVAL((TYPE*)syp->alpha), A10, syp->lda, B10, syp->ldb,
                        ONE, C, syp->ldc);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasNoTrans, nR, syp->N, nL,
                        SVVAL((TYPE*)syp->alpha), A10, syp->lda, B, syp->ldb,
                        ONE, C10, syp->ldc);
   }
   else
   {
      A01 = A + (syp->lda SHIFT);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasNoTrans, nL, syp->N, nR,
                        SVVAL((TYPE*)syp->alpha), A01, syp->lda, B10, syp->ldb,
                        ONE, C, syp->ldc);
      Mjoin(PATL,tompgemm)(AtlasConjTrans, AtlasNoTrans, nR, syp->N, nL,
                        SVVAL((TYPE*)syp->alpha), A01, syp->lda, B, syp->ldb,
                        ONE, C10, syp->ldc);
   }
}
static void ATL_omphemmR_rec
   (ATL_TSYMM_t *syp, ATL_CINT Mblks, ATL_CINT mr, ATL_CINT Nblks, ATL_CINT nr,
    const TYPE *A, const TYPE *B, TYPE *C)
{
   const TYPE *A10, *A01, *B01;
   TYPE *C01;
   const int nb = syp->nb;
   ATL_INT nbR, nbL, rR, rL, nL, nR;
   #ifdef TCPLX
      TYPE ONE[2] = {ATL_rone,ATL_rzero};
   #else
      TYPE ONE = ATL_rone;
   #endif

   nbR = Nblks>>1;
   nbL = Nblks - nbR;
/*
 * Stop recursion once we are no longer getting parallelism
 */
   if (nbR*Mblks < ATL_TGEMM_XOVER)
   {
      Mjoin(PATL,hemm)(syp->side, syp->uplo, syp->M, Nblks*nb+nr,
                       SVVAL((TYPE*)syp->alpha), A, syp->lda, B, syp->ldb,
                       SVVAL((TYPE*)syp->beta), C, syp->ldc);
      return;
   }
   rL = (nbR == nbL) ? nr : 0;
   rR = nr - rL;
   nL = nbL*nb + rL;
   nR = nbR*nb + rR;
   B01 = B + (nL*syp->ldb SHIFT);
   C01 = C + (nL*syp->ldc SHIFT);
   ATL_omphemmL_rec(syp, Mblks, mr, nbL, rL, A, B, C);
   ATL_omphemmL_rec(syp, Mblks, mr, nbR, rR, A+(syp->lda+1)*(nL SHIFT), B01, C01);
   if (syp->uplo == AtlasLower)
   {
      A10 = A + (nL SHIFT);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasNoTrans, syp->M, nL, nR,
                        SVVAL((TYPE*)syp->alpha), B01, syp->ldb, A10, syp->lda,
                        ONE, C, syp->ldc);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasConjTrans, syp->M, nR, nL,
                        SVVAL((TYPE*)syp->alpha), B, syp->ldb, A10, syp->lda,
                        ONE, C01, syp->ldc);
   }
   else
   {
      A01 = A + (syp->lda SHIFT);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasConjTrans, syp->M, nL, nR,
                        SVVAL((TYPE*)syp->alpha), B01, syp->ldb, A01, syp->lda,
                        ONE, C, syp->ldc);
      Mjoin(PATL,tompgemm)(AtlasNoTrans, AtlasNoTrans, syp->M, nR, nL,
                        SVVAL((TYPE*)syp->alpha), B, syp->ldb, A01, syp->lda,
                        ONE, C01, syp->ldc);
   }
}

static void ATL_tomphemm_SYsplit
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc,
    ATL_CINT nb)
/*
 * This routine is specialized for the case where we cannnot split the
 * B matrix, and must instead split the symmetric matrix (A).  It calls
 * a recursive GEMM-based BLAS, that gets its parallel performance from
 * calling threaded GEMM.
 */
{
   ATL_TSYMM_t ss;
   ss.side = Side;
   ss.uplo = Uplo;
   ss.M = M;
   ss.N = N;
   ss.nb = nb;
   ss.alpha = SADD alpha;
   ss.beta  = SADD beta;
   ss.lda = lda;
   ss.ldb = ldb;
   ss.ldc = ldc;
   if (Side == AtlasLeft)
      ATL_omphemmL_rec(&ss, M/nb, M%nb, N/nb, N%nb, A, B, C);
   else
      ATL_omphemmR_rec(&ss, M/nb, M%nb, N/nb, N%nb, A, B, C);

}

/*
 * The XOVER is the min # of blocks required to do parallel operation
 */
#ifndef ATL_TSYMM_XOVER
   #ifdef ATL_TGEMM_XOVER
      #define ATL_TSYMM_XOVER ATL_TGEMM_XOVER
   #else
      #define ATL_TSYMM_XOVER 4  /* want 4 blocks for parallel execution */
   #endif
#endif
/*
 * Once you have achieved enough blocks to do computation, minimum number
 * of blocks to give each processor
 */
#ifndef ATL_TSYMM_ADDP
   #ifdef ATL_TGEMM_ADDP
      #define ATL_TSYMM_ADDP  ATL_TGEMM_ADDP
   #else
      #define ATL_TSYMM_ADDP  1  /* want 1 blocks to add proc to workers */
   #endif
#endif
void Mjoin(PATL,tomphemm)
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc)
{
   ATL_INT n, nblks, tblks, nr, minblks, extrablks, p, i, j;
   ATL_thread_t tp[ATL_NTHREADS];
   ATL_TSYMM_t symms[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   const TYPE *b;
   TYPE *c;
   static int nb=0;

   int iam;
   ATL_LAUNCHSTRUCT_t *lp;

   if (M < 1 || N < 1)
      return;
   if (SCALAR_IS_ZERO(alpha))
   {
      if (!SCALAR_IS_ONE(beta))
         Mjoin(PATL,gescal)(M, N, beta, C, ldc);
      return;
   }
   if (!nb) nb = Mjoin(PATL,GetNB());
   if (Side == AtlasLeft)
   {
      nblks = N / nb;
      nr = N - nblks*nb;
      tblks = ((double)(M*N)) / ( (double)nb * nb );
      p = (nblks+ATL_TSYMM_ADDP-1)/ATL_TSYMM_ADDP;
      if (p < ATL_NTHREADS)  /* N not big enough to give blk to each proc */
      {
/*
 *       If I can't split N, and M is the dominant cost, use recursion to
 *       decompose symmetric matrix; parallelism will come from TGEMM calls
 */
         if (M > (N<<(ATL_NTHRPOW2+2)))
         {
            ATL_tomphemm_SYsplit(Side, Uplo, M, N, alpha, A, lda, B, ldb,
                              beta, C, ldc, nb);
            return;
         }
      }
      else
         p = ATL_NTHREADS;

      if (p < 2)
         goto SERIAL;
/*
 *    Distribute N over the processors
 */
      b = B;
      c = C;
      minblks = nblks / p;
      extrablks = nblks - minblks*p;
      for (i=0; i < p; i++)
      {
         if (i < extrablks)
            n = (minblks+1)*nb;
         else if (i == extrablks)
            n = minblks*nb + nr;
         else
            n = minblks*nb;
         j = ATL_launchorder[i];
         symms[j].A = A;
         symms[j].B = b;
         symms[j].alpha = SADD alpha;
         symms[j].beta = SADD beta;
         symms[j].C = c;
         symms[j].M = M;
         symms[j].N = n;
         symms[j].lda = lda;
         symms[j].ldb = ldb;
         symms[j].ldc = ldc;
         symms[j].side = Side;
         symms[j].uplo = Uplo;
         b = MindxT(b, ATL_MulBySize(ldb)*n);
         c = MindxT(c, ATL_MulBySize(ldc)*n);
      }
      for (; i < ATL_NTHREADS; i++)  /* flag rest of struct as uninitialized */
         symms[ATL_launchorder[i]].M = 0;
   }
   else  /* Side == AtlasRight */
   {
      nblks = M / nb;
      nr = M - nblks*nb;
      tblks = ((double)(M*N)) / ( (double)nb * nb );
      p = (nblks+ATL_TSYMM_ADDP-1)/ATL_TSYMM_ADDP;
      if (p < ATL_NTHREADS)  /* N not big enough to give blk to each proc */
      {
/*
 *       If I can't split M, and N is the dominant cost, use recursion to
 *       decompose symmetric matrix; parallelism will come from TGEMM calls
 */
         if (N > (M<<(ATL_NTHRPOW2+2)))
         {
            ATL_tomphemm_SYsplit(Side, Uplo, M, N, alpha, A, lda, B, ldb,
                              beta, C, ldc, nb);
            return;
         }
      }
      else
         p = ATL_NTHREADS;
      if (p < 2)
         goto SERIAL;
/*
 *    Distribute M over the processors
 */
      b = B;
      c = C;
      minblks = nblks / p;
      extrablks = nblks - minblks*p;
      for (i=0; i < p; i++)
      {
         if (i < extrablks)
            n = (minblks+1)*nb;
         else if (i == extrablks)
            n = minblks*nb + nr;
         else
            n = minblks*nb;
         j = ATL_launchorder[i];
         symms[j].A = A;
         symms[j].B = b;
         symms[j].alpha = SADD alpha;
         symms[j].beta = SADD beta;
         symms[j].C = c;
         symms[j].M = n;
         symms[j].N = N;
         symms[j].lda = lda;
         symms[j].ldb = ldb;
         symms[j].ldc = ldc;
         symms[j].side = Side;
         symms[j].uplo = Uplo;
         b = MindxT(b, ATL_MulBySize(n));
         c = MindxT(c, ATL_MulBySize(n));
      }
      for (; i < ATL_NTHREADS; i++)  /* flag rest of struct as uninitialized */
         symms[ATL_launchorder[i]].M = 0;
   }
   if (p < 2)
   {
SERIAL:
      Mjoin(PATL,hemm)(Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
      return;
   }
   ls.opstruct = (char*) symms;
   ls.opstructstride = (int) ( ((char*)(symms+1)) - (char*)(symms) );
   ls.CombineOpStructs = NULL;
   ls.OpStructIsInit = Mjoin(PATL,StructIsInitHEMM);
   ls.DoWork = Mjoin(PATL,DoWorkHEMM);
   ls.rank2thr = tp;
   /*
    *modifed from the orginal code to thread using omp instead of pthreads
    */
   #pragma omp parallel private(lp,iam)
   {
      iam = omp_get_thread_num();
      tp[iam].rank = iam;
      tp[iam].vp = &ls;
      lp = tp[iam].vp;

      if(lp->OpStructIsInit(lp->opstruct+lp->opstructstride*iam))
	      lp->DoWork(lp,lp->opstruct+lp->opstructstride*iam);
   }
}

