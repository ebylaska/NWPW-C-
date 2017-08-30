#include "atlas_misc.h"
#define ATL_LAUNCHORDER         /* we want static ATL_launchorder array */
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
#include "atlas_omplvl3.h"

/*
 *prototypes of functions found in ATL_ttrmm.c and ATL_ttrsm.c
 */
int Mjoin(PATL,StructIsInitTRMM)(void *vp);
void Mjoin(PATL,DoWorkTRMM)(ATL_LAUNCHSTRUCT_t *lp, void *vp);

#ifndef ATL_TTRSM_XOVER
   #define ATL_TTRSM_XOVER 4   /* want 4 total blocks before adding proc */
#endif
void Mjoin(PATL,tomptrmm)
   (const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
    const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb)
{
   ATL_thread_t tp[ATL_NTHREADS];
   ATL_TTRSM_t trsms[ATL_NTHREADS];
   ATL_LAUNCHSTRUCT_t ls;
   TYPE *b;
   ATL_INT n, nblks, minblks;
   double tblks;
   int nr, p, i, j, extrablks;
   static int nb=0;

   int iam;
   ATL_LAUNCHSTRUCT_t *lp;

   if (M < 1 || N < 1)
      return;
   if (SCALAR_IS_ZERO(alpha))
   {
      Mjoin(PATL,gezero)(M, N, B, ldb);
      return;
   }
/*
 * Distribute RHS over the processors
 */
   if (!nb) nb = Mjoin(PATL,GetNB)();
   if (side == AtlasLeft)
   {
      nblks = N/nb;
      nr = N - nblks*nb;
      tblks = ((double)(M*N)) / ( (double)nb * nb );
      p = (tblks+ATL_TTRSM_XOVER-1)/ATL_TTRSM_XOVER;
      p = Mmin(p, ATL_NTHREADS);
      p = p ? p : 1;

      b = B;
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
         trsms[j].A = A;
         trsms[j].M = M;
         trsms[j].N = n;
         trsms[j].lda = lda;
         trsms[j].ldb = ldb;
         trsms[j].B = b;
         trsms[j].alpha = SADD alpha;
         trsms[j].side = side;
         trsms[j].uplo = uplo;
         trsms[j].TA   = TA;
         trsms[j].diag = diag;
         n *= (ldb << Mjoin(PATL,shift));
         b = MindxT(b, n);
      }
      for (; i < ATL_NTHREADS; i++)  /* flag rest of struct as uninitialized */
         trsms[ATL_launchorder[i]].B = NULL;
   }
   else /* Side == AtlasRight */
   {
      nblks = M/nb;
      nr = M - nblks*nb;
      tblks = (N/nb)*nblks;
      p = (tblks+ATL_TTRSM_XOVER-1)/ATL_TTRSM_XOVER;
      p = Mmin(p, ATL_NTHREADS);
      p = p ? p : 1;

      b = B;
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
         trsms[j].A = A;
         trsms[j].M = n;
         trsms[j].N = N;
         trsms[j].lda = lda;
         trsms[j].ldb = ldb;
         trsms[j].B = b;
         trsms[j].alpha = SADD alpha;
         trsms[j].side = side;
         trsms[j].uplo = uplo;
         trsms[j].TA   = TA;
         trsms[j].diag = diag;
         n <<= Mjoin(PATL,shift);
         b = MindxT(b, n);
      }
   }
   if (p < 2)
   {
      Mjoin(PATL,trmm)(side, uplo, TA, diag, M, N, alpha, A, lda, B, ldb);
      return;
   }
   for (; i < ATL_NTHREADS; i++)  /* flag rest of struct as uninitialized */
      trsms[ATL_launchorder[i]].B = NULL;
   ls.opstruct = (char*) trsms;
   ls.opstructstride = (int) ( ((char*)(trsms+1)) - (char*)(trsms) );
   ls.CombineOpStructs = NULL;
   ls.OpStructIsInit = Mjoin(PATL,StructIsInitTRMM);
   ls.DoWork = Mjoin(PATL,DoWorkTRMM);
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
}

