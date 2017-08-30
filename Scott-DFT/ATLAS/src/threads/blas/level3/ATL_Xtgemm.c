#include "atlas_misc.h"
#define ATL_LAUNCHORDER         /* we want static ATL_launchorder array */
#include "atlas_threads.h"
#include "atlas_tlvl3.h"
/*
 * =========================================================================
 * This file contains support routines for TGEMM that are not type-dependent
 * =========================================================================
 */
int ATL_thrdecompMM_rMN
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC)
/*
 * This routine recursively splits the M & N dimensions over P processors
 */
{
   int pR, pL, rL, rR, np, eltsh;
   ATL_INT nblksL, nblksR, i, j;
   double d;

/*
 * Choose to split either M or N.  Want M < N always, so require
 * N to be twice as big before splitting (or be out of M blocks)
 */
   if (P > 1 && Nblks > 1 && (Mblks < 2 || Nblks >= Mblks+Mblks))
   {
      eltsh = ptmms[indx].eltsh;
      pR = P>>1;    /* on right, take P/2 threads */
      pL = P - pR;  /* on left, take remaining threads */
      d = (pR == pL) ? 0.5 : ((double)pL)/((double)P);    /* percent on left */
      #ifdef DEBUG
         fprintf(stderr, "Cut N\n");
      #endif
      nblksL = (d * Nblks);
      nblksR = Nblks - nblksL;
      if (nblksR < nblksL)
      {
         rL = 0;
         rR = nr;
      }
      else
      {
         rL = nr;
         rR = 0;
      }
      i = (nblksL*ptmms[indx].nb+rL) << eltsh;
      np = ATL_thrdecompMM_rMN(ptmms, TA, TB, Mblks, mr, nblksL, rL, Kblks, kr,
                               A, lda, B, ldb, C, ldc, pL, indx, COPYC);
      np += ATL_thrdecompMM_rMN(ptmms, TA, TB, Mblks, mr, nblksR, rR, Kblks, kr,
                                A,  lda, (TB == AtlasNoTrans) ?
                                MindxT(B,i*ldb) : MindxT(B,i), ldb,
                                MindxT(C,i*ldc), ldc, pR, indx+pL, COPYC);
      return(np);
   }
/*
 * If we have failed to split N, split M if possible
 */
   if (P > 1 && Mblks > 1)
   {
      eltsh = ptmms[indx].eltsh;
      pR = P>>1;    /* on right, take P/2 threads */
      pL = P - pR;  /* on left, take remaining threads */
      d = (pR == pL) ? 0.5 : ((double)pL)/((double)P);    /* percent on left */
      #ifdef DEBUG
         fprintf(stderr, "Cut M\n");
      #endif
      nblksL = (d * Mblks);
      nblksR = Mblks - nblksL;
      if (nblksR < nblksL)
      {
         rL = 0;
         rR = mr;
      }
      else
      {
         rL = mr;
         rR = 0;
      }
      i = (nblksL*ptmms[indx].mb+rL) << eltsh;
      np = ATL_thrdecompMM_rMN(ptmms, TA, TB, nblksL, rL, Nblks, nr, Kblks, kr,
                               A, lda, B, ldb, C, ldc, pL, indx, COPYC);
      np += ATL_thrdecompMM_rMN(ptmms, TA, TB, nblksR, rR, Nblks, nr, Kblks, kr,
                                (TA==AtlasNoTrans)?MindxT(A,i):MindxT(A,i*lda),
                                lda, B, ldb, MindxT(C,i), ldc, pR, indx+pL,
                                COPYC);
      return(np);
   }
/*
 * If no desirable splitting possible, stop recursion
 */
   ptmms[indx].A = A;
   ptmms[indx].B = B;
   ptmms[indx].C = (void*)C;
   ptmms[indx].lda = lda;
   ptmms[indx].ldb = ldb;
   ptmms[indx].ldc = ldc;
   ptmms[indx].M = ptmms[indx].mb*Mblks + mr;
   ptmms[indx].N = ptmms[indx].nb*Nblks + nr;
   ptmms[indx].K = ptmms[indx].kb*Kblks + kr;
   ptmms[indx].ldcw = ptmms[indx].nCw = 0;
   ptmms[indx].nCp = ptmms[indx].ownC = 1;
   ptmms[indx].Cinfp[ATL_NTHREADS-1] = ptmms+indx;
   ptmms[indx].Cw = NULL;
   #ifdef DEBUG
      fprintf(stderr, "%d: M=%d, N=%d, K=%d, ownC=%d, nCp=%d, nCw=%d\n",
              indx, ptmms[indx].M, ptmms[indx].N, ptmms[indx].K,
              ptmms[indx].ownC, ptmms[indx].nCp, ptmms[indx].nCw);
   #endif
   return(1);
}

int ATL_thrdecompMM_rMNK
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC)
/*
 * This routine decomposes the GEMM over P processors by splitting any of
 * the dimensions.  We only call this routine when K is very large or
 * M and N are very small (and thus splitting K, with its associated
 * extra workspace and flops, makes sense).
 */
{
   int pR, pL, rL, rR, np, eltsh;
   ATL_INT nblksL, nblksR, i, j;
   double d;

   eltsh = ptmms[indx].eltsh;
#ifdef DEBUG
   ATL_assert(P > 0);
#endif
   if (P <= 1 || (Mblks <= 1 && Nblks <= 1 && Kblks <= 1))
      goto STOP_REC;
   pR = P>>1;    /* on right, take P/2 threads */
   pL = P - pR;  /* on left, take remaining threads */
   d = (pR == pL) ? 0.5 : ((double)pL)/((double)P);    /* percent on left */
/*
 * Do not consider cutting K unless we have some K blocks, and we have either
 * already done so, or if we are sure that we are within our workspace limit
 */
   if (Kblks > 1 && (COPYC ||
       ((Mblks*ptmms[indx].mb+mr) * (Nblks*ptmms[indx].nb+nr) *
        (ATL_NTHREADS<<eltsh) < ATL_PTMAXMALLOC)))
   {
/*
 *    Before splitting K, ask that we are out of M and N blocks, or that
 *    our K is 4 times M and twice N
 */
      if ( (Mblks < 2 && Nblks < 2) ||
           (Kblks > (Mblks<<2) && Kblks > (Nblks+Nblks)) )
      {
         #ifdef DEBUG
            fprintf(stderr, "Cut K\n");
         #endif
         nblksL = (d * Kblks);
         nblksR = Kblks - nblksL;
         if (nblksR < nblksL)
         {
            rL = 0;
            rR = kr;
         }
         else
         {
            rL = kr;
            rR = 0;
         }
         i = (nblksL*ptmms[indx].kb + rL)<<eltsh;
         np = ATL_thrdecompMM_rMNK(ptmms, TA, TB, Mblks, mr, Nblks, nr,
                                   nblksL, rL, A, lda, B, ldb, C, ldc, pL,
                                   indx, COPYC);
         np += ATL_thrdecompMM_rMNK(ptmms, TA, TB, Mblks, mr, Nblks, nr,
                                    nblksR, rR, (TA==AtlasNoTrans)?
                                    MindxT(A,lda*i):MindxT(A,i), lda,
                                    (TB == AtlasNoTrans)?MindxT(B,i):
                                    MindxT(B,i*ldb), ldb, C, ldc, pR,
                                    indx+pL, 1);
         return(np);
      }
   }
/*
 * Now choose to split either M or N.  Want M < N always, so require
 * N to be twice as big before splitting
 */
   if (Nblks > 1 && (Mblks < 2 || Nblks >= Mblks+Mblks))
   {
      #ifdef DEBUG
         fprintf(stderr, "Cut N\n");
      #endif
      nblksL = (d * Nblks);
      nblksR = Nblks - nblksL;
      if (nblksR < nblksL)
      {
         rL = 0;
         rR = nr;
      }
      else
      {
         rL = nr;
         rR = 0;
      }
      i = (nblksL*ptmms[indx].nb+rL) << eltsh;
      np = ATL_thrdecompMM_rMNK(ptmms, TA, TB, Mblks, mr, nblksL, rL, Kblks, kr,
                                A, lda, B, ldb, C, ldc, pL, indx, COPYC);
      np += ATL_thrdecompMM_rMNK(ptmms, TA, TB, Mblks, mr, nblksR, rR,
                                 Kblks, kr, A,  lda, (TB == AtlasNoTrans) ?
                                 MindxT(B,i*ldb) : MindxT(B,i), ldb,
                                 MindxT(C,i*ldc), ldc, pR, indx+pL, COPYC);
      return(np);
   }
/*
 * If we have failed to split N or K, split M if possible
 */
   if (Mblks > 1)
   {
      #ifdef DEBUG
         fprintf(stderr, "Cut M\n");
      #endif
      nblksL = (d * Mblks);
      nblksR = Mblks - nblksL;
      if (nblksR < nblksL)
      {
         rL = 0;
         rR = mr;
      }
      else
      {
         rL = mr;
         rR = 0;
      }
      i = (nblksL*ptmms[indx].mb+rL) << eltsh;
      np = ATL_thrdecompMM_rMNK(ptmms, TA, TB, nblksL, rL, Nblks, nr, Kblks, kr,
                                A, lda, B, ldb, C, ldc, pL, indx, COPYC);
      np += ATL_thrdecompMM_rMNK(ptmms, TA, TB, nblksR, rR, Nblks, nr,
                                 Kblks, kr,
                                 (TA==AtlasNoTrans)?MindxT(A,i):MindxT(A,i*lda),
                                 lda, B, ldb, MindxT(C,i), ldc, pR, indx+pL,
                                 COPYC);
      return(np);
   }
/*
 * If no desirable splitting possible, stop recursion
 */
STOP_REC:
   ptmms[indx].A = A;
   ptmms[indx].B = B;
   ptmms[indx].C = (void*)C;
   ptmms[indx].lda = lda;
   ptmms[indx].ldb = ldb;
   ptmms[indx].ldc = ldc;
   ptmms[indx].M = ptmms[indx].mb*Mblks + mr;
   ptmms[indx].N = ptmms[indx].nb*Nblks + nr;
   ptmms[indx].K = ptmms[indx].kb*Kblks + kr;
   if (COPYC)
   {
      ptmms[indx].nCw = 1;
      ptmms[indx].nCp = ptmms[indx].ownC = 0;
      ptmms[indx].Cinfp[0] = ptmms+indx;
/*
 *    Make ldcw a multiple of 4 that is not a power of 2
 */
      i = ((ptmms[indx].M + 3)>>2)<<2;
      if (!(i & (i-1)))
         i += 4;
      ptmms[indx].ldcw = i;
   }
   else
   {
      ptmms[indx].ldcw = ptmms[indx].nCw = 0;
      ptmms[indx].nCp = ptmms[indx].ownC = 1;
      ptmms[indx].Cinfp[ATL_NTHREADS-1] = ptmms+indx;
   }
   ptmms[indx].Cw = NULL;
   #ifdef DEBUG
      fprintf(stderr, "%d: M=%d, N=%d, K=%d, ownC=%d, nCp=%d, nCw=%d\n",
              indx, ptmms[indx].M, ptmms[indx].N, ptmms[indx].K,
              ptmms[indx].ownC, ptmms[indx].nCp, ptmms[indx].nCw);
   #endif
   return(1);

}

int ATL_thrdecompMM_M
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC)
{
   int j, i, m, p;
   const char *a=A, *c=C;
   const int eltsh = ptmms[0].eltsh, mb = ptmms[0].mb, n = ptmms[0].nb*Nblks+nr,
             k = ptmms[0].kb*Kblks+kr, minblks = Mblks / P,
             extrablks = Mblks - minblks*P;

   for (p=i=0; i < P; i++)
   {
      m = minblks * mb;
      if (i < extrablks)
         m = (minblks + 1)*mb;
      else if (i == extrablks)
         m = minblks*mb + mr;
      else
         m = minblks*mb;
     if (m)
        p++;

      j = ATL_launchorder[i];  /* use log2-launch order */
      ptmms[j].A = a;
      ptmms[j].B = B;
      ptmms[j].C = (void*)c;
      ptmms[j].lda = lda;
      ptmms[j].ldb = ldb;
      ptmms[j].ldc = ldc;
      ptmms[j].M = m;
      ptmms[j].N = n;
      ptmms[j].K = (m) ? k : 0;
      ptmms[j].ownC = 1;
      ptmms[j].nCp = ptmms[j].nCw = 0;
      ptmms[j].Cw = NULL;
      ptmms[j].ldcw = 0;
      m <<= eltsh;
      a = (TA == AtlasNoTrans) ? MindxT(a,m) : MindxT(a,m*lda);
      c = MindxT(c,m);
   }
   return(p);
}

int ATL_thrdecompMM_N
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC)
{
   int j, i, n, p;
   const char *b=B, *c=C;
   const int eltsh = ptmms[0].eltsh, nb = ptmms[0].nb, m = ptmms[0].mb*Mblks+mr,
             k = ptmms[0].kb*Kblks+kr, minblks = Nblks / P,
             extrablks = Nblks - minblks*P;

   for (p=i=0; i < P; i++)
   {
      n = minblks * nb;
      if (i < extrablks)
         n = (minblks + 1)*nb;
      else if (i == extrablks)
         n = minblks*nb + nr;
      else
         n = minblks*nb;
      if (n)
         p++;

      j = ATL_launchorder[i];  /* use log2-launch order */
      ptmms[j].A = A;
      ptmms[j].B = b;
      ptmms[j].C = (void*)c;
      ptmms[j].lda = lda;
      ptmms[j].ldb = ldb;
      ptmms[j].ldc = ldc;
      ptmms[j].M = m;
      ptmms[j].N = n;
      ptmms[j].K = (n) ? k : 0;
      ptmms[j].ownC = 1;
      ptmms[j].nCp = ptmms[j].nCw = 0;
      ptmms[j].Cw = NULL;
      ptmms[j].ldcw = 0;
      n <<= eltsh;
      b = (TB == AtlasNoTrans) ? MindxT(b,n*ldb) : MindxT(b,n);
      c = MindxT(c,n*ldc);
   }
   return(p);
}
int ATL_thrdecompMM_K
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT Mblks, const int mr, ATL_CINT Nblks, const int nr, ATL_CINT Kblks,
    const int kr, const void *A, ATL_INT lda, const void *B, const ATL_INT ldb,
    const void *C, ATL_CINT ldc, const int P, const int indx, const int COPYC)
{
   int j, i, k, p, ldw;
   const char *a=A, *b=B;
   const int eltsh = ptmms[0].eltsh, kb = ptmms[0].kb, m = ptmms[0].mb*Mblks+mr,
             n = ptmms[0].nb*Nblks+nr, minblks = Kblks / P,
             extrablks = Kblks - minblks*P;

   for (p=i=0; i < P; i++)
   {
      k = minblks * kb;
      if (i < extrablks)
         k = (minblks + 1)*kb;
      else if (i == extrablks)
         k = minblks*kb + kr;
      else
         k = minblks*kb;
      if (n)
         p++;

      j = ATL_launchorder[i];  /* use log2-launch order */
      ptmms[j].A = a;
      ptmms[j].B = b;
      ptmms[j].C = (void*)C;
      ptmms[j].lda = lda;
      ptmms[j].ldb = ldb;
      ptmms[j].ldc = ldc;
      ptmms[j].M = m;
      ptmms[j].N = n;
      ptmms[j].K = k;
      if (i)
      {
         ptmms[j].nCw = 1;
         ptmms[j].nCp = ptmms[j].ownC = 0;
         ldw = ((m + 3)>>2)<<2;  /* make ldw mul of 4 */
         if (!(i & (i-1)))
            ldw += 4;            /* make sure ldw not power of 2 */
         ptmms[j].ldcw = ldw;
         ptmms[j].Cinfp[0] = ptmms+j;
      }
      else
      {
         ptmms[j].ldcw = 0;
         ptmms[j].nCp = ptmms[j].ownC = 1;
         ptmms[j].nCw = 0;
         ptmms[j].Cinfp[ATL_NTHREADS-1] = ptmms+j;
      }
      ptmms[j].Cw = NULL;
      k <<= eltsh;
      a = (TA == AtlasNoTrans) ? MindxT(a,lda*k) : MindxT(a,k);
      b = (TB == AtlasNoTrans) ? MindxT(b,k) : MindxT(b,k*ldb);
   }
   return(p);
}

void ATL_EnforceNonPwr2LO(ATL_TMMNODE_t *ptmms, const int P)
/*
 * If threads aren't a power of 2, then the recursive decomposition will
 * fill in the ptmms array in a different order than the log2 spawn.
 * As long as P >= NTHREADS, all entries are filled in, so the worst that
 * happens is that thread 3 does the work that we expect to be done by 4.
 * However, if P < NTHREADS, then threads that launch starts won't have
 * work, so we must make sure that the 1st P elts in launchorder have
 * work to do.
 */
{
   int i, j, k, kk, h;

   if (P >= ATL_NTHREADS)
      return;
   for (i=0; i < P; i++)
   {
      j = ATL_launchorder[i];
      if (!ptmms[j].K)  /* we have found an empty required entry! */
      {
/*
 *      Search for a filled in entry that will not be used in launch
 */
        for (kk=ATL_NTHREADS-1; kk >= P; kk--)
        {
           k = ATL_launchorder[kk];
           if (ptmms[k].K)  /* found a non-empty entry */
           {
              ptmms[j].A = ptmms[k].A;
              ptmms[j].B = ptmms[k].B;
              ptmms[j].C = ptmms[k].C;
              ptmms[j].lda = ptmms[k].lda;
              ptmms[j].ldb = ptmms[k].ldb;
              ptmms[j].ldc = ptmms[k].ldc;
              ptmms[j].M = ptmms[k].M;
              ptmms[j].N = ptmms[k].N;
              ptmms[j].K = ptmms[k].K;
              ptmms[j].ownC = ptmms[k].ownC;
              ptmms[j].nCp = ptmms[k].nCp;
              ptmms[j].nCw = ptmms[k].nCw;
              ptmms[j].Cw = ptmms[k].Cw;
              ptmms[j].ldcw = ptmms[k].ldcw;
              for (h=0; h < ptmms[j].nCw; h++)
                 ptmms[j].Cinfp[h] = (ptmms[k].Cinfp[h] != ptmms+k) ?
                                     ptmms[k].Cinfp[h] : ptmms+j;
              for (h=ATL_NTHREADS-1; h >= ATL_NTHREADS-ptmms[j].nCp; h--)
                 ptmms[j].Cinfp[h] = (ptmms[k].Cinfp[h] != ptmms+k) ?
                                     ptmms[k].Cinfp[h] : ptmms+j;
              ptmms[k].K = 0;                /* entry k now empty */
              break;
           }
        }
        ATL_assert(kk >= P);
      }
   }
}

int ATL_thrdecompMM
   (ATL_TMMNODE_t *ptmms, const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    ATL_CINT M, ATL_CINT N, ATL_CINT K, const void *A, ATL_INT lda,
    const void *B, const ATL_INT ldb, const void *C, ATL_CINT ldc, const int P,
    int *DivideK)
{
   int np, i;
   ATL_CINT Mblks = M/ptmms[0].mb, mr = M-Mblks*ptmms[0].mb;
   ATL_CINT Nblks = N/ptmms[0].nb, nr = N-Nblks*ptmms[0].nb;
   ATL_CINT Kblks = K/ptmms[0].kb, kr = K-Kblks*ptmms[0].kb;
   ATL_CINT mnblks = ((Nblks) ? Nblks : 1) * ((Mblks) ? Mblks : 1);

  *DivideK = 0;
/*
 * First, consider cutting K, which we only do if the number of Kblks
 * dominates the number of blocks we can find in cutting both M & N,
 */
   if (mnblks < P || Kblks > P*mnblks)
   {
      np = ATL_thrdecompMM_rMNK(ptmms, TA, TB, Mblks, mr, Nblks, nr, Kblks, kr,
                                A, lda, B, ldb, C, ldc, P, 0, 0);
      for (i=0; i < np; i++)
      {
         if (ptmms[i].K > 0 && ptmms[i].K < K)
         {
            *DivideK = 1;
            break;
         }
      }
      if (np < ATL_NTHREADS)
         ATL_EnforceNonPwr2LO(ptmms, np);
      return(np);
   }
/*
 * Divide only the M-dimension to cut down on JIK workspace & improve CE
 * efficiency if we have enough M blocks to make it worthwhile;
 * We ask that we can give each thread at least 4 blocks, and that
 * the N diminsion doesn't dominate
 */
   if ((Mblks >= (P<<2) && Nblks < P*Mblks))
      return(ATL_thrdecompMM_M(ptmms, TA, TB, Mblks, mr, Nblks, nr, Kblks, kr,
                               A, lda, B, ldb, C, ldc, P, 0, 0));
/*
 * If none of these special cases are triggered, recursively divide up C
 */
   np = ATL_thrdecompMM_rMN(ptmms, TA, TB, Mblks, mr, Nblks, nr, Kblks, kr,
                            A, lda, B, ldb, C, ldc, P, 0, 0);
   if (np < ATL_NTHREADS)
      ATL_EnforceNonPwr2LO(ptmms, np);
   return(np);
}

int ATL_StructIsInitMM(void *vp)
{
   return(((ATL_TMMNODE_t*)vp)->K);
}

void ATL_DoWorkMM(ATL_LAUNCHSTRUCT_t *lp, void *vp)
/*
 * Current implementation doesn't need lp, but if we had an error queue or
 * something similar we would need it, so keep it around
 */
{
   ATL_TMMNODE_t *mmp = vp;
/*
 * Allocate space if needed, do operation
 */
   if (mmp->nCw)
   {
/*
 *    If malloc fails, we'll do the operation during the combine
 */
      mmp->Cw = malloc(((mmp->ldcw)<<mmp->eltsh)*mmp->N+ATL_Cachelen);
      if (mmp->Cw)
      {
         mmp->gemmK(mmp->M, mmp->N, mmp->K, mmp->alpha, mmp->A, mmp->lda,
                    mmp->B, mmp->ldb, mmp->zero,
                    ATL_AlignPtr(mmp->Cw), mmp->ldcw);
      }
#ifdef DEBUG
      else
         fprintf(stderr, "%d: unable to allocate C(%dx%d)!!\n",
                 mmp->rank, mmp->M, mmp->N);
#endif
   }
   else  /* do GEMM directly into original C; no possibility of failure! */
      mmp->gemmK(mmp->M, mmp->N, mmp->K, mmp->alpha, mmp->A, mmp->lda,
                 mmp->B, mmp->ldb, mmp->beta, mmp->C, mmp->ldc);
}
