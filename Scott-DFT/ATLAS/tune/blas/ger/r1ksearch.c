#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "atlas_misc.h"
#include "atlas_r1parse.h"
#include "atlas_r1testtime.h"


ATL_r1node_t *TimeAllKernels
(
   int L1CacheElts,             /* size of L1 cache in elements */
   int imf,                     /* index into bp->mflop  */
   char pre,                    /* type/precision prefix */
   int  M, int N, int lda,      /* problem dimensions */
   int Fflops,                  /* what to set Force flops to (see below) */
   ATL_r1node_t *kq             /* queue of kernels */
)
/*
 * This routine times all kernels in kq, putting results in kq->mflop[imf]
 * Fflops: for very large problems, set to 0, else the # of flops to force
 * RETURNS: ptr to fastest kernel
 */
{
   ATL_r1node_t *kp, *bestp=NULL;
   int percL1;
   double mf, mfmax=0.0;

   for (kp=kq; kp; kp = kp->next)
   {
      mf = ((double)kp->CacheElts) / ((double)L1CacheElts)*100.0;
      percL1 = mf;
      mf = TimeR1Kernel(1, 0, kp, pre, M, N, lda, percL1, Fflops, -1);
      kp->mflop[imf] = mf;
      if (mf > mfmax)
      {
         bestp = kp;
         mfmax = mf;
      }
   }
   return(bestp);
}

ATL_r1node_t *TimeAllKernelsForContext
(
   int L1CacheElts,     /* size of L1 cache in elements */
   int imf,             /* index into mflop & see below */
   char pre,            /* type/precision prefix */
   ATL_r1node_t *bp     /* queue of kernels */
)
/*
 * This routine naively assumes that 4*L1CacheSize is a decent marker for
 * the L2 size to use in initial L2-timings.
 *
 * imf: parameter describing type of timings to perform:
 *    0: large out-of-cache, blocked for 85% of L1CacheSize
 *    1: large out-of-cache, blocked for MIN(128K,4*L1CacheSize) (L2-blocked)
 *    2: large out-of-cache, no blocking
 *    3: in-L2 problem, no blocking
 *    4: in-L1 problem, no blocking
 * RETURNS: best-performing kernel in context
 */
{
   ATL_r1node_t *r1p, *r1max=NULL;
   double mf, mfmax=0.0;
   int M, N, lda, percL1, cflush=(-1);

   if (imf == 0 || imf == 1 || imf == 2)
   {
      M = (pre == 's') ? 2000 : 1000;
      N = (pre == 'z') ? 500 : 1000;
      cflush = -1;
      percL1 = imf ? 400 : 85;
      percL1 = (imf == 2) ? 0 : percL1;
   }
   else   /* Time in-cache data with no blocking */
      cflush = percL1 = 0;
   for (r1p=bp; r1p; r1p = r1p->next)
   {
      if (imf == 3)  /* L2-contained data */
      {
         N = ((16+r1p->YU-1)/r1p->YU)*r1p->YU;
         M = 128*1024/pre2size(pre);
         if (M > 4*L1CacheElts)
            M = 4*L1CacheElts;
         M /= N;
      }
      else if (imf == 4)  /* L1-contained data */
      {
         N = ((8+r1p->YU-1)/r1p->YU)*r1p->YU;
         M = (85*L1CacheElts)/(N*100);
      }
      lda = M+8;
      mf = TimeR1Kernel(0, 0, r1p, pre, M, N, lda, percL1, -1, cflush);
      if (mf > mfmax)
      {
         mfmax = mf;
         r1max = r1p;
      }
      r1p->mflop[imf] = mf;
   }
   return(r1max);
}

static ATL_r1node_t *DelBadTestKernels(char pre, ATL_r1node_t *bp)
/*
 * Deletes all kernels that can't pass basic usage test
 */
{
   int die;
   ATL_r1node_t *p, *prev;
   int m, n, lda, i, j;
   fprintf(stdout, "\nBEGIN BASIC KERNEL TESTS:\n");

   prev = p = bp;
   while(p)
   {
      m = n = lda = 1000;
      if (FLAG_IS_SET(p->flag, R1F_FYU))
      {
         i = p->YU;
         n = ((n+i-1)/i)*i;
      }
      if (p->ldamul)
      {
         j = pre2size(pre);
         i = p->ldamul / j;
         assert(p->ldamul == i*j);
         lda = ((lda+i-1)/i)*i;
      }
      if (R1KernelFailsTest(0, pre, m, n, lda, p))
      {
         fprintf(stdout, "   NUKING bad kernel %s(%d)\n", p->rout, p->ID);
         if (p == bp)
            bp = p = KillR1Node(p);
         else
            prev->next = p = KillR1Node(p);
      }
      else
      {
         fprintf(stdout, "   Kernel %s(%d) passes basic test\n",
                 p->rout, p->ID);
         prev = p;
         p = p->next;
      }
   }
   fprintf(stdout, "DONE BASIC KERNEL TESTS:\n\n");
   return(bp);
}

ATL_r1node_t *ChooseKernelBlocking
(ATL_r1node_t *L1,      /* L1 blocked kernel, scope imf=0 */
 ATL_r1node_t *L2,      /* L2 blocked kernel, scope imf=1 */
 ATL_r1node_t *NOB      /* no-blocking kernel, scope imf=2 */
)
/*
 * This routine compares 3 different blocking strategy using 1-3 kernels
 * (i.e., they may all be the same kernel).
 * It is possible that the data we use may stay in a very large L3 cache,
 * so only accept no-blocking if it is significantly faster than doing
 * blocking.  L1-blocking will tend to minimize the number of kernels required,
 * so stress it very slightly more than L2 blocking.
 * RETURNS: cloned node of best blocking/kernel
 */
{
   ATL_r1node_t *best;
   double mf1, mf2, mf3;

   mf1 = L1->mflop[0] * 1.05;  /* give small adv to safest option, L1 blk */
   mf2 = L2->mflop[1] * 1.03;  /* give small adv to blocking over not */
   mf3 = NOB->mflop[2];        /* no blocking loss may vary by size, penalize */
   if (mf1 > mf2 && mf1 > mf3)
      best = L1;
   else if (mf2 > mf1 && mf2 > mf3)
      best = L2;
   else
      best = NOB;
   return(best);
}

double ExhCESrch
/*
 * RETURNS: best mflop found
 */
(
   ATL_r1node_t *r1p,           /* kernel to search with */
   char pre,                    /* type/precision prefix */
   int M, int N, int lda,
   int stride,                  /* stride to search with, real pL = pL*stride */
   int pLL,                     /* lower percL1 (mul by stride for real val) */
   int pLH,                     /* higher percL1 */
   double mfL,                  /* mflops achieved by lower */
   double mfH,                  /* mflops achieved by higher */
   int *pLB                     /* the best percL1 found */
)
{
   int plm, plb;
   double mf;

   plm = (pLH-pLL)>>1;
   if (plm < 1)
   {
      if (mfL < mfH)
      {
         *pLB = pLH;
         return(mfH);
      }
      *pLB = pLL;
      return(mfL);
   }
   plm += pLL;
   mf = TimeR1Kernel(0, 0, r1p, pre, M, N, lda, plm*stride, 0, -1);
   fprintf(stdout, "%6d  %6d  %6d  %6d  %9.2f\n", M, N, lda, plm*stride, mf);
   mfL = ExhCESrch(r1p, pre, M, N, lda, stride, pLL, plm, mfL, mf, pLB);
   mfH = ExhCESrch(r1p, pre, M, N, lda, stride, plm, pLH, mf, mfH, &plb);
   if (mfH > mfL)
   {
      mfL = mfH;
      *pLB = plb;
   }
   return(mfL);
}


void ExhaustiveCESrch
/*
 * Performs an exhaustive search on entire range using recursive halving
 * And modifies CE to be best % of L1 size found.
 */
(
   ATL_r1node_t *r1p,           /* kernel to search with */
   int imf,                     /* set r1p->mflop[imf] to best perf */
   char pre,                    /* type/precision prefix */
   int M, int N, int lda,       /* prob size to tune with */
   int stride,                  /* stride to search with, real pL = pL*stride */
   int pLL,                     /* lower percL1 (mul by stride for real val) */
   int pLH                      /* higher percL1 */
)
{
   double mfH, mfL;
   int percL1;
   fprintf(stdout, "     M       N     lda  percL1       mflop\n");
   fprintf(stdout, "======  ======  ======  ======  ==========\n");

   mfL = TimeR1Kernel(0, 0, r1p, pre, M, N, lda, pLL*stride, 0, -1);
   fprintf(stdout, "%6d  %6d  %6d  %6d  %9.2f\n", M, N, lda, pLL*stride, mfL);
   mfH = TimeR1Kernel(0, 0, r1p, pre, M, N, lda, pLH*stride, 0, -1);
   fprintf(stdout, "%6d  %6d  %6d  %6d  %9.2f\n", M, N, lda, pLH*stride, mfH);
   mfL = ExhCESrch(r1p, pre, M, N, lda, stride, pLL, pLH, mfL, mfH, &percL1);
   fprintf(stdout, "\nBEST CASE %d percent of L1, MFLOP=%.2f\n\n",
           percL1*stride, mfL);
   r1p->mflop[imf] = mfL;
   r1p->CacheElts = percL1*stride;
}

static int GetMaxID(ATL_r1node_t *r1b)
{
   ATL_r1node_t *r1p;
   int maxID=0;

   for (r1p=r1b; r1p; r1p = r1p->next)
      maxID = Mmax(maxID, r1p->ID);

   return(maxID);
}
#ifdef ATL_SSE3
void FillInGenNode
(
   ATL_r1node_t *r1p,   /* data structure to fill in */
   int nmu,             /* unrolling on mu unrolled inner (X) loop */
   int mu,              /* register blk to apply to X */
                        /* total X unrolling is nmu*mu! */
   int nu,              /* unroll&jam on outer (Y) loop */
   int evenlda,         /* assume X&Y have same (mis)align, lda is even */
   int allalign16,      /* assume X,Y,A all aligned to 16 bytes */
   int aptrs            /* use ptrs rather than lda for column indexing */
)
{
   char ln[2048];
   r1p->NXU = nmu;
//   HERE HERE HERE
}
ATL_r1node_t *SrchSSEGen
(
   char pre,  /* precision prefix indicating type */
   int maxID  /* IDs > maxID are safe to use in generation */
)
/*
 * This search finds the best kernel provided by the ATLAS SSE GER generator
 * RETURNS: a list of kernels to be added to the multiple implementation srch
 */
{
   return(NULL);
}
#endif

ATL_r1node_t *SortRestricted
/*
 * Sorts queue of kernels into a queue of unrestricted kernels (can be
 * always be used) and restricted (only used under certain conditions,
 * such as lda*size a multiple 16).  Destroys R1B in process.
 */
(
   ATL_r1node_t *R1B,   /* original queue containing restricted & unrest */
   ATL_r1node_t **R1R   /* queue of only restricted kernels */
)
{
   ATL_r1node_t *r1B=NULL, *r1R=NULL, *r1b=NULL, *r1r=NULL, *r1p;

   for (r1p=R1B; r1p; r1p = r1p->next)
   {
      if (r1p->ldamul > R1flag2size(r1p->flag))
      {
         if (r1R)
            r1r->next = r1p;
         else
            r1R = r1p;
         r1r = r1p;
      }
      else
      {
         if (r1B)
            r1b->next = r1p;
         else
            r1B = r1p;
         r1b = r1p;
      }
   }
   if (r1R)
      r1r->next = NULL;
   assert(r1B);
   r1b->next = NULL;
   *R1R = r1R;
   return(r1B);
}

void PrintUsage(char *name, int ierr, char *flag)
{
   if (ierr > 0)
      fprintf(stderr, "Bad argument #%d: '%s'\n",
              ierr, flag ? flag : "Not enough arguments");
   else if (ierr < 0)
      fprintf(stderr, "ERROR: %s\n", flag);
   fprintf(stderr, "USAGE: %s [flags]:\n", name);
   fprintf(stderr, "   -p [s,d,c,z]: set precision prefix \n");
   exit(ierr ? ierr : -1);
}

void GetFlags(int nargs, char **args, char *pre)
{
   int i;
   char ch;

   *pre = 'd';
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'p':  /* -p <pre> */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);

         ch = tolower(args[i][0]);
         assert(ch == 's' || ch == 'd' || ch == 'c' || ch == 'z');
         *pre = ch;
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }
}

void WriteR1SummFile(char pre, ATL_r1node_t *r1b)
{
   char ln[256];
   FILE *fpout;

   sprintf(ln, "res/%cR1SUMM", pre);
   fpout = fopen(ln, "w");
   assert(fpout);
   fprintf(fpout, "#\n#MFLOP array has following meaning by index:\n");
   fprintf(fpout, "#   0 : Out-of-cache, L1 blocked\n");
   fprintf(fpout, "#   1 : Out-of-cache, L2 blocked\n");
   fprintf(fpout, "#   2 : Out-of-cache, no blocking\n");
   fprintf(fpout, "#   3 : Problem preloaded to L2, no blocking\n");
   fprintf(fpout, "#   4 : Problem preloaded to L1, no blocking\n#\n");
   fprintf(fpout, "#\n#Each kernel context has two kernels:\n");
   fprintf(fpout,
           "#   1st has a restriction and so can't be called all the time\n");
   fprintf(fpout, "#   2nd is used whenever restricted kernel can't be\n");
   fprintf(fpout, "#   -> If kernels are same, no restricted kernel needed\n");
   fprintf(fpout, "#\n");
   fprintf(fpout,
   "# --------------------------------------------------------------\n");
   fprintf(fpout,
   "# Next two lines are GER kernels to use for out-of-cache timings\n");
   fprintf(fpout,
   "# --------------------------------------------------------------\n");
   PrintR1Line(fpout, r1b);
   PrintR1Line(fpout, r1b->next);
   r1b = r1b->next->next;
   fprintf(fpout, "# ---------------------------------------------------------------------------\n");
   fprintf(fpout, "# The next two GER kernels are for use when ops are preloaded to the L2 cache\n");
   fprintf(fpout, "# ---------------------------------------------------------------------------\n");
   PrintR1Line(fpout, r1b);
   PrintR1Line(fpout, r1b->next);
   r1b = r1b->next->next;
   fprintf(fpout, "# ---------------------------------------------------------------------------\n");
   fprintf(fpout, "# The next two GER kernels are for use when ops are preloaded to the L1 cache\n");
   fprintf(fpout, "# ---------------------------------------------------------------------------\n");
   PrintR1Line(fpout, r1b);
   PrintR1Line(fpout, r1b->next);
   r1b = r1b->next->next;
   fprintf(fpout,
           "# -----------------------------------------------------------\n");
   fprintf(fpout,
           "# Last two lines are L1-blocked kernels for building SYR/SYR2\n");
   fprintf(fpout,
           "# -----------------------------------------------------------\n");
   PrintR1Line(fpout, r1b);
   PrintR1Line(fpout, r1b->next);
   fclose(fpout);
}

main(int nargs, char **args)
{
   ATL_r1node_t *r1b, *r1p, *r1r, *r1bestL2b;
   ATL_r1node_t *r1bestOC, *r1bestOCr, *r1bestL1b, *r1bestL1br;
   ATL_r1node_t *r1bestL1, *r1bestL1r, *r1bestL2, *r1bestL2r;
   FILE *fpin, *fpout;
   double mf, percL1;
   int i, L1CacheElts, CE1, CE2, maxID, bigN;
   char ln[128];
   char pre;

   GetFlags(nargs, args, &pre);
   if (pre == 'd' || pre == 'c')
      bigN = 2300;
   else if (pre == 's')
      bigN = 3000;
   else
      bigN = 1500;
   sprintf(ln, "res/%cR1SUMM", pre);
   L1CacheElts = (GetL1CacheSize()*1024) / pre2size(pre);
   r1b = ReadR1File(ln);
   SetAllR1TypeFlags(pre, r1b);
   maxID = GetMaxID(r1b);
   if (r1b)
   {
      if (r1b->next->next->mflop[3] <= 0.0)
      {
/*
 *       Retime out-of-cache kernel
 */
         r1p = r1b;
         r1r = r1b->next->next;
         r1b->next->next = NULL;
         if (r1p->next->ID == r1p->ID)
         {
            TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1p->next);
            r1p->mflop[0] = r1p->next->mflop[0];
         }
         else
            TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1p);
         r1b->next->next = r1r;
/*
 *       Retime in-L2 timings
 */
         r1p = r1b->next->next;
         r1r = r1p->next->next;
         r1p->next->next = NULL;
         if (r1p->ID == r1p->next->ID)
         {
            TimeAllKernelsForContext(L1CacheElts, 3, pre, r1p->next);
            r1p->mflop[3] = r1p->next->mflop[3];
         }
         else
            TimeAllKernelsForContext(L1CacheElts, 3, pre, r1p);
         r1p->next->next = r1r;
/*
 *       Retime in-L1 timings
 */
         r1p = r1p->next->next;
         r1r = r1p->next->next;
         r1p->next->next = NULL;
         if (r1p->ID == r1p->next->ID)
         {
            TimeAllKernelsForContext(L1CacheElts, 4, pre, r1p->next);
            r1p->mflop[4] = r1p->next->mflop[4];
         }
         else
            TimeAllKernelsForContext(L1CacheElts, 4, pre, r1p);
         r1p->next->next = r1r;
/*
 *       Retime L1-blocked out-of-cache timings
 */
         r1p = r1p->next->next;
         if (r1p->next->ID == r1p->ID)
         {
            TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1p->next);
            r1p->mflop[0] = r1p->next->mflop[0];
         }
         else
            TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1p);
         WriteR1SummFile(pre, r1b);
      }
      PrintR1Nodes(stdout, r1b);
      KillAllR1Nodes(r1b);
      exit(0);
   }
   sprintf(ln, "CASES/%cr1cases.idx", pre);
   r1b = ReadR1File(ln);
   fprintf(stdout, "\nCases read in:\n");
   WriteR1File("stdout", r1b);
   r1b = DelBadArchR1Kernels(r1b);
   r1b = DelBadTestKernels(pre, r1b);
   fprintf(stdout, "\nSurviving cases:\n");
   WriteR1File("stdout", r1b);

   #ifdef ATL_SSE3   /* current generator requires SSE3 */
      r1p = SrchSSEGen(pre, maxID);
      if (r1p)
      {
         while (r1p->next)
            r1p = r1p->next;
         r1p->next = r1b;
         r1b = r1p;
      }
   #endif
   r1b = SortRestricted(r1b, &r1r);     /* sort into general & rest kernels */
/*
 * Find best general kernel for L1-blocked out-of-cache behavior; use
 * cache elements of 85% of L1 size
 */
   for (r1p=r1b; r1p; r1p = r1p->next)
      r1p->CacheElts = 0.85 * L1CacheElts;
   fprintf(stdout, "\nBEGIN L1-BLOCKED TUNING\n");
   r1p = TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1b);
   r1bestL1b = CloneR1Node(r1p);
   r1bestL1b->next = NULL;
   fprintf(stdout, "DONE L1-BLOCKED TUNING, CHOSE '%s' (%.2f)\n",
           r1bestL1b->rout, r1bestL1b->mflop[0]);
/*
 * Find best L1 blocking % of cache, convert back to Elts
 */
   ExhaustiveCESrch(r1bestL1b, 0, pre, bigN, bigN, bigN, 2, 25, 50);
   r1bestL1b->CacheElts = L1CacheElts*0.01*r1bestL1b->CacheElts;
/*
 * Find best kernel for L2-blocked out-of-cache context; use min of 4*L1
 * and 128K as good effective L2 estimate, convert back to elts
 */
   CE1 = 4 * L1CacheElts;
   CE1 = Mmin(CE1, 128*1024/pre2size(pre));
   for (r1p=r1b; r1p; r1p = r1p->next)
      r1p->CacheElts = CE1;
   fprintf(stdout, "\nBEGIN L2-BLOCKED TUNING\n");
   r1p = TimeAllKernels(L1CacheElts, 1, pre, bigN, bigN, bigN, 0, r1b);
   r1bestL2b = CloneR1Node(r1p);
   r1bestL2b->next = NULL;
   fprintf(stdout, "DONE L2-BLOCKED TUNING, CHOSE '%s' (%.2f)\n",
           r1bestL2b->rout, r1bestL2b->mflop[1]);
   ExhaustiveCESrch(r1bestL2b, 0, pre, bigN, bigN, bigN,
                    50, 3, 16);
   r1bestL2b->CacheElts = L1CacheElts*0.01*r1bestL2b->CacheElts;

   printf("BEST L1-blocked kernel:\n");
   WriteR1File("stdout", r1bestL1b);
   printf("BEST L2-blocked kernel:\n");
   WriteR1File("stdout", r1bestL2b);
/*
 * Find best kernel for in-L2 and in-L1 usage
 */
   r1bestL1 = TimeAllKernelsForContext(L1CacheElts, 4, pre, r1b);
   r1bestL2 = TimeAllKernelsForContext(L1CacheElts, 3, pre, r1b);
   r1bestOC = TimeAllKernelsForContext(L1CacheElts, 2, pre, r1b);
   r1bestL1 = CloneR1Node(r1bestL1);
   r1bestL2 = CloneR1Node(r1bestL2);
   r1bestL2->CacheElts = r1bestL2b->CacheElts;
   r1bestOC = CloneR1Node(r1bestOC);
   r1bestOC->CacheElts = 0;
   r1bestL1->next = r1bestL2->next = r1bestOC->next = NULL;
   r1bestL1->CacheElts = r1bestL1b->CacheElts;
/*
 * Figure out what type of blocking to use for out-of-cache context
 */
   r1p = ChooseKernelBlocking(r1bestL1, r1bestL2, r1bestOC);
   if (r1p != r1bestOC)
   {
      KillR1Node(r1bestOC);
      r1bestOC = CloneR1Node(r1p);
   }
/*
 * When timing restricted kernels, we will use best CE found by unresticted,
 * and we can tune only the contexts we use: r1bestOC, r1bestL1[b], r1bestL2.
 */
   if (r1r)
   {
/*
 *    Find best restricted out-of-cache, L1blocked, in-L1 & in-L2 kernels
 */
      fprintf(stdout, "\nBEGIN RESTRICTED OUT-OF-CACHE TUNING\n");
      for (r1p=r1r; r1p; r1p = r1p->next)
         r1p->CacheElts = r1bestOC->CacheElts;
      if (r1bestOC->CacheElts > L1CacheElts)
         i = 1;
      else if (r1bestOC->CacheElts)
         i = 0;
      else
         i = 2;
      r1p = TimeAllKernels(L1CacheElts, i, pre, bigN, bigN, bigN, 0, r1r);
      r1bestOCr = CloneR1Node(r1p);
      PrintR1Nodes(stdout, r1bestOCr);
      fprintf(stdout, "\nBEGIN RESTRICTED L1-BLOCKED TUNING TUNING\n");
      for (r1p=r1r; r1p; r1p = r1p->next)
         r1p->CacheElts = r1bestL1b->CacheElts;
      r1p = TimeAllKernels(L1CacheElts, 0, pre, bigN, bigN, bigN, 0, r1r);
      r1bestL1br = CloneR1Node(r1p);
      PrintR1Nodes(stdout, r1bestL1br);
      fprintf(stdout, "\nBEGIN RESTRICTED in-L1 TUNING TUNING\n");
      r1bestL1r = TimeAllKernelsForContext(L1CacheElts, 4, pre, r1r);
      r1bestL1r = CloneR1Node(r1bestL1r);
      PrintR1Nodes(stdout, r1bestL1r);
      fprintf(stdout, "\nBEGIN RESTRICTED in-L2 TUNING TUNING\n");
      for (r1p=r1r; r1p; r1p = r1p->next)
         r1p->CacheElts = r1bestL2b->CacheElts;
      r1bestL2r = TimeAllKernelsForContext(L1CacheElts, 3, pre, r1r);
      r1bestL2r = CloneR1Node(r1bestL2r);
      PrintR1Nodes(stdout, r1bestL2r);
/*
 *    If restricted kern not 5% better than unrestricted, nuke it
 */
      r1bestL1br->next = r1bestL1r->next = r1bestL2r->next =
                         r1bestOCr->next = NULL;
      if (r1bestL2r->mflop[3] < r1bestL2->mflop[3]*1.05)
      {
         KillR1Node(r1bestL2r);
         r1bestL2r = NULL;
      }
      if (r1bestL1r->mflop[4] < r1bestL1->mflop[4]*1.05)
      {
         KillR1Node(r1bestL1r);
         r1bestL1r = NULL;
      }
      if (r1bestL1br->mflop[0] < r1bestL1b->mflop[0]*1.05)
      {
         KillR1Node(r1bestL1br);
         r1bestL1br = NULL;
      }
      if (r1bestOC->CacheElts > L1CacheElts)
         i = 1;
      else if (r1bestOC->CacheElts > 0)
         i = 0;
      else
         i = 2;
      if (r1bestOCr->mflop[i] < r1bestOC->mflop[i]*1.05)
         r1bestOCr = NULL;
   }
   else
      r1bestOCr = r1bestL1r = r1bestL2r = r1bestL1br = NULL;
/*
 * Link restricted with unrestricted, duping unrest if rest NULL
 */
   KillAllR1Nodes(r1b);
   KillAllR1Nodes(r1r);
   r1b = (r1bestOCr) ? r1bestOCr : CloneR1Node(r1bestOC);
   r1p = r1b->next = r1bestOC;
   r1p->next = (r1bestL2r) ? r1bestL2r : CloneR1Node(r1bestL2);
   r1p = r1p->next->next = r1bestL2;
   r1p->next = (r1bestL1r) ? r1bestL1r : CloneR1Node(r1bestL1);
   r1p = r1p->next->next = r1bestL1;
   r1p->next = (r1bestL1br) ? r1bestL1br : CloneR1Node(r1bestL1b);
   r1p = r1p->next->next = r1bestL1b;
   r1p->next = NULL;
   WriteR1SummFile(pre, r1b);
   KillAllR1Nodes(r1b);
   return(0);
}

