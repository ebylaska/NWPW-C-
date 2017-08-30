#include "atlas_misc.h"
#include "atlas_mmtesttime.h"

enum SEARCH {Linear=0, GreedyLinear, GreedySquare};
typedef struct flAgnode ATL_flagnode_t;
struct flAgnode
{
   char *flags;
   ATL_flagnode_t *next;
};

ATL_flagnode_t *NewFlagNode(char *flags)
{
    ATL_flagnode_t *fp;
    int n, i;

    fp = malloc(sizeof(ATL_flagnode_t));
    assert(fp);
    if (flags)
    {
       n = strlen(flags)+1;
       fp->flags = malloc(n*sizeof(char));
       assert(fp->flags);
       for (i=0; i < n; i++)
          fp->flags[i] = flags[i];
    }
    else
       fp->flags = NULL;
    fp->next = NULL;
    return(fp);
}

ATL_flagnode_t *KillFlagNode(ATL_flagnode_t *die)
{
   ATL_flagnode_t *next=NULL;
   if (die)
   {
      next = die->next;
      if (die->flags)
         free(die->flags);
      free(die);
   }
   return(next);
}

void KillAllFlagNodes(ATL_flagnode_t *bp)
{
   while(bp)
      bp = KillFlagNode(bp);
}

int CountFlagNodes(ATL_flagnode_t *bp)
{
   int i;
   for (i=0; bp; bp = bp->next, i++);
   return(i);
}

ATL_flagnode_t *GetGccOptFlags(void)
{
   char *gccflags[] =
{"-O3",
 "-fschedule-insns -fschedule-insns2",
 "-fschedule-insns -fno-schedule-insns2",
 "-fno-schedule-insns -fschedule-insns2",
 "-fno-schedule-insns -fno-schedule-insns2",
 "-fexpensive-optimizations", "-fno-expensive-optimizations",
 "-O1", "-O2", "-Os", "-fno-cprop-registers", "-fcprop-registers",
 "-fcrossjumping", "-fno-crossjumping",
 "-fmodulo-sched", "-fno-branch-count-reg", "-fno-move-loop-invariants",
 "-fno-peephole -fno-peephole2", "-fno-peephole", "-fno-peephole2",
 "-funsafe-loop-optimizations", "-fomit-frame-pointer",
 "-foptimize-register-move", "-fno-optimize-register-move",
 "-fprefetch-loop-arrays", "-fno-prefetch-loop-arrays",
 "-frename-registers", "-fno-rename-registers",
 "-fno-reorder-blocks", "-fno-sched-interblock", "-fno-sched-spec",
 "-fsched-spec-load", "-fsched2-use-superblocks", "-fsched2-use-traces",
 "-fsee", "-freschedule-modulo-scheduled-loops", "-fsingle-precision-constant",
 "-fstrict-aliasing", "-funroll-all-loops", "-funroll-loops",
 "-fno-split-ivs-in-unroller",
 "-fvariable-expansion-in-unroller", "-fno-variable-expansion-in-unroller",
 "-fno-tree-pre", "-fno-tree-fre", "-fno-tree-loop-optimize",
 "-ftree-loop-linear", "-fno-tree-loop-linear", "-fno-tree-loop-im",
 "-ftree-loop-ivcanon", "-fno-tree-loop-ivcanon", "-fivopts", "-fno-ivopts",
 "-fno-tree-dominator-opts", "-ftree-dse", "-fno-tree-dse",
 "-fno-tree-copyrename", "-fno-tree-sink", "-fno-tree-ch", "-fno-tree-ter",
 "-fno-tree-lrs", "-ftree-vectorize -fno-tree-vect-loop-version",
 "-ftree-vectorize -ftree-vect-loop-version", "-fno-tree-salias",
 "-fweb", "-fno-web", "-fno-rerun-loop-opt", "-frerun-loop-opt",
 NULL
};
   int i;
   ATL_flagnode_t *bp, *pf;

   pf = bp = NewFlagNode(gccflags[0]);
   for (i=1; gccflags[i]; i++)
   {
      pf->next = NewFlagNode(gccflags[i]);
      pf = pf->next;
   }
   return(bp);
}

ATL_flagnode_t *ReadFlags(FILE *fpin)
{
   int j, n;
   char ln[1024];
   ATL_flagnode_t *bf=NULL, *pf;  /* base and pointer to flags */

   while (fgets(ln, 1024, fpin))
   {
      n = strlen(ln);
      for (j=n-1; j >= 0 && isspace(ln[j]); j--)
         ln[j] = '\0';  /* strip trailing whitespace */
//      fprintf(stderr, "line = '%s'\n", ln);
      if (bf)
      {
         pf->next = NewFlagNode(ln);
         pf = pf->next;
      }
      else
         pf = bf = NewFlagNode(ln);
   }
   if (fpin != stderr && fpin != stdin)
      fclose(fpin);
   return(bf);
}

char *JoinStrings(char *str0, char *str1)
/*
 * Returns join of two strings (sep by space) in a new string
 */
{
   char *jstr;
   int n, n1, n2, j;

/*
 * Handle degenerate cases
 */
   if (!str0 && !str1)
      return(NULL);
   else if (!str0 || !str1)
   {
      jstr = DupString(str0 ? str0 : str1);
      n = strlen(jstr);
      for (j=n-1; j >= 0 && isspace(jstr[j]); j--)
         jstr[j] = '\0';  /* strip trailing whitespace */
   }
/*
 * Both strings exist, so join them
 */
   else
   {
      n1 = strlen(str0);
      n2 = strlen(str1);
      jstr = malloc(sizeof(char)*(n1+n2+2));
      assert(jstr);

      strcpy(jstr, str0);
      for (j=n1-1; j >= 0 && isspace(jstr[j]); j--)
         jstr[j] = '\0';  /* strip trailing whitespace */
      if (j > 0)
      {
         jstr[j+1] = ' ';
         n1 = j+2;
      }
      else
      {
         n1 = 0;
      }
      strcpy(jstr+n1, str1);
      for (j=n1+n2-1; j >= 0 && isspace(jstr[j]); j--)
         jstr[j] = '\0';  /* strip trailing whitespace */
   }
   return(jstr);
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
   fprintf(stderr,
           "   -1 [0/1]: first line is not/is added to all flags tried\n");
   fprintf(stderr, "   -v # : higher numbers print out more\n");
   fprintf(stderr, "   -f <file> : get flags from file (default stdin)\n");
   fprintf(stderr, "      if <file> is gcc, create gccflags.txt, which\n");
   fprintf(stderr, "      contains most machine independent gcc flags\n");
   fprintf(stderr, "      using baseflags as given in next mandatory arg\n");
   fprintf(stderr,
      "   -m <file> : get mmcase from file (default res/<pre>gMMRES.sum\n");
   fprintf(stderr,
   "   -S [g/G/l] : do O(N) or O(N^2) Gready or O(N) linear search\n");
   exit(-1);
}

char GetFlags(int nargs, char **args, int *verb, int *FirstLine,
              enum SEARCH *srch, FILE **FPFLAG, char **MMNAME)
{
   FILE *fpflag=stdin;
   char *mmname = NULL;
   char pre = 'd';
   int i;

   *srch = GreedyLinear;
   *verb = 1;
   *FirstLine = 0;
   for (i=1; i < nargs; i++)
   {
      if (args[i][0] != '-')
         PrintUsage(args[0], i, args[i]);
      switch(args[i][1])
      {
      case 'v':  /* verbosity */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *verb = atoi(args[i]);
         break;
      case 'p':  /* precision */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         pre = tolower(args[i][0]);
         assert(pre == 's' || pre == 'd' || pre == 'c' || pre == 'z');
         break;
      case '1':  /* 1st line is base flags */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         *FirstLine = atoi(args[i]);
         break;
      case 'f':  /* flag file */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (!strcmp(args[i], "gcc"))
         {
            FILE *fpout;
            ATL_flagnode_t *bp, *fp;
            fpout = fopen("gccflags.txt", "w");
            assert(fpout);
            bp = GetGccOptFlags();
            for (fp=bp; fp; fp = fp->next)
               fprintf(fpout, "%s\n", fp->flags);
            fclose(fpout);
            KillAllFlagNodes(bp);
            exit(0);
         }
         else
         {
            fpflag = fopen(args[i], "r");
            assert(fpflag);
         }
         break;
      case 'S':  /* search algorithm to employ */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         if (args[i][0] == 'G')
            *srch = GreedySquare;
         else if (args[i][0] == 'g')
            *srch = GreedyLinear;
         else if (args[i][0] == 'l')
            *srch = Linear;
         break;
      case 'm':  /* matmul kernel descriptor file */
         if (++i >= nargs)
            PrintUsage(args[0], i, NULL);
         mmname = DupString(args[i]);
         break;
      default:
         PrintUsage(args[0], i, args[i]);
      }
   }

   if (!mmname)
   {
      mmname = malloc(sizeof(char)*16);
      assert(mmname);
      sprintf(mmname, "res/%cgMMRES.sum", pre);
   }

   *MMNAME = mmname;
   *FPFLAG = fpflag;
   return(pre);
}

ATL_flagnode_t *LinearFlagSearch
(
   char pre,                    /* type/precision prefix */
   int verb,                    /* verbosity */
   char *spc,                   /* indent spaces to print */
   double mulB,                 /* muliplier put on base case */
   ATL_mmnode_t *mmp,           /* gemm kernel to time */
   ATL_flagnode_t *baseflags,   /* NULL, or flags to add to all runs */
   ATL_flagnode_t *bp           /* flags to search */
)
/*
 * Tries all candidate flags in bp, and prints the best found.  If mulB != 1,
 * then the present case is advantaged (>1) to avoid adding flags that
 * really don't help performance.
 * RETURNS: pointer to flag that produce the best performance
 */
{
   double mfB=0.0, mf;
   ATL_flagnode_t *fp, *fpB;
   char *jstr = (baseflags) ? baseflags->flags : NULL;
   char *flags, *flagsB=NULL;
   int i;

   if (!spc)
      spc = "";
   printf("%sFINDING BEST FLAG SETTINGS FOR THIS MATMUL KERNEL:\n", spc);
   if (jstr)
   {
      printf("%s   All cases using flags: '%s'\n", spc, jstr);
      fp = NewFlagNode(jstr);  /* always test base flags alone! */
      fp->next = bp;           /* so we're sure new flags are an improvement */
      bp = fp;
   }
   for (i=0,fp=bp; fp; i++,fp = fp->next)
   {
      flags = JoinStrings(jstr, fp->flags);
      if (mmp->cflags)
         free(mmp->cflags);
      mmp->cflags = DupString(flags);
      mf = TimeMMKernel(verb, 1, mmp, pre, mmp->mbB, mmp->nbB, mmp->kbB,
                        0, 0, 0, 1, -1, -1);
      if (verb > 0)
         printf("%s%4d. mf=%.2f, newflags='%s'\n", spc, i, mf, fp->flags);
      if (mf > mfB*mulB)
      {
         mfB = mf;
         fpB = fp;
      }
      free(flags);
   }
   mmp->cflags = NULL;
   printf("%sBEST FLAGS (%.2f) ARE:\n", spc, mfB);
   printf("%s   '%s'\n\n", spc, fpB->flags);
   KillFlagNode(bp);
   return((bp == fpB) ? NULL : fpB);
}

char *GreedyLinearFlagSearch
(
   char pre,                    /* type/precision prefix */
   int verb,                    /* verbosity */
   char *spc,                   /* indent spaces to print */
   double mulB,                 /* muliplier put on base case */
   ATL_mmnode_t *mmp,           /* gemm kernel to time */
   ATL_flagnode_t *baseflags,   /* NULL, or flags to add to all runs */
   ATL_flagnode_t *bp           /* flags to search */
)
/*
 * Tries all the candidate flags in bp, and if anyone of them improves
 * performance, it is added to the default case.
 * RETURNS: string of new base case
 */
{
   double mfB=0.0, mf, mf0;
   ATL_flagnode_t *fp, *fpB;
   char *jstr = (baseflags) ? baseflags->flags : NULL;
   char *flags, *flagsB=NULL, *sp;
   int i;

   if (!spc)
      spc = "";
   printf("%sFINDING BEST FLAG SETTINGS FOR THIS MATMUL KERNEL:\n", spc);
   if (jstr)
   {
      printf("%s   All cases using flags: '%s'\n", spc, jstr);
      fp = NewFlagNode(jstr);  /* always test base flags alone! */
      fp->next = bp;           /* so we're sure new flags are an improvement */
      bp = fp;
   }
   for (i=0,fp=bp; fp; i++,fp = fp->next)
   {
      flags = JoinStrings((fp == bp) ? NULL:jstr, fp->flags);
      if (mmp->cflags)
         free(mmp->cflags);
      mmp->cflags = DupString(flags);
      mf = TimeMMKernel(verb, 1, mmp, pre, mmp->mbB, mmp->nbB, mmp->kbB,
                        0, 0, 0, 1, -1, -1);
      if (bp == fp)
         mf0 = mf;
      if (verb > 0)
         printf("%s%4d. mf=%.2f, flags='%s'\n", spc, i, mf, fp->flags);
      if (mf > mfB*mulB)
      {
         if (verb > 0)
             printf("%s    ---> Adding flag '%s'!\n", spc, fp->flags);
         mfB = mf;
         if (jstr)
            free(jstr);
         jstr = flags;
      }
      else
         free(flags);
   }
   mmp->cflags = NULL;
   printf("%sBEST FLAGS GIVE %.2f%% IMPROVEMENT FROM %.2f to %.2f:\n",
           spc, ((mfB/mf0)-1)*100.0, mf0, mfB);
   printf("%s   '%s'\n\n", spc, jstr);
   KillAllFlagNodes(bp);
   if (baseflags)
      baseflags->flags = jstr;
   return(jstr);
}

void GreedyFlagSearch
(
   char pre,                    /* type/precision prefix */
   int verb,                    /* verbosity */
   ATL_mmnode_t *mmp,           /* gemm kernel to time */
   ATL_flagnode_t *baseflags,   /* NULL, or flags to add to all runs */
   ATL_flagnode_t *bp           /* flags to search */
)
/*
 * For each flag in bp, see if it improves performance over present
 * baseline; if so, add to baseline, else reject.  This algorithm is
 * O(N^2), and will not catch two flags that only help when thrown
 * together unless they are given together in bp!
 */
{
   ATL_flagnode_t *fp, *fpp, *fpB;
   char *sp;
   int N, i;

   N = CountFlagNodes(bp);
   if (verb > 0)
      printf("DOING GREEDY SEARCH USING %d COMBINATIONS\n", N);
   i = 0;
   while (fpB = LinearFlagSearch(pre, verb, "      ", 1.02, mmp, baseflags, bp))
   {
      if (verb > 0)
         printf("   Adding '%s' to list of flags!\n", fpB->flags);
/*
 *    Add new flags to baseflags
 */
      if (!baseflags)
         baseflags = NewFlagNode(fpB->flags);
      else
      {
         sp = JoinStrings(baseflags->flags, fpB->flags);
         free(baseflags->flags);
         baseflags->flags = sp;
      }
/*
 *    Remove selected flags from queue of searched flags
 */
      fpp = NULL;
      for (fp=bp; fp != fpB; fp = fp->next)
         fpp = fp;
      if (fpp)
      {
         fpp->next = KillFlagNode(fpp->next);
      }
      else
         bp = KillFlagNode(bp);
      i++;
      if (verb > 0)
         printf("   GREEDY PASS %d of %d\n", i, N);
   }
}

int main(int nargs, char **args)
{
   ATL_mmnode_t *mmp;
   ATL_flagnode_t *jp=NULL, *bp;
   FILE *fpflag;
   char *mmfile;
   enum SEARCH srch;
   int verb, FirstLine;
   char pre;

   pre = GetFlags(nargs, args, &verb, &FirstLine, &srch, &fpflag, &mmfile);
   mmp = ReadMMFile(mmfile);
   if (!mmp)
   {
      fprintf(stderr, "CANNOT READ FILE '%s'!\n", mmfile);
      exit(-1);
   }
   if (verb > 0)
   {
      printf("\nFINDING BEST FLAGS USING MATMUL KERNEL:\n");
      PrintMMLine(stdout, mmp);
   }
   bp = ReadFlags(fpflag);
   if (FirstLine)
   {
      jp = bp;
      bp = bp->next;
      jp->next = NULL;
   }
   if (srch == Linear)
      LinearFlagSearch(pre, verb, NULL, 1.0, mmp, jp, bp);
   else if (srch == GreedyLinear)
   {
      GreedyLinearFlagSearch(pre, verb, NULL, 1.01, mmp, jp, bp);
      KillAllFlagNodes(jp);
      exit(0);
   }
   else
      GreedyFlagSearch(pre, verb, mmp, jp, bp);
   KillAllFlagNodes(jp);
   KillAllFlagNodes(bp);
   return(0);
}
