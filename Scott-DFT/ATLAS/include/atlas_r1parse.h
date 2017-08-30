#ifndef ATLAS_R1PARSE_H
   #define ATLAS_R1PARSE_H

#include "atlas_genparse.h"


#define R1F_INCACHE     0  /* consider kernel for in-cache gemv */
#define R1F_OUTCACHE    1  /* consider kernel for out-of-cache gemv */
#define R1F_ALLALIGNXY  2  /* X&Y are copied into all legal alignments */
#define R1F_ALIGNX2A    3  /* X forced to same alignment as A */
#define R1F_SINGLE      4  /* single precision */
#define R1F_COMPLEX     5  /* complex arithmetic */
#define R1F_APTRS       6  /* use ptrs rather than lda for column indexing */
#define R1F_X87         7  /* requires the Intel x87 unit */
#define R1F_FYU         8  /* N must be a multiple of YU */
#define R1F_NFLAGS      9
char *R1F_exp[R1F_NFLAGS] =
{
"Consider kernel for in-cache use only",
"Consider kernel for out-of-cache use only",
"X&Y are copied into all legal alignments",
"X forced to have same alignment as A",
"Data uses single precision",
"Data is of complex type",
"use ptrs rather than lda for column indexing",
"Kernel requires the x87 unit for correct operation",
"N must be a multiple of YU"
};

#define R1F_DEFAULT ((1<<R1F_INCACHE) | (1<<R1F_OUTCACHE))
typedef struct R1NODE ATL_r1node_t;
struct R1NODE
{
   double mflop[8];
   ATL_r1node_t *next;
   char *rout, *auth, *comp, *cflags;
   char *str;                   /* tmp string used in generation */
   char *genstr;                /* system(genstr) will generate gened kernel */
   int alignA, alignX, alignY;  /* required alignments */
   int ldamul;                  /* lda must be a multiple of ldamul */
   int ID, YU, XU;              /* unrolling for Y & X vectors */
   int NXU;                     /* # of repetitions of XU */
   int minY, minX;              /* min veclen to call the rout with */
   int SSE;                     /* 0: no SSE, 1: SSE1 req, 2: SSE2 req, etc */
   int asmbits;                 /* valid assemblies in this file */
   int CacheElts;               /* # of cache elts to assume for blocking */
   int flag;
};

static ATL_r1node_t *GetR1Node(void)
{
   ATL_r1node_t *p;
   p = calloc(1, sizeof(ATL_r1node_t));
   assert(p);
   p->flag = R1F_DEFAULT;
   return(p);
}

static ATL_r1node_t *CloneR1Node(ATL_r1node_t *dup)
{
   ATL_r1node_t *p;
   p = malloc(sizeof(ATL_r1node_t));
   assert(p);
   memcpy(p, dup, sizeof(ATL_r1node_t));
   if (dup->rout)
      p->rout = DupString(dup->rout);
   if (dup->auth)
      p->auth = DupString(dup->auth);
   if (dup->comp)
      p->comp = DupString(dup->comp);
   if (dup->cflags)
      p->cflags = DupString(dup->cflags);
   if (dup->str)
      p->str = DupString(dup->str);
   if (dup->genstr)
      p->genstr = DupString(dup->genstr);
   p->next = NULL;
   return(p);
}

static ATL_r1node_t *KillR1Node(ATL_r1node_t *die)
{
   ATL_r1node_t *p=NULL;
   if (die)
   {
      p = die->next;
      if (die->rout)
         free(die->rout);
      if (die->auth)
         free(die->auth);
      if (die->comp)
         free(die->comp);
      if (die->cflags)
         free(die->cflags);
      if (die->str)
         free(die->str);
      if (die->genstr)
         free(die->genstr);
      free(die);
   }
   return(p);
}

static void KillAllR1Nodes(ATL_r1node_t *die)
{
   while (die)
      die = KillR1Node(die);
}

static ATL_r1node_t *ParseR1Line(char *ln)
/*
 * Given a line from a r1 index file (with multiple lines pasted together
 * into one line (ln), return a structure describing that line.
 */
{
   ATL_r1node_t *p;
   char *sp;
   int itmp;
   char ch;

   p = GetR1Node();

   sp = strstr(ln, "LDAMUL=");
   if (sp)
      p->ldamul = atoi(sp+6+1);
   else
      p->ldamul = 0;

   sp = strstr(ln, "CacheElts=");
   if (sp)
      p->CacheElts = atoi(sp+9+1);
   else
      p->CacheElts = 0;

   sp = strstr(ln, "SSE=");
   if (sp)
      p->SSE = atoi(sp+3+1);
   else
      p->SSE = 0;

   sp = strstr(ln, "alignA=");
   if (sp)
      p->alignA = atoi(sp+6+1);
   else
      p->alignA = 0;

   sp = strstr(ln, "alignY=");
   if (sp)
      p->alignY = atoi(sp+6+1);
   else
      p->alignY = 0;

   sp = strstr(ln, "alignX=");
   if (sp)
      p->alignX = atoi(sp+6+1);
   else
      p->alignX = 0;

   sp = strstr(ln, "minX=");
   if (sp)
      p->minX = atoi(sp+4+1);
   else
      p->minX = 0;

   sp = strstr(ln, "minY=");
   if (sp)
      p->minY = atoi(sp+4+1);
   else
      p->minY = 0;

   sp = strstr(ln, "YU=");
   if (sp)
      p->YU = atoi(sp+2+1);
   else
      p->YU = 0;

   sp = strstr(ln, "XU=");
   if (sp)
      p->XU = atoi(sp+2+1);
   else
      p->XU = 0;

   sp = strstr(ln, "ID=");
   if (sp)
      p->ID = atoi(sp+2+1);
   else
      p->ID = 0;

   sp = strstr(ln, "ALIGNX2A=");
   if (sp)
   {
      if (atoi(sp+8+1))
         p->flag |= (1<<R1F_ALIGNX2A);
      else
         p->flag &= ~(1<<R1F_ALIGNX2A);
   }
   sp = strstr(ln, "FYU=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<R1F_FYU);
      else
         p->flag &= ~(1<<R1F_FYU);
   }
   sp = strstr(ln, "ALLALIGNXY=");
   if (sp)
   {
      if (atoi(sp+10+1))
         p->flag |= (1<<R1F_ALLALIGNXY);
      else
         p->flag &= ~(1<<R1F_ALLALIGNXY);
   }
   sp = strstr(ln, "X87=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<R1F_X87);
      else
         p->flag &= ~(1<<R1F_X87);
   }

   sp = strstr(ln, "MFLOP=");
   if (sp)
      GetDoubleArr(sp+6, 8, p->mflop);

   sp = strstr(ln, "ASM=");
   if (sp)
      p->asmbits = asmNames2bitfield(sp+4);



   sp = strstr(ln, "CFLAGS='");
   if (sp)
      p->cflags = GetSingleQuoteString(sp+6+1);
   else
      p->cflags = NULL;

   sp = strstr(ln, "COMP='");
   if (sp)
      p->comp = GetSingleQuoteString(sp+4+1);
   else
      p->comp = NULL;

   sp = strstr(ln, "AUTH='");
   if (sp)
      p->auth = GetSingleQuoteString(sp+4+1);
   else
      p->auth = NULL;

   sp = strstr(ln, "ROUT='");
   if (sp)
      p->rout = GetSingleQuoteString(sp+4+1);
   else
      p->rout = NULL;

   return(p);
}

static void PrintR1Line(FILE *fpout, ATL_r1node_t *np)
{
   int i, j, k;
   char ta, tb;

   if (!np)
      return;
   if (!np->rout)
      np->ID = 0;
   fprintf(fpout, "ID=%d ROUT='%s' AUTH='%s' \\\n",
           np->ID, np->rout ? np->rout : "generated",
           np->auth ? np->auth : "R. Clint Whaley");
   fprintf(fpout, "   ");
   i = 3;
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "CacheElts=%d ", np->CacheElts);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "SSE=%d ", np->SSE);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "alignA=%d ", np->alignA);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "alignY=%d ", np->alignY);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "alignX=%d ", np->alignX);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "minX=%d ", np->minX);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "minY=%d ", np->minY);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "YU=%d ", np->YU);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "XU=%d ", np->XU);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDAMUL=%d ", np->ldamul);

   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "FYU=%d ", FLAG_IS_SET(np->flag, R1F_FYU));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "X87=%d ", FLAG_IS_SET(np->flag, R1F_X87));

   if (np->mflop[0]+np->mflop[1]+np->mflop[2]+np->mflop[3]+np->mflop[4]+
       np->mflop[5]+np->mflop[6] != 0.0)
   {
      if (i > 3) { fprintf(fpout, "\\\n   "); i = 3; }
      i += fprintf(fpout, "MFLOP=%le", np->mflop[0]);
      for (j=7; j && np->mflop[j] == 0.0; j--);
      for (k=1; k <= j; k++)
         i += fprintf(fpout, ",%le", np->mflop[k]);
   }
   if (np->asmbits)
   {
      if (i > 40) { fprintf(fpout, "\\\n   "); i = 3; }
      for (j=0; !(np->asmbits & (1<<j)); j++);
      assert(j < NASMD);
      i += fprintf(fpout, "  ASM=%s", ASMNAM[j]);
      for (j++; j < NASMD; j++)
         if (np->asmbits & (1<<i))
            i += fprintf(fpout, ",%s", ASMNAM[j]);
   }
   if (np->cflags)
   {
      if (i+strlen(np->cflags) > 70) { fprintf(fpout, "\\\n   "); i = 3; }
      i += fprintf(fpout, "  CFLAGS='%s'", np->cflags);
   }
   if (np->comp)
   {
      if (i+strlen(np->comp) > 70) { fprintf(fpout, "\\\n   "); i = 3; }
      i += fprintf(fpout, "  COMP='%s'", np->comp);
   }
   if (i)
      fprintf(fpout, "\n");
}

void PrintR1Nodes(FILE *fpout, ATL_r1node_t *bp)
{
   while (bp)
   {
      PrintR1Line(fpout, bp);
      bp = bp->next;
   }
}

static void WriteR1File(char *file, ATL_r1node_t *nq)
{
   FILE *fpout;

   if (!file || !strcmp(file, "stdout"))
      fpout = stdout;
   else if (!strcmp(file, "stderr"))
      fpout = stderr;
   else
   {
      fpout = fopen(file, "w");
      assert(fpout);
   }
   PrintR1Nodes(fpout, nq);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}

static void WriteR1FileWithPath
   (char pre, char *path, char *file, ATL_r1node_t *nq)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   WriteR1File(ln, nq);
}

static ATL_r1node_t *ReadR1File(char *file)
/*
 * Reads in a standard ATLAS parsable R1 index file, and returns a
 * list of all the kernels defined there.
 */
{
   ATL_r1node_t *nq=NULL, *p;
   FILE *fpin;
   char *ln, *sp;
   int i, j, KeepOn, len;

   if (!file || !strcmp(file, "stdin"))
      fpin = stdin;
   else
      fpin = fopen(file, "r");
   if (!fpin)
      return(NULL);
   nq = p = GetR1Node();
   while (ln = GetJoinedLines(fpin))
   {
      if (ln[0] != '#')
      {
         p->next = ParseR1Line(ln);
         p = p->next;
      }
   }
   fclose(fpin);
   return(KillR1Node(nq));
}

static ATL_r1node_t *ReadR1FileWithPath
   (char pre, char *path, char *file)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   return(ReadR1File(ln));
}
static ATL_r1node_t *DelBadArchR1Kernels(ATL_r1node_t *bp)
/*
 * Weeds out kernels that require SSE/assembly that we haven't got
 */
{
   int asmb=0, die;
   ATL_r1node_t *p, *prev;
   #ifdef ATL_GAS_MIPS
      asmb |= (1<<6);
   #endif
   #ifdef ATL_GAS_PARISC
      asmb |= (1<<5);
   #endif
   #ifdef ATL_GAS_PPC
      asmb |= (1<<4);
   #endif
   #ifdef ATL_GAS_SPARC
      asmb |= (1<<3);
   #endif
   #ifdef ATL_GAS_x8664
      asmb |= (1<<2);
   #endif
   #ifdef ATL_GAS_x8632
      asmb |= (1<<1);
   #endif

   prev = p = bp;
   while (p)
   {
      die = (p->asmbits) ? !(asmb & p->asmbits) : 0;
      #ifndef ATL_SSE3
         if (p->SSE)
         {
            die |= (p->SSE >= 3);
            #ifndef ATL_SSE2
               die |= (p->SSE >= 2);
            #endif
            #ifndef ATL_SSE1
               die |= (p->SSE >= 1);
            #endif
         }
      #endif
      if (die)
      {
         if (p == bp)
            bp = p = KillR1Node(p);
         else
            prev->next = p = KillR1Node(p);
      }
      else
      {
         prev = p;
         p = p->next;
      }
   }
   return(bp);
}

static ATL_r1node_t *FindFastestKernel
(  char pre,             /* precision prefix */
   ATL_r1node_t *bp,  /* kernel queue */
   int imf,              /* which mflop entry to sort by */
   int RESTRICTOK        /* consider restricted kernel? */
)
/*
 * A RESTRICTed kernel is one that requires something that can't be fixed
 * by loop peeling or the like.  Examples include forcing lda to a given
 * multiple, or 16-byte alignment for double complex (can't peel 1/2 of
 * a complex number to make 8-byte aligned array 16).
 * RETURNS: pointer to node in bp that is fastest in context imf wt RESTRCT
 */
{
   double mf;
   ATL_r1node_t *kp, *kmax=bp;
   int size, usize, RKERN;

   if (bp)
   {
      usize = (pre == 'c' || pre == 's') ? 4 : 8;
      if (pre == 'c' || pre == 'd') size = 8;
      else if (pre == 's') size = 4;
      else size = 16;
      mf = bp->mflop[imf];
      for (kp=bp->next; kp; kp = kp->next)
      {
         if (kp->mflop[imf] > mf)
         {
            RKERN = (pre == 'z' || pre == 'c') ? (kp->alignA > usize) : 0;
            RKERN = RKERN | (kp->ldamul > size);
            if (RESTRICTOK | !RKERN)
            {
               mf = kp->mflop[imf];
               kmax = kp;
            }
         }
      }
   }
   return(kmax);
}
static int R1flag2size(int flag)
/*
 * RETURNS: size of type using precision/type bits in flag
 */
{
   int size;

   size = FLAG_IS_SET(flag, R1F_SINGLE) ? 4 : 8;
   size *= FLAG_IS_SET(flag, R1F_COMPLEX) ? 2 : 1;
   return(size);
}

static char R1flag2pre(int flag)
/*
 * RETURNS: correct precision/type prefix based on flag
 */
{
   char pre = 'd';
   if (FLAG_IS_SET(flag, R1F_SINGLE))
      return(FLAG_IS_SET(flag, R1F_COMPLEX) ? 'c' : 's');
   return(FLAG_IS_SET(flag, R1F_COMPLEX) ? 'z' : 'd');
}
static int pre2R1flag(char pre, int flag)
/*
 * RETURNS: flag modified to reflect type/precision indicated by pre
 */
{
   SET_FLAG(flag, R1F_COMPLEX, (pre == 'c' || pre == 'z'));
   SET_FLAG(flag, R1F_SINGLE, (pre == 'c' || pre == 's'));
   return(flag);
}

static void SetAllR1TypeFlags(char pre, ATL_r1node_t *bp)
{
   ATL_r1node_t *p;
   for (p=bp; p; p = p->next)
      p->flag = pre2R1flag(pre, p->flag);
}

void PutKernNameInStr(ATL_r1node_t *r1B)
/*
 * Fills in the proper name for all kernels in r1->str
 */
{
   char ln[32] = {"ATL_dgerk_L0_restrict"};
   char pre;
   const int ipre=4, iL=11, irest=12;

   pre = R1flag2pre(r1B->flag);
   ln[ipre] = pre;
   r1B->str = DupString(ln);
   ln[irest] = '\0';
   r1B->next->str = DupString(ln);
   ln[irest] = '_';
   r1B = r1B->next->next;

   ln[iL] = '2';
   r1B->str = DupString(ln);
   ln[irest] = '\0';
   r1B->next->str = DupString(ln);
   ln[irest] = '_';
   r1B = r1B->next->next;

   ln[iL] = '1';
   r1B->str = DupString(ln);
   ln[irest] = '\0';
   r1B->next->str = DupString(ln);
   ln[irest] = '_';
   r1B = r1B->next->next;

   sprintf(ln, "ATL_%cgerk_L1b_restrict", pre);
   r1B->str = DupString(ln);
   ln[irest+1] = '\0';
   r1B->next->str = DupString(ln);
}

ATL_r1node_t *GetSortedUniqueR1Kerns
   (char pre, ATL_r1node_t *r1kerns, char **aliases)
/*
 * Takes the 8-length queue of rank-1 update kernels:
 *    First 2 are restricted and normal GER out-of-cache kernels
 *    next 2 are restricted & normal in-L2 GER kernels
 *    next 2 are restristed & normal in-L1 GERM kernels
 *    next 2 are restristed & normal out-of-cache, L1-blocked kernels
 *
 * ALIASES: a null-terminated list of string pointers, where pairs
 *          of strings give the correct aliasing: 1st entry is the
 *          routine to be aliased, 2nd is what it should be aliased to.
 *          Every kernel that uses the same actual routine as another
 *          in the 8-length queue is deleted, and a #define is used
 *          to call the appropriate kernel, to avoid unnecessary
 *          code size expansion.
 *          If ALIASES is NULL, then ALIASES is not accessed.
 *          ALIASES must be at least 15 pointers long.
 *
 * RETURNS: new queue with only the unique kernels left (unrestricted
 *          kernels appear first in list), and the p->str entry having
 *          the correct routine/file name.
 * NOTE   : Leaves the original queue intact.
 */
{
   ATL_r1node_t *r1b, *r1p, *r1k, *r1prev;
   char *kern = "gerk";
   int i, ialias=0;
   char *suff[8] = {"_L0", "_L0_restrict", "_L2", "_L2_restrict",
                    "_L1", "_L1_restrict", "_L1b", "_L1b_restrict"};

/*
 * Make sure all routines are present, and there are no extra
 */
  for (i=0, r1p = r1kerns; i < 8; i++, r1p=r1p->next)
     assert(r1p);
  assert(!r1p);
/*
 * Build new queue with the "normal" kernels first
 */
   r1b = CloneR1Node(r1kerns->next);
   r1b->next = r1p = CloneR1Node(r1kerns);
   r1p->next = CloneR1Node(r1kerns->next->next->next);
   r1p->next->next =  CloneR1Node(r1kerns->next->next);
   r1p = r1p->next->next;
   r1p->next = CloneR1Node(r1kerns->next->next->next->next->next);
   r1p->next->next = CloneR1Node(r1kerns->next->next->next->next);
   r1p = r1p->next->next;
   r1p->next = CloneR1Node(r1kerns->next->next->next->next->next->next->next);
   r1p->next->next = CloneR1Node(r1kerns->next->next->next->next->next->next);
   r1p->next->next->next = NULL;
/*
 * Label queue entries with proper kernel names
 */
   for (i=0,r1p = r1b; i < 8; i++, r1p = r1p->next)
   {
       r1p->str = malloc(32*sizeof(char));
       assert(r1p->str);
       sprintf(r1p->str, "ATL_%c%s%s", pre, kern, suff[i]);
   }
/*
 * Add duplicated kernels to alias array, and then get rid of them from Q
 */
   r1prev = r1b;
   r1p = r1b->next;
   while (r1p)
   {
      for (r1k=r1b; r1k != r1p; r1k = r1k->next)
         if (r1k->ID == r1p->ID) break;
      if (r1k != r1p)  /* got duplicate */
      {
         if (aliases)
         {
            aliases[ialias++] = r1p->str;
            r1p->str = NULL;
            aliases[ialias++] = DupString(r1k->str);
         }
         r1prev->next = r1p = KillR1Node(r1p);
      }
      else
      {
        r1prev = r1p;
        r1p = r1p->next;
      }
   }
   if (aliases)
      aliases[ialias] = NULL;
   return(r1b);
}
#endif  /* end atlas_r1parse.h guard */
