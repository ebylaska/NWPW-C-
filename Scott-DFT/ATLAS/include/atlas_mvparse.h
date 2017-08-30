#ifndef ATLAS_MVPARSE_H
   #define ATLAS_MVPARSE_H

#include "atlas_genparse.h"

#define MVF_INCACHE     0  /* consider kernel for in-cache gemv */
#define MVF_OUTCACHE    1  /* consider kernel for out-of-cache gemv */
#define MVF_ALLALIGNXY  2  /* X&Y are copied into all legal alignments */
#define MVF_AXPYBASED   3  /* 0:ddot based, 1: axpy-based */
#define MVF_GEMMBASED   4  /* gemm-based */
#define MVF_LDYTOP      5  /* 0: load Y value after dot product */
#define MVF_CONJDEF     6  /* 1: conj(A) if Conj_ is defined */
#define MVF_X87         7  /* requires the Intel x87 unit */
#define MVF_FYU         8  /* Length of Y must be a multiple of YU */
#define MVF_SINGLE      9  /* 1: single precision, else double */
#define MVF_COMPLEX    10  /* 1: complex type, else real */

#define MVF_DEFAULT ((1<<MVF_INCACHE) | (1<<MVF_OUTCACHE))
typedef struct MVNODE ATL_mvnode_t;
struct MVNODE
{
   double mflop;
   ATL_mvnode_t *next;
   char *rout, *auth, *comp, *cflags;
   char *str;                   /* tmp string used in generation */
   char *genstr;                /* system(genstr) will generate gened kernel */
   int alignA, alignX, alignY;  /* required alignments */
   int ID, YU, XU;              /* unrolling for Y & X vectors */
   int minY, minX;              /* min veclen to call the rout with */
   int SSE;                     /* 0: no SSE, 1: SSE1 req, 2: SSE2 req, etc */
   int asmbits;                 /* valid assemblies in this file */
   enum ATLAS_TRANS TA;
};

static ATL_mvnode_t *GetMVNode(void)
{
   ATL_mvnode_t *p;
   p = calloc(1, sizeof(ATL_mvnode_t));
   assert(p);
   p->TA = AtlasNoTrans;
   p->flag = MVF_DEFAULT;
   return(p);
}

static ATL_mvnode_t *CloneMVNode(ATL_mvnode_t *dup)
{
   ATL_mvnode_t *p;
   p = malloc(sizeof(ATL_mvnode_t));
   assert(p);
   memcpy(p, dup, sizeof(ATL_mvnode_t));
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

static ATL_mvnode_t *KillMVNode(ATL_mvnode_t *die)
{
   ATL_mvnode_t *p=NULL;
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

static void KillAllMVNodes(ATL_mvnode_t *die)
{
   while (die)
      die = KillMVNode(die);
}

static ATL_mvnode_t *ParseMVLine(char *ln)
/*
 * Given a line from a mv index file (with multiple lines pasted together
 * into one line (ln), return a structure describing that line.
 */
{
   ATL_mvnode_t *p;
   char *sp;
   int itmp;
   char ch;

   p = GetMVNode();

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
         p->flag |= (1<<MVF_ALIGNX2A);
      else
         p->flag &= ~(1<<MVF_ALIGNX2A);
   }
   sp = strstr(ln, "FYU=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<MVF_FYU);
      else
         p->flag &= ~(1<<MVF_FYU);
   }
   sp = strstr(ln, "CONJDEF=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MVF_CONJDEF);
      else
         p->flag &= ~(1<<MVF_CONJDEF);
   }
   sp = strstr(ln, "GEMMBASED=");
   if (sp)
   {
      if (atoi(sp+9+1))
         p->flag |= (1<<MVF_GEMMBASED);
      else
         p->flag &= ~(1<<MVF_GEMMBASED);
   }
   sp = strstr(ln, "AXPYBASED=");
   if (sp)
   {
      if (atoi(sp+9+1))
         p->flag |= (1<<MVF_AXPYBASED);
      else
         p->flag &= ~(1<<MVF_AXPYBASED);
   }
   sp = strstr(ln, "ALLALIGNXY=");
   if (sp)
   {
      if (atoi(sp+10+1))
         p->flag |= (1<<MVF_ALLALIGNXY);
      else
         p->flag &= ~(1<<MVF_ALLALIGNXY);
   }
   sp = strstr(ln, "LDYTOP=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MVF_LDYTOP);
      else
         p->flag &= ~(1<<MVF_LDYTOP);
   }
   sp = strstr(ln, "X87=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<MVF_X87);
      else
         p->flag &= ~(1<<MVF_X87);
   }

   sp = strstr(ln, "MFLOP=");
   if (sp)
      GetDoubleArr(sp+6, 8, p->mflop);

   sp = strstr(ln, "ASM=");
   if (sp)
      p->asmbits = asmNames2bitfield(sp+4);


   sp = strstr(ln, "TA='");
   if (sp)
   {
      ch = tolower(sp[4]);
      if (ch == 'n')
         p->TA = AtlasNoTrans;
      else if (ch == 'c')
         p->TA = AtlasConjTrans;
      else if (ch == 't')
         p->TA = AtlasTrans;
      else
         assert(0);
   }
   sp = strstr(ln, "TB='");
   if (sp)
   {
      ch = tolower(sp[4]);
      if (ch == 'n')
         p->TB = AtlasNoTrans;
      else if (ch == 'c')
         p->TB = AtlasConjTrans;
      else if (ch == 't')
         p->TB = AtlasTrans;
      else
         assert(0);
   }
   sp = strstr(ln, "TA='");
   if (sp)
   {
      ch = tolower(sp[4]);
      if (ch == 'n')
         p->TA = AtlasNoTrans;
      else if (ch == 'c')
         p->TA = AtlasConjTrans;
      else if (ch == 't')
         p->TA = AtlasTrans;
      else
         assert(0);
   }

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

static void PrintMVLine(FILE *fpout, ATL_mvnode_t *np)
{
   int i, j, k;
   char ta, tb;

   if (!np)
      return;
   if (!np->rout)
      np->ID = 0;
   if (np->TA == AtlasConjTrans) ta = 'C';
   else if (np->TA == AtlasTrans) ta = 'T';
   else ta = 'N';
   fprintf(fpout, "ID=%d ROUT='%s' AUTH='%s' TA='%c' \\\n",
           np->ID, np->rout ? np->rout : "generated",
           np->auth ? np->auth : "R. Clint Whaley", ta);
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
   i += fprintf(fpout, "GEMMBASED=%d ", FLAG_IS_SET(np->flag, MVF_GEMMBASED));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "AXPYBASED=%d ", FLAG_IS_SET(np->flag, MVF_AXPYBASED));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "ALLALIGNXY=%d ", FLAG_IS_SET(np->flag, MVF_ALLALIGNXY));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDYTOP=%d ", FLAG_IS_SET(np->flag, MVF_LDYTOP));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "FYU=%d ", FLAG_IS_SET(np->flag, MVF_FYU));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "X87=%d ", FLAG_IS_SET(np->flag, MVF_X87));

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

void PrintMVNodes(FILE *fpout, ATL_mvnode_t *bp)
{
   while (bp)
   {
      PrintMVLine(fpout, bp);
      bp = bp->next;
   }
}

static void WriteMVFile(char *file, ATL_mvnode_t *nq)
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
   PrintMVNodes(fpout, nq);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}

static void WriteMVFileWithPath
   (char pre, char *path, char *file, ATL_mvnode_t *nq)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   WriteMVFile(ln, nq);
}

static ATL_mvnode_t *ReadMVFile(char *file)
/*
 * Reads in a standard ATLAS parsable MV index file, and returns a
 * list of all the kernels defined there.
 */
{
   ATL_mvnode_t *nq=NULL, *p;
   FILE *fpin;
   char *ln, *sp;
   int i, j, KeepOn, len;

   if (!file || !strcmp(file, "stdin"))
      fpin = stdin;
   else
      fpin = fopen(file, "r");
   if (!fpin)
      return(NULL);
   nq = p = GetMVNode();
   while (ln = GetJoinedLines(fpin))
   {
      if (ln[0] != '#')
      {
         p->next = ParseMVLine(ln);
         p = p->next;
      }
   }
   fclose(fpin);
   return(KillMVNode(nq));
}

static ATL_mvnode_t *ReadMVFileWithPath
   (char pre, char *path, char *file)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   return(ReadMVFile(ln));
}
static ATL_mvnode_t *DelBadArchMVKernels(ATL_mvnode_t *bp)
/*
 * Weeds out kernels that require SSE/assembly that we haven't got
 */
{
   int asmb=0, die;
   ATL_mvnode_t *p, *prev;
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
            bp = p = KillMVNode(p);
         else
            prev->next = p = KillMVNode(p);
      }
      else
      {
         prev = p;
         p = p->next;
      }
   }
   return(bp);
}

static ATL_mvnode_t *FindFastestKernel
(  char pre,             /* precision prefix */
   ATL_mvnode_t *bp,  /* kernel queue */
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
   ATL_mvnode_t *kp, *kmax=bp;
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
static int MVflag2size(int flag)
/*
 * RETURNS: size of type using precision/type bits in flag
 */
{
   int size;

   size = FLAG_IS_SET(flag, MVF_SINGLE) ? 4 : 8;
   size *= FLAG_IS_SET(flag, MVF_COMPLEX) ? 2 : 1;
   return(size);
}

static char MVflag2pre(int flag)
/*
 * RETURNS: correct precision/type prefix based on flag
 */
{
   char pre = 'd';
   if (FLAG_IS_SET(flag, MVF_SINGLE))
      return(FLAG_IS_SET(flag, MVF_COMPLEX) ? 'c' : 's');
   return(FLAG_IS_SET(flag, MVF_COMPLEX) ? 'z' : 'd');
}
static int pre2MVflag(char pre, int flag)
/*
 * RETURNS: flag modified to reflect type/precision indicated by pre
 */
{
   SET_FLAG(flag, MVF_COMPLEX, (pre == 'c' || pre == 'z'));
   SET_FLAG(flag, MVF_SINGLE, (pre == 'c' || pre == 's'));
   return(flag);
}

static void SetAllMVTypeFlags(char pre, ATL_mvnode_t *bp)
{
   ATL_mvnode_t *p;
   for (p=bp; p; p = p->next)
      p->flag = pre2MVflag(pre, p->flag);
}
#endif  /* end atlas_mvparse.h guard */
void *SortByTrans
   (ATL_mvnode_t *bp,    /* original kernels wt mixture of trans cases */
    ATL_mvnode_t **bN0,  /* No trans cases */
    ATL_mvnode_t **bT0,  /* trans cases */
    ATL_mvnode_t **bNC0, /* ConjNotrans cases */
    ATL_mvnode_t **bTC0, /* Conjtrans cases */
/*
 * Sorts bp into the separate transpose queues, destroying bp in the process.
 * If a bp entry has the CONJDEF property, then its entry is duplicated to
 * put it on both queues (it can be used for normal and conjugate cases).
 */
{
   ATL_mvnode_t *bN=NULL, *bT=NULL, *bNC=NULL, *bTC=NULL,
                   *p, *next, *new;

   for (p=bp; p; p = next)
   {
      next = p->next;
      if (p->TA == AtlasNoTrans)
      {
         p->next = bN;
         bN = p;
         #ifdef TCPLX
            if (p->flag & (1<<MVF_CONJDEF))
            {
               new = CloneMVNode(p);
               new->TA = AtlasConj;
               new->next = bNC;
               bNC = new;
            }
         #endif
      }
      else if (p->TA == AtlasTrans)
      {
         p->next = bT;
         bT = p;
         #ifdef TCPLX
            if (p->flag & (1<<MVF_CONJDEF))
            {
               new = CloneMVNode(p);
               new->TA = AtlasConjTrans;
               new->next = bTC;
               bTC = new;
            }
         #endif
      }
   #ifdef TCPLX
      else if (p->TA == AtlasConjTrans)
      {
         p->next = bTC;
         bTC = p;
      }
      else /* TA == AtlasConj */
      {
         p->next = bNC;
         bNC = p;
      }
   #endif
   }
   *bN0 = bN;
   *bT0 = bT;
   *bNC0 = bNC;
   *bTC0 = bTC;
}
