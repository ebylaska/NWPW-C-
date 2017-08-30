#ifndef ATLAS_MMPARSE_H
   #define ATLAS_MMPARSE_H

#include "atlas_genparse.h"
#include "atlas_enum.h"

#define MMF_LDCTOP      0  /* 1: load C before K-loop, 0: ld after */
#define MMF_X87         1  /* 1: requires the Intel x87 unit */
#define MMF_MRUNTIME    2  /* 1: M dim is run-time variable */
#define MMF_NRUNTIME    3  /* 1: M dim is run-time variable */
#define MMF_KRUNTIME    4  /* 1: M dim is run-time variable */
#define MMF_KUISKB      5  /* 1: KU == KB */
#define MMF_LDISKB      6  /* 1: lda=ldb=KB */
#define MMF_BETAN1      7  /* 1: kernel has special support for BETA = -1 */
#define MMF_LDAB        8  /* 1: lda = ldb */
#define MMF_AOUTER      9  /* 1: MNK loop order, 0: NMK loop order */
#define MMF_LDFLOAT    10  /* 1: use single prec load for double */
#define MMF_STFLOAT    11  /* 1: use single prec load for double */
#define MMF_PFACOLS    12  /* 1: prefetch next mu cols of A */
#define MMF_PFABLK     13  /* 1: prefetch next KBxNB block of A */
#define MMF_PFBCOLS    14  /* 1: prefetch next nu cols of B */
#define MMF_PFCELTS    15  /* 1: pf elts of C at top of loop, load at bottom */
#define MMF_SINGLE     16  /* 1: single precision, else double */
#define MMF_COMPLEX    17  /* 1: complex type, else real */
#define MMF_L14NB      18  /* 1: need to fit all 3 matrices+nextA in L1 */

#define MMF_DEFAULT ( (1<<MMF_LDISKB) | (1<<MMF_LDAB) )
#ifndef  FLAG_IS_SET
   #define FLAG_IS_SET(field_, bit_) ( ((field_) & (1<<(bit_))) != 0 )
#endif

typedef struct MMNode ATL_mmnode_t;
struct MMNode
{
   double mflop[8];     /* 1st entry perf using mbB, nbB, kbB */
   int ID, mu, nu, ku;  /* ID, and unrolling on each loop */
   int kbmin, kbmax;    /* min/max KB this kernel can handle */
   int SSE;             /* 0: no SSE, 1: SSE1 req, 2: SSE2 req, etc */
   int lat, muladd, pref, clean, fftch, iftch, nftch; /* used for gened codes */
   int mbB, nbB, kbB;  /* best blocking dims found by search */
   enum ATLAS_TRANS TA, TB;
   int asmbits;   /* bitfield indicating which assembly(ies) is required */
   char *rout, *auth, *comp, *cflags;
   char *str;                   /* tmp string used in generation */
   char *genstr;                /* system(genstr) will generate gened kernel */
   int flag;
   ATL_mmnode_t *next;
};

static ATL_mmnode_t *GetMMNode(void)
{
   ATL_mmnode_t *p;
   p = calloc(1, sizeof(ATL_mmnode_t));
   assert(p);
   p->TA = AtlasTrans; p->TB = AtlasNoTrans;
   p->flag = MMF_DEFAULT;
   return(p);
}

static ATL_mmnode_t *CloneMMNode(ATL_mmnode_t *dup)
{
   ATL_mmnode_t *p;
   p = malloc(sizeof(ATL_mmnode_t));
   assert(p);
   memcpy(p, dup, sizeof(ATL_mmnode_t));
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

static ATL_mmnode_t *KillMMNode(ATL_mmnode_t *die)
{
   ATL_mmnode_t *p=NULL;
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

static void KillAllMMNodes(ATL_mmnode_t *die)
{
   while (die)
      die = KillMMNode(die);
}

static ATL_mmnode_t *ParseMMLine(char *ln)
/*
 * Given a line from a mm index file (with multiple lines pasted together
 * into one line (ln), return a structure describing that line.
 */
{
   ATL_mmnode_t *p;
   char *sp;
   int itmp;
   char ch;

   p = GetMMNode();

   sp = strstr(ln, "MULADD=");
   if (sp)
      p->muladd = atoi(sp+6+1);
   else
      p->muladd = 0;
   sp = strstr(ln, "PREF=");
   if (sp)
      p->pref = atoi(sp+4+1);
   else
      p->pref = 0;
   sp = strstr(ln, "LAT=");
   if (sp)
      p->lat = atoi(sp+3+1);
   else
      p->lat = 0;
   sp = strstr(ln, "NFTCH=");
   if (sp)
      p->nftch = atoi(sp+5+1);
   else
      p->nftch = 0;
   sp = strstr(ln, "IFTCH=");
   if (sp)
      p->iftch = atoi(sp+5+1);
   else
      p->iftch = 0;
   sp = strstr(ln, "FFTCH=");
   if (sp)
      p->fftch = atoi(sp+5+1);
   else
      p->fftch = 0;
   sp = strstr(ln, "KBMAX=");
   if (sp)
      p->kbmax = atoi(sp+5+1);
   else
      p->kbmax = 0;
   sp = strstr(ln, "KBMIN=");
   if (sp)
      p->kbmin = atoi(sp+5+1);
   else
      p->kbmin = 0;
   sp = strstr(ln, "KU=");
   if (sp)
      p->ku = atoi(sp+2+1);
   else
      p->ku = 0;
   sp = strstr(ln, "NU=");
   if (sp)
      p->nu = atoi(sp+2+1);
   else
      p->nu = 0;
   sp = strstr(ln, "MU=");
   if (sp)
      p->mu = atoi(sp+2+1);
   else
      p->mu = 0;
   sp = strstr(ln, "MB=");
   if (sp)
      p->mbB = atoi(sp+2+1);
   else
      p->mbB = 0;
   sp = strstr(ln, "NB=");
   if (sp)
      p->nbB = atoi(sp+2+1);
   else
      p->nbB = 0;
   sp = strstr(ln, "KB=");
   if (sp)
      p->kbB = atoi(sp+2+1);
   else
      p->kbB = 0;
   sp = strstr(ln, "SSE=");
   if (sp)
      p->SSE = atoi(sp+3+1);
   else
      p->SSE = 0;

   sp = strstr(ln, "ID=");
   if (sp)
      p->ID = atoi(sp+2+1);
   else
      p->ID = 0;


   sp = strstr(ln, "L14NB=");
   if (sp)
   {
      if (atoi(sp+5+1))
         p->flag |= (1<<MMF_L14NB);
      else
         p->flag &= ~(1<<MMF_L14NB);
   }
   sp = strstr(ln, "PFCELTS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFCELTS);
      else
         p->flag &= ~(1<<MMF_PFCELTS);
   }
   sp = strstr(ln, "PFBCOLS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFBCOLS);
      else
         p->flag &= ~(1<<MMF_PFBCOLS);
   }
   sp = strstr(ln, "PFABLK=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_PFABLK);
      else
         p->flag &= ~(1<<MMF_PFABLK);
   }
   sp = strstr(ln, "PFACOLS=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_PFACOLS);
      else
         p->flag &= ~(1<<MMF_PFACOLS);
   }
   sp = strstr(ln, "STFLOAT=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_STFLOAT);
      else
         p->flag &= ~(1<<MMF_STFLOAT);
   }
   sp = strstr(ln, "LDFLOAT=");
   if (sp)
   {
      if (atoi(sp+7+1))
         p->flag |= (1<<MMF_LDFLOAT);
      else
         p->flag &= ~(1<<MMF_LDFLOAT);
   }
   sp = strstr(ln, "AOUTER=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_AOUTER);
      else
         p->flag &= ~(1<<MMF_AOUTER);
   }
   sp = strstr(ln, "LDAB=");
   if (sp)
   {
      if (atoi(sp+4+1))
         p->flag |= (1<<MMF_LDAB);
      else
         p->flag &= ~(1<<MMF_LDAB);
   }
   sp = strstr(ln, "BETAN1=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_BETAN1);
      else
         p->flag &= ~(1<<MMF_BETAN1);
   }
   sp = strstr(ln, "LDISKB=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_LDISKB);
      else
         p->flag &= ~(1<<MMF_LDISKB);
   }
   sp = strstr(ln, "KUISKB=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_KUISKB);
      else
         p->flag &= ~(1<<MMF_KUISKB);
   }
   sp = strstr(ln, "KRUNTIME=");
   if (sp)
   {
      if (atoi(sp+8+1))
         p->flag |= (1<<MMF_KRUNTIME);
      else
         p->flag &= ~(1<<MMF_KRUNTIME);
   }
   sp = strstr(ln, "NRUNTIME=");
   if (sp)
   {
      if (atoi(sp+8+1))
         p->flag |= (1<<MMF_NRUNTIME);
      else
         p->flag &= ~(1<<MMF_NRUNTIME);
   }
   sp = strstr(ln, "MRUNTIME=");
   if (sp)
   {
      if (atoi(sp+8+1))
         p->flag |= (1<<MMF_MRUNTIME);
      else
         p->flag &= ~(1<<MMF_MRUNTIME);
   }
   sp = strstr(ln, "LDCTOP=");
   if (sp)
   {
      if (atoi(sp+6+1))
         p->flag |= (1<<MMF_LDCTOP);
      else
         p->flag &= ~(1<<MMF_LDCTOP);
   }
   sp = strstr(ln, "X87=");
   if (sp)
   {
      if (atoi(sp+3+1))
         p->flag |= (1<<MMF_X87);
      else
         p->flag &= ~(1<<MMF_X87);
   }

   sp = strstr(ln, "MFLOP=");
   if (sp)
      GetDoubleArr(sp+6, 8, p->mflop);

   sp = strstr(ln, "ASM=");
   if (sp)
      p->asmbits = asmNames2bitfield(sp+4);


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

static void PrintMMLine(FILE *fpout, ATL_mmnode_t *np)
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
   if (np->TB == AtlasConjTrans) tb = 'C';
   else if (np->TB == AtlasTrans) tb = 'T';
   else tb = 'N';
   fprintf(fpout, "ID=%d ROUT='%s' AUTH='%s' TA='%c' TB='%c' \\\n",
           np->ID, np->rout ? np->rout : "generated",
           np->auth ? np->auth : "R. Clint Whaley", ta, tb);

   fprintf(fpout, "   ");
   i = 3;
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "MULADD=%d ", np->muladd);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "PREF=%d ", np->pref);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LAT=%d ", np->lat);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "NFTCH=%d ", np->nftch);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "IFTCH=%d ", np->iftch);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "FFTCH=%d ", np->fftch);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "KBMAX=%d ", np->kbmax);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "KBMIN=%d ", np->kbmin);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "KU=%d ", np->ku);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "NU=%d ", np->nu);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "MU=%d ", np->mu);

   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   if (np->mbB != 0)
      i += fprintf(fpout, "MB=%d ", np->mbB);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   if (np->nbB != 0)
      i += fprintf(fpout, "NB=%d ", np->nbB);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   if (np->kbB != 0)
      i += fprintf(fpout, "KB=%d ", np->kbB);
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "L14NB=%d ", FLAG_IS_SET(np->flag, MMF_L14NB));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "PFBCOLS=%d ", FLAG_IS_SET(np->flag, MMF_PFBCOLS));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "PFABLK=%d ", FLAG_IS_SET(np->flag, MMF_PFABLK));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "PFACOLS=%d ", FLAG_IS_SET(np->flag, MMF_PFACOLS));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "STFLOAT=%d ", FLAG_IS_SET(np->flag, MMF_STFLOAT));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDFLOAT=%d ", FLAG_IS_SET(np->flag, MMF_LDFLOAT));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "AOUTER=%d ", FLAG_IS_SET(np->flag, MMF_AOUTER));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDAB=%d ", FLAG_IS_SET(np->flag, MMF_LDAB));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "BETAN1=%d ", FLAG_IS_SET(np->flag, MMF_BETAN1));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDISKB=%d ", FLAG_IS_SET(np->flag, MMF_LDISKB));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "KUISKB=%d ", FLAG_IS_SET(np->flag, MMF_KUISKB));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "KRUNTIME=%d ", FLAG_IS_SET(np->flag, MMF_KRUNTIME));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "NRUNTIME=%d ", FLAG_IS_SET(np->flag, MMF_NRUNTIME));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "MRUNTIME=%d ", FLAG_IS_SET(np->flag, MMF_MRUNTIME));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "LDCTOP=%d ", FLAG_IS_SET(np->flag, MMF_LDCTOP));
   if (i > 70) { fprintf(fpout, "\\\n   "); i = 3; }
   i += fprintf(fpout, "X87=%d ", FLAG_IS_SET(np->flag, MMF_X87));

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

void PrintMMNodes(FILE *fpout, ATL_mmnode_t *bp)
{
   while (bp)
   {
      PrintMMLine(fpout, bp);
      bp = bp->next;
   }
}

static void WriteMMFile(char *file, ATL_mmnode_t *nq)
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
   PrintMMNodes(fpout, nq);
   if (fpout != stdout && fpout != stderr)
      fclose(fpout);
}

static void WriteMMFileWithPath
   (char pre, char *path, char *file, ATL_mmnode_t *nq)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   WriteMMFile(ln, nq);
}

static ATL_mmnode_t *ReadMMFile(char *file)
/*
 * Reads in a standard ATLAS parsable MM index file, and returns a
 * list of all the kernels defined there.
 */
{
   ATL_mmnode_t *nq=NULL, *p;
   FILE *fpin;
   char *ln, *sp;
   int i, j, KeepOn, len;

   if (!file || !strcmp(file, "stdin"))
      fpin = stdin;
   else
      fpin = fopen(file, "r");
   if (!fpin)
      return(NULL);
   nq = p = GetMMNode();
   while (ln = GetJoinedLines(fpin))
   {
      if (ln[0] != '#')
      {
         p->next = ParseMMLine(ln);
         p = p->next;
      }
   }
   fclose(fpin);
   return(KillMMNode(nq));
}

static ATL_mmnode_t *ReadMMFileWithPath
   (char pre, char *path, char *file)
{
   char ln[2048];
   sprintf(ln, "%s/%c%s", path, pre, file);
   return(ReadMMFile(ln));
}
static ATL_mmnode_t *DelBadArchMMKernels(ATL_mmnode_t *bp)
/*
 * Weeds out kernels that require SSE/assembly that we haven't got
 */
{
   int asmb=0, die;
   ATL_mmnode_t *p, *prev;
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
            bp = p = KillMMNode(p);
         else
            prev->next = p = KillMMNode(p);
      }
      else
      {
         prev = p;
         p = p->next;
      }
   }
   return(bp);
}

#endif  /* end atlas_mmparse.h guard */
