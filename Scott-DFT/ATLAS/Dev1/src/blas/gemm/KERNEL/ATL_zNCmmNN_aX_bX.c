#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_zJIK36x36x36NN0x0x0_aX_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=N, MB=36, NB=36, KB=36, 
 * lda=0, ldb=0, ldc=0, mu=2, nu=2, ku=36, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const double *stM = A + 72;
   const double *stN = B + (((ldb) << 6)+((ldb) << 3));
   const double BetaAlpha = beta / alpha;
   const int incAk = (((lda) << 1));
   const int incAm = 4 - (((lda) << 6)+((lda) << 3)), incAn = -72;
   #define incBk 72
   const int incBm = -72, incBn = (((ldb) << 2));
   #define incAk0 incAk
   const int incBk0 = ((incBk) / 36);
   #define incCm 4
   const int incCn = (((ldc) << 2)) - 72;
   double *pC0=C, *pC1=pC0+(((ldc) << 1));
   const double *pA0=A;
   const double *pB0=B, *pB1=pB0+(((ldb) << 1));
   register int k;
   register double rA0, rA1;
   register double rB0, rB1;
   register double m0;
   register double rC0_0, rC1_0, rC0_1, rC1_1;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[2];
         rC1_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[2];
         rC1_1 *= rA0;
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         m0 = rA0 * rB0;
         rB1 = *pB1;

/*
 *       Completely unrolled K-loop
 */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[2];
         rA1 = pA0[2];
         rB1 = pB1[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[4];
         rA1 = pA0[2];
         rB1 = pB1[4];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[6];
         rA1 = pA0[2];
         rB1 = pB1[6];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[8];
         rA1 = pA0[2];
         rB1 = pB1[8];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[10];
         rA1 = pA0[2];
         rB1 = pB1[10];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[12];
         rA1 = pA0[2];
         rB1 = pB1[12];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[14];
         rA1 = pA0[2];
         rB1 = pB1[14];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[16];
         rA1 = pA0[2];
         rB1 = pB1[16];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[18];
         rA1 = pA0[2];
         rB1 = pB1[18];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[20];
         rA1 = pA0[2];
         rB1 = pB1[20];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[22];
         rA1 = pA0[2];
         rB1 = pB1[22];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[24];
         rA1 = pA0[2];
         rB1 = pB1[24];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[26];
         rA1 = pA0[2];
         rB1 = pB1[26];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[28];
         rA1 = pA0[2];
         rB1 = pB1[28];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[30];
         rA1 = pA0[2];
         rB1 = pB1[30];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[32];
         rA1 = pA0[2];
         rB1 = pB1[32];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[34];
         rA1 = pA0[2];
         rB1 = pB1[34];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[36];
         rA1 = pA0[2];
         rB1 = pB1[36];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[38];
         rA1 = pA0[2];
         rB1 = pB1[38];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[40];
         rA1 = pA0[2];
         rB1 = pB1[40];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[42];
         rA1 = pA0[2];
         rB1 = pB1[42];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[44];
         rA1 = pA0[2];
         rB1 = pB1[44];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[46];
         rA1 = pA0[2];
         rB1 = pB1[46];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[48];
         rA1 = pA0[2];
         rB1 = pB1[48];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[50];
         rA1 = pA0[2];
         rB1 = pB1[50];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[52];
         rA1 = pA0[2];
         rB1 = pB1[52];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[54];
         rA1 = pA0[2];
         rB1 = pB1[54];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[56];
         rA1 = pA0[2];
         rB1 = pB1[56];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[58];
         rA1 = pA0[2];
         rB1 = pB1[58];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[60];
         rA1 = pA0[2];
         rB1 = pB1[60];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[62];
         rA1 = pA0[2];
         rB1 = pB1[62];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[64];
         rA1 = pA0[2];
         rB1 = pB1[64];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[66];
         rA1 = pA0[2];
         rB1 = pB1[66];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[68];
         rA1 = pA0[2];
         rB1 = pB1[68];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[70];
         rA1 = pA0[2];
         rB1 = pB1[70];
         rC1_1 += m0;
         m0 = rA0 * rB0;
/*
 *       Drain pipe on last iteration of K-loop
 */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rC1_1 += m0;
         pB0 += incBk;
         pB1 += incBk;
         rB0 = alpha;
         rC0_0 *= rB0;
         rC0_1 *= rB0;
         rC1_0 *= rB0;
         rC1_1 *= rB0;
         *pC0 = rC0_0;
         pC0[2] = rC1_0;
         *pC1 = rC0_1;
         pC1[2] = rC1_1;
         pC0 += incCm;
         pC1 += incCm;
         pA0 += incAm;
         pB0 += incBm;
         pB1 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pA0 += incAn;
      pB0 += incBn;
      pB1 += incBn;
   }
   while(pB0 != stN);
}
#ifdef incAm
   #undef incAm
#endif
#ifdef incAn
   #undef incAn
#endif
#ifdef incAk
   #undef incAk
#endif
#ifdef incBm
   #undef incBm
#endif
#ifdef incBn
   #undef incBn
#endif
#ifdef incBk
   #undef incBk
#endif
#ifdef incCm
   #undef incCm
#endif
#ifdef incCn
   #undef incCn
#endif
#ifdef incCk
   #undef incCk
#endif
#ifdef Mb
   #undef Mb
#endif
#ifdef Nb
   #undef Nb
#endif
#ifdef Kb
   #undef Kb
#endif
