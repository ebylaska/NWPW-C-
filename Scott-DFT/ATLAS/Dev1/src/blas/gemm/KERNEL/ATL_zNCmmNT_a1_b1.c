#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_zJIK36x36x36NT0x0x0_a1_b1
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=T, MB=36, NB=36, KB=36, 
 * lda=0, ldb=0, ldc=0, mu=2, nu=2, ku=36, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const double *stM = A + 72;
   const double *stN = B + 72;
   const int incAk = (((lda) << 1));
   const int incAm = 4 - (((lda) << 6)+((lda) << 3)), incAn = -72;
   const int incBk = (((ldb) << 1)), incBm = -(((ldb) << 6)+((ldb) << 3));
   #define incBn 4
   #define incAk0 incAk
   #define incBk0 incBk
   #define incCm 4
   const int incCn = (((ldc) << 2)) - 72;
   double *pC0=C, *pC1=pC0+(((ldc) << 1));
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1;
   register double rB0, rB1;
   register double m0;
   register double rC0_0, rC1_0, rC0_1, rC1_1;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rC0_0 = *pC0;
         rC1_0 = pC0[2];
         rC0_1 = *pC1;
         rC1_1 = pC1[2];
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         m0 = rA0 * rB0;
         rB1 = pB0[2];

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
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[2];
         rB1 = pB0[2];
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
         pB0 += incBk;
         rC1_1 += m0;
         *pC0 = rC0_0;
         pC0[2] = rC1_0;
         *pC1 = rC0_1;
         pC1[2] = rC1_1;
         pC0 += incCm;
         pC1 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pA0 += incAn;
      pB0 += incBn;
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
