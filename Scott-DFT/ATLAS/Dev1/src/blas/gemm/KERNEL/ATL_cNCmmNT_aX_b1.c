#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_cJIK44x44x44NT0x0x0_aX_b1
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=T, MB=44, NB=44, KB=44, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=44, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + 88;
   const float *stN = B + 88;
   const float BetaAlpha = beta / alpha;
   const int incAk = (((lda) << 1));
   const int incAm = 8 - (88*(lda)), incAn = -88;
   const int incBk = (((ldb) << 1)), incBm = -(88*(ldb));
   #define incBn 4
   #define incCm 8
   const int incCn = (((ldc) << 2)) - 88;
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3;
   register float rB0, rB1;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC0_1, rC1_1, rC2_1, rC3_1;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[2];
         rC1_0 *= rA0;
         rC2_0 = pC0[4];
         rC2_0 *= rA0;
         rC3_0 = pC0[6];
         rC3_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[2];
         rC1_1 *= rA0;
         rC2_1 = pC1[4];
         rC2_1 *= rA0;
         rC3_1 = pC1[6];
         rC3_1 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB0[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         rB0 = alpha;
         rC0_0 *= rB0;
         rC0_1 *= rB0;
         rC1_0 *= rB0;
         rC1_1 *= rB0;
         rC2_0 *= rB0;
         rC2_1 *= rB0;
         rC3_0 *= rB0;
         rC3_1 *= rB0;
         *pC0 = rC0_0;
         pC0[2] = rC1_0;
         pC0[4] = rC2_0;
         pC0[6] = rC3_0;
         *pC1 = rC0_1;
         pC1[2] = rC1_1;
         pC1[4] = rC2_1;
         pC1[6] = rC3_1;
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
