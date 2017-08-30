#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_sJIK64x64x64NT0x0x0_a1_b1
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=T, MB=64, NB=64, KB=64, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=64, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + 64;
   const float *stN = B + 64;
   const int incAk = (lda);
   const int incAm = 4 - (((lda) << 6)), incAn = -64;
   const int incBk = (ldb), incBm = -(((ldb) << 6));
   #define incBn 2
   #define incCm 4
   const int incCn = (((ldc) << 1)) - 64;
   float *pC0=C, *pC1=pC0+(ldc);
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
         rC0_0 = *pC0;
         rC1_0 = pC0[1];
         rC2_0 = pC0[2];
         rC3_0 = pC0[3];
         rC0_1 = *pC1;
         rC1_1 = pC1[1];
         rC2_1 = pC1[2];
         rC3_1 = pC1[3];
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
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
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB0[1];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         pC0[2] = rC2_0;
         pC0[3] = rC3_0;
         *pC1 = rC0_1;
         pC1[1] = rC1_1;
         pC1[2] = rC2_1;
         pC1[3] = rC3_1;
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