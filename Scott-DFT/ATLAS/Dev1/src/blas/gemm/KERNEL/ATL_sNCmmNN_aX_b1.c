#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_sJIK64x64x64NN0x0x0_aX_b1
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=N, MB=64, NB=64, KB=64, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=64, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + 64;
   const float *stN = B + (((ldb) << 6));
   const float BetaAlpha = beta / alpha;
   const int incAk = (lda);
   const int incAm = 4 - (((lda) << 6)), incAn = -64;
   #define incBk 64
   const int incBm = -64, incBn = (((ldb) << 1));
   #define incCm 4
   const int incCn = (((ldc) << 1)) - 64;
   float *pC0=C, *pC1=pC0+(ldc);
   const float *pA0=A;
   const float *pB0=B, *pB1=pB0+(ldb);
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
         rC1_0 = pC0[1];
         rC1_0 *= rA0;
         rC2_0 = pC0[2];
         rC2_0 *= rA0;
         rC3_0 = pC0[3];
         rC3_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[1];
         rC1_1 *= rA0;
         rC2_1 = pC1[2];
         rC2_1 *= rA0;
         rC3_1 = pC1[3];
         rC3_1 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = *pB1;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[1];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[1];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[2];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[2];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[3];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[3];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[4];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[4];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[5];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[5];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[6];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[6];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[7];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[7];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[8];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[8];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[9];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[9];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[10];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[10];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[11];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[11];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[12];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[12];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[13];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[13];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[14];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[14];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[15];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[15];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[16];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[16];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[17];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[17];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[18];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[18];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[19];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[19];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[20];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[20];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[21];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[21];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[22];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[22];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[23];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[23];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[24];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[24];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[25];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[25];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[26];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[26];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[27];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[27];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[28];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[28];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[29];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[29];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[30];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[30];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[31];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[31];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[32];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[32];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[33];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[33];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[34];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[34];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[35];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[35];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[36];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[36];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[37];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[37];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[38];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[38];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[39];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[39];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[40];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[40];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[41];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[41];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[42];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[42];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[43];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[43];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[44];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[44];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[45];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[45];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[46];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[46];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[47];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[47];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[48];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[48];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[49];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[49];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[50];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[50];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[51];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[51];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[52];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[52];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[53];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[53];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[54];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[54];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[55];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[55];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[56];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[56];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[57];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[57];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[58];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[58];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[59];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[59];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[60];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[60];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[61];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[61];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[62];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[62];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[63];
         rA1 = pA0[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA0[2];
         rA3 = pA0[3];
         rB1 = pB1[63];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         pB1 += incBk;
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
