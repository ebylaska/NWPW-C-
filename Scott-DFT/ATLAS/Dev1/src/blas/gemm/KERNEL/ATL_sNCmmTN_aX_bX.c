#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_sJIK64x64x64TN0x0x0_aX_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=64, NB=64, KB=64, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=64, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + (((lda) << 6));
   const float *stN = B + (((ldb) << 6));
   const float BetaAlpha = beta / alpha;
   #define incAk 64
   const int incAm = ((((lda) << 2)) - 64), incAn = -(((lda) << 6));
   #define incBk 64
   const int incBm = -64, incBn = (((ldb) << 1));
   #define incCm 4
   const int incCn = (((ldc) << 1)) - 64;
   float *pC0=C, *pC1=pC0+(ldc);
   const float *pA0=A, *pA1=pA0+(lda), *pA2=pA1+(lda), *pA3=pA2+(lda);
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
         rA1 = *pA1;
         rC0_0 += rA0 * rB0;
         rA2 = *pA2;
         rA3 = *pA3;
         rB1 = *pB1;
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rA1 = pA1[1];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[1];
         rA3 = pA3[1];
         rB1 = pB1[1];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA1[2];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[2];
         rA3 = pA3[2];
         rB1 = pB1[2];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rA1 = pA1[3];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[3];
         rA3 = pA3[3];
         rB1 = pB1[3];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA1[4];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[4];
         rA3 = pA3[4];
         rB1 = pB1[4];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rA1 = pA1[5];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[5];
         rA3 = pA3[5];
         rB1 = pB1[5];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA1[6];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[6];
         rA3 = pA3[6];
         rB1 = pB1[6];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rA1 = pA1[7];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[7];
         rA3 = pA3[7];
         rB1 = pB1[7];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA1[8];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[8];
         rA3 = pA3[8];
         rB1 = pB1[8];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rA1 = pA1[9];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[9];
         rA3 = pA3[9];
         rB1 = pB1[9];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA1[10];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[10];
         rA3 = pA3[10];
         rB1 = pB1[10];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rA1 = pA1[11];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[11];
         rA3 = pA3[11];
         rB1 = pB1[11];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA1[12];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[12];
         rA3 = pA3[12];
         rB1 = pB1[12];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rA1 = pA1[13];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[13];
         rA3 = pA3[13];
         rB1 = pB1[13];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA1[14];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[14];
         rA3 = pA3[14];
         rB1 = pB1[14];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rA1 = pA1[15];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[15];
         rA3 = pA3[15];
         rB1 = pB1[15];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA1[16];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[16];
         rA3 = pA3[16];
         rB1 = pB1[16];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rA1 = pA1[17];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[17];
         rA3 = pA3[17];
         rB1 = pB1[17];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA1[18];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[18];
         rA3 = pA3[18];
         rB1 = pB1[18];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rA1 = pA1[19];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[19];
         rA3 = pA3[19];
         rB1 = pB1[19];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA1[20];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[20];
         rA3 = pA3[20];
         rB1 = pB1[20];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rA1 = pA1[21];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[21];
         rA3 = pA3[21];
         rB1 = pB1[21];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA1[22];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[22];
         rA3 = pA3[22];
         rB1 = pB1[22];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rA1 = pA1[23];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[23];
         rA3 = pA3[23];
         rB1 = pB1[23];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA1[24];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[24];
         rA3 = pA3[24];
         rB1 = pB1[24];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rA1 = pA1[25];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[25];
         rA3 = pA3[25];
         rB1 = pB1[25];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA1[26];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[26];
         rA3 = pA3[26];
         rB1 = pB1[26];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rA1 = pA1[27];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[27];
         rA3 = pA3[27];
         rB1 = pB1[27];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA1[28];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[28];
         rA3 = pA3[28];
         rB1 = pB1[28];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[29];
         rB0 = pB0[29];
         rA1 = pA1[29];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[29];
         rA3 = pA3[29];
         rB1 = pB1[29];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA1[30];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[30];
         rA3 = pA3[30];
         rB1 = pB1[30];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[31];
         rB0 = pB0[31];
         rA1 = pA1[31];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[31];
         rA3 = pA3[31];
         rB1 = pB1[31];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA1[32];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[32];
         rA3 = pA3[32];
         rB1 = pB1[32];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[33];
         rB0 = pB0[33];
         rA1 = pA1[33];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[33];
         rA3 = pA3[33];
         rB1 = pB1[33];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA1[34];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[34];
         rA3 = pA3[34];
         rB1 = pB1[34];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[35];
         rB0 = pB0[35];
         rA1 = pA1[35];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[35];
         rA3 = pA3[35];
         rB1 = pB1[35];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA1[36];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[36];
         rA3 = pA3[36];
         rB1 = pB1[36];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[37];
         rB0 = pB0[37];
         rA1 = pA1[37];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[37];
         rA3 = pA3[37];
         rB1 = pB1[37];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA1[38];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[38];
         rA3 = pA3[38];
         rB1 = pB1[38];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[39];
         rB0 = pB0[39];
         rA1 = pA1[39];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[39];
         rA3 = pA3[39];
         rB1 = pB1[39];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA1[40];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[40];
         rA3 = pA3[40];
         rB1 = pB1[40];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[41];
         rB0 = pB0[41];
         rA1 = pA1[41];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[41];
         rA3 = pA3[41];
         rB1 = pB1[41];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA1[42];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[42];
         rA3 = pA3[42];
         rB1 = pB1[42];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[43];
         rB0 = pB0[43];
         rA1 = pA1[43];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[43];
         rA3 = pA3[43];
         rB1 = pB1[43];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA1[44];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[44];
         rA3 = pA3[44];
         rB1 = pB1[44];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[45];
         rB0 = pB0[45];
         rA1 = pA1[45];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[45];
         rA3 = pA3[45];
         rB1 = pB1[45];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA1[46];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[46];
         rA3 = pA3[46];
         rB1 = pB1[46];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[47];
         rB0 = pB0[47];
         rA1 = pA1[47];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[47];
         rA3 = pA3[47];
         rB1 = pB1[47];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[48];
         rB0 = pB0[48];
         rA1 = pA1[48];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[48];
         rA3 = pA3[48];
         rB1 = pB1[48];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[49];
         rB0 = pB0[49];
         rA1 = pA1[49];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[49];
         rA3 = pA3[49];
         rB1 = pB1[49];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[50];
         rB0 = pB0[50];
         rA1 = pA1[50];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[50];
         rA3 = pA3[50];
         rB1 = pB1[50];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[51];
         rB0 = pB0[51];
         rA1 = pA1[51];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[51];
         rA3 = pA3[51];
         rB1 = pB1[51];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[52];
         rB0 = pB0[52];
         rA1 = pA1[52];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[52];
         rA3 = pA3[52];
         rB1 = pB1[52];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[53];
         rB0 = pB0[53];
         rA1 = pA1[53];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[53];
         rA3 = pA3[53];
         rB1 = pB1[53];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[54];
         rB0 = pB0[54];
         rA1 = pA1[54];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[54];
         rA3 = pA3[54];
         rB1 = pB1[54];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[55];
         rB0 = pB0[55];
         rA1 = pA1[55];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[55];
         rA3 = pA3[55];
         rB1 = pB1[55];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[56];
         rB0 = pB0[56];
         rA1 = pA1[56];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[56];
         rA3 = pA3[56];
         rB1 = pB1[56];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[57];
         rB0 = pB0[57];
         rA1 = pA1[57];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[57];
         rA3 = pA3[57];
         rB1 = pB1[57];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[58];
         rB0 = pB0[58];
         rA1 = pA1[58];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[58];
         rA3 = pA3[58];
         rB1 = pB1[58];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[59];
         rB0 = pB0[59];
         rA1 = pA1[59];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[59];
         rA3 = pA3[59];
         rB1 = pB1[59];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[60];
         rB0 = pB0[60];
         rA1 = pA1[60];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[60];
         rA3 = pA3[60];
         rB1 = pB1[60];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[61];
         rB0 = pB0[61];
         rA1 = pA1[61];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[61];
         rA3 = pA3[61];
         rB1 = pB1[61];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[62];
         rB0 = pB0[62];
         rA1 = pA1[62];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[62];
         rA3 = pA3[62];
         rB1 = pB1[62];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[63];
         rB0 = pB0[63];
         rA1 = pA1[63];
         rC0_0 += rA0 * rB0;
         rA2 = pA2[63];
         rA3 = pA3[63];
         rB1 = pB1[63];
         rC1_0 += rA1 * rB0;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pA1 += incAk;
         pA2 += incAk;
         pA3 += incAk;
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
         pA1 += incAm;
         pA2 += incAm;
         pA3 += incAm;
         pB0 += incBm;
         pB1 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pA0 += incAn;
      pA1 += incAn;
      pA2 += incAn;
      pA3 += incAn;
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
