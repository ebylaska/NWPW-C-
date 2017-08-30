#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_cJIK44x44x44TN0x0x0_aX_b1
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=44, NB=44, KB=44, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=44, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + (88*(lda));
   const float *stN = B + (88*(ldb));
   const float BetaAlpha = beta / alpha;
   #define incAk 88
   const int incAm = ((((((lda) << 2)) - 44) << 1)), incAn = -(88*(lda));
   #define incBk 88
   const int incBm = -88, incBn = (((ldb) << 2));
   #define incCm 8
   const int incCn = (((ldc) << 2)) - 88;
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
   const float *pA0=A, *pA1=pA0+(((lda) << 1)), *pA2=pA1+(((lda) << 1)), *pA3=pA2+(((lda) << 1));
   const float *pB0=B, *pB1=pB0+(((ldb) << 1));
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
         rA1 = *pA1;
         rA2 = *pA2;
         rC1_0 += rA1 * rB0;
         rA3 = *pA3;
         rB1 = *pB1;
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[2];
         rA2 = pA2[2];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[2];
         rB1 = pB1[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[4];
         rA2 = pA2[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[4];
         rB1 = pB1[4];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[6];
         rA2 = pA2[6];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[6];
         rB1 = pB1[6];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[8];
         rA2 = pA2[8];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[8];
         rB1 = pB1[8];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[10];
         rA2 = pA2[10];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[10];
         rB1 = pB1[10];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[12];
         rA2 = pA2[12];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[12];
         rB1 = pB1[12];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[14];
         rA2 = pA2[14];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[14];
         rB1 = pB1[14];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[16];
         rA2 = pA2[16];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[16];
         rB1 = pB1[16];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[18];
         rA2 = pA2[18];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[18];
         rB1 = pB1[18];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[20];
         rA2 = pA2[20];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[20];
         rB1 = pB1[20];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[22];
         rA2 = pA2[22];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[22];
         rB1 = pB1[22];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[24];
         rA2 = pA2[24];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[24];
         rB1 = pB1[24];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[26];
         rA2 = pA2[26];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[26];
         rB1 = pB1[26];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[28];
         rA2 = pA2[28];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[28];
         rB1 = pB1[28];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[30];
         rA2 = pA2[30];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[30];
         rB1 = pB1[30];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[32];
         rA2 = pA2[32];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[32];
         rB1 = pB1[32];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[34];
         rA2 = pA2[34];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[34];
         rB1 = pB1[34];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[36];
         rA2 = pA2[36];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[36];
         rB1 = pB1[36];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[38];
         rA2 = pA2[38];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[38];
         rB1 = pB1[38];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[40];
         rA2 = pA2[40];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[40];
         rB1 = pB1[40];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[42];
         rA2 = pA2[42];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[42];
         rB1 = pB1[42];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[44];
         rA2 = pA2[44];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[44];
         rB1 = pB1[44];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[46];
         rA2 = pA2[46];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[46];
         rB1 = pB1[46];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[48];
         rB0 = pB0[48];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[48];
         rA2 = pA2[48];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[48];
         rB1 = pB1[48];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[50];
         rB0 = pB0[50];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[50];
         rA2 = pA2[50];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[50];
         rB1 = pB1[50];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[52];
         rB0 = pB0[52];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[52];
         rA2 = pA2[52];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[52];
         rB1 = pB1[52];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[54];
         rB0 = pB0[54];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[54];
         rA2 = pA2[54];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[54];
         rB1 = pB1[54];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[56];
         rB0 = pB0[56];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[56];
         rA2 = pA2[56];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[56];
         rB1 = pB1[56];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[58];
         rB0 = pB0[58];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[58];
         rA2 = pA2[58];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[58];
         rB1 = pB1[58];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[60];
         rB0 = pB0[60];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[60];
         rA2 = pA2[60];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[60];
         rB1 = pB1[60];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[62];
         rB0 = pB0[62];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[62];
         rA2 = pA2[62];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[62];
         rB1 = pB1[62];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[64];
         rB0 = pB0[64];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[64];
         rA2 = pA2[64];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[64];
         rB1 = pB1[64];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[66];
         rB0 = pB0[66];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[66];
         rA2 = pA2[66];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[66];
         rB1 = pB1[66];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[68];
         rB0 = pB0[68];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[68];
         rA2 = pA2[68];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[68];
         rB1 = pB1[68];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[70];
         rB0 = pB0[70];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[70];
         rA2 = pA2[70];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[70];
         rB1 = pB1[70];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[72];
         rB0 = pB0[72];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[72];
         rA2 = pA2[72];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[72];
         rB1 = pB1[72];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[74];
         rB0 = pB0[74];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[74];
         rA2 = pA2[74];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[74];
         rB1 = pB1[74];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[76];
         rB0 = pB0[76];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[76];
         rA2 = pA2[76];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[76];
         rB1 = pB1[76];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[78];
         rB0 = pB0[78];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[78];
         rA2 = pA2[78];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[78];
         rB1 = pB1[78];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[80];
         rB0 = pB0[80];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[80];
         rA2 = pA2[80];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[80];
         rB1 = pB1[80];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[82];
         rB0 = pB0[82];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[82];
         rA2 = pA2[82];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[82];
         rB1 = pB1[82];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[84];
         rB0 = pB0[84];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[84];
         rA2 = pA2[84];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[84];
         rB1 = pB1[84];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         rA0 = pA0[86];
         rB0 = pB0[86];
         rC0_0 += rA0 * rB0;
         rA1 = pA1[86];
         rA2 = pA2[86];
         rC1_0 += rA1 * rB0;
         rA3 = pA3[86];
         rB1 = pB1[86];
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
