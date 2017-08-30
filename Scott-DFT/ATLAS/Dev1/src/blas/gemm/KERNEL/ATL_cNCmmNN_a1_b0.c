#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_cJIK44x44x44NN0x0x0_a1_b0
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=N, MB=44, NB=44, KB=44, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=44, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const float *stM = A + 88;
   const float *stN = B + (88*(ldb));
   const int incAk = (((lda) << 1));
   const int incAm = 8 - (88*(lda)), incAn = -88;
   #define incBk 88
   const int incBm = -88, incBn = (((ldb) << 2));
   #define incCm 8
   const int incCn = (((ldc) << 2)) - 88;
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
   const float *pA0=A;
   const float *pB0=B, *pB1=pB0+(((ldb) << 1));
   register int k;
   register float rA0, rA1, rA2, rA3;
   register float rB0, rB1;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC0_1, rC1_1, rC2_1, rC3_1;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 = rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 = rA1 * rB0;
         rA3 = pA0[6];
         rB1 = *pB1;
         rC2_0 = rA2 * rB0;
         rC3_0 = rA3 * rB0;
         rC0_1 = rA0 * rB1;
         rC1_1 = rA1 * rB1;
         rC2_1 = rA2 * rB1;
         rC3_1 = rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[2];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[2];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[4];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[4];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[6];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[6];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[8];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[8];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[10];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[10];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[12];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[12];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[14];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[14];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[16];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[16];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[18];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[18];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[20];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[20];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[22];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[22];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[24];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[24];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[26];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[26];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[28];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[28];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[30];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[30];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[32];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[32];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[34];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[34];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[36];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[36];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[38];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[38];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[40];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[40];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[42];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[42];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[44];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[44];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[46];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[46];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[48];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[48];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[50];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[50];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[52];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[52];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[54];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[54];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[56];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[56];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[58];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[58];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[60];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[60];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[62];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[62];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[64];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[64];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[66];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[66];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[68];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[68];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[70];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[70];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[72];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[72];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[74];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[74];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[76];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[76];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[78];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[78];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[80];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[80];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[82];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[82];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[84];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[84];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[86];
         rC0_0 += rA0 * rB0;
         rA1 = pA0[2];
         rA2 = pA0[4];
         rC1_0 += rA1 * rB0;
         rA3 = pA0[6];
         rB1 = pB1[86];
         rC2_0 += rA2 * rB0;
         rC3_0 += rA3 * rB0;
         rC0_1 += rA0 * rB1;
         rC1_1 += rA1 * rB1;
         rC2_1 += rA2 * rB1;
         rC3_1 += rA3 * rB1;
         pA0 += incAk;
         pB0 += incBk;
         pB1 += incBk;
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