#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_dJIK48x48x48TN0x0x0_a1_b0
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=48, NB=48, KB=48, 
 * lda=0, ldb=0, ldc=0, mu=2, nu=2, ku=48, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const double *stM = A + (((lda) << 5)+((lda) << 4));
   const double *stN = B + (((ldb) << 5)+((ldb) << 4));
   #define incAk 48
   const int incAm = ((((lda) << 1)) - 48), incAn = -(((lda) << 5)+((lda) << 4));
   #define incBk 48
   const int incBm = -48, incBn = (((ldb) << 1));
   const int incAk0 = ((incAk) / 48), incBk0 = ((incBk) / 48);
   #define incCm 2
   const int incCn = (((ldc) << 1)) - 48;
   double *pC0=C, *pC1=pC0+(ldc);
   const double *pA0=A, *pA1=pA0+(lda);
   const double *pB0=B, *pB1=pB0+(ldb);
   register int k;
   register double rA0, rA1;
   register double rB0, rB1;
   register double m0;
   register double rC0_0, rC1_0, rC0_1, rC1_1;
   do /* N-loop */
   {
      do /* M-loop */
      {
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = *pA1;
         m0 = rA0 * rB0;
         rB1 = *pB1;

/*
 *       Completely unrolled K-loop
 */
         rC0_0 = m0;
         m0 = rA1 * rB0;
         rC1_0 = m0;
         m0 = rA0 * rB1;
         rC0_1 = m0;
         m0 = rA1 * rB1;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rA1 = pA1[1];
         rB1 = pB1[1];
         rC1_1 = m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rA1 = pA1[2];
         rB1 = pB1[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rA1 = pA1[3];
         rB1 = pB1[3];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rA1 = pA1[4];
         rB1 = pB1[4];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rA1 = pA1[5];
         rB1 = pB1[5];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rA1 = pA1[6];
         rB1 = pB1[6];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rA1 = pA1[7];
         rB1 = pB1[7];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rA1 = pA1[8];
         rB1 = pB1[8];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rA1 = pA1[9];
         rB1 = pB1[9];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rA1 = pA1[10];
         rB1 = pB1[10];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rA1 = pA1[11];
         rB1 = pB1[11];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rA1 = pA1[12];
         rB1 = pB1[12];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rA1 = pA1[13];
         rB1 = pB1[13];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rA1 = pA1[14];
         rB1 = pB1[14];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rA1 = pA1[15];
         rB1 = pB1[15];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rA1 = pA1[16];
         rB1 = pB1[16];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rA1 = pA1[17];
         rB1 = pB1[17];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rA1 = pA1[18];
         rB1 = pB1[18];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rA1 = pA1[19];
         rB1 = pB1[19];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rA1 = pA1[20];
         rB1 = pB1[20];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rA1 = pA1[21];
         rB1 = pB1[21];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rA1 = pA1[22];
         rB1 = pB1[22];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rA1 = pA1[23];
         rB1 = pB1[23];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rA1 = pA1[24];
         rB1 = pB1[24];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rA1 = pA1[25];
         rB1 = pB1[25];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rA1 = pA1[26];
         rB1 = pB1[26];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rA1 = pA1[27];
         rB1 = pB1[27];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rA1 = pA1[28];
         rB1 = pB1[28];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[29];
         rB0 = pB0[29];
         rA1 = pA1[29];
         rB1 = pB1[29];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[30];
         rB0 = pB0[30];
         rA1 = pA1[30];
         rB1 = pB1[30];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[31];
         rB0 = pB0[31];
         rA1 = pA1[31];
         rB1 = pB1[31];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[32];
         rB0 = pB0[32];
         rA1 = pA1[32];
         rB1 = pB1[32];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[33];
         rB0 = pB0[33];
         rA1 = pA1[33];
         rB1 = pB1[33];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[34];
         rB0 = pB0[34];
         rA1 = pA1[34];
         rB1 = pB1[34];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[35];
         rB0 = pB0[35];
         rA1 = pA1[35];
         rB1 = pB1[35];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[36];
         rB0 = pB0[36];
         rA1 = pA1[36];
         rB1 = pB1[36];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[37];
         rB0 = pB0[37];
         rA1 = pA1[37];
         rB1 = pB1[37];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[38];
         rB0 = pB0[38];
         rA1 = pA1[38];
         rB1 = pB1[38];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[39];
         rB0 = pB0[39];
         rA1 = pA1[39];
         rB1 = pB1[39];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[40];
         rB0 = pB0[40];
         rA1 = pA1[40];
         rB1 = pB1[40];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[41];
         rB0 = pB0[41];
         rA1 = pA1[41];
         rB1 = pB1[41];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[42];
         rB0 = pB0[42];
         rA1 = pA1[42];
         rB1 = pB1[42];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[43];
         rB0 = pB0[43];
         rA1 = pA1[43];
         rB1 = pB1[43];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[44];
         rB0 = pB0[44];
         rA1 = pA1[44];
         rB1 = pB1[44];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[45];
         rB0 = pB0[45];
         rA1 = pA1[45];
         rB1 = pB1[45];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[46];
         rB0 = pB0[46];
         rA1 = pA1[46];
         rB1 = pB1[46];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         rA0 = pA0[47];
         rB0 = pB0[47];
         rA1 = pA1[47];
         rB1 = pB1[47];
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
         rC1_1 += m0;
         pA0 += incAk;
         pA1 += incAk;
         pB0 += incBk;
         pB1 += incBk;
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         *pC1 = rC0_1;
         pC1[1] = rC1_1;
         pC0 += incCm;
         pC1 += incCm;
         pA0 += incAm;
         pA1 += incAm;
         pB0 += incBm;
         pB1 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
      pC1 += incCn;
      pA0 += incAn;
      pA1 += incAn;
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
