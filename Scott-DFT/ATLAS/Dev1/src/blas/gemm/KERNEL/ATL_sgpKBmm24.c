#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_sJIK0x0x24TN1x1x24_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=24, ldb=24, ldc=0, mu=1, nu=1, ku=24, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const float *stM = A + (((Mb) << 4)+((Mb) << 3));
   const float *stN = B + (((Nb) << 4)+((Nb) << 3));
   #define incAk 24
   const int incAm = 0, incAn = -(((Mb) << 4)+((Mb) << 3));
   #define incBk 24
   const int incBm = -24, incBn = 24;
   #define incCm 1
   const int incCn = (ldc) - (Mb);
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0;
   register float rB0;
   register float rC0_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rA0 = *pA0;
         rB0 = *pB0;
         rC0_0 += rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rC0_0 += rA0 * rB0;
         pA0 += incAk;
         pB0 += incBk;
         *pC0 = rC0_0;
         pC0 += incCm;
         pA0 += incAm;
         pB0 += incBm;
      }
      while(pA0 != stM);
      pC0 += incCn;
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
static void ATL_sJIK0x0x24TN4x1x24_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=24, ldb=24, ldc=0, mu=4, nu=1, ku=24, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   #define Nb N
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+((Mb) << 3));
   const float *stN = B + (((Nb) << 4)+((Nb) << 3));
   #define incAk 24
   const int incAm = 72, incAn = -(((Mb) << 4)+((Mb) << 3));
   #define incBk 24
   const int incBm = -24, incBn = 24;
   #define incCm 4
   const int incCn = (ldc) - (Mb);
   float *pC0=C;
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3;
   register float rB0;
   register float rC0_0, rC1_0, rC2_0, rC3_0;
   if (pA0 != stM)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
            rA0 = beta;
            rC0_0 = *pC0;
            rC0_0 *= rA0;
            rC1_0 = pC0[1];
            rC1_0 *= rA0;
            rC2_0 = pC0[2];
            rC2_0 *= rA0;
            rC3_0 = pC0[3];
            rC3_0 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[24];
            rA2 = pA0[48];
            rA3 = pA0[72];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[25];
            rA2 = pA0[49];
            rA3 = pA0[73];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[26];
            rA2 = pA0[50];
            rA3 = pA0[74];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[27];
            rA2 = pA0[51];
            rA3 = pA0[75];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[28];
            rA2 = pA0[52];
            rA3 = pA0[76];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[29];
            rA2 = pA0[53];
            rA3 = pA0[77];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[30];
            rA2 = pA0[54];
            rA3 = pA0[78];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[31];
            rA2 = pA0[55];
            rA3 = pA0[79];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[32];
            rA2 = pA0[56];
            rA3 = pA0[80];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[33];
            rA2 = pA0[57];
            rA3 = pA0[81];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[34];
            rA2 = pA0[58];
            rA3 = pA0[82];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[35];
            rA2 = pA0[59];
            rA3 = pA0[83];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[36];
            rA2 = pA0[60];
            rA3 = pA0[84];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[37];
            rA2 = pA0[61];
            rA3 = pA0[85];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[38];
            rA2 = pA0[62];
            rA3 = pA0[86];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[39];
            rA2 = pA0[63];
            rA3 = pA0[87];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[40];
            rA2 = pA0[64];
            rA3 = pA0[88];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[41];
            rA2 = pA0[65];
            rA3 = pA0[89];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[42];
            rA2 = pA0[66];
            rA3 = pA0[90];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA0[43];
            rA2 = pA0[67];
            rA3 = pA0[91];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA0[44];
            rA2 = pA0[68];
            rA3 = pA0[92];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA0[45];
            rA2 = pA0[69];
            rA3 = pA0[93];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA0[46];
            rA2 = pA0[70];
            rA3 = pA0[94];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA0[47];
            rA2 = pA0[71];
            rA3 = pA0[95];
            rC0_0 += rA0 * rB0;
            rC1_0 += rA1 * rB0;
            rC2_0 += rA2 * rB0;
            rC3_0 += rA3 * rB0;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            pC0[1] = rC1_0;
            pC0[2] = rC2_0;
            pC0[3] = rC3_0;
            pC0 += incCm;
            pA0 += incAm;
            pB0 += incBm;
         }
         while(pA0 != stM);
         pC0 += incCn;
         pA0 += incAn;
         pB0 += incBn;
      }
      while(pB0 != stN);
   }
   if (k=M-Mb)
      ATL_sJIK0x0x24TN1x1x24_a1_bX(k, N, 24, alpha, ca + (((Mb) << 4)+((Mb) << 3)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
static void ATL_sJIK0x0x24TN1x2x24_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=24, ldb=24, ldc=0, mu=1, nu=2, ku=24, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+((Mb) << 3));
   const float *stN = B + (((Nb) << 4)+((Nb) << 3));
   #define incAk 24
   const int incAm = 0, incAn = -(((Mb) << 4)+((Mb) << 3));
   #define incBk 24
   const int incBm = -24, incBn = 48;
   #define incCm 1
   const int incCn = (((ldc) << 1)) - (Mb);
   float *pC0=C, *pC1=pC0+(ldc);
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0;
   register float rB0, rB1;
   register float rC0_0, rC0_1;
   if (pB0 != stN)
   {
      do /* N-loop */
      {
         do /* M-loop */
         {
            rA0 = beta;
            rC0_0 = *pC0;
            rC0_0 *= rA0;
            rC0_1 = *pC1;
            rC0_1 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rB1 = pB0[24];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rB1 = pB0[25];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rB1 = pB0[26];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rB1 = pB0[27];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rB1 = pB0[28];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rB1 = pB0[29];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rB1 = pB0[30];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rB1 = pB0[31];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rB1 = pB0[32];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rB1 = pB0[33];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rB1 = pB0[34];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rB1 = pB0[35];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rB1 = pB0[36];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rB1 = pB0[37];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rB1 = pB0[38];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rB1 = pB0[39];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rB1 = pB0[40];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rB1 = pB0[41];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rB1 = pB0[42];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rB1 = pB0[43];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rB1 = pB0[44];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rB1 = pB0[45];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rB1 = pB0[46];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rB1 = pB0[47];
            rC0_0 += rA0 * rB0;
            rC0_1 += rA0 * rB1;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            *pC1 = rC0_1;
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
   if (k=N-Nb)
      ATL_sJIK0x0x24TN1x1x24_a1_bX(M, k, 24, alpha, ca, lda, cb + (((Nb) << 4)+((Nb) << 3)), ldb, beta, cc + (Nb*ldc), ldc);
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
void ATL_sJIK0x0x24TN24x24x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=24, 
 * lda=24, ldb=24, ldc=0, mu=4, nu=2, ku=24, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+((Mb) << 3));
   const float *stN = B + (((Nb) << 4)+((Nb) << 3));
   #define incAk 24
   const int incAm = 72, incAn = -(((Mb) << 4)+((Mb) << 3));
   #define incBk 24
   const int incBm = -24, incBn = 48;
   #define incCm 4
   const int incCn = (((ldc) << 1)) - (Mb);
   float *pC0=C, *pC1=pC0+(ldc);
   const float *pA0=A;
   const float *pB0=B;
   register int k;
   register float rA0, rA1, rA2, rA3;
   register float rB0, rB1;
   register float rC0_0, rC1_0, rC2_0, rC3_0, rC0_1, rC1_1, rC2_1, rC3_1;
   if (pA0 != stM)
   {
      if (pB0 != stN)
      {
         do /* N-loop */
         {
            do /* M-loop */
            {
               rA0 = beta;
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
               rA1 = pA0[24];
               rA2 = pA0[48];
               rA3 = pA0[72];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[24];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[1];
               rB0 = pB0[1];
               rA1 = pA0[25];
               rA2 = pA0[49];
               rA3 = pA0[73];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[25];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[2];
               rB0 = pB0[2];
               rA1 = pA0[26];
               rA2 = pA0[50];
               rA3 = pA0[74];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[26];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[3];
               rB0 = pB0[3];
               rA1 = pA0[27];
               rA2 = pA0[51];
               rA3 = pA0[75];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[27];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[4];
               rB0 = pB0[4];
               rA1 = pA0[28];
               rA2 = pA0[52];
               rA3 = pA0[76];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[28];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[5];
               rB0 = pB0[5];
               rA1 = pA0[29];
               rA2 = pA0[53];
               rA3 = pA0[77];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[29];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[6];
               rB0 = pB0[6];
               rA1 = pA0[30];
               rA2 = pA0[54];
               rA3 = pA0[78];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[30];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[7];
               rB0 = pB0[7];
               rA1 = pA0[31];
               rA2 = pA0[55];
               rA3 = pA0[79];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[31];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[8];
               rB0 = pB0[8];
               rA1 = pA0[32];
               rA2 = pA0[56];
               rA3 = pA0[80];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[32];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[9];
               rB0 = pB0[9];
               rA1 = pA0[33];
               rA2 = pA0[57];
               rA3 = pA0[81];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[33];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[10];
               rB0 = pB0[10];
               rA1 = pA0[34];
               rA2 = pA0[58];
               rA3 = pA0[82];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[34];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[11];
               rB0 = pB0[11];
               rA1 = pA0[35];
               rA2 = pA0[59];
               rA3 = pA0[83];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[35];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[12];
               rB0 = pB0[12];
               rA1 = pA0[36];
               rA2 = pA0[60];
               rA3 = pA0[84];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[36];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[13];
               rB0 = pB0[13];
               rA1 = pA0[37];
               rA2 = pA0[61];
               rA3 = pA0[85];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[37];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[14];
               rB0 = pB0[14];
               rA1 = pA0[38];
               rA2 = pA0[62];
               rA3 = pA0[86];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[38];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[15];
               rB0 = pB0[15];
               rA1 = pA0[39];
               rA2 = pA0[63];
               rA3 = pA0[87];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[39];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[16];
               rB0 = pB0[16];
               rA1 = pA0[40];
               rA2 = pA0[64];
               rA3 = pA0[88];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[40];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[17];
               rB0 = pB0[17];
               rA1 = pA0[41];
               rA2 = pA0[65];
               rA3 = pA0[89];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[41];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[18];
               rB0 = pB0[18];
               rA1 = pA0[42];
               rA2 = pA0[66];
               rA3 = pA0[90];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[42];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[19];
               rB0 = pB0[19];
               rA1 = pA0[43];
               rA2 = pA0[67];
               rA3 = pA0[91];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[43];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[20];
               rB0 = pB0[20];
               rA1 = pA0[44];
               rA2 = pA0[68];
               rA3 = pA0[92];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[44];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[21];
               rB0 = pB0[21];
               rA1 = pA0[45];
               rA2 = pA0[69];
               rA3 = pA0[93];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[45];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[22];
               rB0 = pB0[22];
               rA1 = pA0[46];
               rA2 = pA0[70];
               rA3 = pA0[94];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[46];
               rC1_0 += rA1 * rB0;
               rC2_0 += rA2 * rB0;
               rC3_0 += rA3 * rB0;
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[23];
               rB0 = pB0[23];
               rA1 = pA0[47];
               rA2 = pA0[71];
               rA3 = pA0[95];
               rC0_0 += rA0 * rB0;
               rB1 = pB0[47];
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
   }
   if (k=N-Nb)
      ATL_sJIK0x0x24TN4x1x24_a1_bX(M, k, 24, alpha, ca, lda, cb + (((Nb) << 4)+((Nb) << 3)), ldb, beta, cc + (Nb*ldc), ldc);
   if (Nb && (k=M-Mb))
      ATL_sJIK0x0x24TN1x2x24_a1_bX(k, Nb, 24, alpha, ca + (((Mb) << 4)+((Mb) << 3)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
