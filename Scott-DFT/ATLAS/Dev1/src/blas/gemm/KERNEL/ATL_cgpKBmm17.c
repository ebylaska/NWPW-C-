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
static void ATL_cJIK0x0x17TN1x1x17_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=17, 
 * lda=17, ldb=17, ldc=0, mu=1, nu=1, ku=17, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const float *stM = A + (((Mb) << 4)+Mb);
   const float *stN = B + (((Nb) << 4)+Nb);
   #define incAk 17
   const int incAm = 0, incAn = -(((Mb) << 4)+Mb);
   #define incBk 17
   const int incBm = -17, incBn = 17;
   #define incCm 2
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
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
static void ATL_cJIK0x0x17TN4x1x17_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=17, 
 * lda=17, ldb=17, ldc=0, mu=4, nu=1, ku=17, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   #define Nb N
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+Mb);
   const float *stN = B + (((Nb) << 4)+Nb);
   #define incAk 17
   const int incAm = 51, incAn = -(((Mb) << 4)+Mb);
   #define incBk 17
   const int incBm = -17, incBn = 17;
   #define incCm 8
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
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
            rC1_0 = pC0[2];
            rC1_0 *= rA0;
            rC2_0 = pC0[4];
            rC2_0 *= rA0;
            rC3_0 = pC0[6];
            rC3_0 *= rA0;
            rA0 = *pA0;
            rB0 = *pB0;
            rC0_0 += rA0 * rB0;
            rA1 = pA0[17];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[34];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[51];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[18];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[35];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[52];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[19];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[36];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[53];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[20];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[37];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[54];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[21];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[38];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[55];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[22];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[39];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[56];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[23];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[40];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[57];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[24];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[41];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[58];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[25];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[42];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[59];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[26];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[43];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[60];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[27];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[44];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[61];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[28];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[45];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[62];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[29];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[46];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[63];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[30];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[47];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[64];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[31];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[48];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[65];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[32];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[49];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[66];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[33];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[50];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[67];
            rC3_0 += rA3 * rB0;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            pC0[2] = rC1_0;
            pC0[4] = rC2_0;
            pC0[6] = rC3_0;
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
      ATL_cJIK0x0x17TN1x1x17_a1_bX(k, N, 17, alpha, ca + (((Mb) << 4)+Mb), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
static void ATL_cJIK0x0x17TN1x2x17_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=17, 
 * lda=17, ldb=17, ldc=0, mu=1, nu=2, ku=17, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+Mb);
   const float *stN = B + (((Nb) << 4)+Nb);
   #define incAk 17
   const int incAm = 0, incAn = -(((Mb) << 4)+Mb);
   #define incBk 17
   const int incBm = -17, incBn = 34;
   #define incCm 2
   const int incCn = (((ldc) << 2)) - (((Mb) << 1));
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
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
            rC0_0 += rA0 * rB0;
            rB1 = pB0[17];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[18];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[19];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[20];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[21];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[22];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[23];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[24];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[25];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[26];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[27];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[28];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[29];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[30];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[31];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[32];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[33];
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
      ATL_cJIK0x0x17TN1x1x17_a1_bX(M, k, 17, alpha, ca, lda, cb + (((Nb) << 4)+Nb), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
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
void ATL_cJIK0x0x17TN17x17x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=17, 
 * lda=17, ldb=17, ldc=0, mu=4, nu=2, ku=17, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (((Mb) << 4)+Mb);
   const float *stN = B + (((Nb) << 4)+Nb);
   #define incAk 17
   const int incAm = 51, incAn = -(((Mb) << 4)+Mb);
   #define incBk 17
   const int incBm = -17, incBn = 34;
   #define incCm 8
   const int incCn = (((ldc) << 2)) - (((Mb) << 1));
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
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
               rA1 = pA0[17];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[34];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[51];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[17];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[1];
               rB0 = pB0[1];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[18];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[35];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[52];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[18];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[2];
               rB0 = pB0[2];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[19];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[36];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[53];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[19];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[3];
               rB0 = pB0[3];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[20];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[37];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[54];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[20];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[4];
               rB0 = pB0[4];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[21];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[38];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[55];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[21];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[5];
               rB0 = pB0[5];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[22];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[39];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[56];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[22];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[6];
               rB0 = pB0[6];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[23];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[40];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[57];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[23];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[7];
               rB0 = pB0[7];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[24];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[41];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[58];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[24];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[8];
               rB0 = pB0[8];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[25];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[42];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[59];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[25];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[9];
               rB0 = pB0[9];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[26];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[43];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[60];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[26];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[10];
               rB0 = pB0[10];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[27];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[44];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[61];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[27];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[11];
               rB0 = pB0[11];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[28];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[45];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[62];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[28];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[12];
               rB0 = pB0[12];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[29];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[46];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[63];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[29];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[13];
               rB0 = pB0[13];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[30];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[47];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[64];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[30];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[14];
               rB0 = pB0[14];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[31];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[48];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[65];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[31];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[15];
               rB0 = pB0[15];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[32];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[49];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[66];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[32];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[16];
               rB0 = pB0[16];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[33];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[50];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[67];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[33];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               pA0 += incAk;
               pB0 += incBk;
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
   }
   if (k=N-Nb)
      ATL_cJIK0x0x17TN4x1x17_a1_bX(M, k, 17, alpha, ca, lda, cb + (((Nb) << 4)+Nb), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
   if (Nb && (k=M-Mb))
      ATL_cJIK0x0x17TN1x2x17_a1_bX(k, Nb, 17, alpha, ca + (((Mb) << 4)+Mb), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
