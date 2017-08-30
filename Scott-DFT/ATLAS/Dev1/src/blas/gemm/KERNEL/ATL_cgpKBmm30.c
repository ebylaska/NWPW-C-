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
static void ATL_cJIK0x0x30TN1x1x30_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=1, nu=1, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const float *stM = A + (30*(Mb));
   const float *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 0, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 30;
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
         rA0 = pA0[24];
         rB0 = pB0[24];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rC0_0 += rA0 * rB0;
         rA0 = pA0[29];
         rB0 = pB0[29];
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
static void ATL_cJIK0x0x30TN4x1x30_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=4, nu=1, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   #define Nb N
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (30*(Mb));
   const float *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 90, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 30;
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
            rA1 = pA0[30];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[60];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[90];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[31];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[61];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[91];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[32];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[62];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[92];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[33];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[63];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[93];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[34];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[64];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[94];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[35];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[65];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[95];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[36];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[66];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[96];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[37];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[67];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[97];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[38];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[68];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[98];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[39];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[69];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[99];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[40];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[70];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[100];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[41];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[71];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[101];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[42];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[72];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[102];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[43];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[73];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[103];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[44];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[74];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[104];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[45];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[75];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[105];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[46];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[76];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[106];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[47];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[77];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[107];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[48];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[78];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[108];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[49];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[79];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[109];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[50];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[80];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[110];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[51];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[81];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[111];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[52];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[82];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[112];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[53];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[83];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[113];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[54];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[84];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[114];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[55];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[85];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[115];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[56];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[86];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[116];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[57];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[87];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[117];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[58];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[88];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[118];
            rC3_0 += rA3 * rB0;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rC0_0 += rA0 * rB0;
            rA1 = pA0[59];
            rC1_0 += rA1 * rB0;
            rA2 = pA0[89];
            rC2_0 += rA2 * rB0;
            rA3 = pA0[119];
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
      ATL_cJIK0x0x30TN1x1x30_a1_bX(k, N, 30, alpha, ca + (30*(Mb)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
static void ATL_cJIK0x0x30TN1x2x30_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=1, nu=2, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (30*(Mb));
   const float *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 0, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 60;
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
            rB1 = pB0[30];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[31];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[32];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[33];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[34];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[35];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[36];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[37];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[38];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[39];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[40];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[41];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[42];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[43];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[44];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[45];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[46];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[47];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[48];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[49];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[50];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[51];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[52];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[53];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[54];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[55];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[56];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[57];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[58];
            rC0_1 += rA0 * rB1;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rC0_0 += rA0 * rB0;
            rB1 = pB0[59];
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
      ATL_cJIK0x0x30TN1x1x30_a1_bX(M, k, 30, alpha, ca, lda, cb + (30*(Nb)), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
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
void ATL_cJIK0x0x30TN30x30x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=4, nu=2, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   const int Nb = (N>>1)<<1;
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (30*(Mb));
   const float *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 90, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 60;
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
               rA1 = pA0[30];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[60];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[90];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[30];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[1];
               rB0 = pB0[1];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[31];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[61];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[91];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[31];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[2];
               rB0 = pB0[2];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[32];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[62];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[92];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[32];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[3];
               rB0 = pB0[3];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[33];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[63];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[93];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[33];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[4];
               rB0 = pB0[4];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[34];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[64];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[94];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[34];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[5];
               rB0 = pB0[5];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[35];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[65];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[95];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[35];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[6];
               rB0 = pB0[6];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[36];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[66];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[96];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[36];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[7];
               rB0 = pB0[7];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[37];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[67];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[97];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[37];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[8];
               rB0 = pB0[8];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[38];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[68];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[98];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[38];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[9];
               rB0 = pB0[9];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[39];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[69];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[99];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[39];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[10];
               rB0 = pB0[10];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[40];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[70];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[100];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[40];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[11];
               rB0 = pB0[11];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[41];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[71];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[101];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[41];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[12];
               rB0 = pB0[12];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[42];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[72];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[102];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[42];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[13];
               rB0 = pB0[13];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[43];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[73];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[103];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[43];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[14];
               rB0 = pB0[14];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[44];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[74];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[104];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[44];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[15];
               rB0 = pB0[15];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[45];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[75];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[105];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[45];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[16];
               rB0 = pB0[16];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[46];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[76];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[106];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[46];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[17];
               rB0 = pB0[17];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[47];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[77];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[107];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[47];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[18];
               rB0 = pB0[18];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[48];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[78];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[108];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[48];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[19];
               rB0 = pB0[19];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[49];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[79];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[109];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[49];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[20];
               rB0 = pB0[20];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[50];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[80];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[110];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[50];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[21];
               rB0 = pB0[21];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[51];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[81];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[111];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[51];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[22];
               rB0 = pB0[22];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[52];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[82];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[112];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[52];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[23];
               rB0 = pB0[23];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[53];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[83];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[113];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[53];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[24];
               rB0 = pB0[24];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[54];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[84];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[114];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[54];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[25];
               rB0 = pB0[25];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[55];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[85];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[115];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[55];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[26];
               rB0 = pB0[26];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[56];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[86];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[116];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[56];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[27];
               rB0 = pB0[27];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[57];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[87];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[117];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[57];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[28];
               rB0 = pB0[28];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[58];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[88];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[118];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[58];
               rC0_1 += rA0 * rB1;
               rC1_1 += rA1 * rB1;
               rC2_1 += rA2 * rB1;
               rC3_1 += rA3 * rB1;
               rA0 = pA0[29];
               rB0 = pB0[29];
               rC0_0 += rA0 * rB0;
               rA1 = pA0[59];
               rC1_0 += rA1 * rB0;
               rA2 = pA0[89];
               rC2_0 += rA2 * rB0;
               rA3 = pA0[119];
               rC3_0 += rA3 * rB0;
               rB1 = pB0[59];
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
      ATL_cJIK0x0x30TN4x1x30_a1_bX(M, k, 30, alpha, ca, lda, cb + (30*(Nb)), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
   if (Nb && (k=M-Mb))
      ATL_cJIK0x0x30TN1x2x30_a1_bX(k, Nb, 30, alpha, ca + (30*(Mb)), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
