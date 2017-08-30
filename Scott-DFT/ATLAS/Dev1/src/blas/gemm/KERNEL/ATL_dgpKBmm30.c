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
static void ATL_dJIK0x0x30TN1x1x30_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=1, nu=1, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (30*(Mb));
   const double *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 0, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 30;
   const int incAk0 = ((incAk) / 30), incBk0 = ((incBk) / 30);
   #define incCm 1
   const int incCn = (ldc) - (Mb);
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0;
   register double m0;
   register double rC0_0;
   do /* N-loop */
   {
      do /* M-loop */
      {
         rA0 = beta;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         m0 = rA0 * rB0;
         rA0 = pA0[1];
         rB0 = pB0[1];

/*
 *       Completely unrolled K-loop
 */
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[2];
         rB0 = pB0[2];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[3];
         rB0 = pB0[3];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[4];
         rB0 = pB0[4];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[5];
         rB0 = pB0[5];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[6];
         rB0 = pB0[6];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[7];
         rB0 = pB0[7];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[8];
         rB0 = pB0[8];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[9];
         rB0 = pB0[9];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[10];
         rB0 = pB0[10];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[11];
         rB0 = pB0[11];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[12];
         rB0 = pB0[12];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[13];
         rB0 = pB0[13];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[14];
         rB0 = pB0[14];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[15];
         rB0 = pB0[15];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[16];
         rB0 = pB0[16];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[17];
         rB0 = pB0[17];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[18];
         rB0 = pB0[18];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[19];
         rB0 = pB0[19];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[20];
         rB0 = pB0[20];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[21];
         rB0 = pB0[21];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[22];
         rB0 = pB0[22];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[23];
         rB0 = pB0[23];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[24];
         rB0 = pB0[24];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[25];
         rB0 = pB0[25];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[26];
         rB0 = pB0[26];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[27];
         rB0 = pB0[27];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[28];
         rB0 = pB0[28];
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rA0 = pA0[29];
         rB0 = pB0[29];
/*
 *       Drain pipe on last iteration of K-loop
 */
         rC0_0 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
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
static void ATL_dJIK0x0x30TN2x1x30_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=2, nu=1, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (30*(Mb));
   const double *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 30, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 30;
   const int incAk0 = ((incAk) / 30), incBk0 = ((incBk) / 30);
   #define incCm 2
   const int incCn = (ldc) - (Mb);
   double *pC0=C;
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1;
   register double rB0;
   register double m0;
   register double rC0_0, rC1_0;
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
/*
 *          Start pipeline
 */
            rA0 = *pA0;
            rB0 = *pB0;
            rA1 = pA0[30];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[31];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[32];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[33];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[34];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[35];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[36];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[37];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[38];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[39];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[40];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[41];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[42];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[43];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[44];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[45];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[46];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[47];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[48];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rA1 = pA0[49];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rA1 = pA0[50];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rA1 = pA0[51];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rA1 = pA0[52];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rA1 = pA0[53];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rA1 = pA0[54];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rA1 = pA0[55];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rA1 = pA0[56];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rA1 = pA0[57];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rA1 = pA0[58];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rA1 = pA0[59];
            rC1_0 += m0;
            m0 = rA0 * rB0;
/*
 *          Drain pipe on last iteration of K-loop
 */
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            pA0 += incAk;
            pB0 += incBk;
            *pC0 = rC0_0;
            pC0[1] = rC1_0;
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
      ATL_dJIK0x0x30TN1x1x30_a1_bX(k, N, 30, alpha, ca + (30*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
static void ATL_dJIK0x0x30TN1x2x30_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=1, nu=2, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (30*(Mb));
   const double *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 0, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 60;
   const int incAk0 = ((incAk) / 30), incBk0 = ((incBk) / 30);
   #define incCm 1
   const int incCn = (((ldc) << 1)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0;
   register double rB0, rB1;
   register double m0;
   register double rC0_0, rC0_1;
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
/*
 *          Start pipeline
 */
            rA0 = *pA0;
            rB0 = *pB0;
            rB1 = pB0[30];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rB1 = pB0[31];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rB1 = pB0[32];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rB1 = pB0[33];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rB1 = pB0[34];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rB1 = pB0[35];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rB1 = pB0[36];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rB1 = pB0[37];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rB1 = pB0[38];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rB1 = pB0[39];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rB1 = pB0[40];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rB1 = pB0[41];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rB1 = pB0[42];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rB1 = pB0[43];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rB1 = pB0[44];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rB1 = pB0[45];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rB1 = pB0[46];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rB1 = pB0[47];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rB1 = pB0[48];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[19];
            rB0 = pB0[19];
            rB1 = pB0[49];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[20];
            rB0 = pB0[20];
            rB1 = pB0[50];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[21];
            rB0 = pB0[21];
            rB1 = pB0[51];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[22];
            rB0 = pB0[22];
            rB1 = pB0[52];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[23];
            rB0 = pB0[23];
            rB1 = pB0[53];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[24];
            rB0 = pB0[24];
            rB1 = pB0[54];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[25];
            rB0 = pB0[25];
            rB1 = pB0[55];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[26];
            rB0 = pB0[26];
            rB1 = pB0[56];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[27];
            rB0 = pB0[27];
            rB1 = pB0[57];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[28];
            rB0 = pB0[28];
            rB1 = pB0[58];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[29];
            rB0 = pB0[29];
            rB1 = pB0[59];
            rC0_1 += m0;
            m0 = rA0 * rB0;
/*
 *          Drain pipe on last iteration of K-loop
 */
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
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
      ATL_dJIK0x0x30TN1x1x30_a1_bX(M, k, 30, alpha, ca, lda, cb + (30*(Nb)), ldb, beta, cc + (Nb*ldc), ldc);
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
void ATL_dJIK0x0x30TN30x30x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=30, 
 * lda=30, ldb=30, ldc=0, mu=2, nu=2, ku=30, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (30*(Mb));
   const double *stN = B + (30*(Nb));
   #define incAk 30
   const int incAm = 30, incAn = -(30*(Mb));
   #define incBk 30
   const int incBm = -30, incBn = 60;
   const int incAk0 = ((incAk) / 30), incBk0 = ((incBk) / 30);
   #define incCm 2
   const int incCn = (((ldc) << 1)) - (Mb);
   double *pC0=C, *pC1=pC0+(ldc);
   const double *pA0=A;
   const double *pB0=B;
   register int k;
   register double rA0, rA1;
   register double rB0, rB1;
   register double m0;
   register double rC0_0, rC1_0, rC0_1, rC1_1;
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
               rC0_1 = *pC1;
               rC0_1 *= rA0;
               rC1_1 = pC1[1];
               rC1_1 *= rA0;
/*
 *             Start pipeline
 */
               rA0 = *pA0;
               rB0 = *pB0;
               rA1 = pA0[30];
               m0 = rA0 * rB0;
               rB1 = pB0[30];

/*
 *             Completely unrolled K-loop
 */
               rC0_0 += m0;
               m0 = rA1 * rB0;
               rC1_0 += m0;
               m0 = rA0 * rB1;
               rC0_1 += m0;
               m0 = rA1 * rB1;
               rA0 = pA0[1];
               rB0 = pB0[1];
               rA1 = pA0[31];
               rB1 = pB0[31];
               rC1_1 += m0;
               m0 = rA0 * rB0;
               rC0_0 += m0;
               m0 = rA1 * rB0;
               rC1_0 += m0;
               m0 = rA0 * rB1;
               rC0_1 += m0;
               m0 = rA1 * rB1;
               rA0 = pA0[2];
               rB0 = pB0[2];
               rA1 = pA0[32];
               rB1 = pB0[32];
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
               rA1 = pA0[33];
               rB1 = pB0[33];
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
               rA1 = pA0[34];
               rB1 = pB0[34];
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
               rA1 = pA0[35];
               rB1 = pB0[35];
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
               rA1 = pA0[36];
               rB1 = pB0[36];
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
               rA1 = pA0[37];
               rB1 = pB0[37];
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
               rA1 = pA0[38];
               rB1 = pB0[38];
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
               rA1 = pA0[39];
               rB1 = pB0[39];
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
               rA1 = pA0[40];
               rB1 = pB0[40];
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
               rA1 = pA0[41];
               rB1 = pB0[41];
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
               rA1 = pA0[42];
               rB1 = pB0[42];
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
               rA1 = pA0[43];
               rB1 = pB0[43];
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
               rA1 = pA0[44];
               rB1 = pB0[44];
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
               rA1 = pA0[45];
               rB1 = pB0[45];
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
               rA1 = pA0[46];
               rB1 = pB0[46];
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
               rA1 = pA0[47];
               rB1 = pB0[47];
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
               rA1 = pA0[48];
               rB1 = pB0[48];
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
               rA1 = pA0[49];
               rB1 = pB0[49];
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
               rA1 = pA0[50];
               rB1 = pB0[50];
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
               rA1 = pA0[51];
               rB1 = pB0[51];
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
               rA1 = pA0[52];
               rB1 = pB0[52];
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
               rA1 = pA0[53];
               rB1 = pB0[53];
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
               rA1 = pA0[54];
               rB1 = pB0[54];
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
               rA1 = pA0[55];
               rB1 = pB0[55];
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
               rA1 = pA0[56];
               rB1 = pB0[56];
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
               rA1 = pA0[57];
               rB1 = pB0[57];
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
               rA1 = pA0[58];
               rB1 = pB0[58];
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
               rA1 = pA0[59];
               rB1 = pB0[59];
               rC1_1 += m0;
               m0 = rA0 * rB0;
/*
 *             Drain pipe on last iteration of K-loop
 */
               rC0_0 += m0;
               m0 = rA1 * rB0;
               rC1_0 += m0;
               m0 = rA0 * rB1;
               rC0_1 += m0;
               m0 = rA1 * rB1;
               rC1_1 += m0;
               pA0 += incAk;
               pB0 += incBk;
               *pC0 = rC0_0;
               pC0[1] = rC1_0;
               *pC1 = rC0_1;
               pC1[1] = rC1_1;
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
      ATL_dJIK0x0x30TN2x1x30_a1_bX(M, k, 30, alpha, ca, lda, cb + (30*(Nb)), ldb, beta, cc + (Nb*ldc), ldc);
   if (Nb && (k=M-Mb))
      ATL_dJIK0x0x30TN1x2x30_a1_bX(k, Nb, 30, alpha, ca + (30*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
