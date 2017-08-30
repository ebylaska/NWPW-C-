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
static void ATL_dJIK0x0x19TN1x1x19_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=19, 
 * lda=19, ldb=19, ldc=0, mu=1, nu=1, ku=19, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (19*(Mb));
   const double *stN = B + (19*(Nb));
   #define incAk 19
   const int incAm = 0, incAn = -(19*(Mb));
   #define incBk 19
   const int incBm = -19, incBn = 19;
   const int incAk0 = ((incAk) / 19), incBk0 = ((incBk) / 19);
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
static void ATL_dJIK0x0x19TN2x1x19_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=19, 
 * lda=19, ldb=19, ldc=0, mu=2, nu=1, ku=19, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (19*(Mb));
   const double *stN = B + (19*(Nb));
   #define incAk 19
   const int incAm = 19, incAn = -(19*(Mb));
   #define incBk 19
   const int incBm = -19, incBn = 19;
   const int incAk0 = ((incAk) / 19), incBk0 = ((incBk) / 19);
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
            rA1 = pA0[19];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[20];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[21];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[22];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[23];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[24];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[25];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[26];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rA1 = pA0[27];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rA1 = pA0[28];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rA1 = pA0[29];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rA1 = pA0[30];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rA1 = pA0[31];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rA1 = pA0[32];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rA1 = pA0[33];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rA1 = pA0[34];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rA1 = pA0[35];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rA1 = pA0[36];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rA1 = pA0[37];
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
      ATL_dJIK0x0x19TN1x1x19_a1_bX(k, N, 19, alpha, ca + (19*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
static void ATL_dJIK0x0x19TN1x2x19_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=19, 
 * lda=19, ldb=19, ldc=0, mu=1, nu=2, ku=19, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (19*(Mb));
   const double *stN = B + (19*(Nb));
   #define incAk 19
   const int incAm = 0, incAn = -(19*(Mb));
   #define incBk 19
   const int incBm = -19, incBn = 38;
   const int incAk0 = ((incAk) / 19), incBk0 = ((incBk) / 19);
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
            rB1 = pB0[19];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rB1 = pB0[20];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rB1 = pB0[21];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rB1 = pB0[22];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rB1 = pB0[23];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rB1 = pB0[24];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rB1 = pB0[25];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rB1 = pB0[26];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[8];
            rB0 = pB0[8];
            rB1 = pB0[27];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[9];
            rB0 = pB0[9];
            rB1 = pB0[28];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[10];
            rB0 = pB0[10];
            rB1 = pB0[29];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[11];
            rB0 = pB0[11];
            rB1 = pB0[30];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[12];
            rB0 = pB0[12];
            rB1 = pB0[31];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[13];
            rB0 = pB0[13];
            rB1 = pB0[32];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[14];
            rB0 = pB0[14];
            rB1 = pB0[33];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[15];
            rB0 = pB0[15];
            rB1 = pB0[34];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[16];
            rB0 = pB0[16];
            rB1 = pB0[35];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[17];
            rB0 = pB0[17];
            rB1 = pB0[36];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[18];
            rB0 = pB0[18];
            rB1 = pB0[37];
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
      ATL_dJIK0x0x19TN1x1x19_a1_bX(M, k, 19, alpha, ca, lda, cb + (19*(Nb)), ldb, beta, cc + (Nb*ldc), ldc);
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
void ATL_dJIK0x0x19TN19x19x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=19, 
 * lda=19, ldb=19, ldc=0, mu=2, nu=2, ku=19, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (19*(Mb));
   const double *stN = B + (19*(Nb));
   #define incAk 19
   const int incAm = 19, incAn = -(19*(Mb));
   #define incBk 19
   const int incBm = -19, incBn = 38;
   const int incAk0 = ((incAk) / 19), incBk0 = ((incBk) / 19);
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
               rA1 = pA0[19];
               m0 = rA0 * rB0;
               rB1 = pB0[19];

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
               rA1 = pA0[20];
               rB1 = pB0[20];
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
               rA1 = pA0[21];
               rB1 = pB0[21];
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
               rA1 = pA0[22];
               rB1 = pB0[22];
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
               rA1 = pA0[23];
               rB1 = pB0[23];
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
               rA1 = pA0[24];
               rB1 = pB0[24];
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
               rA1 = pA0[25];
               rB1 = pB0[25];
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
               rA1 = pA0[26];
               rB1 = pB0[26];
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
               rA1 = pA0[27];
               rB1 = pB0[27];
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
               rA1 = pA0[28];
               rB1 = pB0[28];
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
               rA1 = pA0[29];
               rB1 = pB0[29];
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
               rA1 = pA0[30];
               rB1 = pB0[30];
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
               rA0 = pA0[13];
               rB0 = pB0[13];
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
               rA0 = pA0[14];
               rB0 = pB0[14];
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
               rA0 = pA0[15];
               rB0 = pB0[15];
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
               rA0 = pA0[16];
               rB0 = pB0[16];
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
               rA0 = pA0[17];
               rB0 = pB0[17];
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
               rA0 = pA0[18];
               rB0 = pB0[18];
               rA1 = pA0[37];
               rB1 = pB0[37];
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
      ATL_dJIK0x0x19TN2x1x19_a1_bX(M, k, 19, alpha, ca, lda, cb + (19*(Nb)), ldb, beta, cc + (Nb*ldc), ldc);
   if (Nb && (k=M-Mb))
      ATL_dJIK0x0x19TN1x2x19_a1_bX(k, Nb, 19, alpha, ca + (19*(Mb)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
