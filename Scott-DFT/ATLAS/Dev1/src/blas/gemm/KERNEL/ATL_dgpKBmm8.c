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
static void ATL_dJIK0x0x8TN1x1x8_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=8, 
 * lda=8, ldb=8, ldc=0, mu=1, nu=1, ku=8, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   const double *stM = A + (((Mb) << 3));
   const double *stN = B + (((Nb) << 3));
   #define incAk 8
   const int incAm = 0, incAn = -(((Mb) << 3));
   #define incBk 8
   const int incBm = -8, incBn = 8;
   const int incAk0 = ((incAk) >> 3), incBk0 = ((incBk) >> 3);
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
static void ATL_dJIK0x0x8TN2x1x8_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=8, 
 * lda=8, ldb=8, ldc=0, mu=2, nu=1, ku=8, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   #define Nb N
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (((Mb) << 3));
   const double *stN = B + (((Nb) << 3));
   #define incAk 8
   const int incAm = 8, incAn = -(((Mb) << 3));
   #define incBk 8
   const int incBm = -8, incBn = 8;
   const int incAk0 = ((incAk) >> 3), incBk0 = ((incBk) >> 3);
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
            rA1 = pA0[8];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rA1 = pA0[9];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rA1 = pA0[10];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rA1 = pA0[11];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rA1 = pA0[12];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rA1 = pA0[13];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rA1 = pA0[14];
            rC1_0 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rA1 = pA0[15];
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
      ATL_dJIK0x0x8TN1x1x8_a1_bX(k, N, 8, alpha, ca + (((Mb) << 3)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
static void ATL_dJIK0x0x8TN1x2x8_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=8, 
 * lda=8, ldb=8, ldc=0, mu=1, nu=2, ku=8, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (((Mb) << 3));
   const double *stN = B + (((Nb) << 3));
   #define incAk 8
   const int incAm = 0, incAn = -(((Mb) << 3));
   #define incBk 8
   const int incBm = -8, incBn = 16;
   const int incAk0 = ((incAk) >> 3), incBk0 = ((incBk) >> 3);
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
            rB1 = pB0[8];
            m0 = rA0 * rB0;

/*
 *          Completely unrolled K-loop
 */
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[1];
            rB0 = pB0[1];
            rB1 = pB0[9];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[2];
            rB0 = pB0[2];
            rB1 = pB0[10];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[3];
            rB0 = pB0[3];
            rB1 = pB0[11];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[4];
            rB0 = pB0[4];
            rB1 = pB0[12];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[5];
            rB0 = pB0[5];
            rB1 = pB0[13];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[6];
            rB0 = pB0[6];
            rB1 = pB0[14];
            rC0_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA0 * rB1;
            rA0 = pA0[7];
            rB0 = pB0[7];
            rB1 = pB0[15];
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
      ATL_dJIK0x0x8TN1x1x8_a1_bX(M, k, 8, alpha, ca, lda, cb + (((Nb) << 3)), ldb, beta, cc + (Nb*ldc), ldc);
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
void ATL_dJIK0x0x8TN8x8x0_a1_bX
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=8, 
 * lda=8, ldb=8, ldc=0, mu=2, nu=2, ku=8, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>1)<<1;
   const int Nb = (N>>1)<<1;
   const double *ca=A, *cb=B;
   double *cc=C;
   const double *stM = A + (((Mb) << 3));
   const double *stN = B + (((Nb) << 3));
   #define incAk 8
   const int incAm = 8, incAn = -(((Mb) << 3));
   #define incBk 8
   const int incBm = -8, incBn = 16;
   const int incAk0 = ((incAk) >> 3), incBk0 = ((incBk) >> 3);
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
               rA1 = pA0[8];
               m0 = rA0 * rB0;
               rB1 = pB0[8];

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
               rA1 = pA0[9];
               rB1 = pB0[9];
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
               rA1 = pA0[10];
               rB1 = pB0[10];
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
               rA1 = pA0[11];
               rB1 = pB0[11];
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
               rA1 = pA0[12];
               rB1 = pB0[12];
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
               rA1 = pA0[13];
               rB1 = pB0[13];
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
               rA1 = pA0[14];
               rB1 = pB0[14];
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
               rA1 = pA0[15];
               rB1 = pB0[15];
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
      ATL_dJIK0x0x8TN2x1x8_a1_bX(M, k, 8, alpha, ca, lda, cb + (((Nb) << 3)), ldb, beta, cc + (Nb*ldc), ldc);
   if (Nb && (k=M-Mb))
      ATL_dJIK0x0x8TN1x2x8_a1_bX(k, Nb, 8, alpha, ca + (((Mb) << 3)), lda, cb, ldb, beta, cc + (Mb), ldc);
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
