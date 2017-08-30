#ifndef ATL_RESTRICT
#if defined(__STDC_VERSION__) && (__STDC_VERSION__/100 >= 1999)
   #define ATL_RESTRICT restrict
#else
   #define ATL_RESTRICT
#endif
#endif
void ATL_dJIK48x48x48NN0x0x0_aX_b1
   (const int M, const int N, const int K, const double alpha, const double * ATL_RESTRICT A, const int lda, const double * ATL_RESTRICT B, const int ldb, const double beta, double * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=N, TB=N, MB=48, NB=48, KB=48, 
 * lda=0, ldb=0, ldc=0, mu=2, nu=2, ku=24, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const double *stM = A + 48;
   const double *stN = B + (((ldb) << 5)+((ldb) << 4));
   const double BetaAlpha = beta / alpha;
   const int incAk = (lda);
   const int incAm = 2 - (((lda) << 5)+((lda) << 4)), incAn = -48;
   #define incBk 24
   const int incBm = -48, incBn = (((ldb) << 1));
   #define incAk0 incAk
   const int incBk0 = ((incBk) / 24);
   #define incCm 2
   const int incCn = (((ldc) << 1)) - 48;
   double *pC0=C, *pC1=pC0+(ldc);
   const double *pA0=A;
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
         rA0 = BetaAlpha;
         rC0_0 = *pC0;
         rC0_0 *= rA0;
         rC1_0 = pC0[1];
         rC1_0 *= rA0;
         rC0_1 = *pC1;
         rC0_1 *= rA0;
         rC1_1 = pC1[1];
         rC1_1 *= rA0;
/*
 *       Start pipeline
 */
         rA0 = *pA0;
         rB0 = *pB0;
         rA1 = pA0[1];
         m0 = rA0 * rB0;
         rB1 = *pB1;

         for (k=0; k < 24; k += 24) /* easy loop to unroll */
         {
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[1];
            rA1 = pA0[1];
            rB1 = pB1[1];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[2];
            rA1 = pA0[1];
            rB1 = pB1[2];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[3];
            rA1 = pA0[1];
            rB1 = pB1[3];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[4];
            rA1 = pA0[1];
            rB1 = pB1[4];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[5];
            rA1 = pA0[1];
            rB1 = pB1[5];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[6];
            rA1 = pA0[1];
            rB1 = pB1[6];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[7];
            rA1 = pA0[1];
            rB1 = pB1[7];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[8];
            rA1 = pA0[1];
            rB1 = pB1[8];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[9];
            rA1 = pA0[1];
            rB1 = pB1[9];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[10];
            rA1 = pA0[1];
            rB1 = pB1[10];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[11];
            rA1 = pA0[1];
            rB1 = pB1[11];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[12];
            rA1 = pA0[1];
            rB1 = pB1[12];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[13];
            rA1 = pA0[1];
            rB1 = pB1[13];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[14];
            rA1 = pA0[1];
            rB1 = pB1[14];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[15];
            rA1 = pA0[1];
            rB1 = pB1[15];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[16];
            rA1 = pA0[1];
            rB1 = pB1[16];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[17];
            rA1 = pA0[1];
            rB1 = pB1[17];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[18];
            rA1 = pA0[1];
            rB1 = pB1[18];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[19];
            rA1 = pA0[1];
            rB1 = pB1[19];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[20];
            rA1 = pA0[1];
            rB1 = pB1[20];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[21];
            rA1 = pA0[1];
            rB1 = pB1[21];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[22];
            rA1 = pA0[1];
            rB1 = pB1[22];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[23];
            rA1 = pA0[1];
            rB1 = pB1[23];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            rC0_0 += m0;
            m0 = rA1 * rB0;
            rC1_0 += m0;
            m0 = rA0 * rB1;
            rC0_1 += m0;
            m0 = rA1 * rB1;
            pA0 += incAk;
            rA0 = *pA0;
            rB0 = pB0[24];
            rA1 = pA0[1];
            rB1 = pB1[24];
            rC1_1 += m0;
            m0 = rA0 * rB0;
            pB0 += incBk;
            pB1 += incBk;
         } /* end K-loop */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[1];
         rA1 = pA0[1];
         rB1 = pB1[1];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[2];
         rA1 = pA0[1];
         rB1 = pB1[2];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[3];
         rA1 = pA0[1];
         rB1 = pB1[3];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[4];
         rA1 = pA0[1];
         rB1 = pB1[4];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[5];
         rA1 = pA0[1];
         rB1 = pB1[5];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[6];
         rA1 = pA0[1];
         rB1 = pB1[6];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[7];
         rA1 = pA0[1];
         rB1 = pB1[7];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[8];
         rA1 = pA0[1];
         rB1 = pB1[8];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[9];
         rA1 = pA0[1];
         rB1 = pB1[9];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[10];
         rA1 = pA0[1];
         rB1 = pB1[10];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[11];
         rA1 = pA0[1];
         rB1 = pB1[11];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[12];
         rA1 = pA0[1];
         rB1 = pB1[12];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[13];
         rA1 = pA0[1];
         rB1 = pB1[13];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[14];
         rA1 = pA0[1];
         rB1 = pB1[14];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[15];
         rA1 = pA0[1];
         rB1 = pB1[15];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[16];
         rA1 = pA0[1];
         rB1 = pB1[16];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[17];
         rA1 = pA0[1];
         rB1 = pB1[17];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[18];
         rA1 = pA0[1];
         rB1 = pB1[18];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[19];
         rA1 = pA0[1];
         rB1 = pB1[19];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[20];
         rA1 = pA0[1];
         rB1 = pB1[20];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[21];
         rA1 = pA0[1];
         rB1 = pB1[21];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[22];
         rA1 = pA0[1];
         rB1 = pB1[22];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk;
         rA0 = *pA0;
         rB0 = pB0[23];
         rA1 = pA0[1];
         rB1 = pB1[23];
         rC1_1 += m0;
         m0 = rA0 * rB0;
         pB0 += (23*(incBk0));
         pB1 += (23*(incBk0));
/*
 *       Drain pipe on last iteration of K-loop
 */
         rC0_0 += m0;
         m0 = rA1 * rB0;
         rC1_0 += m0;
         m0 = rA0 * rB1;
         rC0_1 += m0;
         m0 = rA1 * rB1;
         pA0 += incAk0;
         rC1_1 += m0;
         pB0 += incBk0;
         pB1 += incBk0;
         rB0 = alpha;
         rC0_0 *= rB0;
         rC0_1 *= rB0;
         rC1_0 *= rB0;
         rC1_1 *= rB0;
         *pC0 = rC0_0;
         pC0[1] = rC1_0;
         *pC1 = rC0_1;
         pC1[1] = rC1_1;
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
