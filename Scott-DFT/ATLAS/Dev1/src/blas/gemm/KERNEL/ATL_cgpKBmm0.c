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
static void ATL_cJIK0x0x0TN1x1x1_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=0, 
 * lda=0, ldb=0, ldc=0, mu=1, nu=1, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   #define Nb N
   #define Kb K
   const float *stM = A + (lda*Mb);
   const float *stN = B + (ldb*Nb);
   #define incAk 1
   const int incAm = ((lda) - Kb), incAn = -(Mb*lda);
   #define incBk 1
   const int incBm = -(Kb), incBn = (ldb);
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
         for (k=0; k < K; k++) /* easy loop to unroll */
         {
            rA0 = *pA0;
            rB0 = *pB0;
            rC0_0 += rA0 * rB0;
            pA0 += incAk;
            pB0 += incBk;
         }
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
static void ATL_cJIK0x0x0TN4x1x1_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=0, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=1, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   #define Nb N
   #define Kb K
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (lda*Mb);
   const float *stN = B + (ldb*Nb);
   #define incAk 1
   const int incAm = ((((lda) << 2)) - Kb), incAn = -(Mb*lda);
   #define incBk 1
   const int incBm = -(Kb), incBn = (ldb);
   #define incCm 8
   const int incCn = (((ldc) << 1)) - (((Mb) << 1));
   float *pC0=C;
   const float *pA0=A, *pA1=pA0+(lda), *pA2=pA1+(lda), *pA3=pA2+(lda);
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
            for (k=0; k < K; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rC0_0 += rA0 * rB0;
               rA1 = *pA1;
               rC1_0 += rA1 * rB0;
               rA2 = *pA2;
               rC2_0 += rA2 * rB0;
               rA3 = *pA3;
               rC3_0 += rA3 * rB0;
               pA0 += incAk;
               pA1 += incAk;
               pA2 += incAk;
               pA3 += incAk;
               pB0 += incBk;
            }
            *pC0 = rC0_0;
            pC0[2] = rC1_0;
            pC0[4] = rC2_0;
            pC0[6] = rC3_0;
            pC0 += incCm;
            pA0 += incAm;
            pA1 += incAm;
            pA2 += incAm;
            pA3 += incAm;
            pB0 += incBm;
         }
         while(pA0 != stM);
         pC0 += incCn;
         pA0 += incAn;
         pA1 += incAn;
         pA2 += incAn;
         pA3 += incAn;
         pB0 += incBn;
      }
      while(pB0 != stN);
   }
   if (k=M-Mb)
      ATL_cJIK0x0x0TN1x1x1_a1_bX(k, N, K, alpha, ca + (Mb*lda), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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
static void ATL_cJIK0x0x0TN1x2x1_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=0, 
 * lda=0, ldb=0, ldc=0, mu=1, nu=2, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   #define Mb M
   const int Nb = (N>>1)<<1;
   #define Kb K
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (lda*Mb);
   const float *stN = B + (ldb*Nb);
   #define incAk 1
   const int incAm = ((lda) - Kb), incAn = -(Mb*lda);
   #define incBk 1
   const int incBm = -(Kb), incBn = (((ldb) << 1));
   #define incCm 2
   const int incCn = (((ldc) << 2)) - (((Mb) << 1));
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
   const float *pA0=A;
   const float *pB0=B, *pB1=pB0+(ldb);
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
            for (k=0; k < K; k++) /* easy loop to unroll */
            {
               rA0 = *pA0;
               rB0 = *pB0;
               rC0_0 += rA0 * rB0;
               rB1 = *pB1;
               rC0_1 += rA0 * rB1;
               pA0 += incAk;
               pB0 += incBk;
               pB1 += incBk;
            }
            *pC0 = rC0_0;
            *pC1 = rC0_1;
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
   if (k=N-Nb)
      ATL_cJIK0x0x0TN1x1x1_a1_bX(M, k, K, alpha, ca, lda, cb + (Nb*ldb), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
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
void ATL_cJIK0x0x0TN0x0x0_a1_bX
   (const int M, const int N, const int K, const float alpha, const float * ATL_RESTRICT A, const int lda, const float * ATL_RESTRICT B, const int ldb, const float beta, float * ATL_RESTRICT C, const int ldc)
/*
 * matmul with TA=T, TB=N, MB=0, NB=0, KB=0, 
 * lda=0, ldb=0, ldc=0, mu=4, nu=2, ku=1, pf=0
 * Generated by ATLAS/tune/blas/gemm/emit_mm.c (3.9.23)
 */
{
   const int Mb = (M>>2)<<2;
   const int Nb = (N>>1)<<1;
   #define Kb K
   const float *ca=A, *cb=B;
   float *cc=C;
   const float *stM = A + (lda*Mb);
   const float *stN = B + (ldb*Nb);
   #define incAk 1
   const int incAm = ((((lda) << 2)) - Kb), incAn = -(Mb*lda);
   #define incBk 1
   const int incBm = -(Kb), incBn = (((ldb) << 1));
   #define incCm 8
   const int incCn = (((ldc) << 2)) - (((Mb) << 1));
   float *pC0=C, *pC1=pC0+(((ldc) << 1));
   const float *pA0=A, *pA1=pA0+(lda), *pA2=pA1+(lda), *pA3=pA2+(lda);
   const float *pB0=B, *pB1=pB0+(ldb);
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
               for (k=0; k < K; k++) /* easy loop to unroll */
               {
                  rA0 = *pA0;
                  rB0 = *pB0;
                  rC0_0 += rA0 * rB0;
                  rA1 = *pA1;
                  rC1_0 += rA1 * rB0;
                  rA2 = *pA2;
                  rC2_0 += rA2 * rB0;
                  rA3 = *pA3;
                  rC3_0 += rA3 * rB0;
                  rB1 = *pB1;
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
               }
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
   }
   if (k=N-Nb)
      ATL_cJIK0x0x0TN4x1x1_a1_bX(M, k, K, alpha, ca, lda, cb + (Nb*ldb), ldb, beta, cc + (((Nb*ldc) << 1)), ldc);
   if (Nb && (k=M-Mb))
      ATL_cJIK0x0x0TN1x2x1_a1_bX(k, Nb, K, alpha, ca + (Mb*lda), lda, cb, ldb, beta, cc + (((Mb) << 1)), ldc);
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