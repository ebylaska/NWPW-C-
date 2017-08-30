#ifndef ATL_UCLEANK
   #define ATL_zgpKBmm ATL_zpKBmm
#endif

void ATL_zJIK0x0x47TN47x47x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x48TN48x48x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x49TN49x49x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x50TN50x50x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x51TN51x51x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x52TN52x52x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
void ATL_zJIK0x0x0TN0x0x0_a1_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);
typedef void (*MMfunc)(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc);

void ATL_zgpKBmm(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   static MMfunc mmfunc[  6] = {
                         ATL_zJIK0x0x47TN47x47x0_a1_bX,
                         ATL_zJIK0x0x48TN48x48x0_a1_bX,
                         ATL_zJIK0x0x49TN49x49x0_a1_bX,
                         ATL_zJIK0x0x50TN50x50x0_a1_bX,
                         ATL_zJIK0x0x51TN51x51x0_a1_bX,
                         ATL_zJIK0x0x52TN52x52x0_a1_bX,
                        };
   MMfunc mm;

   if (K <= 46)
   {
      ATL_zJIK0x0x0TN0x0x0_a1_bX(M, N, K, alpha, A, lda, B, ldb, -beta, C, ldc);
      ATL_zJIK0x0x0TN0x0x0_a1_bX(M, N, K, alpha, A, lda, B+N*ldb, ldb, beta, C+1, ldc);
      ATL_zJIK0x0x0TN0x0x0_a1_bX(M, N, K, alpha, A+M*lda, lda, B+N*ldb, ldb, -1.0, C, ldc);
      ATL_zJIK0x0x0TN0x0x0_a1_bX(M, N, K, alpha, A+M*lda, lda, B, ldb, 1.0, C+1, ldc);
   }
   else
   {
      mm = mmfunc[K-47];
      mm(M, N, K, alpha, A, lda, B, ldb, -beta, C, ldc);
      mm(M, N, K, alpha, A, lda, B+N*ldb, ldb, beta, C+1, ldc);
      mm(M, N, K, alpha, A+M*lda, lda, B+N*ldb, ldb, -1.0, C, ldc);
      mm(M, N, K, alpha, A+M*lda, lda, B, ldb, 1.0, C+1, ldc);
   }
}
#ifndef ATL_UCLEANK
void ATL_zpKBmm_b0(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   ATL_zgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
void ATL_zpKBmm_b1(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   ATL_zgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
void ATL_zpKBmm_bX(const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
   ATL_zgpKBmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
#endif
