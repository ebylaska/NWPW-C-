#define DREAL
#include "atlas_misc.h"
void ATL_dupMBmm0_2_0_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dgpMBmm_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);

void ATL_dpMBmm_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc)
{

   if (M == (((((M) >> 1)) << 1)))
   {
      ATL_dupMBmm0_2_0_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
   else ATL_dgpMBmm_bX(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
