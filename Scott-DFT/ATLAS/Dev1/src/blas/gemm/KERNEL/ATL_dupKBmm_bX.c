#define DREAL
#include "atlas_misc.h"
void ATL_dupKBmm2_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm4_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm6_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm8_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm10_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm12_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm14_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm16_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm18_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm20_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm22_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm24_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm26_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm28_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm30_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm32_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm34_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm36_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm38_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm40_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm42_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm44_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dupKBmm46_2_1_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dpKBmm_bX
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);
void ATL_dgpKBmm
   (const int M, const int N, const int K, const double alpha,
    const double *A, const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc);

typedef void (*MMfunc)(const int, const int, const int, const TYPE,
                       const TYPE *, const int, const TYPE *, const int,
                       const TYPE, TYPE *, const int);

void ATL_dpKBmm_bX
   (const int M, const int N, const int K, const TYPE alpha,
    const TYPE *A, const int lda, const TYPE *B, const int ldb,
    const TYPE beta, TYPE *C, const int ldc)
{

   static MMfunc mmfunc[48] = 
   {
      NULL,
      ATL_dgpKBmm,
      ATL_dupKBmm2_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm4_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm6_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm8_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm10_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm12_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm14_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm16_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm18_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm20_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm22_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm24_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm26_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm28_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm30_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm32_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm34_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm36_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm38_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm40_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm42_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm44_2_1_bX,
      ATL_dgpKBmm,
      ATL_dupKBmm46_2_1_bX,
      ATL_dgpKBmm
   };

   mmfunc[K](M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
