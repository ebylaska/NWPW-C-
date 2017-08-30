
void Mjoin(PATL,tompgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const SCALAR alpha,
                       const TYPE *A, ATL_CINT lda, const TYPE *B, ATL_CINT ldb,
                       const SCALAR beta, TYPE *C, ATL_CINT ldc);
void Mjoin(PATL,tvompgemm)(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
                       ATL_CINT M, ATL_CINT N, ATL_CINT K, const void *alpha,
                       const void *A, ATL_CINT lda, const void *B, ATL_CINT ldb,
                       const void *beta, void *C, ATL_CINT ldc);
void Mjoin(PATL,tompsymm)
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,tomphemm)
   (const enum ATLAS_SIDE Side, const enum ATLAS_UPLO Uplo,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#endif
void Mjoin(PATL,tompsyr2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,tompher2k)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans,
    ATL_CINT N, ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const TYPE *B, ATL_CINT ldb, const TYPE beta, TYPE *C, ATL_CINT ldc);
#endif
void Mjoin(PATL,tomptrmm)
   (const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
    const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb);
void Mjoin(PATL,tomptrsm)
   (const enum ATLAS_SIDE side, const enum ATLAS_UPLO uplo,
    const enum ATLAS_TRANS TA, const enum ATLAS_DIAG diag,
    ATL_CINT M, ATL_CINT N, const SCALAR alpha,
    const TYPE *A, ATL_CINT lda, TYPE *B, ATL_CINT ldb);
void Mjoin(PATL,tompsyrk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const SCALAR alpha, const TYPE *A, ATL_CINT lda,
    const SCALAR beta, TYPE *C, ATL_CINT ldc);
#ifdef TCPLX
void Mjoin(PATL,tompherk)
   (const enum ATLAS_UPLO Uplo, const enum ATLAS_TRANS Trans, ATL_CINT N,
    ATL_CINT K, const TYPE alpha, const TYPE *A, ATL_CINT lda,
    const TYPE beta, TYPE *C, ATL_CINT ldc);
#endif


