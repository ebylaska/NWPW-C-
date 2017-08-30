/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef ATLAS_LAPACK_H
   #define ATLAS_LAPACK_H

#include "atlas_misc.h"
#include "cblas.h"

/*
 * Define LAPACK flag arguments as powers of two so we can | them together
 * for calls to ILAENV
 */
enum ATL_LAFLG
   {LAUpper=1, LALower=2, LARight=4, LALeft=8, LAUnit=16, LANonunit=32,
    LASreal=(1<<27), LADreal=(1<<28), LAScplx=(1<<29), LADcplx=(1<<30)};
/*
 * We can overload QR names to one by filling in LAFLG, giving Side for
 * the trapezoidal matrix, and type of triangle:
 *    LAormqr + LARight+LAUpper --> ORMQR
 *    LAormqr + LARight+LALower --> ORMQL
 *    LAormqr + LALeft +LAUpper --> ORMRQ
 *    LAormqr + LALeft +LALower --> ORMLQ
 *    LAgeqrf + LARight+LAUpper --> GEQRF
 *    LAgeqrf + LARight+LALower --> GEQLF
 *    LAgeqrf + LALeft +LAUpper --> GERQF
 *    LAgeqrf + LALeft +LALower --> GELQF
 */
enum ATL_LAROUT
   {LAunknown=0, LAgetrf=1, LAgeqrf=(1<<1),  /* handles ge[qr,rq,lq,ql]f */
    LAormqr=(1<<2),   /* handles all [[d,s]orm,[c,z]unm][qr,ql,rq,lq] */
    LArorgen=(1<<3), /* general [D,S]OR* routine needing constrained NB */
    LAcungen=(1<<4), /* general [Z,C]UN* routine needing constrained NB */
    LAgehrd=(1<<5), LAgebrd=(1<<6), LAgetri=(1<<7),
    LApotrf=(1<<8), LAsytrf=(1<<9), LAsytrd=(1<<10), LAhetrf=(1<<11),
    LAhetrd=(1<<12), LAhegst=(1<<13), LAhbgst=(1<<14), LAhpgst=(1<<15),
    LAspgst=(1<<16), LAsbgst=(1<<17), LAsygst=(1<<18), LAstebz=(1<<19),
    LAgbtrf=(1<<20), LApbtrf=(1<21), LAtrtri=(1<<22), LAlauum=(1<<23)
    };

enum ATL_ISPEC
   {LAIS_OPT_NB=1, LAIS_MIN_NB=2, LAIS_NBXOVER=3, LAIS_NEIGSHFT=4,
    LAIS_MINCSZ=5, LAIS_SVDXOVER=6, LAIS_NPROC=7, LAIS_MSQRXOVER=8,
    LAIS_MAXDCSPSZ=9, LAIS_NTNAN=10, LAIS_NTINF=11};
/*
 * Comments from lapack's ILAENV
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*/

enum ATL_LADIRECT
   {LAForward=1, LABackward=2 };
enum ATL_LASTOREV
   {LARowStore=1, LAColumnStore=2 };

#ifdef PATL

#include "atlas_cblastypealias.h"
/*
 * predefined type macro names
 */
#define ATL_printMat     Mjoin(PATL,printMat)
#define ATL_larfp        Mjoin(PATL,larfp)
#define ATL_lacgv        Mjoin(PATL,lacgv)
#define ATL_lapy3        Mjoin(PATL,lapy3)
#define ATL_lapy2        Mjoin(PATL,lapy2)
#define ATL_larft_block  Mjoin(PATL,larft_block)
#define ATL_ladiv        Mjoin(PATL,ladiv)
#define ATL_larft        Mjoin(PATL,larft)
#define ATL_larfg        Mjoin(PATL,larfg)
#define ATL_larf         Mjoin(PATL,larf)
#define ATL_larfb        Mjoin(PATL,larfb)
#define ATL_gelqr        Mjoin(PATL,gelqr)
#define ATL_gelqf        Mjoin(PATL,gelqf)
#define ATL_gelq2        Mjoin(PATL,gelq2)
#define ATL_geqlr        Mjoin(PATL,geqlr)
#define ATL_geqlf        Mjoin(PATL,geqlf)
#define ATL_geql2        Mjoin(PATL,geql2)
#define ATL_gerqr        Mjoin(PATL,gerqr)
#define ATL_gerqf        Mjoin(PATL,gerqf)
#define ATL_gerq2        Mjoin(PATL,gerq2)
#define ATL_geqrr        Mjoin(PATL,geqrr)
#define ATL_geqrf        Mjoin(PATL,geqrf)
#define ATL_geqr2        Mjoin(PATL,geqr2)
#define ATL_getriR       Mjoin(PATL,getriR)
#define ATL_getriC       Mjoin(PATL,getriC)
#define ATL_getri        Mjoin(PATL,getri)
#define ATL_lauumRL      Mjoin(PATL,lauumRL)
#define ATL_lauumRU      Mjoin(PATL,lauumRU)
#define ATL_lauumCL      Mjoin(PATL,lauumCL)
#define ATL_lauumCU      Mjoin(PATL,lauumCU)
#define ATL_lauum        Mjoin(PATL,lauum)
#define ATL_trtriRL      Mjoin(PATL,trtriRL)
#define ATL_trtriRU      Mjoin(PATL,trtriRU)
#define ATL_trtriCL      Mjoin(PATL,trtriCL)
#define ATL_trtriCU      Mjoin(PATL,trtriCU)
#define ATL_trtri        Mjoin(PATL,trtri)
#define ATL_potrfU       Mjoin(PATL,potrfU)
#define ATL_potrfL       Mjoin(PATL,potrfL)
#define ATL_potrs        Mjoin(PATL,potrs)
#define ATL_potrf        Mjoin(PATL,potrf)
#define ATL_getrfR       Mjoin(PATL,getrfR)
#define ATL_getrfC       Mjoin(PATL,getrfC)
#define ATL_getrs        Mjoin(PATL,getrs)
#define ATL_getrf        Mjoin(PATL,getrf)
#define ATL_laswp        Mjoin(PATL,laswp)
#define ATL_lamch        Mjoin(Mjoin(ATL_,UPR),lamch)

#endif

#ifdef ATL_USEPTHREADS
   #define ATL_ilaenv ATL_itlaenv
#endif
int clapack_ilaenv(enum ATL_ISPEC ISPEC, enum ATL_LAROUT ROUT,
                   unsigned int OPTS, int N1, int N2, int N3, int N4);
int ATL_ilaenv(enum ATL_ISPEC ISPEC, enum ATL_LAROUT ROUT, unsigned int OPTS,
               int N1, int N2, int N3, int N4);
int ATL_sgetri(const enum CBLAS_ORDER Order, const int N, float *A, const int lda,
               const int *ipiv, float *wrk, int *lwrk);
int ATL_sgetriR(const int N, float *A, const int lda, const int *ipiv,
                float *wrk, const int lwrk);
int ATL_sgetriC(const int N, float *A, const int lda, const int *ipiv,
                float *wrk, const int lwrk);
void ATL_slauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, float *A, const int lda);
int ATL_spotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, float *A, const int lda);
void ATL_spotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const float *A, const int lda,
                float *B, const int ldb);
int ATL_sgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               float *A, const int lda, int *ipiv);
void ATL_sgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const float *A, const int lda,
                const int *ipiv, float *B, const int ldb);
void ATL_slaswp(const int N, float *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_sgetrfC(const int M, const int N, float *A, const int lda,
                int *ipiv);
int ATL_sgetrfR(const int M, const int N, float *A, const int lda,
                int *ipiv);
void ATL_slauumRU(const int N, float *A, const int lda);
void ATL_slauumRL(const int N, float *A, const int lda);
void ATL_slauumCU(const int N, float *A, const int lda);
void ATL_slauumCL(const int N, float *A, const int lda);
int ATL_sgeqr2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_sgerq2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_sgeql2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_sgelq2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_sgeqrf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_sgerqf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_sgeqlf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_sgelqf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_sgeqrr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_sgeqlr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_sgelqr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_sgerqr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
void ATL_slarfb(const enum CBLAS_SIDE SIDE, const enum CBLAS_TRANSPOSE TRANS,
                const enum ATL_LADIRECT  DIRECT, const enum ATL_LASTOREV STOREV,
                int M, int N, int K, float *V, int LDV, float *T, int LDT,
                float *C, int LDC, float  *WORK, int LDWORK);
void ATL_slarfg( const int N, float *ALPHA, float *X, int INCX, float *TAU);
void ATL_slarfp( const int N, float *ALPHA, float *X, int INCX, float *TAU);
void ATL_slarft_block(const enum ATL_LADIRECT  DIRECT,
                      const enum ATL_LASTOREV STOREV, int N, int K,
                      int left, int right, float *V, int LDV,
                      float *T, int LDT );
void ATL_slarft(const enum ATL_LADIRECT  DIRECT,
                const enum ATL_LASTOREV STOREV, int N,
                int K, float *V, int LDV, float  *TAU, float *T, int LDT );
void ATL_sprintMat
   (char *mat, const int M, const int N, float *A, const int lda0);
float  ATL_slapy2(float X, float Y);
void ATL_slarf(const enum CBLAS_SIDE SIDE, const int M, const int N,
               float *V, int INCV, float TAU, float *C,
               int LDC, float *WORK);
float ATL_slamch(char);
int ATL_spotrfU(const int N, float *A, const int lda);
int ATL_spotrfL(const int N, float *A, const int lda);
int ATL_strtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               float *A, const int lda);
int ATL_strtriRU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriRL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriCU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_strtriCL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);

int ATL_dgetri(const enum CBLAS_ORDER Order, const int N, double *A, const int lda,
               const int *ipiv, double *wrk, int *lwrk);
int ATL_dgetriR(const int N, double *A, const int lda, const int *ipiv,
                double *wrk, const int lwrk);
int ATL_dgetriC(const int N, double *A, const int lda, const int *ipiv,
                double *wrk, const int lwrk);
void ATL_dlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, double *A, const int lda);
int ATL_dpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, double *A, const int lda);
void ATL_dpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const double *A, const int lda,
                double *B, const int ldb);
int ATL_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               double *A, const int lda, int *ipiv);
void ATL_dgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const double *A, const int lda,
                const int *ipiv, double *B, const int ldb);
void ATL_dlaswp(const int N, double *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_dgetrfC(const int M, const int N, double *A, const int lda,
                int *ipiv);
int ATL_dgetrfR(const int M, const int N, double *A, const int lda,
                int *ipiv);
void ATL_dlauumRU(const int N, double *A, const int lda);
void ATL_dlauumRL(const int N, double *A, const int lda);
void ATL_dlauumCU(const int N, double *A, const int lda);
void ATL_dlauumCL(const int N, double *A, const int lda);
int ATL_dgeqr2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_dgerq2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_dgeql2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_dgelq2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_dgeqrf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_dgerqf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_dgeqlf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_dgelqf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_dgeqrr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_dgeqlr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_dgelqr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_dgerqr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
void ATL_dlarfb(const enum CBLAS_SIDE SIDE, const enum CBLAS_TRANSPOSE TRANS,
                const enum ATL_LADIRECT  DIRECT, const enum ATL_LASTOREV STOREV,
                int M, int N, int K, double *V, int LDV, double *T, int LDT,
                double *C, int LDC, double  *WORK, int LDWORK);
void ATL_dlarfg( const int N, double *ALPHA, double *X, int INCX, double *TAU);
void ATL_dlarfp( const int N, double *ALPHA, double *X, int INCX, double *TAU);
void ATL_dlarft_block(const enum ATL_LADIRECT  DIRECT,
                      const enum ATL_LASTOREV STOREV, int N, int K,
                      int left, int right, double *V, int LDV,
                      double *T, int LDT );
void ATL_dlarft(const enum ATL_LADIRECT  DIRECT,
                const enum ATL_LASTOREV STOREV, int N,
                int K, double *V, int LDV, double  *TAU, double *T, int LDT );
void ATL_dprintMat
   (char *mat, const int M, const int N, double *A, const int lda0);
double  ATL_dlapy2(double X, double Y);
void ATL_dlarf(const enum CBLAS_SIDE SIDE, const int M, const int N,
               double *V, int INCV, double TAU, double *C,
               int LDC, double *WORK);
double ATL_dlamch(char);
int ATL_dpotrfU(const int N, double *A, const int lda);
int ATL_dpotrfL(const int N, double *A, const int lda);
int ATL_dtrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               double *A, const int lda);
int ATL_dtrtriRU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriRL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriCU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_dtrtriCL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);

int ATL_cgetri(const enum CBLAS_ORDER Order, const int N, float *A, const int lda,
               const int *ipiv, float *wrk, int *lwrk);
int ATL_cgetriR(const int N, float *A, const int lda, const int *ipiv,
                float *wrk, const int lwrk);
int ATL_cgetriC(const int N, float *A, const int lda, const int *ipiv,
                float *wrk, const int lwrk);
void ATL_clauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, float *A, const int lda);
int ATL_cpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, float *A, const int lda);
void ATL_cpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const float *A, const int lda,
                float *B, const int ldb);
int ATL_cgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               float *A, const int lda, int *ipiv);
void ATL_cgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const float *A, const int lda,
                const int *ipiv, float *B, const int ldb);
void ATL_claswp(const int N, float *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_cgetrfC(const int M, const int N, float *A, const int lda,
                int *ipiv);
int ATL_cgetrfR(const int M, const int N, float *A, const int lda,
                int *ipiv);
void ATL_clauumRU(const int N, float *A, const int lda);
void ATL_clauumRL(const int N, float *A, const int lda);
void ATL_clauumCU(const int N, float *A, const int lda);
void ATL_clauumCL(const int N, float *A, const int lda);
int ATL_cgeqr2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_cgerq2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_cgeql2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_cgelq2( const int M, const int N, float *A, int LDA, float  *TAU, float *WORK);
int ATL_cgeqrf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_cgerqf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_cgeqlf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_cgelqf(int M, int N, float  *A, int LDA, float  *TAU,
               float *WORK, int LWORK);
int ATL_cgeqrr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_cgeqlr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_cgelqr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
int ATL_cgerqr(int M, int N, float *A, int LDA, float  *TAU,
               float *ws_QR2, float *ws_T, int LDT,
               float *WORKM, int buildT);
void ATL_clarfb(const enum CBLAS_SIDE SIDE, const enum CBLAS_TRANSPOSE TRANS,
                const enum ATL_LADIRECT  DIRECT, const enum ATL_LASTOREV STOREV,
                int M, int N, int K, float *V, int LDV, float *T, int LDT,
                float *C, int LDC, float  *WORK, int LDWORK);
void ATL_clarfg( const int N, float *ALPHA, float *X, int INCX, float *TAU);
void ATL_clarfp( const int N, float *ALPHA, float *X, int INCX, float *TAU);
void ATL_clarft_block(const enum ATL_LADIRECT  DIRECT,
                      const enum ATL_LASTOREV STOREV, int N, int K,
                      int left, int right, float *V, int LDV,
                      float *T, int LDT );
void ATL_clarft(const enum ATL_LADIRECT  DIRECT,
                const enum ATL_LASTOREV STOREV, int N,
                int K, float *V, int LDV, float  *TAU, float *T, int LDT );
void ATL_cprintMat
   (char *mat, const int M, const int N, float *A, const int lda0);
float  ATL_clapy2(float X, float Y);
int ATL_cpotrfRU(const int N, float *A, const int lda);
int ATL_cpotrfRL(const int N, float *A, const int lda);
void ATL_cladiv( float *X, float *Y, float  *Z);
void  ATL_clacgv( const int N, float *X, int INCX);
float  ATL_clapy3(float X, float Y, float Z);
void ATL_clarf(const enum CBLAS_SIDE  SIDE, const int M, const int N,
               float *V, int INCV, float *TAU, float *C, int LDC,
               float *WORK);
int ATL_cpotrfU(const int N, float *A, const int lda);
int ATL_cpotrfL(const int N, float *A, const int lda);
int ATL_ctrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               float *A, const int lda);
int ATL_ctrtriRU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriRL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriCU(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);
int ATL_ctrtriCL(const enum CBLAS_DIAG Diag, const int N, float *A,
                 const int lda);

int ATL_zgetri(const enum CBLAS_ORDER Order, const int N, double *A, const int lda,
               const int *ipiv, double *wrk, int *lwrk);
int ATL_zgetriR(const int N, double *A, const int lda, const int *ipiv,
                double *wrk, const int lwrk);
int ATL_zgetriC(const int N, double *A, const int lda, const int *ipiv,
                double *wrk, const int lwrk);
void ATL_zlauum(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, double *A, const int lda);
int ATL_zpotrf(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const int N, double *A, const int lda);
void ATL_zpotrs(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const int NRHS, const double *A, const int lda,
                double *B, const int ldb);
int ATL_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
               double *A, const int lda, int *ipiv);
void ATL_zgetrs(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
                const int N, const int NRHS, const double *A, const int lda,
                const int *ipiv, double *B, const int ldb);
void ATL_zlaswp(const int N, double *A, const int lda0, const int K1,
               const int K2, const int *ipiv, const int inci);
int ATL_zgetrfC(const int M, const int N, double *A, const int lda,
                int *ipiv);
int ATL_zgetrfR(const int M, const int N, double *A, const int lda,
                int *ipiv);
void ATL_zlauumRU(const int N, double *A, const int lda);
void ATL_zlauumRL(const int N, double *A, const int lda);
void ATL_zlauumCU(const int N, double *A, const int lda);
void ATL_zlauumCL(const int N, double *A, const int lda);
int ATL_zgeqr2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_zgerq2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_zgeql2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_zgelq2( const int M, const int N, double *A, int LDA, double  *TAU, double *WORK);
int ATL_zgeqrf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_zgerqf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_zgeqlf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_zgelqf(int M, int N, double  *A, int LDA, double  *TAU,
               double *WORK, int LWORK);
int ATL_zgeqrr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_zgeqlr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_zgelqr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
int ATL_zgerqr(int M, int N, double *A, int LDA, double  *TAU,
               double *ws_QR2, double *ws_T, int LDT,
               double *WORKM, int buildT);
void ATL_zlarfb(const enum CBLAS_SIDE SIDE, const enum CBLAS_TRANSPOSE TRANS,
                const enum ATL_LADIRECT  DIRECT, const enum ATL_LASTOREV STOREV,
                int M, int N, int K, double *V, int LDV, double *T, int LDT,
                double *C, int LDC, double  *WORK, int LDWORK);
void ATL_zlarfg( const int N, double *ALPHA, double *X, int INCX, double *TAU);
void ATL_zlarfp( const int N, double *ALPHA, double *X, int INCX, double *TAU);
void ATL_zlarft_block(const enum ATL_LADIRECT  DIRECT,
                      const enum ATL_LASTOREV STOREV, int N, int K,
                      int left, int right, double *V, int LDV,
                      double *T, int LDT );
void ATL_zlarft(const enum ATL_LADIRECT  DIRECT,
                const enum ATL_LASTOREV STOREV, int N,
                int K, double *V, int LDV, double  *TAU, double *T, int LDT );
void ATL_zprintMat
   (char *mat, const int M, const int N, double *A, const int lda0);
double  ATL_zlapy2(double X, double Y);
int ATL_zpotrfRU(const int N, double *A, const int lda);
int ATL_zpotrfRL(const int N, double *A, const int lda);
void ATL_zladiv( double *X, double *Y, double  *Z);
void  ATL_zlacgv( const int N, double *X, int INCX);
double  ATL_zlapy3(double X, double Y, double Z);
void ATL_zlarf(const enum CBLAS_SIDE  SIDE, const int M, const int N,
               double *V, int INCV, double *TAU, double *C, int LDC,
               double *WORK);
int ATL_zpotrfU(const int N, double *A, const int lda);
int ATL_zpotrfL(const int N, double *A, const int lda);
int ATL_ztrtri(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
               const enum CBLAS_DIAG Diag, const int N,
               double *A, const int lda);
int ATL_ztrtriRU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriRL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriCU(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);
int ATL_ztrtriCL(const enum CBLAS_DIAG Diag, const int N, double *A,
                 const int lda);

#endif
