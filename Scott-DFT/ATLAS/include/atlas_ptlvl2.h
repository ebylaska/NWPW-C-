/* ---------------------------------------------------------------------
 *
 * -- Automatically Tuned Linear Algebra Software (ATLAS)
 *    (C) Copyright 2000 All Rights Reserved
 *
 * -- ATLAS routine -- Version 3.2 -- December 25, 2000
 *
 * Author         : Antoine P. Petitet
 * Originally developed at the University of Tennessee,
 * Innovative Computing Laboratory, Knoxville TN, 37996-1301, USA.
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the ATLAS group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------
 */
#ifndef ATLAS_PTLVL2_H
#define ATLAS_PTLVL2_H
#include "atlas_ptmisc.h"
/*
 * =====================================================================
 * Include files
 * =====================================================================
 */
#include "atlas_ptmisc.h"

#include "atlas_level2.h"
/*
 * =====================================================================
 * macro constants
 * =====================================================================
 */

/*
 * =====================================================================
 * macro functions
 * =====================================================================
 */

/*
 * =====================================================================
 * typedef definitions
 * =====================================================================
 */
typedef struct
{
   const void                 * a, * al, * x, * be;
   void                       * y;
   enum ATLAS_TRANS           ta;
   int                        incx, incy, kl, ku, la, m, n;
} PT_GBMV_T;

typedef struct
{
   const void                 * a, * al, * x, * be;
   void                       * y;
   enum ATLAS_TRANS           ta;
   int                        incx, incy, la, m, n;
} PT_GEMV_T;

typedef struct
{
   const void                 * al, * x, * y;
   void                       * a;
   int                        incx, incy, la, m, n;
} PT_GER_T;

typedef struct
{
   const void                 * al, * x, * y;
   void                       * a;
   enum ATLAS_UPLO            up;
   int                        incx, incy, la, n;
} PT_SPR2_T;

typedef struct
{
   const void                 * al, * x, * y;
   void                       * a;
   enum ATLAS_UPLO            up;
   int                        incx, incy, la, n;
} PT_SYR2_T;

typedef struct
{
   const void                 * al, * x;
   void                       * a;
   enum ATLAS_UPLO            up;
   int                        incx, la, n;
} PT_SPR_T;

typedef struct
{
   const void                 * al, * x;
   void                       * a;
   enum ATLAS_UPLO            up;
   int                        incx, la, n;
} PT_SYR_T;

typedef struct
{
   const void                 * a, * al, * x, * be;
   void                       * y;
   enum ATLAS_UPLO            up;
   int                        incx, incy, k, la, n;
} PT_SBMV_T;

typedef struct
{
   const void                 * a, * al, * x, * be;
   void                       * y;
   enum ATLAS_UPLO            up;
   int                        incx, incy, la, n;
} PT_SPMV_T;

typedef struct
{
   const void                 * a, * al, * x, * be;
   void                       * y;
   enum ATLAS_UPLO            up;
   int                        incx, incy, la, n;
} PT_SYMV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, k, la, n;
} PT_TBMV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, la, n;
} PT_TPMV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, la, n;
} PT_TRMV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, k, la, n;
} PT_TBSV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, la, n;
} PT_TPSV_T;

typedef struct
{
   const void                 * a;
   void                       * x;
   enum ATLAS_UPLO            up;
   enum ATLAS_TRANS           ta;
   enum ATLAS_DIAG            di;
   int                        incx, la, n;
} PT_TRSV_T;

/*
 * =====================================================================
 * Function prototypes
 * =====================================================================
 */
PT_FUN_ARG_T      Mjoin( PATL, ptgbmv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptgemv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttbmv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttpmv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttrmv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttbsv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttpsv0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pttrsv0   )          ( PT_FUN_ARG_T );

#ifdef TREAL

PT_FUN_ARG_T      Mjoin( PATL, ptger0    )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptsbmv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptspmv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptsymv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptspr0    )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptsyr0    )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptspr20   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptsyr20   )          ( PT_FUN_ARG_T );

#else

PT_FUN_ARG_T      Mjoin( PATL, ptgerc0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptgeru0   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pthbmv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pthpmv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pthemv20  )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pthpr0    )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pther0    )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pthpr20   )          ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, pther20   )          ( PT_FUN_ARG_T );

#endif
/*
 * =====================================================================
 * Prototypes for the Level 2 multi-threaded ATLAS BLAS routines
 * =====================================================================
 */
void              Mjoin( PATL, ptgbmv )
(  const enum ATLAS_TRANS,           const int,       const int,
   const int,       const int,       const SCALAR,    const TYPE *,
   const int,       const TYPE *,    const int,       const SCALAR,
   TYPE *,          const int );
void              Mjoin( PATL, ptgemv )
(  const enum ATLAS_TRANS,           const int,       const int,
   const SCALAR,    const TYPE *,    const int,       const TYPE *,
   const int,       const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, pttbmv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const TYPE *,    const int,       TYPE *,          const int );
void              Mjoin( PATL, pttpmv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const TYPE *,
   TYPE *,          const int );
void              Mjoin( PATL, pttrmv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const TYPE *,
   const int,       TYPE *,          const int );
void              Mjoin( PATL, pttbsv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const int,
   const TYPE *,    const int,       TYPE *,          const int );
void              Mjoin( PATL, pttpsv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const TYPE *,
   TYPE *,          const int );
void              Mjoin( PATL, pttrsv )
(  const enum ATLAS_UPLO,            const enum ATLAS_TRANS,
   const enum ATLAS_DIAG,            const int,       const TYPE *,
   const int,       TYPE *,          const int );

#ifdef TREAL

void              Mjoin( PATL, ptger  )
(  const int,       const int,       const SCALAR,    const TYPE *,
   const int,       const TYPE *,    const int,       TYPE *,
   const int );
void              Mjoin( PATL, ptsbmv )
(  const enum ATLAS_UPLO,            const int,       const int,
   const SCALAR,    const TYPE *,    const int,       const TYPE *,
   const int,       const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, ptspmv )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const TYPE *,    const int,       const SCALAR,
   TYPE *,          const int );
void              Mjoin( PATL, ptsymv )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, ptspr  )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       TYPE * );
void              Mjoin( PATL, ptsyr  )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       TYPE *,          const int );
void              Mjoin( PATL, ptspr2 )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   TYPE * );
void              Mjoin( PATL, ptsyr2 )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   TYPE *,          const int );

#else

void              Mjoin( PATL, ptgerc )
(  const int,       const int,       const SCALAR,    const TYPE *,
   const int,       const TYPE *,    const int,       TYPE *,
   const int );
void              Mjoin( PATL, ptgeru )
(  const int,       const int,       const SCALAR,    const TYPE *,
   const int,       const TYPE *,    const int,       TYPE *,
   const int );
void              Mjoin( PATL, pthbmv )
(  const enum ATLAS_UPLO,            const int,       const int,
   const SCALAR,    const TYPE *,    const int,       const TYPE *,
   const int,       const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, pthpmv )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const TYPE *,    const int,       const SCALAR,
   TYPE *,          const int );
void              Mjoin( PATL, pthemv )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, pthpr  )
(  const enum ATLAS_UPLO,            const int,       const TYPE,
   const TYPE *,    const int,       TYPE * );
void              Mjoin( PATL, pther  )
(  const enum ATLAS_UPLO,            const int,       const TYPE,
   const TYPE *,    const int,       TYPE *,          const int );
void              Mjoin( PATL, pthpr2 )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   TYPE * );
void              Mjoin( PATL, pther2 )
(  const enum ATLAS_UPLO,            const int,       const SCALAR,
   const TYPE *,    const int,       const TYPE *,    const int,
   TYPE *,          const int );

#endif

#endif
/*
 * End of atlas_ptlvl2.h
 */
