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
#ifndef ATL_PTLVL1_H
#define ATL_PTLVL1_H
/*
 * =====================================================================
 * Include files
 * =====================================================================
 */
#include "atlas_ptmisc.h"

#include "atlas_level1.h"
/*
 * =====================================================================
 * Cross-over points (default)
 * =====================================================================
 */
#define    ATL_XOVER_L1_DEFAULT         128
/*
 * =====================================================================
 * typedef definitions
 * =====================================================================
 */
typedef struct
{
   const void                 * x;
   void                       * scale, * ssq;
   int                        incx, n;
} PT_NRM2_T;

typedef struct
{
   const void                 * x;
   void                       * sum;
   int                        incx, n;
} PT_ASUM_T;

typedef struct
{
   const void                 * x;
   void                       * amax;
   int                        incx, indx, n;
} PT_AMAX_T;

typedef struct
{
   const void                 * alpha;
   void                       * x;
   int                        incx, n;
} PT_SCAL_T;

typedef struct
{
   const void                 * alpha, * x;
   void                       * y;
   int                        incx, incy, n;
} PT_AXPY_T;

typedef struct
{
   const void                 * x;
   void                       * y;
   int                        incx, incy, n;
} PT_COPY_T;

typedef struct
{
   void                       * x, * y;
   int                        incx, incy, n;
} PT_SWAP_T;

typedef struct
{
   const void                 * c, * s;
   void                       * x, * y;
   int                        incx, incy, n;
} PT_ROT_T;

typedef struct
{
   const void                 * param;
   void                       * x, * y;
   int                        incx, incy, n;
} PT_ROTM_T;

typedef struct
{
   const void                 * x, * y;
   void                       * dot;
   int                        incx, incy, n;
} PT_DOT_T;
/*
 * =====================================================================
 * Function prototypes
 * =====================================================================
 */
PT_FUN_ARG_T      Mjoin( PATL, ptnrm20     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptasum0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptamax0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptscal0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptaxpy0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptcopy0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptswap0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptrot0      )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptrotm0     )        ( PT_FUN_ARG_T );
#ifdef TREAL
PT_FUN_ARG_T      Mjoin( PATL, ptdot0      )        ( PT_FUN_ARG_T );
#else
PT_FUN_ARG_T      Mjoin( PATL, ptdotc0     )        ( PT_FUN_ARG_T );
PT_FUN_ARG_T      Mjoin( PATL, ptdotu0     )        ( PT_FUN_ARG_T );
#endif
PT_FUN_ARG_T      ATL_dsptdot0                      ( PT_FUN_ARG_T );
/*
 * =====================================================================
 * Prototypes for the Level 1 multi-threaded ATLAS BLAS routines
 * =====================================================================
 */
double            ATL_dsptdot
(  const int,       const float *,   const int,       const float *,
   const int );
float             ATL_sdsptdot
(  const int,       const float,     const float *,   const int,
   const float *,   const int );

#ifdef TREAL

void              Mjoin( PATL, ptrotg )
(  TYPE *,          TYPE *,          TYPE *,          TYPE * );
void              Mjoin( PATL, ptrotmg )
(  TYPE *,          TYPE *,          TYPE *,          const TYPE,
   TYPE * );
TYPE              Mjoin( PATL, ptnrm2 )
(  const int,       const TYPE *,    const int );
TYPE              Mjoin( PATL, ptasum )
(  const int,       const TYPE *,    const int );
void              Mjoin( PATL, ptrot  )
(  const int,       TYPE *,          const int,       TYPE *,
   const int,       const SCALAR,    const SCALAR );
void              Mjoin( PATL, ptrotm )
(  const int,       TYPE *,          const int,       TYPE *,
   const int,       const TYPE * );
TYPE              Mjoin( PATL, ptdot  )
(  const int,       const TYPE *,    const int,       const TYPE *,
   const int );

#else

void              Mjoin( PATL, ptrotg )
(  TYPE  *,         const TYPE *,    TYPE *,          TYPE * );
TYPE              Mjoin( Mjoin( PATLU, PRE ), ptnrm2 )
(  const int,       const TYPE *,    const int );
TYPE              Mjoin( Mjoin( PATLU, PRE ), ptasum )
(  const int,       const TYPE *,    const int );
void              Mjoin( Mjoin( PATL,  UPR ), ptscal )
(  const int,       const TYPE,      TYPE *,          const int );
void              Mjoin( Mjoin( PATL,  UPR ), ptrot  )
(  const int,       TYPE *,          const int,       TYPE *,
   const int,       const TYPE,      const TYPE );
void              Mjoin( PATL, ptdotc_sub  )
(  const int,       const TYPE *,    const int,       const TYPE *,
   const int,       SCALAR );
void              Mjoin( PATL, ptdotu_sub  )
(  const int,       const TYPE *,    const int,       const TYPE *,
   const int,       SCALAR );

#endif

int               Mjoin( Mjoin( Mjoin( ATL_, i ), PRE ), ptamax )
(  const int,       const TYPE *,    const int );
void              Mjoin( PATL, ptscal )
(  const int,       const SCALAR,    TYPE *,          const int );
void              Mjoin( PATL, ptaxpy )
(  const int,       const SCALAR,    const TYPE *,    const int,
   TYPE *,          const int );
void              Mjoin( PATL, ptcopy )
(  const int,       const TYPE *,    const int,       TYPE *,
   const int );
void              Mjoin( PATL, ptswap )
(  const int,       TYPE *,          const int,       TYPE *,
   const int );

#endif
/*
 * End of atlas_ptlvl1.h
 */
