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
/*
 * Include files
 */
#include "atlas_ptmisc.h"
#include "atlas_ptlvl1.h"
#include "atlas_ptlevel1.h"

#if   defined( SREAL )
#ifndef ATL_XOVER_SAXPY
#define    ATL_XOVER_SAXPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_DAXPY
#define    ATL_XOVER_DAXPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( SCPLX )
#ifndef ATL_XOVER_CAXPY
#define    ATL_XOVER_CAXPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DCPLX )
#ifndef ATL_XOVER_ZAXPY
#define    ATL_XOVER_ZAXPY     ATL_XOVER_L1_DEFAULT
#endif
#endif

PT_FUN_ARG_T Mjoin( PATL, ptaxpy0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_AXPY_T                  * arg = (PT_AXPY_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
#ifdef TREAL
   Mjoin( PATL, axpy )( arg->n, (SCALAR)(*((TYPE *)(arg->alpha))),
      (TYPE *)(arg->x), arg->incx, (TYPE *)(arg->y), arg->incy );
#else
   Mjoin( PATL, axpy )( arg->n, (SCALAR)(arg->alpha), (TYPE *)(arg->x),
                        arg->incx, (TYPE *)(arg->y), arg->incy );
#endif
   return( NULL );
/*
 * End of Mjoin( PATL, ptaxpy0 )
 */
}

void Mjoin( PATL, ptaxpy )
(
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * X,
   const int                  INCX,
   TYPE                       * Y,
   const int                  INCY
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptaxpy ) performs the following operation:
 *
 *    y := y + alpha * x,
 *
 * where alpha is a scalar and x and y are two n-vectors.
 *
 * For a  more  detailed description of  the arguments of this function,
 * see the reference implementation in the  ATLAS/src/blas/reference di-
 * rectory.
 *
 * ---------------------------------------------------------------------
 */
/*
 * .. Local Variables ..
 */
   pthread_t                  pid[ATL_NTHREADS];
   PT_AXPY_T                  arg[ATL_NTHREADS];
   pthread_attr_t             attr;
#ifdef TREAL
   TYPE                       alpha0 = (TYPE)(ALPHA);
#define  incx2                INCX
#define  incy2                INCY
#endif
   TYPE                       * alpha;
#ifdef TCPLX
   int                        incx2 = 2 * INCX, incy2 = 2 * INCY;
#endif
   int                        i, ib, nb, nn, nthreads;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N <= 0 ) || SCALAR_IS_ZERO( ALPHA ) ) { return; }
#if   defined( SREAL )
   else if( ( N <= ATL_XOVER_SAXPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, axpy )( N, ALPHA, X, INCX, Y, INCY ); return; }

   alpha = (void *)(&alpha0);
   nthreads = ( ATL_XOVER_SAXPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SAXPY - 1 ) / ATL_XOVER_SAXPY );
#elif defined( DREAL )
   else if( ( N <= ATL_XOVER_DAXPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, axpy )( N, ALPHA, X, INCX, Y, INCY ); return; }

   alpha = (void *)(&alpha0);
   nthreads = ( ATL_XOVER_DAXPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DAXPY - 1 ) / ATL_XOVER_DAXPY );
#elif defined( SCPLX )
   else if( ( N <= ATL_XOVER_CAXPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, axpy )( N, ALPHA, X, INCX, Y, INCY ); return; }

   alpha = (void *)(ALPHA);
   nthreads = ( ATL_XOVER_CAXPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_CAXPY - 1 ) / ATL_XOVER_CAXPY );
#elif defined( DCPLX )
   else if( ( N <= ATL_XOVER_ZAXPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, axpy )( N, ALPHA, X, INCX, Y, INCY ); return; }

   alpha = (void *)(ALPHA);
   nthreads = ( ATL_XOVER_ZAXPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_ZAXPY - 1 ) / ATL_XOVER_ZAXPY );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      arg[i].n = ( ib = Mmin( nb, nn ) ); arg[i].alpha = alpha;
      arg[i].x = (void *)(X); arg[i].incx = INCX;
      arg[i].y = (void *)(Y); arg[i].incy = INCY;
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptaxpy0 ),
                                  (void *)(&arg[i]) ));
      X += ib * incx2; Y += ib * incy2; nn -= ib;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * End of Mjoin( PATL, ptaxpy )
 */
}
