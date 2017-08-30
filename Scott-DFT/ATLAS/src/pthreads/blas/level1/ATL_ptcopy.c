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
#ifndef ATL_XOVER_SCOPY
#define    ATL_XOVER_SCOPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_DCOPY
#define    ATL_XOVER_DCOPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( SCPLX )
#ifndef ATL_XOVER_CCOPY
#define    ATL_XOVER_CCOPY     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DCPLX )
#ifndef ATL_XOVER_ZCOPY
#define    ATL_XOVER_ZCOPY     ATL_XOVER_L1_DEFAULT
#endif
#endif

PT_FUN_ARG_T Mjoin( PATL, ptcopy0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_COPY_T                  * arg = (PT_COPY_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
   Mjoin( PATL, copy )( arg->n, (TYPE *)(arg->x), arg->incx,
                        (TYPE *)(arg->y), arg->incy );
   return( NULL );
/*
 * End of Mjoin( PATL, ptcopy0 )
 */
}

void Mjoin( PATL, ptcopy )
(
   const int                  N,
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
 * Mjoin( PATL, ptcopy ) copies an n-vector x into an n-vector y.
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
   PT_COPY_T                  arg[ATL_NTHREADS];
   pthread_attr_t             attr;
#ifdef TREAL
#define  incx2                INCX
#define  incy2                INCY
#else
   int                        incx2 = 2 * INCX, incy2 = 2 * INCY;
#endif
   int                        i, ib, nb, nn, nthreads;
/* ..
 * .. Executable Statements ..
 *
 */
   if( N <= 0 ) { return; }
#if   defined( SREAL )
   else if( ( N <= ATL_XOVER_SCOPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, copy )( N, X, INCX, Y, INCY ); return; }

   nthreads = ( ATL_XOVER_SCOPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SCOPY - 1 ) / ATL_XOVER_SCOPY );
#elif defined( DREAL )
   else if( ( N <= ATL_XOVER_DCOPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, copy )( N, X, INCX, Y, INCY ); return; }

   nthreads = ( ATL_XOVER_DCOPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DCOPY - 1 ) / ATL_XOVER_DCOPY );
#elif defined( SCPLX )
   else if( ( N <= ATL_XOVER_CCOPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, copy )( N, X, INCX, Y, INCY ); return; }

   nthreads = ( ATL_XOVER_CCOPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_CCOPY - 1 ) / ATL_XOVER_CCOPY );
#elif defined( DCPLX )
   else if( ( N <= ATL_XOVER_ZCOPY ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, copy )( N, X, INCX, Y, INCY ); return; }

   nthreads = ( ATL_XOVER_ZCOPY <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_ZCOPY - 1 ) / ATL_XOVER_ZCOPY );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      arg[i].n = ( ib = Mmin( nb, nn ) );
      arg[i].x = (void *)(X); arg[i].incx = INCX;
      arg[i].y = (void *)(Y); arg[i].incy = INCY;
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptcopy0 ),
                                  (void *)(&arg[i]) ));
      X += ib * incx2; Y += ib * incy2; nn -= ib;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * End of Mjoin( PATL, ptcopy )
 */
}
