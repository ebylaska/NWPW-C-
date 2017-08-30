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
#ifndef ATL_XOVER_SROTM
#define    ATL_XOVER_SROTM     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_DROTM
#define    ATL_XOVER_DROTM     ATL_XOVER_L1_DEFAULT
#endif
#endif

PT_FUN_ARG_T Mjoin( PATL, ptrotm0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_ROTM_T                  * arg = (PT_ROTM_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
   Mjoin( PATL, rotm )( arg->n, (TYPE *)(arg->x), arg->incx, (TYPE *)(arg->y),
                        arg->incy, (TYPE *)(arg->param) );
   return( NULL );
/*
 * End of Mjoin( PATL, ptrotm0 )
 */
}

void Mjoin( PATL, ptrotm )
(
   const int                  N,
   TYPE                       * X,
   const int                  INCX,
   TYPE                       * Y,
   const int                  INCY,
   const TYPE                 * PARAM
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptrotm ) applies a modified-Givens rotation to the two n-vectors
 * x and y.
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
   PT_ROTM_T                  arg[ATL_NTHREADS];
   pthread_attr_t             attr;
#define  incx2                INCX
#define  incy2                INCY
   int                        i, ib, nb, nn, nthreads;
/* ..
 * .. Executable Statements ..
 *
 */
#if   defined( SREAL )
   if( ( N <= 0 ) || ( PARAM[0] == 2.0f ) ) { return; }
   else if( ( N <= ATL_XOVER_SROTM ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, rotm )( N, X, INCX, Y, INCY, PARAM ); return; }

   nthreads = ( ATL_XOVER_SROTM <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SROTM - 1 ) / ATL_XOVER_SROTM );
#elif defined( DREAL )
   if( ( N <= 0 ) || ( PARAM[0] == 2.0 ) ) { return; }
   else if( ( N <= ATL_XOVER_DROTM ) || ( ATL_NTHREADS <= 1 ) )
   { Mjoin( PATL, rotm )( N, X, INCX, Y, INCY, PARAM ); return; }

   nthreads = ( ATL_XOVER_DROTM <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DROTM - 1 ) / ATL_XOVER_DROTM );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      arg[i].n = ( ib = Mmin( nb, nn ) );
      arg[i].x = (void *)(X); arg[i].incx = INCX;
      arg[i].y = (void *)(Y); arg[i].incy = INCY;
      arg[i].param = (void *)(PARAM);
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptrotm0 ),
                                  (void *)(&arg[i]) ));
      X += ib * incx2; Y += ib * incy2; nn -= ib;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * End of Mjoin( PATL, ptrotm )
 */
}
