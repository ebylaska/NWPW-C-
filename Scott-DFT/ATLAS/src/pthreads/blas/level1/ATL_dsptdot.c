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

#ifndef ATL_XOVER_DSDOT
#define    ATL_XOVER_DSDOT     ATL_XOVER_L1_DEFAULT
#endif

PT_FUN_ARG_T ATL_dsptdot0( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_DOT_T                     * arg = (PT_DOT_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
   *((double *)(arg->dot)) = ATL_dsdot( arg->n, (float *)(arg->x), arg->incx,
                                        (float *)(arg->y), arg->incy );
   return( NULL );
/*
 * End of ATL_dsptdot0
 */
}

double ATL_dsptdot
(
   const int                  N,
   const float                * X,
   const int                  INCX,
   const float                * Y,
   const int                  INCY
)
{
/*
 * Purpose
 * =======
 *
 * ATL_dsptdot  returns  the dot product x^T * y of two n-vectors x and
 * y.  The result is internally computed using double precision arithme-
 * tic.
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
   PT_DOT_T                   arg[ATL_NTHREADS];
   double                     dot[ATL_NTHREADS];
   pthread_attr_t             attr;
   double                     d;
   int                        i, ib, nb, nn, nthreads;
#define  incx2                INCX
#define  incy2                INCY
/* ..
 * .. Executable Statements ..
 *
 */
   if( N < 1 ) { return( 0.0 ); }
   else if( ( N <= ATL_XOVER_DSDOT  ) || ( ATL_NTHREADS <= 1 ) )
   { return( ATL_dsdot( N, X, INCX, Y, INCY ) ); }

   nthreads = ( ATL_XOVER_DSDOT  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DSDOT  - 1 ) / ATL_XOVER_DSDOT  );
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      dot[i]     = 0.0; arg[i].n = ( ib = Mmin( nb, nn ) );
      arg[i].x   = (void *)(X); arg[i].incx = INCX;
      arg[i].y   = (void *)(Y); arg[i].incy = INCY;
      arg[i].dot = (void *)(&dot[i]);
      ATL_assert(!pthread_create( &pid[i], &attr, ATL_dsptdot0,
                                  (void *)(&arg[i]) ));
      X += ib * incx2; Y += ib * incy2; nn -= ib;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * Combine the results
 */
   d = dot[0];
   for( i = 1; i < nthreads; i++ ) { d += dot[i]; }
   return( d );
/*
 * End of ATL_dsptdot
 */
}
