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
#ifndef ATL_XOVER_ISAMAX
#define    ATL_XOVER_ISAMAX    ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_IDAMAX
#define    ATL_XOVER_IDAMAX    ATL_XOVER_L1_DEFAULT
#endif
#elif defined( SCPLX )
#ifndef ATL_XOVER_ICAMAX
#define    ATL_XOVER_ICAMAX    ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DCPLX )
#ifndef ATL_XOVER_IZAMAX
#define    ATL_XOVER_IZAMAX    ATL_XOVER_L1_DEFAULT
#endif
#endif

PT_FUN_ARG_T Mjoin( PATL, ptamax0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_AMAX_T                  * arg = (PT_AMAX_T *)(ARGS);
   int                        indx;
/* ..
 * .. Executable Statements ..
 *
 */
   arg->indx = Mjoin( Mjoin( ATL_i, PRE ), amax )( arg->n, (TYPE *)(arg->x),
                                                   arg->incx );
   if( arg->n == 0 ) { *((TYPE *)(arg->amax)) = ATL_rzero; }
   else
   {
#ifdef TREAL
      indx = (arg->indx) * (arg->incx);
      *((TYPE *)(arg->amax)) = Mabs( ((TYPE *)(arg->x))[indx] );
#else
      indx = (arg->indx) * (arg->incx) * 2;
      *((TYPE *)(arg->amax)) = Mabs( ((TYPE *)(arg->x))[indx]   ) +
                               Mabs( ((TYPE *)(arg->x))[indx+1] );
#endif
   }

   return( NULL );
/*
 * End of Mjoin( PATL, ptamax0 )
 */
}

int Mjoin( Mjoin( Mjoin( ATL_, i ), PRE ), ptamax )
(
   const int                  N,
   const TYPE                 * X,
   const int                  INCX
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptamax )  returns the index in an n-vector x of the first ele-
 * ment having maximum absolute value.
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
   PT_AMAX_T                  arg[ATL_NTHREADS];
   TYPE                       amax[ATL_NTHREADS];
   pthread_attr_t             attr;
   TYPE                       smax;
   int                        i, ib, il, imax, nb, nn, nthreads;
#ifdef TREAL
#define  incx2                INCX
#else
   int                        incx2 = 2 * INCX;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if( N < 1 ) { return( 0 ); }
#if   defined( SREAL )
   else if( ( N <= ATL_XOVER_ISAMAX  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( ATL_i, PRE ), amax )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_ISAMAX <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_ISAMAX - 1 ) / ATL_XOVER_ISAMAX );
#elif defined( DREAL )
   else if( ( N <= ATL_XOVER_IDAMAX  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( ATL_i, PRE ), amax )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_IDAMAX <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_IDAMAX - 1 ) / ATL_XOVER_IDAMAX );
#elif defined( SCPLX )
   else if( ( N <= ATL_XOVER_ICAMAX  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( ATL_i, PRE ), amax )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_ICAMAX <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_ICAMAX - 1 ) / ATL_XOVER_ICAMAX );
#elif defined( DCPLX )
   else if( ( N <= ATL_XOVER_IZAMAX  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( ATL_i, PRE ), amax )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_IZAMAX <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_IZAMAX - 1 ) / ATL_XOVER_IZAMAX );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      amax[i] = ATL_rzero; arg[i].n = ( ib = Mmin( nb, nn ) );
      arg[i].x = X; arg[i].incx = INCX; arg[i].amax = (void *)(&amax[i]);
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptamax0 ),
                                  (void *)(&arg[i]) ));
      X += ib * incx2; nn -= nb;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * Combine the results
 */
   imax = 0; smax = ATL_rzero;

   for( i = 0, il = 0, nn = N; i < nthreads; i++ )
   {
      ib = Mmin( nn, nb );
      if( amax[i] > smax ) { imax = il + arg[i].indx; smax = amax[i]; }
      nn -= ib; il += ib;
   }

   return( imax );
/*
 * End of Mjoin( PATL, ptamax )
 */
}
