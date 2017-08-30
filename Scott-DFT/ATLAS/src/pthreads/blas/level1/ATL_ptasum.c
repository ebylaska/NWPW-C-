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
#ifndef ATL_XOVER_SASUM
#define    ATL_XOVER_SASUM     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_DASUM
#define    ATL_XOVER_DASUM     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( SCPLX )
#ifndef ATL_XOVER_SCASUM
#define    ATL_XOVER_SCASUM    ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DCPLX )
#ifndef ATL_XOVER_DZASUM
#define    ATL_XOVER_DZASUM    ATL_XOVER_L1_DEFAULT
#endif
#endif

PT_FUN_ARG_T Mjoin( PATL, ptasum0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_ASUM_T                  * arg = (PT_ASUM_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
#ifdef TREAL
   *((TYPE *)(arg->sum)) = Mjoin( PATL, asum )( arg->n, (TYPE *)(arg->x),
                                                arg->incx );
#else
   *((TYPE *)(arg->sum)) = Mjoin( Mjoin( PATLU, PRE ), asum )( arg->n,
                              (TYPE *)(arg->x), arg->incx );
#endif
   return( NULL );
/*
 * End of Mjoin( PATL, ptasum0 )
 */
}

#ifdef TREAL
TYPE Mjoin( PATL, ptasum )
#else
TYPE Mjoin( Mjoin( PATLU, PRE ), ptasum )
#endif
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
 * Mjoin( PATL,  returns the sum of absolute values of the entries of an
 * n-vector x.
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
   PT_ASUM_T                  arg[ATL_NTHREADS];
   TYPE                       sum[ATL_NTHREADS];
   pthread_attr_t             attr;
   TYPE                       s;
   int                        i, ib, nb, nn, nthreads;
#ifdef TREAL
#define  incx2                INCX
#else
   int                        incx2 = 2 * INCX;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if(    ( N  < 1 ) || ( INCX < 1 ) ) { return( ATL_rzero ); }
#if   defined( SREAL )
   else if( ( N <= ATL_XOVER_SASUM  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( PATL, asum )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_SASUM  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SASUM  - 1 ) / ATL_XOVER_SASUM  );
#elif defined( DREAL )
   else if( ( N <= ATL_XOVER_DASUM  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( PATL, asum )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_DASUM  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DASUM  - 1 ) / ATL_XOVER_DASUM  );
#elif defined( SCPLX )
   else if( ( N <= ATL_XOVER_SCASUM ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( PATLU, PRE ), asum )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_SCASUM <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SCASUM - 1 ) / ATL_XOVER_SCASUM );
#elif defined( DCPLX )
   else if( ( N <= ATL_XOVER_DZASUM  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( PATLU, PRE ), asum )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_DZASUM  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DZASUM  - 1 ) / ATL_XOVER_DZASUM  );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0; i < nthreads; i++ )
   {
      sum[i]      = ATL_rzero;          arg[i].n    = ( ib = Mmin( nb, nn ) );
      arg[i].x    = (void *)(X);        arg[i].incx = INCX;
      arg[i].sum  = (void *)(&sum[i]);
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptasum0 ),
                                  (void *)(&arg[i]) ));
      X += ib * incx2; nn -= ib;
   }

   for( i = 0; i < nthreads; i++ )
   {
      ATL_assert(!pthread_join( pid[i], NULL ));
   }

   ATL_thread_free( &attr );
/*
 * Combine the results
 */
   s = sum[0]; for( i = 1; i < nthreads; i++ ) { s += sum[i]; } return( s );
/*
 * End of Mjoin( PATL, ptasum )
 */
}
