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
#ifndef ATL_XOVER_SNRM2
#define    ATL_XOVER_SNRM2     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DREAL )
#ifndef ATL_XOVER_DNRM2
#define    ATL_XOVER_DNRM2     ATL_XOVER_L1_DEFAULT
#endif
#elif defined( SCPLX )
#ifndef ATL_XOVER_SCNRM2
#define    ATL_XOVER_SCNRM2    ATL_XOVER_L1_DEFAULT
#endif
#elif defined( DCPLX )
#ifndef ATL_XOVER_DZNRM2
#define    ATL_XOVER_DZNRM2    ATL_XOVER_L1_DEFAULT
#endif
#endif

#if 0           /* temporarily disabled by RCW */
PT_FUN_ARG_T Mjoin( PATL, ptnrm20 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_NRM2_T                  * arg = (PT_NRM2_T *)(ARGS);
/* ..
 * .. Executable Statements ..
 *
 */
   Mjoin( PATL, ssq )( arg->n, (TYPE *)(arg->x), arg->incx,
                       (TYPE *)(arg->scale), (TYPE *)(arg->ssq) );
   return( NULL );
/*
 * End of Mjoin( PATL, ptnrm20 )
 */
}
#endif

#ifdef TREAL
TYPE Mjoin( PATL, ptnrm2 )
#else
TYPE Mjoin( Mjoin( PATLU, PRE ), ptnrm2 )
#endif
(
   const int                  N,
   const TYPE                 * X,
   const int                  INCX
)
{
#if 1  /* threading in nrm2 temporarily disabled by RCW on 5/21/08 */
   #ifdef TREAL
      return( Mjoin( PATL, nrm2 )( N, X, INCX ) );
   #else
      return( Mjoin( Mjoin( PATLU, PRE ), nrm2 )( N, X, INCX ) );
   #endif
#else
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptnrm2 ) computes the 2-norm of an n-vector x.
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
   PT_NRM2_T                  arg[ATL_NTHREADS];
   TYPE                       scs[ATL_NTHREADS << 1];
   pthread_attr_t             attr;
   TYPE                       scale, ssq, temp;
   int                        i, i2, ib, nb, nn, nthreads;
#ifdef TREAL
#define  incx2                INCX
#else
   int                        incx2 = 2 * INCX;
#endif
/* ..
 * .. Executable Statements ..
 *
 */
   if(    ( N  < 1 ) || ( INCX < 1 ) ) { return( ATL_rzero    ); }
#if   defined( SREAL )
   else if( N == 1 )                   { return( Mabs( X[0] ) ); }
   else if( ( N <= ATL_XOVER_SNRM2  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( PATL, nrm2 )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_SNRM2  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SNRM2  - 1 ) / ATL_XOVER_SNRM2  );
#elif defined( DREAL )
   else if( N == 1 )                   { return( Mabs( X[0] ) ); }
   else if( ( N <= ATL_XOVER_DNRM2  ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( PATL, nrm2 )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_DNRM2  <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DNRM2  - 1 ) / ATL_XOVER_DNRM2  );
#elif defined( SCPLX )
   else if( ( N <= ATL_XOVER_SCNRM2 ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( PATLU, PRE ), nrm2 )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_SCNRM2 <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_SCNRM2 - 1 ) / ATL_XOVER_SCNRM2 );
#elif defined( DCPLX )
   else if( ( N <= ATL_XOVER_DZNRM2 ) || ( ATL_NTHREADS <= 1 ) )
   { return( Mjoin( Mjoin( PATLU, PRE ), nrm2 )( N, X, INCX ) ); }

   nthreads = ( ATL_XOVER_DZNRM2 <= 0 ? ATL_NTHREADS :
                ( N + ATL_XOVER_DZNRM2 - 1 ) / ATL_XOVER_DZNRM2 );
#endif
   if( nthreads > ATL_NTHREADS ) { nthreads = ATL_NTHREADS; }

   nb = ( ( nn = N ) + nthreads - 1 ) / nthreads;

   ATL_thread_init( &attr );

   for( i = 0, i2 = 0; i < nthreads; i++, i2 += 2 )
   {
      scs[i2]  = ATL_rzero; scs[i2+1] = ATL_rone;
      arg[i].n = ( ib = Mmin( nb, nn ) );
      arg[i].x = (void *)(X); arg[i].incx = INCX;
      arg[i].scale = (void *)(&scs[i2]); arg[i].ssq = (void *)(&scs[i2+1]);
      ATL_assert(!pthread_create( &pid[i], &attr, Mjoin( PATL, ptnrm20 ),
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
   scale = scs[0]; ssq = scs[1];

   for( i = 1, i2 = 2; i < nthreads; i++, i2 += 2 )
   {
      if( scale >= scs[i2] )
      {
         if( scale != ATL_rzero )
         { temp = scs[i2] / scale; ssq += ( temp * temp ) * scs[i2+1]; }
      }
      else
      {
         temp = scale / scs[i2];
         ssq  = scs[i2+1] + ( temp * temp ) * ssq; scale = scs[i2];
      }
   }
#if defined( SREAL ) || defined( SCPLX )
   return( scale * (TYPE)(sqrt( (double)(ssq) )) );
#else
   return( scale * sqrt( ssq ) );
#endif
#endif
/*
 * End of Mjoin( PATL, ptnrm2 )
 */
}
