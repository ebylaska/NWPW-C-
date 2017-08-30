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
#include "atlas_kernel2.h"
#include "atlas_reflvl2.h"
#include "atlas_lvl2.h"
#include "atlas_ptlevel1.h"
#include "atlas_ptlvl2.h"
#include "atlas_ptlevel2.h"

PT_FUN_ARG_T Mjoin( PATL, ptgemv0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GEMV_T                  * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_GEMV_T *)(root->arg);
#ifdef TREAL
   Mjoin( PATL, gemv )( arg->ta, arg->m, arg->n, (SCALAR)(*((TYPE *)(arg->al))),
      (TYPE *)(arg->a), arg->la, (TYPE *)(arg->x), arg->incx,
      (SCALAR)(*((TYPE *)(arg->be))), (TYPE *)(arg->y), arg->incy );
#else
   Mjoin( PATL, gemv )( arg->ta, arg->m, arg->n, (SCALAR)(arg->al),
      (TYPE *)(arg->a), arg->la, (TYPE *)(arg->x), arg->incx,
      (SCALAR)(arg->be), (TYPE *)(arg->y), arg->incy );
#endif
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptgemv0 )
 */
}

void Mjoin( PATL, ptgemv )
(
   const enum ATLAS_TRANS     TRANS,
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * A,
   const int                  LDA,
   const TYPE                 * X,
   const int                  INCX,
   const SCALAR               BETA,
   TYPE                       * Y,
   const int                  INCY
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgemv ) performs one of the matrix-vector operations
 *
 *    y := alpha * op( A ) * x + beta * y,
 *
 * where op( X ) is one of
 *
 *    op( X ) = X   or   op( X ) = conjg( X  )   or
 *
 *    op( X ) = X'  or   op( X ) = conjg( X' ).
 *
 * where  alpha and beta are scalars, x and y are vectors and op( A ) is
 * an m by n matrix.
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
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) ||
       ( ( SCALAR_IS_ZERO( ALPHA ) ) && ( SCALAR_IS_ONE( BETA ) ) ) ) return;

   if( SCALAR_IS_ZERO( ALPHA ) )
   {
      if( !( SCALAR_IS_ONE( BETA ) ) )
         Mjoin( PATL, ptscal )( M, BETA, Y, INCY );
      return;
   }

   Mjoin( PATL, gemv )( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY );
/*
   ATL_thread_init( &attr );
   root = ATL_Sgemv( &type, 0, nthreads, nb, TRANS, M, N, (void *)(alpha),
                     (void *)(A), LDA, (void *)(X), INCX, (void *)(beta),
                     (void *)(Y), INCY );
   ATL_thread_tree( root, &attr );
   ATL_join_tree  ( root  );
   ATL_free_tree  ( root  );
   ATL_thread_free( &attr );
*/
/*
 * End of Mjoin( PATL, ptgemv )
 */
}
