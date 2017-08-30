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

PT_FUN_ARG_T Mjoin( PATL, ptsyr0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_SYR_T                   * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_SYR_T *)(root->arg);
#ifdef TREAL
   Mjoin( PATL, syr )( arg->up, arg->n, (SCALAR)(*((TYPE *)(arg->al))),
      (TYPE *)(arg->x), arg->incx, (TYPE *)(arg->a), arg->la );
#else
   Mjoin( PATL, syr )( arg->up, arg->n, (SCALAR)(arg->al), (TYPE *)(arg->x),
      arg->incx, (TYPE *)(arg->a), arg->la );
#endif
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, ptsyr0 )
 */
}

void Mjoin( PATL, ptsyr )
(
   const enum ATLAS_UPLO      UPLO,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * X,
   const int                  INCX,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptsyr ) performs the symmetric rank 1 operation
 *
 *    A := alpha * x * x' + A,
 *
 * where  alpha is a scalar, x is an n-element vector and A is an n by n
 * symmetric matrix.
 *
 * For a  more  detailed description of  the arguments of this function,
 * see the reference implementation in the  ATLAS/src/blas/reference di-
 * rectory.
 *
 * ---------------------------------------------------------------------
 */
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( N <= 0 ) || ( SCALAR_IS_ZERO( ALPHA ) ) ) return;

   Mjoin(PATL,syr)(UPLO, N, ALPHA, X, INCX, A, LDA);
/*
   ATL_thread_init( &attr );
   root = ATL_Ssyr( &type, 0, nthreads, nb, UPLO, N, (void *)(ALPHA),
                    (void *)(X), INCX, (void *)(A), LDA );
   ATL_thread_tree( root, &attr );
   ATL_join_tree  ( root  );
   ATL_free_tree  ( root  );
   ATL_thread_free( &attr );
*/
/*
 * End of Mjoin( PATL, ptsyr )
 */
}
