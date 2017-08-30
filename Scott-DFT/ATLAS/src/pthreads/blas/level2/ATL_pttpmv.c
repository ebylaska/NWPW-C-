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

PT_FUN_ARG_T Mjoin( PATL, pttpmv0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_TPMV_T                  * arg;
/* ..
 * .. Executable Statements ..
 *
 */
   ATL_wait_tree( root );
   arg = (PT_TPMV_T *)(root->arg);
   Mjoin( PATL, tpmv )( arg->up, arg->ta, arg->di, arg->n, (TYPE *)(arg->a),
                        (TYPE *)(arg->x), arg->incx );
   ATL_signal_tree( root );

   return( NULL );
/*
 * End of Mjoin( PATL, pttpmv0 )
 */
}

void Mjoin( PATL, pttpmv )
(
   const enum ATLAS_UPLO      UPLO,
   const enum ATLAS_TRANS     TRANS,
   const enum ATLAS_DIAG      DIAG,
   const int                  N,
   const TYPE                 * A,
   TYPE                       * X,
   const int                  INCX
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, pttpmv ) performs one of the matrix-vector operations
 *
 *    x := A * x,   or   x := conjg( A  ) * x,   or
 *
 *    x := A'* x,   or   x := conjg( A' ) * x,
 *
 * where x is an n-element vector and  A is an n by n unit, or non-unit,
 * upper or lower triangular matrix, supplied in packed form.
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
   if( N == 0 ) return;
   Mjoin(PATL,tpmv)(UPLO, TRANS, DIAG, N, A, X, INCX);
/*
   ATL_thread_init( &attr );
   root = ATL_Stpmv( &type, 0, nthreads, nb, UPLO, TRANS, DIAG, N, (void *)(A),
                     LDA, (void *)(X), INCX );
   ATL_thread_tree( root, &attr );
   ATL_join_tree  ( root  );
   ATL_free_tree  ( root  );
   ATL_thread_free( &attr );
*/
/*
 * End of Mjoin( PATL, pttpmv )
 */
}
