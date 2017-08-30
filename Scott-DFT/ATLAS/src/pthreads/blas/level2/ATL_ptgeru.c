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

PT_FUN_ARG_T Mjoin( PATL, ptgeru0 )( PT_FUN_ARG_T ARGS )
{
/*
 * .. Local Variables ..
 */
   PT_TREE_T                  root = (PT_TREE_T)(ARGS);
   PT_GER_T                   * arg = ARGS;
/* ..
 * .. Executable Statements ..
 *
 */
   Mjoin( PATL, geru )( arg->m, arg->n, (SCALAR)(arg->al), (TYPE *)(arg->x),
      arg->incx, (TYPE *)(arg->y), arg->incy, (TYPE *)(arg->a), arg->la );

   return( NULL );
/*
 * End of @(rname)
 */
}
void Mjoin( PATL, ptgeru )
(
   const int                  M,
   const int                  N,
   const SCALAR               ALPHA,
   const TYPE                 * X,
   const int                  INCX,
   const TYPE                 * Y,
   const int                  INCY,
   TYPE                       * A,
   const int                  LDA
)
{
/*
 * Purpose
 * =======
 *
 * Mjoin( PATL, ptgeru ) performs the rank 1 operation
 *
 *    A := alpha * x * y' + A,
 *
 * where alpha is a scalar,  x is an m-element vector, y is an n-element
 * vector and A is an m by n matrix.
 *
 * For a  more  detailed description of  the arguments of this function,
 * see the reference implementation in the  ATLAS/src/blas/reference di-
 * rectory.
 *
 * ---------------------------------------------------------------------
 */
   pthread_t pid[ATL_NTHREADS];
   PT_GER_T arg[ATL_NTHREADS];
   pthread_attr_t attr;
   int i, nb, nn, nr, J, nthreads;
/* ..
 * .. Executable Statements ..
 *
 */
   if( ( M == 0 ) || ( N == 0 ) || ( SCALAR_IS_ZERO( ALPHA ) ) ) return;

   nb = N / ATL_NTHREADS;
   nthreads = nb ? ATL_NTHREADS : 1;
   if (nthreads < 2)
   {
      Mjoin( PATL, geru )( M, N, ALPHA, X, INCX, Y, INCY, A, LDA );
      return;
   }
   nr = N - nb*nthreads;
   ATL_thread_init(&attr);
   for (J=i=0; i < nthreads; i++)
   {
      nn = (i < nr) ? nb+1 : nb;
      arg[i].al = SADD ALPHA;
      arg[i].x = X;
      arg[i].incx = INCX;
      arg[i].y = Y + INCY*J;
      arg[i].incy = INCY;
      arg[i].a = A + LDA*J;
      arg[i].la = LDA;
      arg[i].m = M;
      arg[i].n = nn;
      if (i != nthreads-1)
         ATL_assert(!pthread_create(&pid[i], &attr, Mjoin( PATL, ptgeru0 ), arg+i));
      J += nn SHIFT;
   }
   Mjoin( PATL, ptgeru0 )(arg+nthreads-1);
   for (i=0; i < nthreads-1; i++)
      ATL_assert(!pthread_join(pid[i], NULL));
   ATL_thread_free(&attr);

/*
   ATL_thread_init( &attr );
   root = ATL_Sger( &type, 0, nthreads, nb, M, N, (void *)(ALPHA), (void *)(X),
                    INCX, (void *)(Y), INCY, (void *)(A), LDA );
   ATL_thread_tree( root, &attr );
   ATL_join_tree  ( root  );
   ATL_free_tree  ( root  );
   ATL_thread_free( &attr );
*/
/*
 * End of Mjoin( PATL, ptgeru )
 */
}
