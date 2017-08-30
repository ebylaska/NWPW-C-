/*
 *             Automatically Tuned Linear Algebra Software v3.9.23
 * Copyright (C) 2009 Siju Samuel
 *
 * Code contributers : Siju Samuel, Anthony M. Castaldo, R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*-----------------------------------------------------------------------------
* This is the C translation of the standard LAPACK Fortran routine:
*      SUBROUTINE C/Z LACGV( N, X, INCX )
*
*  Purpose
*  =======
*
*  ATL_lacgv.c conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) complex array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*          NOTE : complex numbers are stored as,
*          real(single/complex), imaginary(single/complex)
*          in concequtive memory locations.
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
-----------------------------------------------------------------------------*/
#include "atlas_misc.h"
#include "cblas.h"
#include "atlas_lapack.h"

//Compiled only to precisions single complex and double complex.
void  ATL_lacgv( const int N, TYPE *X, int INCX)
{
   int i, ioff;

   if ( INCX == 1)
   {
      //imaginary part = Negetive of imaginary(X[i])
      for ( i = 0; i < N; i++)
      {
          *(X + (i SHIFT) + 1 ) = 0.0 - *(X + (i SHIFT) + 1 ) ;
      }
   }
   else
   {
      //TODO : Test this section on i off
      ioff =0;
      if ( INCX < 0 ){
         ioff = 0 - ( N-1)*INCX;
      }

      for ( i =0; i < N; i=i++)
      {
          //Negetaive of imaginary number is taken
          *(X+ ( ioff SHIFT ) + 1  ) = 0.0 - *(X+ (ioff SHIFT) + 1 );
          ioff = ioff + INCX;                    //increment the offset
      }
   } //else on INCX

   return;
}
