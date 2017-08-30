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
 *     SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
 *
 * ATL_larfp.c :
 * void ATL_larfp( const int N, TYPE *ALPHA, TYPE *X, int INCX, TYPE *TAU)
 *        NOTE : ATL_larfp.c  will get compiled to 4 precisions
 *               single precision real,      double precision real
 *               single precision complex,   double precision complex
 *  Purpose
 *  =======
 *
 *  Real Precision
 *  --------------
 *
 *  ATL_larfp generates a real/complex  elementary reflector H of order n, such
 *  that
 *
 *        H * ( alpha ) = ( beta ),   H' * H = I.
 *            (   x   )   (   0  )
 *
 *  where alpha and beta are scalars, and x is an (n-1)-element real
 *  vector. H is represented in the form
 *
 *        H = I - tau * ( 1 ) * ( 1 v' ) ,                     (Real Precisions)
 *                      ( v )
 *
 *        H = I - tau * ( 1 ) * ( 1 conjugate(v)' ) ,       (Complex Precisions)
 *                      ( v )
 *
 *  where tau is a real/complex scalar and v is a real/complex (n-1)-element
 *  vector.
 *
 *  If the elements of x are all zero, then tau = 0 and H is taken to be
 *  the unit matrix.
 *
 *  Otherwise  1 <= tau <= 2.
 *
 *
 *  Arguments
 *  =========
 *
 *  N       (input) INTEGER
 *          The order of the elementary reflector.
 *
 *  ALPHA   (input/output)
 *          On entry, the value alpha.
 *          On exit, it is overwritten with the value beta.
 *
 *  X       (input/output)   array pointer, dimension
 *                         (1+(N-2)*abs(INCX))
 *          On entry, the vector x.
 *          On exit, it is overwritten with the vector v.
 *
 *  INCX    (input) INTEGER
 *          The increment between elements of X. INCX > 0.
 *
 *  TAU     (output)
 *          value of tau.
 *-----------------------------------------------------------------------------*/
#include "atlas_misc.h"
#include <math.h>
#include "cblas.h"
#include "atlas_lapack.h"
#include "atlas_lamch.h"

/*---------------------------------------------------------------------------*/
/* HighLevel Logic : * Real Precision                                        */
/* --------------------------------------------------------------------------*/
/* On entry from ATL_geqr2(or any other variant lq, ql, rq) ,                */
/* *Alpha points at A[i,i] and *X points  at A[i+1, i].                      */
/* If N==1, cblas_nrm2 returns zero. The norm is actually                    */
/* found in two parts; XNORM is all of the column except for A[i,i], and if  */
/* that is zero, we return TAU of zero (so H = I).                           */
/* Otherwise, we combine XNORM and A[i,i] into BETAp using lapy2. So BETAp   */
/* is the actual norm2 of A[i:m, i].                                         */
/* We set BETA to BETAp but with the opposite sign as A[i,i]. This is done   */
/* to ensure that TAU is in [1,2].                                           */
/* We set TAU to (BETA-A[i,i])/BETA, and scale A[i+1:m,i] by 1/(A[i,i]-BETA).*/
/* Finally, we replace A[i,i] with BETA.                                     */
/* Note that |A[i,i]-BETA| > |A[i,i]|. Treat the column as X=[x1, x2, x3],   */
/* and assume x1 > 0. Then we return X' (prime, not transpose) as:           */
/* x1' = -||X||  (The element actually on the diagonal, and part of R.)      */
/* x2' = x2/(x1+||X||)                                                       */
/* x3' = x3/(x1+||X||)                                                       */
/* TAU = (-||X||-x1)/-||X|| = (x1+||X||)/||X||                               */
/*                                                                           */
/* This is NOT a textbook Householder reflection, because [1 v]^T does not   */
/* have a norm of 1. The purpose seems to be to save the first storage       */
/* location; the '1' is not stored anywhere. The choice of whether to make   */
/* A[i,i] either +/- ||X|| is strictly to control the range of TAU, but has  */
/* the added value of making the Householder matrix unique for any given A.  */
/* Note that |x1| < | (x1+||X||) |, i.e. the magnitude is always increased,  */
/* so the magnitude of xi/(x1+||X||) is always decreased. I'm not sure if    */
/* that is intended to reduce error or not.                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/* Complex Precisions                                                        */
/*---------------------------------------------------------------------------*/
/* The  highlevel logic for Complex precision remains same as that of Real   */
/* precisions. The difference are                                            */
/*    BETAp is calculated using real part of A[i,i], imaginary part of A[i,i]*/
/*    and XNORM using lapy3.                                                 */
/*                                                                           */
/*    ATL_ladiv is called to apply complex number devision before            */
/*    performing  the scaling operation for A[i+1:m,i] by 1/(A[i,i]-BETA).   */
/*                                                                           */
/*   NOTE :                                                                  */
/*    For Real precision and Complex precision, the codes are kept seperately*/
/*    for clarity. ( Many code might be similar)                             */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void ATL_larfp( const int N, TYPE *ALPHA, TYPE *X, int INCX, TYPE *TAU)
{
#ifdef TREAL

   TYPE ONE=1.0, ZERO=0.0, BETA, BETAp, RSAFMN, SAFMAX, XNORM;
   TYPE TWO=2.0;
   int    j, KNT;

   if (N < 0)
   {
      *TAU = ZERO;
      return;
   }

   // Get the norm2 .
   XNORM = cblas_nrm2(N-1, X, INCX);

   if (XNORM == ZERO)
   {
      /*
       *    H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0
       */
      if ( (*ALPHA) > ZERO )
      {
         /* When TAU.eq.ZERO, the vector is special-cased to be               */
         /* all zeros in the application routines.  We do not need            */
         /* to clear it.                                                      */
         *TAU = ZERO;
      }
      else
      {
         /* However, the application routines rely on explicit                */
         /* zero checks when TAU.ne.ZERO, and we must clear X.                */
         *TAU = TWO;
         for ( j = 0; j < (N-1) ; j++)
         {
            X[ 1 + (j-1)*INCX ] = ZERO;
         }
         *ALPHA = 0. - (*ALPHA);
      }
   }
   else                                  /* XNORM NE ZERO                     */
   {
      BETAp = ATL_lapy2((*ALPHA), XNORM);/* Get sqrt(a^2+b^2)                 */
      BETA = BETAp;                      /* Assume ALPHA < 0                  */
      if ((*ALPHA) < 0) BETA = 0.-BETAp; /* Change the sign if alpha is -ve   */

      /* DLAMCH('S') returns 'safe minimum' (sfmin) such that        */
      /* 1/sfmin does not overflow.                                  */
      /* DLAMCH('E') returns epsilon for single/double  precision.   */



      KNT = 0;
      if (BETAp < ATL_laSAFMIN)
      {
         /*
          *       XNORM, BETA may be inaccurate; scale X and recompute them
          */

         RSAFMN = ONE / ATL_laSAFMIN;          /* Set a maximum                     */

         while (BETAp < ATL_laSAFMIN)
         {
            KNT++;
            cblas_scal(N-1, RSAFMN, X, INCX);
            BETA *= RSAFMN;
            BETAp *= RSAFMN;
            *ALPHA *= RSAFMN;
         }
         /*
          *    New BETA is at most 1, at least SAFMIN
          */
         XNORM = cblas_nrm2(N-1, X, INCX);
         BETAp = ATL_lapy2((*ALPHA), XNORM);  /* Will always be positive      */
         BETA = BETAp;                        /* Assume ALPHA < 0             */
         if ((*ALPHA) > 0) BETA =0.0  -BETAp; /* -SIGN(BETA, ALPHA)           */

      }  /* BATAp < SAFMIN        */

      *ALPHA = *ALPHA + BETA;
      if ( BETA < ZERO )
      {
         BETA = 0. - BETA;
         *TAU = (0. - *ALPHA)/BETA;
      }
      else
      {
         *ALPHA = XNORM * (XNORM/(*ALPHA));
         *TAU   = *ALPHA/BETA;
         *ALPHA = 0.0 - *ALPHA ;
      }
      cblas_scal(N-1, ONE / (*ALPHA), X, INCX);

      for (j=0; j<KNT; j++)
      {
          BETA *= ATL_laSAFMIN;
      }

      *ALPHA = BETA;

   } /* else on XNORM !=0 */
   return;
/*
 * End of  Real Precision ATL_larfp   REAL precision
 */
#else
/*
 * Beginning of  Complex  Precision ATL_larfp
 */
   TYPE ONE=1.0, ZERO=0.0, BETA, BETAp, RSAFMN, SAFMAX, XNORM, ALPHI, ALPHR;
   TYPE ONEVAL[2] =  {ATL_rone, ATL_rzero};
   TYPE TWO[2] =  {2.0, ATL_rzero};
   int j, KNT;

   if ( N < 0)
   {
      /*
       *    H  =  I
       */
      *(TAU)  = 0.0;
      *(TAU + 1) = 0.0;
      return;
   }

   //Get the nrm2
   XNORM = cblas_nrm2(N-1, X, INCX);

   ALPHR = *( ALPHA) ;
   ALPHI = *( ALPHA + 1) ;

   if (XNORM == ZERO )
   {
     /*
      *    H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.
      */
      *(TAU)  = 0.0;
      *(TAU + 1) = 0.0;
      if ( ALPHI == ZERO )
      {
         if(ALPHR > ZERO)
         {
            /* When TAU.eq.ZERO, the vector is special-cased to be
             * all zeros in the application routines.  We do not need
             * to clear it.
             */

            *(TAU)  = 0.0;
            *(TAU + 1) = 0.0;
         }
         else
         {
            /* However, the application routines rely on explicit
             * zero checks when TAU.ne.ZERO, and we must clear X.
             */

            *(TAU)  = 2.0;
            *(TAU + 1) = 0.0;

            for ( j = 0; j < ((N-1) SHIFT ) ; j++) /*complex multiply by 2    */
            {
               X[ 2 + (j-1)*INCX ] = ZERO;
            }
            *(ALPHA) = 0. - *(ALPHA);
            *(ALPHA+1) = 0. - *(ALPHA+1);
         }
      }                                     /* ALPHI != ZERO                  */
      else
      {
         /* Only "reflecting" the diagonal entry to be real                   */
         /*            and non-negative.                                      */

         XNORM = ATL_lapy2( ALPHR, ALPHI ); /* Get sqrt(a^2+b^2)              */

         //TAU = DCMPLX( ONE - ALPHR / XNORM, -ALPHI / XNORM )
         *(TAU) = ( ONE-ALPHR / XNORM) ;
         *(TAU + 1) =  (0.0 - ALPHI) /XNORM ;

         for ( j = 0; j < ((N-1) SHIFT ) ; j++) /*complex multiply by 2       */
         {
            X[ 2 + (j-1)*INCX ] = ZERO;
         }
         *(ALPHA) = XNORM;                   /* Real Part of alpha            */
         *(ALPHA + 1) = 0.0;                 /* Set Imaginary part to Zero    */
      }                                      /* ALPHI = 0 ends                */
   }
   else                                      /* XNORM != 0                    */
   {
      /* XNORM != 0                                                           */
      BETAp = ATL_lapy3(ALPHR, ALPHI, XNORM);/* Get sqrt( a^2 + b^2 + c^2)    */
      BETA = BETAp;                          /* Assume ALPHA < 0              */
      if ( (*ALPHA) < 0) BETA = 0. - BETAp;  /*Change the sign if alpha is -ve*/

      /* DLAMCH('S') returns 'safe minimum' (sfmin) such that        */
      /* 1/sfmin does not overflow.                                  */
      /* DLAMCH('E') returns epsilon for single/double  precision.   */


      RSAFMN = ONE / ATL_laSAFMIN ;
      KNT = 0;

      if ( BETAp  <  ATL_laSAFMIN )
      {
         /*
          *     XNORM, BETA may be inaccurate; scale X and recompute them
          *
          */
         while ( BETAp < ATL_laSAFMIN )
         {
            KNT++;
            #ifdef DCPLX
                cblas_zdscal(N-1, RSAFMN, X, INCX);
            #else
                cblas_csscal(N-1, RSAFMN, X, INCX);
            #endif
            BETA *= RSAFMN;
            BETAp *= RSAFMN;
            ALPHI = ALPHI*RSAFMN;
            ALPHR = ALPHR*RSAFMN;
         }

         /*
          *    New BETA is at most 1, at least SAFMIN
          *
          */
         XNORM = cblas_nrm2(N-1, X, INCX);
         *(ALPHA) = ALPHR;
         *(ALPHA + 1) = ALPHI;

         BETAp = ATL_lapy3(ALPHR, ALPHI, XNORM);  /* Will always be positive   */
         BETA = BETAp;
         if (ALPHR > 0) BETA = 0.0 -BETAp;        /* -SIGN(BETA, ALPHR)        */

      } /* BATAp < SAFMIN        */

      /* alpha = alpha+ beta. Imaginary part remains same                  */
      *ALPHA = *ALPHA + BETA;

      if ( BETA < ZERO )
      {
          BETA = 0.0 - BETA;

         /* TAU = -ALPHA/BETA                                              */
         *TAU     = (0. - *ALPHA)/BETA;
         *(TAU+1) = (0. - *(ALPHA+1))/BETA;
      }
      else
      {
         ALPHR = ALPHI * (ALPHI/ (*ALPHA));
         ALPHR = ALPHR + XNORM * (XNORM/(*ALPHA));

         //TAU = DCMPLX( ALPHR/BETA, -ALPHI/BETA )
         *TAU     = ALPHR/BETA;
         *(TAU+1) = (0.0 - ALPHI)/BETA;

         *ALPHA     = 0.0 - ALPHR;
         *(ALPHA+1) = ALPHI;
      }

      // Perform complex division before scaling the X vector
      ATL_ladiv( ONEVAL,  ALPHA  , ALPHA); /*ALPHA will have  the result   */
      cblas_scal(N-1, ALPHA, X, INCX);

      /* If BETA is subnormal, it may lose relative accuracy               */

      for (j=0; j<KNT; j++)
      {
         BETA *= ATL_laSAFMIN;
      }
      *(ALPHA) = BETA;                    /* Real Part of alpha            */
      *(ALPHA + 1) = 0.0;                 /* Set Imaginary part to Zero    */
   }                                         /* XNORM != 0 ends               */
   return;
#endif
} /* END ATL_larfp */
