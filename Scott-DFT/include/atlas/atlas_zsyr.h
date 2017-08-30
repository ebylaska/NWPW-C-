/*
 * This file generated on line 467 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_ZSYR_H
   #define ATLAS_ZSYR_H

#include "atlas_zr1.h"

#define ATL_s1NU 1

#define ATL_NOBLOCK_S1 1
#define ATL_GetPartS1 ATL_GetPartR1

#define ATL_HER1U_nu(A_, lda_, x_, xt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0r=*(x_), x0i=(x_)[1]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]; \
   aa[lda0_+0] += x0r*xt0r-x0i*xt0i; \
   aa[lda0_+1] = 0.0; \
}
#define ATL_HER1L_nu(A_, lda_, x_, xt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0r=*(x_), x0i=(x_)[1]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1]; \
   aa[lda0_+0] += x0r*xt0r-x0i*xt0i; \
   aa[lda0_+1] = 0.0; \
}

#endif
