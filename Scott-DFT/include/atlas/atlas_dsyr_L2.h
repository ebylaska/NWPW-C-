/*
 * This file generated on line 467 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_DSYR_L2_H
   #define ATLAS_DSYR_L2_H

#include "atlas_dr1_L2.h"

#define ATL_s1NU 1

#define ATL_SYR1U_nu(A_, lda_, x_, y_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0_=*(x_); \
   const TYPE y0_=*(y_); \
   aa[lda0_+0] += x0_*y0_; \
}
#define ATL_SYR1L_nu(A_, lda_, x_, y_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0; \
   const TYPE x0_=*(x_); \
   const TYPE y0_=*(y_); \
   aa[lda0_+0] += x0_*y0_; \
}

#endif
