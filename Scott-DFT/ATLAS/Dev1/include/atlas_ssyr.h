/*
 * This file generated on line 467 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_SSYR_H
   #define ATLAS_SSYR_H

#include "atlas_sr1.h"

#define ATL_s1NU 1

#define ATL_NOBLOCK_S1 1
#define ATL_GetPartS1 ATL_GetPartR1

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
