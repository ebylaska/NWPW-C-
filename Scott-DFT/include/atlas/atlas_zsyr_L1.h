/*
 * This file generated on line 467 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_ZSYR_L1_H
   #define ATLAS_ZSYR_L1_H

#include "atlas_zr1_L1.h"

#define ATL_s1NU 4

#define ATL_HER1U_nu(A_, lda_, x_, xt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0, lda1_ = lda0_+(lda_)+(lda_), lda2_ = lda1_+(lda_)+(lda_), lda3_ = lda2_+(lda_)+(lda_); \
   const TYPE x0r=*(x_), x0i=(x_)[1], x1r=(x_)[2], x1i=(x_)[3], x2r=(x_)[4], x2i=(x_)[5], x3r=(x_)[6], x3i=(x_)[7]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1], xt1r=(xt_)[2], xt1i=(xt_)[3], xt2r=(xt_)[4], xt2i=(xt_)[5], xt3r=(xt_)[6], xt3i=(xt_)[7]; \
   aa[lda0_+0] += x0r*xt0r-x0i*xt0i; \
   aa[lda0_+1] = 0.0; \
   aa[lda1_+0] += x0r*xt1r-x0i*xt1i; \
   aa[lda1_+1] += x0r*xt1i+x0i*xt1r; \
   aa[lda1_+2] += x1r*xt1r-x1i*xt1i; \
   aa[lda1_+3] = 0.0; \
   aa[lda2_+0] += x0r*xt2r-x0i*xt2i; \
   aa[lda2_+1] += x0r*xt2i+x0i*xt2r; \
   aa[lda2_+2] += x1r*xt2r-x1i*xt2i; \
   aa[lda2_+3] += x1r*xt2i+x1i*xt2r; \
   aa[lda2_+4] += x2r*xt2r-x2i*xt2i; \
   aa[lda2_+5] = 0.0; \
   aa[lda3_+0] += x0r*xt3r-x0i*xt3i; \
   aa[lda3_+1] += x0r*xt3i+x0i*xt3r; \
   aa[lda3_+2] += x1r*xt3r-x1i*xt3i; \
   aa[lda3_+3] += x1r*xt3i+x1i*xt3r; \
   aa[lda3_+4] += x2r*xt3r-x2i*xt3i; \
   aa[lda3_+5] += x2r*xt3i+x2i*xt3r; \
   aa[lda3_+6] += x3r*xt3r-x3i*xt3i; \
   aa[lda3_+7] = 0.0; \
}
#define ATL_HER1L_nu(A_, lda_, x_, xt_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0, lda1_ = lda0_+(lda_)+(lda_), lda2_ = lda1_+(lda_)+(lda_), lda3_ = lda2_+(lda_)+(lda_); \
   const TYPE x0r=*(x_), x0i=(x_)[1], x1r=(x_)[2], x1i=(x_)[3], x2r=(x_)[4], x2i=(x_)[5], x3r=(x_)[6], x3i=(x_)[7]; \
   const TYPE xt0r=*(xt_), xt0i=(xt_)[1], xt1r=(xt_)[2], xt1i=(xt_)[3], xt2r=(xt_)[4], xt2i=(xt_)[5], xt3r=(xt_)[6], xt3i=(xt_)[7]; \
   aa[lda0_+0] += x0r*xt0r-x0i*xt0i; \
   aa[lda0_+1] = 0.0; \
   aa[lda0_+2] += x1r*xt0r-x1i*xt0i; \
   aa[lda0_+3] += x1r*xt0i+x1i*xt0r; \
   aa[lda0_+4] += x2r*xt0r-x2i*xt0i; \
   aa[lda0_+5] += x2r*xt0i+x2i*xt0r; \
   aa[lda0_+6] += x3r*xt0r-x3i*xt0i; \
   aa[lda0_+7] += x3r*xt0i+x3i*xt0r; \
   aa[lda1_+2] += x1r*xt1r-x1i*xt1i; \
   aa[lda1_+3] = 0.0; \
   aa[lda1_+4] += x2r*xt1r-x2i*xt1i; \
   aa[lda1_+5] += x2r*xt1i+x2i*xt1r; \
   aa[lda1_+6] += x3r*xt1r-x3i*xt1i; \
   aa[lda1_+7] += x3r*xt1i+x3i*xt1r; \
   aa[lda2_+4] += x2r*xt2r-x2i*xt2i; \
   aa[lda2_+5] = 0.0; \
   aa[lda2_+6] += x3r*xt2r-x3i*xt2i; \
   aa[lda2_+7] += x3r*xt2i+x3i*xt2r; \
   aa[lda3_+6] += x3r*xt3r-x3i*xt3i; \
   aa[lda3_+7] = 0.0; \
}

#endif
