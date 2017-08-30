/*
 * This file generated on line 603 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_DSYR2_H
   #define ATLAS_DSYR2_H

#include "atlas_dr1kernels.h"
#define ATL_s2CacheElts 2048
#define ATL_s2MU 8
#define ATL_s2NU 4
#define ATL_R1OC ATL_dgerk_L1b
#define ATL_R1OCr ATL_dgerk_L1b
#define ATL_R1IC ATL_dgerk_L1
#define ATL_R1ICr ATL_dgerk_L1

#define ATL_GetPartS2(A_, lda_, mb_, nb_) { (mb_) = 200; (nb_) = ATL_s2NU; }

#define ATL_SYR2U_nu(A_, lda_, x_, y_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0, lda1_ = lda0_+(lda_), lda2_ = lda1_+(lda_), lda3_ = lda2_+(lda_); \
   const TYPE x0_=*(x_), x1_=(x_)[1], x2_=(x_)[2], x3_=(x_)[3]; \
   const TYPE y0_=*(y_), y1_=(y_)[1], y2_=(y_)[2], y3_=(y_)[3]; \
   aa[lda0_+0] += x0_*y0_ + y0_*x0_; \
   aa[lda1_+0] += x0_*y1_ + y0_*x1_; \
   aa[lda1_+1] += x1_*y1_ + y1_*x1_; \
   aa[lda2_+0] += x0_*y2_ + y0_*x2_; \
   aa[lda2_+1] += x1_*y2_ + y1_*x2_; \
   aa[lda2_+2] += x2_*y2_ + y2_*x2_; \
   aa[lda3_+0] += x0_*y3_ + y0_*x3_; \
   aa[lda3_+1] += x1_*y3_ + y1_*x3_; \
   aa[lda3_+2] += x2_*y3_ + y2_*x3_; \
   aa[lda3_+3] += x3_*y3_ + y3_*x3_; \
}
#define ATL_SYR2L_nu(A_, lda_, x_, y_) \
{ \
   TYPE *aa=(A_); \
   ATL_CINT lda0_ = 0, lda1_ = lda0_+(lda_), lda2_ = lda1_+(lda_), lda3_ = lda2_+(lda_); \
   const TYPE x0_=*(x_), x1_=(x_)[1], x2_=(x_)[2], x3_=(x_)[3]; \
   const TYPE y0_=*(y_), y1_=(y_)[1], y2_=(y_)[2], y3_=(y_)[3]; \
   aa[lda0_+0] += x0_*y0_ + y0_*x0_; \
   aa[lda0_+1] += x1_*y0_ + y1_*x0_; \
   aa[lda0_+2] += x2_*y0_ + y2_*x0_; \
   aa[lda0_+3] += x3_*y0_ + y3_*x0_; \
   aa[lda1_+1] += x1_*y1_ + y1_*x1_; \
   aa[lda1_+2] += x2_*y1_ + y2_*x1_; \
   aa[lda1_+3] += x3_*y1_ + y3_*x1_; \
   aa[lda2_+2] += x2_*y2_ + y2_*x2_; \
   aa[lda2_+3] += x3_*y2_ + y3_*x2_; \
   aa[lda3_+3] += x3_*y3_ + y3_*x3_; \
}

#endif
