/*
 * This file generated on line 811 of /home/bylaska/Codes/Scott-DFT/ATLAS/Dev1/..//tune/blas/ger/r1hgen.c
 */
#ifndef ATLAS_SR1KERNELS_H
   #define ATLAS_SR1KERNELS_H

void ATL_sgerk_L0(ATL_CINT, ATL_CINT, const float, const float*, ATL_CINT, const float*, ATL_CINT, float*, ATL_CINT);
void ATL_sgerk_L1(ATL_CINT, ATL_CINT, const float, const float*, ATL_CINT, const float*, ATL_CINT, float*, ATL_CINT);

#define ATL_sgerk_L0_restrict    ATL_sgerk_L0
#define ATL_sgerk_L2             ATL_sgerk_L0
#define ATL_sgerk_L2_restrict    ATL_sgerk_L0
#define ATL_sgerk_L1_restrict    ATL_sgerk_L1
#define ATL_sgerk_L1b            ATL_sgerk_L0
#define ATL_sgerk_L1b_restrict   ATL_sgerk_L0

#endif /* end guard around atlas_sr1kernels.h */
