#include <iostream>
using namespace std;

#include <cassert>

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"

#include <cufft.h>

#ifndef REAL
#define REAL float
#endif

cufftHandle fftPlanMany_C2R, fftPlanMany_R2C; cufftHandle fftPlan_C2R, fftPlan_R2C;

template <typename Real>
void inline initFFTs(int *dim, int neall) {
  // Create a batched 3D plan
  if(sizeof(Real) == sizeof(float) ) {
    cufftPlanMany(&fftPlanMany_C2R, 3, dim, NULL, 1, 0, NULL, 1, 0, CUFFT_C2R,neall);
    cufftPlanMany(&fftPlanMany_R2C, 3, dim, NULL, 1, 0, NULL, 1, 0, CUFFT_R2C,neall);
    cufftPlan3d(&fftPlan_C2R, dim[0], dim[1], dim[2], CUFFT_C2R);
    cufftPlan3d(&fftPlan_R2C, dim[0], dim[1], dim[2], CUFFT_R2C);
  } else {
    cufftPlanMany(&fftPlanMany_C2R, 3, dim, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D,neall);
    cufftPlanMany(&fftPlanMany_R2C, 3, dim, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z,neall);
    cufftPlan3d(&fftPlan_C2R, dim[0], dim[1], dim[2], CUFFT_Z2D);
    cufftPlan3d(&fftPlan_R2C, dim[0], dim[1], dim[2], CUFFT_D2Z);
  }

#ifdef FFT_ALL
  cufftSetCompatibilityMode(fftPlanMany_C2R, CUFFT_COMPATIBILITY_FFTW_ALL);
  cufftSetCompatibilityMode(fftPlanMany_R2C, CUFFT_COMPATIBILITY_FFTW_ALL);
  cufftSetCompatibilityMode(fftPlan_C2R, CUFFT_COMPATIBILITY_FFTW_ALL);
  cufftSetCompatibilityMode(fftPlan_R2C, CUFFT_COMPATIBILITY_FFTW_ALL); 
#else
  //cufftSetCompatibilityMode(fftPlanMany_C2R, CUFFT_COMPATIBILITY_NATIVE);
  //cufftSetCompatibilityMode(fftPlanMany_R2C, CUFFT_COMPATIBILITY_NATIVE);
  //cufftSetCompatibilityMode(fftPlan_C2R, CUFFT_COMPATIBILITY_NATIVE);
  //cufftSetCompatibilityMode(fftPlan_R2C, CUFFT_COMPATIBILITY_NATIVE); 
#endif
}

inline void _FFTerror(int ret) {
  switch(ret) {
  case CUFFT_SETUP_FAILED: cerr << "SETUP_FAILED" << endl; break;
  case CUFFT_INVALID_PLAN: cerr << "INVALID_PLAN" << endl; break;
  case CUFFT_INVALID_VALUE: cerr << "INVALID_VALUE" << endl; break;
  case CUFFT_EXEC_FAILED: cerr << "EXEC_FAILED" << endl; break;
  default: cerr << "UNKNOWN ret code " << ret << endl;
  }
}

//template specialization to handle different data types (float,double) 
inline void crFFT_(cufftHandle myFFTplan, float* A, float* B ) {
  int ret=cufftExecC2R(myFFTplan, (cufftComplex*)A,(cufftReal*) B);
  
  if(ret != CUFFT_SUCCESS) {
    cerr << "C2R FFT failed! ret code " << ret << endl; _FFTerror(ret); exit(1);
  }
}

inline void crFFT_(cufftHandle myFFTplan, double* A, double* B ) {
  int ret = cufftExecZ2D(myFFTplan, (cufftDoubleComplex*)A,(cufftDoubleReal*) B);
  
  if(ret != CUFFT_SUCCESS) {
    cerr << "Z2D FFT failed! ret code " << ret << endl; _FFTerror(ret); exit(1);
  }
}

inline void rcFFT_(cufftHandle myFFTplan, float* A, float* B ) {
  int ret = cufftExecR2C(myFFTplan, (cufftReal*)A,(cufftComplex*) B);

  if(ret != CUFFT_SUCCESS) {
    cerr << "C R2C FFT failed!" << endl; _FFTerror(ret); exit(1);
  }
}

inline void rcFFT_(cufftHandle myFFTplan, double* A, double* B ) {
  int ret = cufftExecD2Z(myFFTplan, (cufftDoubleReal*)A,(cufftDoubleComplex*) B);

  if(ret != CUFFT_SUCCESS) {
    cerr << "D2Z FFT failed!" << endl; _FFTerror(ret); exit(1);
  }
}

template <typename Real> inline void rcFFT_neall( Real* A, Real* B ) { rcFFT_(fftPlanMany_R2C, A,B); }
template <typename Real> inline void crFFT_neall( Real* A, Real* B ) { crFFT_(fftPlanMany_C2R, A,B); }
template <typename Real> inline void rcFFT_single( Real* A, Real* B ) { rcFFT_(fftPlan_R2C, A,B); } 
template <typename Real> inline void crFFT_single( Real* A, Real* B ) { crFFT_(fftPlan_C2R, A,B); }

template <typename Real> static void test(int *dim, int neall) {
  int nfft = dim[0]*dim[1]*dim[2];
  int n2ft3d = (dim[0]+2)*dim[1]*dim[2];

  // fill the test data
  thrust::host_vector<Real> h_testData(n2ft3d*neall);
  thrust::host_vector<Real> h_testData1(n2ft3d*neall);

  for(int fft=0; fft < neall; fft++)
    for(int k=0; k < n2ft3d; k++)
      h_testData[k+fft*n2ft3d] = k;

  for(int fft=0; fft < neall; fft++)
    for(int k=0; k < dim[2]; k++)
      for(int j=0; j < dim[1]; j++) {
	h_testData[dim[0]    + j*(dim[0]+2) + k*(dim[0]+2)*dim[1]+fft*n2ft3d] = 0;
	h_testData[dim[0]+ 1 + j*(dim[0]+2) + k*(dim[0]+2)*dim[1]+fft*n2ft3d] = 0;
      }
  
  thrust::device_vector<Real> d_testData = h_testData;
  rcFFT_neall(thrust::raw_pointer_cast(&d_testData[0]),
	      thrust::raw_pointer_cast(&d_testData[0]));
  crFFT_neall(thrust::raw_pointer_cast(&d_testData[0]),
	      thrust::raw_pointer_cast(&d_testData[0]));
  h_testData1 = d_testData;

  for(int fft=0; fft < neall; fft++)
    for (int i=0; i<n2ft3d; ++i) {
      h_testData1[i+fft*n2ft3d] /= (Real)nfft;
      assert(h_testData1[i+fft*n2ft3d] == h_testData1[i+fft*n2ft3d]);
      //cout << "test data " << h_testData1[i+fft*n2ft3d] << " should be " 
	   //<< h_testData[i+fft*n2ft3d] << endl;
    }
}

main(int argc, char *argv[])
{
  if(argc < 5) {
    cerr << "Use dim[0] dim[1] dim[2] numberFFT" << endl;
    exit(1);
  }
  
  int dim[] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
  int neall=atoi(argv[4]);

  cerr << "dim[0] = " << dim[0] << endl;
  cerr << "dim[1] = " << dim[1] << endl;
  cerr << "dim[2] = " << dim[2] << endl;
  cerr << "neall = " << neall << endl;
  cerr << "sizeof(REAL) is " << sizeof(REAL) << " bytes" << endl;
  
  initFFTs<REAL>(dim,neall);
  test<REAL>(dim, neall);
}

