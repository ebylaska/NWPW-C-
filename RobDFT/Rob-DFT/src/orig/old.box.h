#ifndef BOX_H
#define BOX_H
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
using namespace std;

#include <stdlib.h> // needed for exit

#include "psiHeader.h"
#include "CommonData.h"

// Thrust headers
#include <thrust/transform_reduce.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <cmath>

#define UNDEFINED_SPACE 0
#define R_SPACE 1
#define K_SPACE 2
#define M_SPACE 3

template <typename Real>
class box {
 protected:
  int spaceType;
  PsiHeader<Real>* pH;
  CommonData<Real>* cd;
  thrust::host_vector<Real> h_psi[2];
  thrust::device_vector<Real> d_psi[2];
  thrust::host_vector<Real> *psi, *psi1, *psi2;

  void errorExit(char* s) { cerr << s << endl; exit(1);}

 public:
  box()
    {
      spaceType=UNDEFINED_SPACE;
      memset(this,0, sizeof(box<Real>) ); //fill the class regardless of size
    }
  ~box()
    {
    }
  char* read(char* buf, PsiHeader<Real>* psiHeader, 
	     CommonData<Real>* myCommonData)
  {
    pH = psiHeader;
    if(!cd) cd = myCommonData;
    int size = (cd->nx+2) * cd->ny * cd->nz;
    
    d_psi[0] = thrust::device_vector<Real>(size);
    d_psi[1] = thrust::device_vector<Real>(size);

    h_psi[0] = thrust::host_vector<Real>(size);
    h_psi[1] = thrust::host_vector<Real>(size);

    // temp to check code
    psi = psi1 = &h_psi[0];
    psi2 = &h_psi[1];
    
    for (int k=0; k < 2*(cd->nfft3d); ++k)
      h_psi[0][k] = (Real) ((double*)buf)[k];
    d_psi[0] = h_psi[0]; // move to GPU
    
    switchPsi();
    return(buf+sizeof(double)*size);
  }
  inline Real* getPsiPtr() {return(&(*psi)[0]); }
  char* write()
  {
    int nBytes=0;
    return(nBytes);
  }
  inline void switchPsi()
    {
      thrust::host_vector<Real> *tmp = psi1;
      psi1 = psi2;
      psi2 = tmp;
    }

#ifdef HOST_SIDE
  inline void kspace_operator()
  {
    // tg * psi
    for(int j=0; j < 2*cd->nfft3d; ++j) {
      register Real psiJ = (*psi)[j];
      register Real Hpsi = psiJ  * cd->h_2tg[j];
      (*psi2)[j] = psiJ  + cd->dte*Hpsi;
    }
  }
#else
  struct kspace_op_param {
    const Real* psi;
    Real* psi1;
    const Real* tg;
    const Real dte;
  };

  struct kspace_op
  {
    const kspace_op_param k_param;
  kspace_op(kspace_op_param _param): k_param(_param) {};
    __host__ __device__
    Real operator()(const int& index) const { 
      register Real psi = k_param.psi[index];
      register Real Hpsi = psi  * k_param.tg[index]; // Hpsi = psiJ * tg[j];
      k_param.psi1[index] =  psi + k_param.dte * Hpsi;
      //cerr << index << " " << psi << " " << k_param.psi[index] << " " << k_param.dte << endl;
      return (Real) index;
    }
  };
  struct kspace_operator_functor
  {
    const Real dte;
  kspace_operator_functor(Real _dte): dte(_dte) {};
    template <typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
      // psi1: get<0>
      // tg  : get<1>
      // psi2: get<2>
      register Real psiJ = thrust::get<0>(t);
      register Real Hpsi = psiJ  * thrust::get<1>(t); // Hpsi = psiJ * tg[j];
      thrust::get<2>(t) =  psiJ + dte * Hpsi;
    }
  };

  inline void kspace_operator()
  {
    int n = 2*cd->nfft3d;
    thrust::plus<Real> binary_op;
    kspace_op_param k_param = {&(*psi)[0], &(*psi1)[0], &(cd->h_2tg[0]), cd->dte};
    kspace_op unary_op(k_param);
    thrust::plus<Real> kspace_bin_op;
#define TRANSFORM_REDUCE
#ifdef TRANSFORM_REDUCE
    Real sum = thrust::transform_reduce(psi->begin(), psi->begin()+n, unary_op, 0., binary_op);
#else
    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(psi1->begin(),
								  cd->h_2tg.begin(),
								  psi2->begin())),
		     thrust::make_zip_iterator(thrust::make_tuple(psi1->begin()+n,
								  cd->h_2tg.end(),
								  psi2->begin()+n) ),
                     kspace_operator_functor(cd->dte));
#endif
  }
#endif
  inline void fftExecC2R_3d()
    {
    }
  inline void fftExecR2C_3d()
    {
    }
};

#endif
