#ifndef STRFAC_H
#define STRFAC_H

#include <cmath>
#include "CommonData.h"
#include "Ion.h"

template <typename Real, typename Real2>
class Strfac
{
 protected:
  CommonData<Real> *cd;
  Ion *myion;
  int nxh,nfft3d;
  Real unita[9];
  thrust::device_vector<Real> d_wx1;
  thrust::device_vector<Real> d_wy1;
  thrust::device_vector<Real> d_wz1;
  thrust::host_vector<Real> wx1;
  thrust::host_vector<Real> wy1;
  thrust::host_vector<Real> wz1;

  thrust::host_vector<int> i_indx;
  thrust::host_vector<int> j_indx;
  thrust::host_vector<int> k_indx;
  thrust::device_vector<int> d_i_indx;
  thrust::device_vector<int> d_j_indx;
  thrust::device_vector<int> d_k_indx;
  
 public:
  
  /* Constructors */
  
  
  /*********************************
   *                               *
   *         Strfac::Strfac        *
   *                               *
   *********************************/
  Strfac() {}
  
  //TODO: pointer to inion can go away and cause a problem!
  Strfac(Ion* inion, CommonData<Real> *incd)
    {
      cd = incd;
      myion  = inion;
      
      nxh = cd->nx/2;
      nfft3d = (nxh+1)*cd->ny*cd->nz;
      
      wx1 = thrust::host_vector<Real>(2*(myion->nion)*cd->nx);
      wy1 = thrust::host_vector<Real>(2*(myion->nion)*cd->nx);
      wz1 = thrust::host_vector<Real>(2*(myion->nion)*cd->nx);
      i_indx = thrust::host_vector<int>(nfft3d);
      j_indx = thrust::host_vector<int>(nfft3d);
      k_indx = thrust::host_vector<int>(nfft3d);
      for (int k=0; k < cd->nz;   ++k)
	for (int j=0; j< cd->ny;   ++j)
	  for (int i=0; i<= nxh; ++i) {
	    int indx = i + j*(nxh+1) + k*(nxh+1)* cd->ny;
	    i_indx[indx] = i;
	    j_indx[indx] = j;
	    k_indx[indx] = k;
	  }
      d_i_indx = i_indx; d_j_indx = j_indx; d_k_indx = k_indx;
    }
  
  
  /*********************************
   *                               *
   *          Strfac::phafac       *
   *                               *
   *********************************/
  void phafac()
  {
    Real a,b,sw1,sw2,sw3;
    Real cw1x,cw2x,cw3x;
    Real cw1y,cw2y,cw3y;

    int nx = cd->nx;
    int ny = cd->ny;
    int nz = cd->nz;
    
    Real pi  = 4.00*atan(1.0);
    
    int nxh = nx/2;
    int nyh = ny/2;
    int nzh = nz/2;
    
    for (int i=0; i<(myion->nion); ++i) {
	sw1 = cd->unitg[0]*myion->rion1[0+3*i]
          + cd->unitg[1]*myion->rion1[1+3*i]
          + cd->unitg[2]*myion->rion1[2+3*i]+pi;
	sw2 = cd->unitg[3]*myion->rion1[0+3*i]
          + cd->unitg[4]*myion->rion1[1+3*i]
          + cd->unitg[5]*myion->rion1[2+3*i]+pi;
	sw3 = cd->unitg[6]*myion->rion1[0+3*i]
          + cd->unitg[7]*myion->rion1[1+3*i]
          + cd->unitg[8]*myion->rion1[2+3*i]+pi;
	
	cw1x=cos(sw1); cw1y=-sin(sw1);
	cw2x=cos(sw2); cw2y=-sin(sw2);
	cw3x=cos(sw3); cw3y=-sin(sw3);
	
	wx1[2*i*nx] = 1.0; wx1[2*i*nx+1] = 0.0;
	wy1[2*i*ny] = 1.0; wy1[2*i*ny+1] = 0.0;
	wz1[2*i*nz] = 1.0; wz1[2*i*nz+1] = 0.0;
	for (int k=1; k<=nxh; ++k) {
	  a = wx1[2*(k-1 + i*nx)];
	  b = wx1[2*(k-1 + i*nx)+1];
	  wx1[2*(k + i*nx)]   = a*cw1x - b*cw1y;
	  wx1[2*(k + i*nx)+1] = a*cw1y + b*cw1x;
	  wx1[2*(nx-k + i*nx)]   =  wx1[2*(k + i*nx)];
	  wx1[2*(nx-k + i*nx)+1] = -wx1[2*(k + i*nx)+1];
	}
	for (int k=1; k<=nyh; ++k) {
	  a = wy1[2*(k-1 + i*ny)];
	  b = wy1[2*(k-1 + i*ny)+1];
	  wy1[2*(k + i*ny)]   = a*cw2x - b*cw2y;
	  wy1[2*(k + i*ny)+1] = a*cw2y + b*cw2x;
	  wy1[2*(ny-k + i*ny)]   =  wy1[2*(k + i*ny)];
	  wy1[2*(ny-k + i*ny)+1] = -wy1[2*(k + i*ny)+1];
	}
	for (int k=1; k<=nzh; ++k) {
	  a = wz1[2*(k-1 + i*nz)];
	  b = wz1[2*(k-1 + i*nz)+1];
	  wz1[2*(k + i*nz)]   = a*cw3x - b*cw3y;
	  wz1[2*(k + i*nz)+1] = a*cw3y + b*cw3x;
	  wz1[2*(nz-k + i*nz)]   =  wz1[2*(k + i*nz)];
	  wz1[2*(nz-k + i*nz)+1] = -wz1[2*(k + i*nz)+1];
	}
	
	wx1[2*(nxh+i*nx)] = 0.0; wx1[2*(nxh+i*nx)+1] = 0.0;
	wy1[2*(nyh+i*ny)] = 0.0; wy1[2*(nyh+i*ny)+1] = 0.0;
	wz1[2*(nzh+i*nz)] = 0.0; wz1[2*(nzh+i*nz)+1] = 0.0;
      }
      d_wx1 = wx1; d_wy1 = wy1; d_wz1 = wz1;
  }
  
  inline void strfac_sub(const int npack,
			 const Real exi[],
			 const Real exj[],
			 const Real exk[],
			 Real strx[])
  {
    Real ai,aj,ak,c,d;
    Real bi,bj,bk;
    for (int i=0; i<npack; ++i) {
      ai = exi[2*i_indx[i]]; bi = exi[2*i_indx[i]+1];
      aj = exj[2*j_indx[i]]; bj = exj[2*j_indx[i]+1];
      ak = exk[2*k_indx[i]]; bk = exk[2*k_indx[i]+1];
      c  = aj*ak - bj*bk;
      d  = aj*bk + ak*bj;
      //      cerr << "strfac_sub i " << i << " npack " << npack << " ai " << ai << " aj " << aj << " ak " << ak << " c " << c << " d " << d << endl;
      strx[2*i]   = (ai*c - bi*d);
      strx[2*i+1] = (ai*d + bi*c);
    }
  }
  
  /*********************************
   *                               *
   *       Strfac::strfac_pack      *
   *                               *
   *********************************/
  inline void strfac_pack(const int ii, Real *ss) {
    strfac_sub(nfft3d,
	       &wx1[2*ii*cd->nx],
	       &wy1[2*ii*cd->ny],
	       &wz1[2*ii*cd->nz],
	       ss);
  }
  inline void strfac_sub(const int npack,
			 const thrust::host_vector<Real> exi, int exiOffset,
			 const thrust::host_vector<Real> exj, int exjOffset,
			 const thrust::host_vector<Real> exk, int exkOffset,
			 thrust::host_vector<Real2> strx)
  {
    Real ai,aj,ak,c,d;
    Real bi,bj,bk;
    for (int i=0; i<npack; ++i) {
      ai = exi[exiOffset+2*i_indx[i]]; bi = exi[exiOffset+2*i_indx[i]+1];
      aj = exj[exjOffset+2*j_indx[i]]; bj = exj[exjOffset+2*j_indx[i]+1];
      ak = exk[exkOffset+2*k_indx[i]]; bk = exk[exkOffset+2*k_indx[i]+1];
      c  = aj*ak - bj*bk;
      d  = aj*bk + ak*bj;
      strx[i].x = (ai*c - bi*d);
      strx[i].y = (ai*d + bi*c);
    }
  }
  
  inline void strfac_pack(const int ii, thrust::host_vector<Real2> &ss) {
    strfac_sub(nfft3d, wx1,(2*ii*cd->nx), wy1,(2*ii*cd->ny),
	       wz1,(2*ii*cd->nz), ss);
  }

  struct strfac_sub_functor
  {
    const Real* exi; const Real* exj; const Real* exk;
    const int* i_indx; const int* j_indx; const int* k_indx;
    Real2* strx;
  strfac_sub_functor(const Real* _exi, const Real* _exj, const Real* _exk,
		     const int* _i_indx, const int* _j_indx, const int* _k_indx,
		     Real2* _strx) :
    exi(_exi), exj(_exj), exk(_exk), i_indx(_i_indx), j_indx(_j_indx), k_indx(_k_indx),
      strx(_strx) {};

    __device__
    void operator()(unsigned int i)
    {
      Real ai = exi[2*i_indx[i]]; Real bi = exi[2*i_indx[i]+1];
      Real aj = exj[2*j_indx[i]]; Real bj = exj[2*j_indx[i]+1];
      Real ak = exk[2*k_indx[i]]; Real bk = exk[2*k_indx[i]+1];
      Real c  = aj*ak - bj*bk;
      Real d  = aj*bk + ak*bj;
      strx[i].x = (ai*c - bi*d);
      strx[i].y = (ai*d + bi*c);
    }
  };
  void strfac_sub(const int npack,
			 const thrust::device_vector<Real> exi, int exiOffset,
			 const thrust::device_vector<Real> exj, int exjOffset,
			 const thrust::device_vector<Real> exk, int exkOffset,
			 thrust::device_vector<Real2> strx)
  {
    strfac_sub_functor unary_op(thrust::raw_pointer_cast(&exi[exiOffset]),
				thrust::raw_pointer_cast(&exj[exjOffset]),
				thrust::raw_pointer_cast(&exk[exkOffset]),
				thrust::raw_pointer_cast(&d_i_indx[0]),
				thrust::raw_pointer_cast(&d_j_indx[0]),
				thrust::raw_pointer_cast(&d_k_indx[0]),
				thrust::raw_pointer_cast(&strx[0]) );
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->nfft3d),
		     unary_op);
  }
  
  inline void strfac_pack(const int ii, thrust::device_vector<Real2> &ss) {
    strfac_sub(nfft3d, d_wx1,(2*ii*cd->nx), d_wy1,(2*ii*cd->ny),
	       d_wz1,(2*ii*cd->nz), ss);
  }
};
#endif
  
  
