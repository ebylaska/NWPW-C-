#ifndef BOXARRAY_H
#define BOXARRAY_H

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
using namespace std;

#include <string.h>
#include <cublas.h>
#include <cufft.h>
#include "psiHeader.h"
#include "rtdb.h"
#include "Ion.h"
#include "Strfac.h"
#include "control.h"
#include "CommonData.h"
#include "PseudopotentialData.h"

// Thrust headers
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/copy.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <cmath>

template <typename Real, typename Real2>
class boxArray {
 protected:
  cufftHandle fftPlanMany_C2R, fftPlanMany_R2C;
  cufftHandle fftPlan_C2R, fftPlan_R2C;
  Ion* myIon;
  // computational variables
  thrust::device_vector<Real2> psiVect[2],Hpsi;
  thrust::device_vector<Real> dn,psi_r,vall;
  thrust::device_vector<Real> dng;
  thrust::device_vector<Real2> *psi1, *psi2;
  Strfac<Real, Real2> h_strfac;
  PseudoPotentialData<Real> h_pseudoData;

  // housekeeping variables
  PsiHeader<Real> psiHeader;
  CommonData<Real> *cd;
  //PseudoPotentialData<Real> *pspData;
  int neall;

  // methods
  void errorExit(char* s) { cerr << s << endl; exit(1);}

  inline void switchPsi()
    {
      thrust::device_vector<Real2> *tmp = psi1;
      psi1 = psi2;
      psi2 = tmp;
    }

 public:
  // stuff for the k-space operator
  struct kspace_op
  {
    const Real2* psi;
    const Real* tg;
    const int tg_size;

  kspace_op( const Real2* _psi,
	     const Real* _tg,
	     const int _tg_size
	     ): psi(_psi), tg(_tg), tg_size(_tg_size) {};
    __host__ __device__
    Real2 operator()(const Real2& psi1) const { 
      int j = &psi1 - psi;
      Real tg_j = tg[j%tg_size]; // Note older GPUs will serialize on tg
      register Real2 ret;
      ret.x = ( psi1.x * tg_j );
      ret.y = ( psi1.y * tg_j );
      return ret;
    }
  };

  inline void kspace_operator()
  {
#ifdef Kspace_HOST_SIDE
    thrust::host_vector<Real2> tmp = *psi1;
    for(int i=0; i < Hpsi.size(); i++) {
      double tg_j = cd->h_tg[i%cd->nfft3d]; 
      tmp[i].x *= tg_j;
      tmp[i].y *= tg_j;
    }
    Hpsi = tmp; // big transfer
#else
    kspace_op unary_op( thrust::raw_pointer_cast(&((*psi1)[0])),
			thrust::raw_pointer_cast(&((cd->d_tg)[0])),
			cd->nfft3d);
    thrust::transform(psi1->begin(), psi1->end(), Hpsi.begin(), unary_op);
#endif
    v_nonlocal(); 
  }
  
  // calculated values from inner_loop
  Real *E[], deltae, deltar, deltac;
  ~boxArray() {
    delete cd;
    delete myIon;
  }
  boxArray(char* rootFilename) {
    int fd;
    struct stat statbuf;
    string s(rootFilename);
    string filename(s + string(".movecs"));
    string dbFilename(s + string(".db"));
    
    if( (fd=open(filename.c_str(), O_RDONLY)) < 0) errorExit((char*)"cannot open file");
    if(fstat(fd,&statbuf) < 0) errorExit((char*)"cannot stat file");
    cerr << "NWChem filename " << filename << " size " << statbuf.st_size << endl;
    char *mmapBuf = (char *) mmap(0,statbuf.st_size,PROT_READ,MAP_SHARED,fd,0);
    RTDB myrtdb( dbFilename.c_str(), "old");
    control_read(myrtdb);
    
    PsiHeader<Real> tmp(mmapBuf);
    psiHeader = tmp;
    char *buf = mmapBuf + psiHeader.get_nBytes();
    
    cd = new CommonData<Real>(psiHeader, myrtdb);
    
    neall = psiHeader.get_ne(0) + psiHeader.get_ne(1);
    
    // allocate and read wavefunctions 
    // bummer! the nx+2 prevents power of two optimizations
    int size = neall * (cd->nx+2) * cd->ny * cd->nz;
    psiVect[0] = thrust::device_vector<Real2>(size);
    Real2 czero={0.,0.};
    thrust::fill(psiVect[0].begin(),psiVect[0].end(), czero);
    psiVect[1] = psiVect[0];

    psi1 = &psiVect[0]; psi2 = &psiVect[1];
    
    Hpsi = thrust::device_vector<Real2>(size);
    psi_r = thrust::device_vector<Real>( 2*size );
    
    {
      thrust::host_vector<Real2> tmp(size);
      int index=0;
      for (int i=0; i < neall*(cd->nfft3d); i++) {
	tmp[i].x = (Real) ((double*)buf)[index++];
	tmp[i].y = (Real) ((double*)buf)[index++];
      }
      (*psi1) = tmp; // move to GPU
    }
    buf += neall*2*cd->nfft3d*sizeof(double);
    switchPsi();
    
    munmap(mmapBuf,statbuf.st_size);
    close(fd);
    
    myIon = new Ion(myrtdb);
    h_strfac = Strfac<Real, Real2>(myIon, cd);
    
    // read in Pseudopotential data
    h_pseudoData = PseudoPotentialData<Real>(rootFilename, myIon, cd->nfft3d);
    cerr << "npsp " << h_pseudoData.npsp << endl;
    cerr << "vnl[0][0] " << h_pseudoData.h_vnl[0][0] << endl;
    
    cublasInit();
    initFFTs();
    int nsize = psiHeader.get_ispin() * cd->nfft3d;
    dn = thrust::device_vector<Real>( cd->n2ft3d );
    dng = thrust::device_vector<Real>( cd->n2ft3d );
    vall = thrust::device_vector<Real>( cd->n2ft3d );
  } 
    // stuff for the Ke operator
  struct Ke_unary_op {
    const Real2* psi;
    const Real* tg;
    const Real* summer;
    const int nfft3d;
    const int neall;
    
  Ke_unary_op(const Real2* _psi, const Real* _tg, const Real* _summer,
	      const int _nfft3d, const int _neall)
  : psi(_psi), tg(_tg), summer(_summer), nfft3d(_nfft3d), neall(_neall) {};
    
    __device__
    Real operator()(unsigned int tid)
    {
      register Real sum=0.;
      register Real tmp = tg[tid] * summer[tid];
      
      for(int i=0; i < neall; i++) {
	register Real2 p = psi[tid + i*nfft3d];
	sum += (p.x * p.x + p.y * p.y) * tmp;
      }
      return(sum);
    }
  };
  
  inline Real kineticEnergy()
  {
    //#define KE_HOST_SIDE
#ifdef KE_HOST_SIDE
    Real sum=0.0;
    thrust::host_vector<Real2> h_psi = *psi1;
    
    for(int i=0; i < h_psi.size(); i++) {
      int j = i%cd->nfft3d;
      sum += (  h_psi[i].x * h_psi[i].x + h_psi[i].y * h_psi[i].y)
	* cd->h_tg[j] * cd->h_summer[1][j];
    }
#else
    Ke_unary_op my_op(thrust::raw_pointer_cast(&(*psi1)[0]),
			 thrust::raw_pointer_cast(&cd->d_tg[0]),
			 thrust::raw_pointer_cast(&cd->d_summer[1][0]),
			 cd->nfft3d,
			 neall);
    Real sum = thrust::transform_reduce(thrust::counting_iterator<unsigned int>(0),
                                        thrust::counting_iterator<unsigned int>(cd->nfft3d),
                                        my_op,
                                        (Real) 0.,
                                        thrust::plus<Real>());
#endif
    if (psiHeader.get_ispin()==1) sum *= 2.0;
    return -sum;
  }

  Real normCheck(thrust::device_vector<Real2> &a)
  {
    Real sum=0.0;
    thrust::host_vector<Real2> tmp(a);

    for(int i=0; i < neall; i++) {
      Real2 *wavef = &tmp[i*cd->nfft3d];
      for(int j=0; j < cd->nfft3d; j++){
	sum += (  wavef[j].x*wavef[j].x + wavef[j].y*wavef[j].y)
	  * cd->h_summer[1][j];
      }
      cerr << "Norm check: cumulative sum " << sum << " " << i << endl;
    }
    return sum;
  }

  Real normCheck()
  {
    Real sum=0.0;
    thrust::host_vector<Real2> tmp(*psi1);

    for(int i=0; i < neall; i++) {
      Real2 *wavef = &tmp[i*cd->nfft3d];
      for(int j=0; j < cd->nfft3d; j++){
	sum += (  wavef[j].x*wavef[j].x + wavef[j].y*wavef[j].y)
	  * cd->h_summer[1][j];
      }
      cerr << "Norm check: cumulative sum " << sum << " " << i << endl;
    }
    return sum;
  }

  Real psiHpsiCheck()
  {
    Real sum=0.0;
    thrust::host_vector<Real2> tmp1 = *psi1;
    thrust::host_vector<Real2> tmp2 = Hpsi;

    for(int i=0; i < neall; i++) {
      Real2 *wavef1 = &tmp1[i*cd->nfft3d];
      Real2 *wavef2 = &tmp2[i*cd->nfft3d];
      for(int j=0; j < cd->nfft3d; j++){
	sum += (  wavef1[j].x*wavef2[j].x + wavef1[j].y*wavef2[j].y)
	  * cd->h_summer[1][j];
      }
      cerr << "<psiHpsi> cumulative sum " << sum << " " << i << endl;
    }
    return sum;
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
  inline void rcFFT_neall( Real* A, Real* B ) { rcFFT_(fftPlanMany_R2C, A,B); }
  inline void crFFT_neall( Real* A, Real* B ) { crFFT_(fftPlanMany_C2R, A,B); }
  inline void rcFFT_single( Real* A, Real* B ) { rcFFT_(fftPlan_R2C, A,B); }
  inline void crFFT_single( Real* A, Real* B ) { crFFT_(fftPlan_C2R, A,B); }

  struct calc_dng_functor
  {
    const Real scal1;
    const int ispin,n2ft3d;
    const Real* dn;
    Real* dng;
  calc_dng_functor(const Real _scal1, const int _ispin, const int _n2ft3d,
		   const Real *_dn, Real *_dng) : scal1(_scal1), ispin(_ispin),
      n2ft3d(_n2ft3d), dn(_dn), dng(_dng) {};

    __device__
    void operator()(unsigned int i)
    {
      //TODO eliminate subtraction by 1
      dng[i] = dn[i] + dn[(ispin-1)*n2ft3d + i];
      dng[i] *= scal1;
    }
  };
  void calc_dng(Real scal1, int ispin, 
		thrust::device_vector<Real> &dn,
		thrust::device_vector<Real> &dng)
  {
    calc_dng_functor unary_op( scal1, ispin, cd->n2ft3d,
			       thrust::raw_pointer_cast(&dn[0]),
			       thrust::raw_pointer_cast(&dng[0]));

    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->n2ft3d),
		     unary_op);
  }

  struct calc_vc_functor 
  {
    const Real* dng;
    const Real* masker;
    Real* vc;
  calc_vc_functor(const Real *_dng,
		   const Real *_masker,
		   Real *_vc) : dng(_dng), masker(_masker), vc(_vc) {};

    __device__
    void operator()(unsigned int i) {
      vc[i] = masker[i>>1] *  dng[i]; 
    }
  };
  void vcolumb(const thrust::device_vector<Real> &d_dng,
	       thrust::device_vector<Real> &d_vc) 
  {
    //for(int i=0; i < cd->n2ft3d; i++) h_vc[i] = cd->h_masker[1][i/2] *  h_dng[i]; 
    calc_vc_functor unary_op(thrust::raw_pointer_cast(&d_dng[0]),
			     thrust::raw_pointer_cast(&cd->d_masker[1][0]),
			     thrust::raw_pointer_cast(&d_vc[0]) );
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->n2ft3d),
		     unary_op);
  }
  struct calc_vall_functor
  {
    const Real scal2;
    const Real* vc;
    const Real* vl;
    Real* vall;
  calc_vall_functor(const Real _scal2, 
		    const Real *_vc, 
		    const Real *_vl,
		    Real* _vall) : scal2(_scal2), vc(_vc), vl(_vl), vall(_vall) {};
    __device__
    void operator()(unsigned int i)
    {
      vall[i]   = scal2 * vl[i] + vc[i];
    }
  };
  void calc_vall(const Real scal2,
		 const thrust::device_vector<Real> &d_vc, 
		 const thrust::device_vector<Real> &d_vl, 
		 thrust::device_vector<Real> &d_vall) {
    calc_vall_functor unary_op(scal2,
			       thrust::raw_pointer_cast(&d_vc[0]),
			       thrust::raw_pointer_cast(&d_vl[0]),
			       thrust::raw_pointer_cast(&d_vall[0]) );
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->n2ft3d),
		     unary_op);
  }

  struct calc_cc_daxpy_functor
  {
    const Real a;
    const Real* vpsi;
    Real2* hpsi;
  calc_cc_daxpy_functor(const Real _a, 
		    const Real *_vpsi, 
		    Real2* _hpsi) : a(_a), vpsi(_vpsi), hpsi(_hpsi) {};
    __device__
    void operator()(unsigned int i)
    {
      hpsi[i].x += a * vpsi[2*i];
      hpsi[i].y += a * vpsi[2*i+1];
    }
  };
  void calc_cc_daxpy(Real scal1,
		     thrust::device_vector<Real> &vpsi,
		     thrust::device_vector<Real2> &hpsi,
		     int indx1)
  {
    //for(int i=0; i < cd->nfft3d; i++) {
      //h_Hpsi[indx1+i].x += (-scal1)* h_vpsi[2*i];
      //h_Hpsi[indx1+i].y += (-scal1)* h_vpsi[2*i+1];
    //}
    calc_cc_daxpy_functor unary_op(-scal1,
			       thrust::raw_pointer_cast(&vpsi[0]),
			       thrust::raw_pointer_cast(&hpsi[indx1]) );
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->nfft3d),
		     unary_op);
  }


  struct calc_update_functor
  {
    const Real2 *Hpsi;
    const Real2 *psi1;
    const Real2 dte;
    Real2 *psi2;
  calc_update_functor(const Real2* _Hpsi,
		      const Real2 _dte, 
		      const Real2* _psi1, 
		      Real2* _psi2) : Hpsi(_Hpsi), dte(_dte), psi1(_psi1), psi2(_psi2) {};
    __device__
    void operator()(unsigned int i)
    {
      psi2[i].x = Hpsi[i].x * dte.x;
      psi2[i].y = Hpsi[i].y * dte.y;
      psi2[i].x += psi1[i].x;
      psi2[i].y += psi1[i].y;
    }
  };
  void calc_update( thrust::device_vector<Real2> &Hpsi,
		    Real dte,
		    thrust::device_vector<Real2> &psi1,
		    thrust::device_vector<Real2> &psi2)
  {
    Real2 dte2 = {dte,dte};
    calc_update_functor unary_op(thrust::raw_pointer_cast(&Hpsi[0]),
				 dte2,
				 thrust::raw_pointer_cast(&psi1[0]),
				 thrust::raw_pointer_cast(&psi2[0])
				 );
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->nfft3d),
		     unary_op);
    //for(int i=0; i < h_Hpsi.size(); i++) {
      // gg_SMul(dte,Hpsi,psi2); 
      //h_psi2[i].x = cd->dte * h_Hpsi[i].x;
      //h_psi2[i].y = cd->dte * h_Hpsi[i].y;
      // gg_Sum2(psi1,psi2);
      //h_psi2[i].x += h_psi1[i].x;
      //h_psi2[i].y += h_psi1[i].y;
    //}
  }


  void inner_loop()
  {
    Real Ke;
    int ispin = psiHeader.get_ispin();

    int it_in = control_loop(0);
    //it_in = 1000;

    h_strfac.phafac();

    // compute Vlocal(G)
    thrust::host_vector<Real> h_vl(cd->n2ft3d);
    thrust::fill(h_vl.begin(),h_vl.end(), 0.); //Note: Eric's code does not initalize vl!
    v_local(h_vl);
    thrust::device_vector<Real> d_vl = h_vl;

    // Debug stuff (copy psi2 -> psi1 and calculate Initial Ke
    *psi1 = *psi2; Ke = kineticEnergy();
    cerr << "Initial Ke is " << Ke << endl;
    Ke = normCheck();
    cerr << "Initial norm is " << Ke << endl;

    cerr << "inner_loop, it_in " << it_in << endl;

    //TODO: locate this allocation better
    thrust::device_vector<Real> d_vc(cd->n2ft3d);

    for(int iter=0; iter < it_in; iter++) {
      Real2 czero={0.,0.};
      thrust::fill(Hpsi.begin(),Hpsi.end(), czero);

      // copy psi1(G,n,ms) = psi2(G,n,ms). Done via pointer switch
      switchPsi();

      // calculates Hpsi(G,n,ms) = KE*psi1(G,n,ms) + VNL*psi1(G,n,ms)
      kspace_operator();

      //DEBUG
      Ke = -2.0*psiHpsiCheck();
      cerr << "Initial <psi|Hpsi> is " << Ke << endl;
      
      // **** add r-space operator to k-space operator ****
      // convert  psi(R,n,ms)= FFT(psi1(G,n,ms)). 
      cerr << "NORM CHECK " << endl;
      normCheck();
/*
      crFFT_neall( (Real*) thrust::raw_pointer_cast(&(*psi1)[0]),
		   (Real*) thrust::raw_pointer_cast(&psi_r[0]) );
      cerr << "NORM CHECK " << endl;
      normCheck();
      
      // compute the density dn(R,ms)= sum(n=1,nelc(ms)) psi(R,n,ms)*psi(R,n,ms)
      const Real scal2 = 1./cd->omega;
      aSumSqr(scal2,psi_r,dn);

#ifdef DEBUG_DNG
      { // DEBUG
	thrust::host_vector<Real> h_dn=dn; // TODO get rid of host allocation
	Real sum=0.0;
	for (int i=0; i< cd->n2ft3d; ++i) sum += h_dn[i];
	cerr << "DEBUG***** dn=" << sum*((30.0*30.0*30.0)/(48.0*48.0*48.0)) << endl;
      }
#endif

      // compute the FFT of dn, rho(G) = fft(dn(R,1)+dn(R,ispin)) R->C
      const Real scal1 = 1./((Real) cd->nx * (Real) cd->ny * (Real) cd->nz);
      
      //TODO: NOTE using h_dn calculated above without a transfer
      //  calculate rrr_Sum(dn,&dn[(ispin-1)*n2ft3d],tmp);
      /// r_SMul(scal1,tmp);
      calc_dng(scal1, ispin, dn, dng); 
      
      rcFFT_single( (Real*) thrust::raw_pointer_cast(&dng[0]),
		    (Real*) thrust::raw_pointer_cast(&dng[0]) );
      
      // compute Vcoulomb(G) = (4*pi/G^2)*.rho(G)
      /// vcoulomb(dng,vc);
      vcolumb(dng,d_vc);

      // apply rspace operators
      // compute vall(R) = FFT(Vcoulomb(G)+vlocal(G) C->R

      //   cc_SMul(0,scal2,vl,vall);
      //   cc_Sum2(0,vc,vall);
      calc_vall(scal2, d_vc, d_vl,vall);

      //mygrid->cr_fft3d(vall);
      //vall = h_vall; //TODO eliminate this transfer
      crFFT_single(thrust::raw_pointer_cast(&vall[0]),
		   thrust::raw_pointer_cast(&vall[0]));


      {
	// compute Vxc(R,ms) ~= alpha*(dn(R,1)+dn(R,ispin))^1/3
	thrust::host_vector<Real> h_xcp(2*ispin*cd->nfft3d);
	thrust::host_vector<Real> h_xce(2*ispin*cd->nfft3d);
	thrust::host_vector<Real> h_dn=dn; // TODO get rid of host allocation

	//TODO convert to GPU-based version
	v_exc(ispin,cd->n2ft3d,(Real*)&h_dn[0],(Real*) &h_xcp[0],(Real*) &h_xce[0]); 
	
	thrust::device_vector<Real> d_xcp = h_xcp;
	thrust::device_vector<Real> d_xc3 = h_xce;

	// Hpsi(G,n,ms) += FFT(Vall(R)+Vxc(R,ms))*psi(R,n,ms)) r->C
	int indx1 = 0;
	int indx2 = 0;
	const Real scal1 = 1./((Real) cd->nx * (Real) cd->ny * (Real) cd->nz);
	thrust::device_vector<Real> d_tmp(cd->n2ft3d);
	thrust::device_vector<Real> d_vpsi(cd->n2ft3d);
	for (int ms=0; ms < ispin; ++ms) {
	  //mygrid->rrr_Sum(vall,&xcp[ms*n2ft3d],tmp);
	  //for(int i=0; i < cd->n2ft3d; i++) h_tmp[i] = h_vall[i] + h_xcp[ms*cd->n2ft3d+i];
	  d_vpsi = d_tmp;
	  int msOffset = ms+cd->n2ft3d;
	  thrust::transform(d_xcp.begin()+msOffset, d_xcp.begin()+msOffset+cd->n2ft3d,
			    d_tmp.begin(), d_tmp.end(),
			    thrust::plus<Real>());

	  //TODO: make this a multiplan fft
	  for (int ne=0; ne< psiHeader.get_ne(ms); ++ne) {
	    //mygrid->rrr_Mul(tmp,&psi_r[indx2],vpsi);
	    //for(int i=0; i < cd->n2ft3d; i++) h_vpsi[i] = h_psi_r[indx2+i] * h_tmp[i];
	    d_vpsi = d_tmp;
	    thrust::transform(psi_r.begin()+indx2, psi_r.begin()+indx2+cd->n2ft3d,
			      d_tmp.begin(), d_tmp.end(),
			      thrust::multiplies<Real>());
	    //mygrid->rc_fft3d(vpsi);
	    Real *pt = (Real*) thrust::raw_pointer_cast(&d_vpsi[0]);
	    rcFFT_single( (Real*) thrust::raw_pointer_cast(&d_vpsi[0]),
			  (Real*) thrust::raw_pointer_cast(&d_vpsi[0]) );
	    //mygrid->c_pack(1,vpsi);
	    //mygrid->cc_daxpy(1,(-scal1),vpsi,&Hpsi[indx1]);
	    calc_cc_daxpy(-scal1,d_vpsi,Hpsi,indx1);

	    indx1 += cd->nfft3d;
	    indx2 += cd->n2ft3d;
	  }
	}
      }
          
      //*** update ****
*/

      // psi2(G,n,ms) = psi1(G,n,ms) + dte*Hpsi(G,n,ms)
      // steepest descents gg_SMul(dte,Hpsi,psi2); gg_Sum2(psi1,psi2);
      calc_update(Hpsi, cd->dte, *psi1, *psi2);
	
      // compute lagrange multipliers, Lambda(n,n,ms) = Function(dte,psi1,psi2)
      // psi2(G,n,ms) += Lambda(n,n',ms) * psi1(G,n',ms)
      check_Vl_Dng(d_vl, dng, cd->d_summer[0]);
      ggm_lambda();

      /* loop finished */
    }

    Ke = kineticEnergy();
    cerr << "END Ke is " << Ke << endl;

    Real2 czero={0.,0.};
    thrust::fill(Hpsi.begin(),Hpsi.end(), czero);
    v_nonlocal();
    Ke = -2.0*psiHpsiCheck();
    cerr << "END nonlocal is " << Ke << endl;

  }
  void check_Vl_Dng(thrust::device_vector<Real> vl,
		    thrust::device_vector<Real> dng,
		    thrust::device_vector<Real> summer)
  {
    thrust::host_vector<Real> h_vl = vl;
    thrust::host_vector<Real> h_dng = dng;
    thrust::host_vector<Real> h_summer = summer;

    cerr << "size vl " << vl.size() << endl;
    cerr << "size dng " << dng.size() << endl;
    cerr << "size summer " << summer.size() << endl;

    double sum=0.;
    for(int i=0; i < vl.size(); i++) sum += h_vl[i];
    cerr << "vl sum " << sum << endl;

    sum=0.;
    for(int i=0; i < dng.size(); i++) sum += h_dng[i];
    cerr << "dng sum " << sum << endl;

    sum=0.;
    for(int i=0; i < dng.size(); i++)
      sum += h_vl[i] * h_dng[i] * h_summer[i/2];
    cerr << "dng.size()\n\n";
    cerr << "dng[0]" << h_dng[0] << endl;
    cerr << "dng[1]" << h_dng[1] << endl;
    cerr << "Vl[0]" << h_vl[0] << endl;
    cerr << "Vl[1]" << h_vl[1] << endl;
    cerr << "Vl[2]" << h_vl[2] << endl;
    cerr << "Vl[3]" << h_vl[3] << endl;
    cerr << "Vl[4]" << h_vl[4] << endl;
    cerr << "Vl[5]" << h_vl[5] << endl;
    cerr << "Vlocal Sum " << sum << endl;
  }

#include "ggmLambda.h"

  void inline initFFTs()
  {
    //TODO: double check variables
    int dim[] = { cd->nx, cd->ny, cd->nz};
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
    
    cufftSetCompatibilityMode(fftPlanMany_C2R, CUFFT_COMPATIBILITY_NATIVE);
    cufftSetCompatibilityMode(fftPlanMany_R2C, CUFFT_COMPATIBILITY_NATIVE);
    cufftSetCompatibilityMode(fftPlan_C2R, CUFFT_COMPATIBILITY_NATIVE);
    cufftSetCompatibilityMode(fftPlan_R2C, CUFFT_COMPATIBILITY_NATIVE);
  }

  void inline hr_aSumSqr(const Real alpha,
			 thrust::host_vector<Real> &psi_r,
			 thrust::host_vector<Real> &dn)
  {
    int n,ms,k,indx0,indx1;
    int ispin = psiHeader.get_ispin();

    thrust::fill(dn.begin(),dn.end(), 0.);

    indx0 = 0;
    indx1 = 0;
    for (ms=0; ms < ispin; ++ms) {
      for (n=0; n < neall; ++n) {
	for (k=0; k < 2*cd->nfft3d; ++k) {
	  dn[indx0+k] += alpha*psi_r[indx1+k]*psi_r[indx1+k];
	}
	indx1 += 2*cd->nfft3d;
      }
      indx0 += 2*cd->nfft3d;
    }
  }

  struct aSumSqr_functor
  {
    const Real* psi_r;
    Real* dn;
    const Real alpha;
    const int ispin, neall, n2ft3d;
  aSumSqr_functor(const Real* _psi_r, Real* _dn, Real _alpha, int _ispin, int _neall, int _n2ft3d) :
    psi_r(_psi_r), dn(_dn), alpha(_alpha),
      ispin(_ispin), neall(_neall), n2ft3d(_n2ft3d) {};

    __device__
    void operator()(unsigned int k)
    {
      int indx0 = 0, indx1 = 0;
      for (int ms=0; ms < ispin; ++ms) {
	for (int n=0; n < neall; ++n) {
	  dn[indx0+k] += alpha*psi_r[indx1+k]*psi_r[indx1+k];
	  indx1 += n2ft3d;
	}
	indx0 += n2ft3d;
      }
    }
  };
  void inline aSumSqr(const Real alpha,
			 thrust::device_vector<Real> &psi_r,
			 thrust::device_vector<Real> &dn)
  {
    int ispin = psiHeader.get_ispin();

    thrust::fill(dn.begin(),dn.end(), 0.);

    aSumSqr_functor unary_op(thrust::raw_pointer_cast(&psi_r[0]),
			     thrust::raw_pointer_cast(&dn[0]),
			     alpha, ispin, neall, cd->n2ft3d);

    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->n2ft3d),
		     unary_op);
  }

  // methods for v_nonlocal
  struct tcc_Mul_functor
  {
    const Real* a;
    const Real2* b;
    Real2* c;
  tcc_Mul_functor(const Real* _a,
		  const Real2* _b,
		  Real2* _c): a(_a), b(_b), c(_c) {};
    __device__
    void operator()(unsigned int tid)
    {
      Real tmp=a[tid];
      c[tid].x   = b[tid].x * tmp;
      c[tid].y   = b[tid].y * tmp;
    }
  };
  void inline tcc_Mul(const thrust::device_vector<Real>  &a,
		      const thrust::device_vector<Real2> &b,
		      thrust::device_vector<Real2> &c, const int offset)
  {
    int n = offset+cd->nfft3d;
    assert(a.begin()+n <= a.end()); assert(c.begin()+n <= c.end()); assert(b.begin()+cd->nfft3d <= b.end()); 

    tcc_Mul_functor unary_op(thrust::raw_pointer_cast(&a[offset]),
			     thrust::raw_pointer_cast(&b[0]),
			     thrust::raw_pointer_cast(&c[offset]));
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->nfft3d),
		     unary_op);
  }
  
  void inline tcc_Mul(thrust::host_vector<Real>  &a,
		      thrust::host_vector<Real2> &b,
		      thrust::host_vector<Real2> &c, int offset)
  {
    int n = offset+cd->nfft3d;
    assert(a.begin()+n <= a.end()); assert(c.begin()+n <= c.end()); assert(b.begin()+cd->nfft3d <= b.end()); 

    for (int i=0; i< cd->nfft3d; ++i) {
      c[i].x   = b[i].x * a[i];
      c[i].y   = b[i].y * a[i];
    }
  }
  struct tcc_iMul_functor
  {
    const Real* a;
    const Real2* b;
    Real2* c;
  tcc_iMul_functor(const Real* _a,
		  const Real2* _b,
		  Real2* _c): a(_a), b(_b), c(_c) {};
    __device__
    void operator()(unsigned int tid)
    {
      Real tmp=a[tid];
      c[tid].x   = -b[tid].x * tmp;
      c[tid].y   = b[tid].y * tmp;
    }
  };
  void inline tcc_iMul(const thrust::device_vector<Real>  &a,
		       const thrust::device_vector<Real2> &b,
		       thrust::device_vector<Real2> &c, const int offset)
  {
    int n = offset+cd->nfft3d;
    assert(a.begin()+n <= a.end()); assert(c.begin()+n <= c.end()); assert(b.begin()+cd->nfft3d <= b.end()); 

    tcc_iMul_functor unary_op(thrust::raw_pointer_cast(&a[offset]),
			      thrust::raw_pointer_cast(&b[0]),
			      thrust::raw_pointer_cast(&c[offset]));
    thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		     thrust::counting_iterator<unsigned int>(cd->nfft3d),
		     unary_op);
  }
  void inline tcc_iMul(thrust::host_vector<Real> &a,
		       thrust::host_vector<Real2> &b,
		       thrust::host_vector<Real2> &c, int offset)
  {
    assert(a.begin()+offset+cd->nfft3d <= a.end()); 
    assert(c.begin()+offset+cd->nfft3d <= c.end()); 
    assert(b.begin()+cd->nfft3d <= b.end()); 

    for (int i=0; i< cd->nfft3d; ++i) {
      c[i].x = -b[i].y * a[i];
      c[i].y =  b[i].x * a[i];
    }
  }

  inline static void myAxpy(int n, const Real alpha, Real *x, Real *y) {
    for(int i=0; i < n; i++) y[i] += alpha * x[i];
  }
  inline static void myAxpy(const int n, const Real alpha, const thrust::host_vector<Real> &x, const int xOffset,
			    thrust::host_vector<Real> &y, const int yOffset) {
    for(int i=0; i < n; i++) y[i+yOffset] += alpha * x[i+xOffset];
  }
    
  //TODO Look at doing this in one call on GPU
  inline void Multiply_Gijl_sw1(const int nn, const int nprj, const int nmax, const int lmax,
				thrust::host_vector<int> &n_prj,
				thrust::host_vector<int> &l_prj,
				thrust::host_vector<int> &m_prj,
				thrust::host_vector<Real> &G, 
				thrust::host_vector<Real> &sw1, 
				thrust::host_vector<Real> &sw2)
  {
    int a,b,na,nb;
    int nmax2 = nmax*nmax;
    //int nnn = nn*nprj;
    int nna = nn;

    // fill sw2 with zeros
    thrust::fill(sw2.begin(),sw2.end(), 0.);

    //TODO: make this one call by creating an na and nb vector. Call the single kernel with those
    //vectors to perform all the vector ops concurrently
    for (b=0; b<nprj; ++b)
      for (a=0; a<nprj; ++a)
	if ((l_prj[a]==l_prj[b]) && (m_prj[a]==m_prj[b])) {
	  na = n_prj[a]-1;
	  nb = n_prj[b]-1;
	  //cublasSaxpy(nna, G[nb + na*nmax + nmax2 * l_prj[a]], &sw1[a*nn],1, &sw2[b*nn],1);
	  myAxpy(nna, G[nb + na*nmax + nmax2 * l_prj[a]], sw1,(a*nn), sw2,(b*nn));
	}
  }
  
  //TODO: make this one GPU call
  struct cc_indot_unary_op {
    const Real2* psi; 
    const Real2* prj;
    const Real* summer;
    
  cc_indot_unary_op(const Real2* _psi, const Real2* _prj, const Real* _summer)
  : psi(_psi), prj(_psi), summer(_summer) {};
    
    __device__
    Real operator()(unsigned int tid)
    {
      register Real2 r2psi = psi[tid];
      register Real2 r2prj = prj[tid];
      register Real mySummer = summer[tid];
      return (r2prj.x * mySummer * r2psi.x) + (r2prj.y * mySummer * r2psi.y);
    }
  };

  void cc_indot(const int neall, const int nfft3d, 
		const thrust::device_vector<Real2> &prj, const int prjOffset,
		const thrust::device_vector<Real> &summer,
		const thrust::device_vector<Real2> &psi1,
		thrust::host_vector<Real> &sw1, const int sw1Offset) 
  {
    assert( (prj.begin()+prjOffset+nfft3d) <= prj.end()); 
    assert( (summer.begin()+nfft3d) <= summer.end()); 
    assert( (psi1.begin()+neall*nfft3d+nfft3d) <= psi1.end()); 
    assert( (sw1.begin()+sw1Offset+neall) <= sw1.end()); 

    for (int n=0; n < neall; ++n) {
      cc_indot_unary_op unary_op(thrust::raw_pointer_cast(&psi1[n*nfft3d]),
				 thrust::raw_pointer_cast(&prj[prjOffset]),
				 thrust::raw_pointer_cast(&cd->d_summer[1][0]));
      
      Real tmp =  thrust::transform_reduce(thrust::counting_iterator<unsigned int>(0),
				      thrust::counting_iterator<unsigned int>(nfft3d),
				      unary_op,
				      (Real) 0.,
				      thrust::plus<Real>());
      sw1[sw1Offset + n] = tmp;
    }
  }
  inline void cc_indot(const int neall, const int nfft3d, 
		       const thrust::host_vector<Real2> &prj, const int prjOffset,
		       const thrust::host_vector<Real> &summer,
		       const thrust::host_vector<Real2> &psi1,
		       thrust::host_vector<Real> &sw1, const int sw1Offset) 
  {
    assert(prj.begin()+prjOffset+nfft3d <= prj.end()); 
    assert(summer.begin()+nfft3d <= summer.end()); 
    assert(psi1.begin()+neall*nfft3d+nfft3d <= psi1.end()); 

    for (int n=0; n < neall; ++n) {
      Real sum=0.;
      for(int i=0; i < cd->nfft3d; i++) {
	sum += (prj[prjOffset + i].x * summer[i] * psi1[i+n*nfft3d].x)
	  + (prj[prjOffset + i].y * summer[i] * psi1[i+n*nfft3d].y);
      }
      cerr << "HOST cc_indot " << sum << endl;
      sw1[sw1Offset + n] = sum;
    }
  }
    
  void v_nonlocal()
  {
    int ii,ia,l,sd_function;
    int neall = psiHeader.get_ne(0) + psiHeader.get_ne(1);

    thrust::host_vector<Real> h_sw1(neall*h_pseudoData.nprj_max);
    thrust::host_vector<Real> h_sw2(neall*h_pseudoData.nprj_max);

    thrust::device_vector<Real> d_sw2(neall*h_pseudoData.nprj_max);
    thrust::device_vector<Real2> d_prjtmp(h_pseudoData.nprj_max*cd->nfft3d);
    thrust::device_vector<Real2> d_exi(cd->nfft3d);
    
    for (ii=0; ii<(h_pseudoData.myIonPtr->nion); ++ii) {
      ia = h_pseudoData.myIonPtr->katm[ii];
      if (h_pseudoData.nprj[ia] > 0) {

	/* structure factor */
	h_strfac.strfac_pack(ii, d_exi);
	
	/* generate sw1 and projectors */
	for (l=0; l<h_pseudoData.nprj[ia]; ++ l) {
	  int prjOffset = l*cd->nfft3d;
	  
	  sd_function = !(h_pseudoData.l_projector[ia][l] & 1);
	  if (sd_function) {
	    tcc_Mul( h_pseudoData.d_vnl[ia], d_exi,d_prjtmp, prjOffset);
	  } else {
	    tcc_iMul( h_pseudoData.d_vnl[ia], d_exi,d_prjtmp, prjOffset); 
	  }
	  //cc_pack_indot(1,nn,psi,prj,&h_sw1[l*nn]);
	  // TODO: Double-check (The code in the method looks wrong)
	  cc_indot(neall, cd->nfft3d, d_prjtmp, prjOffset, cd->d_summer[1], *psi1, h_sw1, l*neall);
	}

	/* sw2 = Gijl*sw1 */
	//Multiply_Gijl_sw1(nn,h_pseudoData.nprj[ia],nmax[ia],lmax[ia],
	//	  n_projector[ia],l_projector[ia],m_projector[ia],
	//	  Gijl[ia],sw1,sw2);
	Multiply_Gijl_sw1(neall, h_pseudoData.nprj[ia], h_pseudoData.nmax[ia], h_pseudoData.lmax[ia],
			  h_pseudoData.n_projector[ia], h_pseudoData.l_projector[ia], h_pseudoData.m_projector[ia],
			  h_pseudoData.h_Gijl[ia], h_sw1, h_sw2);
	
	/* do Kleinman-Bylander Multiplication */
	Real scal = 1.0/cd->omega;
	for(int i=0; i < h_sw2.size(); i++) h_sw2[i] *= scal;
	d_sw2 = h_sw2;
	
        int ntmp = h_pseudoData.nprj[ia];
	// dgemm_("N","T",&nshift,&nn,&ntmp, &rmone, prjtmp,&nshift,
	//        sw2,   &nn, &rone, Hpsi,&nshift);
	blasGemm('N', 'T', cd->n2ft3d, neall, ntmp, -1.,
	      (Real*) thrust::raw_pointer_cast(&d_prjtmp[0]), cd->n2ft3d,
	      (Real*) thrust::raw_pointer_cast(&d_sw2[0]), neall, 1.,
	      (Real*) thrust::raw_pointer_cast(&Hpsi[0]),cd->n2ft3d);
      } /*if nprj>0*/
    } /*ii*/
  }
  inline void v_local(thrust::host_vector<Real> &vout) {
    thrust::host_vector<Real> h_exi(cd->n2ft3d);

    for (int ii=0; ii<(h_pseudoData.myIonPtr->nion); ++ii) {
      int ia = h_pseudoData.myIonPtr->katm[ii];
      h_strfac.strfac_pack(ii, (Real*) &h_exi[0]);
      
      //tcc_MulSum2(1,vl[ia],exi,vout);
      for(int i=0; i < cd->n2ft3d; i++)
	vout[i]   += h_exi[i] * h_pseudoData.vl[ia][i>>1];
    }
  }
  
  /**********************************************
   *                                            *
   *                v_exc                       *              *
   *                                            *
   **********************************************/
  
  /* computes the vosko lda xc energy density and potential
   */
  
  /*---- parameters given by vosko et al -----------------*/
#define ap   3.109070e-02
#define af   1.554530e-02
#define x0p -1.049800e-01
#define x0f -3.250000e-01
#define bp   3.727440e+00
#define bf   7.060420e+00
#define cp   1.293520e+01
#define cf   1.805780e+01
  /*------------------------------------------------------*/
  
  /*     constants calculated from vosko's parameters     */
#define xp   -4.581653e-01
#define xf   -5.772521e-01
#define qp    6.151991e+00
#define qf    4.730927e+00
#define xx0p  1.255491e+01
#define xx0f  1.586879e+01
#define cp1   3.109070e-02
#define cf1   1.554530e-02
#define cp2   9.690228e-04
#define cf2   2.247860e-03
#define cp3   1.049800e-01
#define cf3   3.250000e-01
#define cp4   3.878329e-02
#define cf4   3.878329e-02
#define cp5   3.075995e+00
#define cf5   2.365463e+00
#define cp6   1.863720e+00
#define cf6   3.530210e+00
#define dp1   6.218140e-02
#define df1   3.109060e-02
#define dp2   1.938045e-03
#define df2   4.495720e-03
#define dp3   1.049800e-01
#define df3   3.250000e-01
#define dp4  -3.205972e-02
#define df4  -1.779316e-02
#define dp5  -1.192972e-01
#define df5  -1.241661e-01
#define dp6   1.863720e+00
#define df6   3.530210e+00
#define dp7   9.461748e+00
#define df7   5.595417e+00
#define fc    1.923661e+00
#define fd    2.564881e+00
#define crs   7.876233e-01
  
  /* other constants */
#define one3rd (1.00/3.00)
#define for3rd (4.00/3.00)
#define one6th (1.00/6.00)
#define dncut  1.0e-20
  
  void v_exc(const int ispin, const int n2ft3d, Real *dn, 
	     Real *xcp, Real *xce) {
    thrust::host_vector<Real> x(n2ft3d);
    
    /* local variables */
    Real xx,xx1,rho,zup,zdw,f,df;
    Real *rhoup,*xcpup,*xceup;
    Real *rhodn,*xcpdn,*xcedn;
    
    //pi      = 4.00*atan(1.00);
    
    /* define arrays and such */
    rhoup = dn;   rhodn = &dn[ (ispin-1)*n2ft3d];
    xcpup = xcp;  xcpdn = &xcp[(ispin-1)*n2ft3d];
    xceup = xce;  xcedn = &xce[(ispin-1)*n2ft3d];
    
    /* square root of wigner radius */
    for (int k=0; k<n2ft3d; ++k) {
      rho=rhoup[k]+rhodn[k]+dncut;
      x[k]=crs/powf(rho,one6th);
    }
    
    
    /* paramagnetic correlation energy & potential */
    for (int k=0; k<n2ft3d; ++k) {
      xx=1.00/(x[k]*(x[k]+bp)+cp);
      xx1 = x[k] + cp3;
      xx1 *= xx1;
      xceup[k]= cp1*log(xx*x[k]*x[k])
	+ cp2*log(xx*xx1)
	+ cp4*atan(cp5/(x[k]+cp6));
      
      xx1 = x[k]+dp6;
      xx1 *= xx1;
      xcpup[k]=xceup[k]
	-one6th*x[k]*(
		      dp1/x[k]+dp2/(x[k]+dp3)
		      +dp4*xx*(2.00*x[k]+bp)
		      +dp5/(xx1+dp7) );
    }
    
    /* paramagnetic exchange energy & potential */
    for (int k=0; k<n2ft3d; ++k) {
      xceup[k]=xceup[k]+(xp/(x[k]*x[k]));
      xcpup[k]=xcpup[k]+for3rd*(xp/(x[k]*x[k]));
    }
    
    /* do unrestricted part */
    if (ispin==2) {
      /* ferromagnetic correlation energy & potential */
      for (int k=0; k<n2ft3d; ++k) {
	xx=1.00/(x[k]*(x[k]+bf)+cf);
	xx1 = x[k]+cf6;
	xx1 *= xx1;
	xcedn[k]=cf1*log(xx*x[k]*x[k])
	  +cf2*log(xx*xx1)
	  +cf4*atan(cf5/(x[k]+cf6));
	
	xx1 = x[k]+df6;
	xx1 *= xx1;
	xcpdn[k]=xcedn[k]
	  -one6th*x[k]*(
			df1/x[k]+df2/(x[k]+df3)
			+df4*xx*(2.00*x[k]+bf)
			+df5/(xx1+df7) );
      }
      
      /* ferromagnetic exchange-energy & potential */
      for (int k=0; k<n2ft3d; ++k) {
	xcedn[k]=xcedn[k]+(xf/(x[k]*x[k]));
	xcpdn[k]=xcpdn[k]+for3rd*(xf/(x[k]*x[k]));
      }
      
      /* spin polarized exchange-correlation potential */
      for (int k=0; k<n2ft3d; ++k) {
	rho=rhoup[k]+rhodn[k]+dncut;
	zup=2.00*rhoup[k]/rho;
	zdw=2.00*rhodn[k]/rho;
	f=(zup*(powf(zup,one3rd))+zdw*(powf(zdw,one3rd))-2.00)*fc;
	xcpup[k]=(1.00-f)*xcpup[k]+f*xcpdn[k];
	df=(powf(zup,one3rd)-powf(zdw,one3rd))*(xceup[k]-xcedn[k])*fd;
	xcpdn[k]=xcpup[k]-zup*df;
	xcpup[k]=xcpup[k]+zdw*df;
	xceup[k]=xceup[k]+f*(xcedn[k]-xceup[k]);
      }
    }
  }
  
};

#endif
