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
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
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
  thrust::device_vector<Real2> dng;
  thrust::device_vector<Real2> *psi1, *psi2;
  Strfac<Real> h_strfac;
  PseudoPotentialData<Real> h_pseudoData;

  // housekeeping variables
  PsiHeader<Real> psiHeader;
  CommonData<Real> *cd;
  //PseudoPotentialData<Real> *pspData;
  int nWavefunc;

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

  kspace_op(
	    const Real2* _psi,
	    const Real* _tg,
	    const int _tg_size
	    ): psi(_psi), tg(_tg), tg_size(_tg_size) {};
    __host__ __device__
    Real2 operator()(const Real2& psi1) const { 
      int j = &psi1 - psi;
      double tg_j = tg[j%tg_size]; // Note older GPUs will serialize on tg
      register Real2 ret;
      ret.x = ( psi1.x * tg_j );
      ret.y = ( psi1.y * tg_j );
      return ret;
    }
  };

  inline void kspace_operator()
  {
    kspace_op unary_op( thrust::raw_pointer_cast(&((*psi1)[0])),
			thrust::raw_pointer_cast(&((cd->d_tg)[0])),
			cd->nfft3d);
    thrust::transform(psi1->begin(), psi1->end(), Hpsi.begin(), unary_op);
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
    
    nWavefunc = psiHeader.get_ne(0) + psiHeader.get_ne(1);
    
    // allocate and read wavefunctions 
    // bummer! the nx+2 prevents power of two optimizations
    int size = nWavefunc * (cd->nx+2) * cd->ny * cd->nz;
    psiVect[0] = thrust::device_vector<Real2>(size);
    psiVect[1] = thrust::device_vector<Real2>(size);
    psi1 = &psiVect[0]; psi2 = &psiVect[1];
    
    Hpsi = thrust::device_vector<Real2>(size);
    
    {
      thrust::host_vector<Real2> tmp(size);
      int index=0;
      for (int i=0; i < nWavefunc*(cd->nfft3d); i++) {
	tmp[i].x = (Real) ((double*)buf)[index++];
	tmp[i].y = (Real) ((double*)buf)[index++];
      }
      (*psi1) = tmp; // move to GPU
    }
    buf += nWavefunc*2*cd->nfft3d*sizeof(double);
    switchPsi();
    
    munmap(mmapBuf,statbuf.st_size);
    close(fd);
    
    myIon = new Ion(myrtdb);
    h_strfac = Strfac<Real>(myIon, cd);
    
    // read in Pseudopotential data
    h_pseudoData = PseudoPotentialData<Real>(rootFilename, myIon, cd->nfft3d);
    cerr << "npsp " << h_pseudoData.npsp << endl;
    cerr << "vnl[0][0] " << h_pseudoData.vnl[0][0] << endl;
    
    cublasInit();
    initFFTs();
    int nsize = psiHeader.get_ispin() * cd->nfft3d;
    dn = thrust::device_vector<Real>( 2*cd->nfft3d );
    psi_r = thrust::device_vector<Real>( nWavefunc*2*cd->nfft3d );
    vall = thrust::device_vector<Real>( 2*cd->nfft3d );
  } 
  
  // stuff for the Ke operator
  struct Ke_unary_op : public thrust::unary_function<unsigned int,Real>
    {
      const Real2* psi;
      const Real* tg;
      const Real* summer;
      const int nfft3d;
      
    Ke_unary_op(const Real2* _psi, const Real* _tg,
		const Real* _summer, const int _nfft3d )
      : psi(_psi), tg(_tg), summer(_summer), nfft3d(_nfft3d) {};
      
      __host__ __device__
	Real operator()(unsigned int tid)
      {
	register Real2 p = psi[tid];
	register int i = (tid)% nfft3d;
	return ((p.x * p.x + p.y * p.y) * tg[i] * summer[i]);
      }
    };

  inline Real kineticEnergy()
  {
    //#define KE_HOST_SIDE
#ifdef KE_HOST_SIDE
    Real sum=0.0;
    thrust::host_vector<Real2> tmp(*psi1);
    
    for(int i=0; i < nWavefunc; i++) {
      Real2 *wavef = &tmp[i*cd->nfft3d];
      for(int j=0; j < cd->nfft3d; j++) {
	Real tmp_tg = cd->h_tg[j%cd->nfft3d];
        sum += (  wavef[j].x * wavef[j].x + wavef[j].y * wavef[j].y)
         * cd->h_tg[j%cd->nfft3d]
         * cd->h_summer[1][j%cd->nfft3d];
      }
    }
#else
    int n = nWavefunc * cd->nfft3d;
    Ke_unary_op unary_op(thrust::raw_pointer_cast(&(*psi1)[0]),
			 thrust::raw_pointer_cast(&cd->d_tg[0]),
			 thrust::raw_pointer_cast(&cd->d_summer[1][0]),
			 cd->nfft3d);
    Real sum = thrust::transform_reduce(thrust::counting_iterator<int>(0),
					thrust::counting_iterator<int>(n),
					unary_op,
					(Real) 0.,
					thrust::plus<Real>());
#endif
    if (psiHeader.get_ispin()==1) sum *= 2.0;
    sum = -sum;
    return sum;
  }

  Real normCheck()
  {
    Real sum=0.0;
    thrust::host_vector<Real2> tmp(*psi1);

    for(int i=0; i < nWavefunc; i++) {
      Real2 *wavef = &tmp[i*cd->nfft3d];
      for(int j=0; j < cd->nfft3d; j++){
	sum += (  wavef[j].x*wavef[j].x + wavef[j].y*wavef[j].y)
	  * cd->h_summer[1][j];
      }
      cerr << "cumulative sum " << sum << " " << i << endl;
    }
    return sum;
  }

  Real psiHpsiCheck()
  {
    Real sum=0.0;
    thrust::host_vector<Real2> tmp1(*psi1);
    thrust::host_vector<Real2> tmp2(Hpsi);

    for(int i=0; i < nWavefunc; i++) {
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

  void inner_loop()
  {
    Real Ke;
    int ispin = psiHeader.get_ispin();

    int it_in = control_loop(0);
    //it_in = 1000;

    h_strfac.phafac();

    switchPsi();
    Ke = kineticEnergy();
    cerr << "Initial Ke is " << Ke << endl;
    Ke = normCheck();
    cerr << "Initial norm is " << Ke << endl;
    switchPsi();

    cerr << "inner_loop, it_in " << it_in << endl;

    for(int iter=0; iter < it_in; iter++) {
      //cerr << nWavefunc << " " << normCheck() << endl;

      // copy psi1(G,n,ms) = psi2(G,n,ms). Done via pointer switch
      switchPsi();

      // calculates Hpsi(G,n,ms) = KE*psi1(G,n,ms) + VNL*psi1(G,n,ms)
      Real2 czero={0.,0.};
      thrust::fill(Hpsi.begin(),Hpsi.end(), czero);
      kspace_operator();


      Ke = -2.0*psiHpsiCheck();
      cerr << "Initial <psi|Hpsi> is " << Ke << endl;


      // **** add r-space operator to k-space operator ****
      // convert  psi(R,n,ms)= FFT(psi1(G,n,ms)). NOTE! done in place in psi1
      Real *d_psi1 = (Real*) thrust::raw_pointer_cast(&(*psi1)[0]);
      Real *d_psi_r = (Real*) thrust::raw_pointer_cast(&psi_r[0]);
      // Might be able to do in place
      if(cufftExecC2R(fftPlanMany_C2R, (cufftComplex*)d_psi1,(cufftReal*) d_psi_r) != CUFFT_SUCCESS) {
	cerr << "A C2R FFT failed!" << endl;
	exit(1);
      }
      
      // compute the density dn(R,ms)= sum(n=1,nelc(ms)) psi(R,n,ms)*psi(R,n,ms)
      Real scal2 = 1/cd->omega;
      thrust::host_vector<Real> h_dn(dn.size()); // TODO get rid of transfer
      thrust::host_vector<Real> h_psi_r = psi_r; // TODO get rid of transfer
      hr_aSumSqr(scal2,h_psi_r,h_dn);
      dn = h_dn; // TODO get rid of transfer

      Real sum=0.0;
      for (int i=0; i<2*cd->nfft3d; ++i)
         sum += dn[i];
      cerr << "dn=" << sum*(30.0*30.0*30.0/(48.0*48.0*48.0)) << endl;


      // compute the FFT of dn, rho(G) = fft(dn(R,1)+dn(R,ispin)) R->C
      thrust::host_vector<Real2> h_vc(psi1->size());
      {
	thrust::host_vector<Real2> h_rho = *psi2; // Note psi2 is used as a temp for rho
	int n2ft3d = 2*psi1->size();
	const Real scal1 = 1./((Real) cd->nx * (Real) cd->ny * (Real) cd->nz);

	thrust::host_vector<Real2> h_dng(cd->nfft3d);
	//TODO: note using h_dn calculated above without a transfer
	// calculate rrr_Sum and r_SMul
	for(int i=0; i < h_dng.size(); i++) {
	  h_dng[i].x = h_dn[2*i] + h_dn[(ispin-1)*cd->nfft3d + 2*i];
	  h_dng[i].y = h_dn[2*i+1] + h_dn[(ispin-1)*cd->nfft3d + 2*i+1];

	  h_dng[i].x *= scal1;
	  h_dng[i].y *= scal1;
	}

	dng = h_dng; //TODO eliminate this transfer
	Real *pt = (Real*) thrust::raw_pointer_cast(&dng[0]);
	if(cufftExecR2C(fftPlan_R2C, (cufftReal*)pt,(cufftComplex*) pt) != CUFFT_SUCCESS) {
	  cerr << "C R2C FFT failed!" << endl;
	  exit(1);
	}
	h_dng = dng; //TODO eliminate this transfer

        // compute Vcoulomb(G) = (4*pi/G^2)*.rho(G)
	//vcoulomb(dng,vcout);
	for(int i=0; i < h_dng.size(); i++) {
	  Real myTG=cd->h_tg[i%cd->nfft3d];
	  if(myTG != 0.) {
	    h_vc[i].x = ((4.*M_PI)/cd->h_tg[i%cd->nfft3d]) *  h_dng[i].x; 
	    h_vc[i].y = ((4.*M_PI)/cd->h_tg[i%cd->nfft3d]) *  h_dng[i].y; 
	  } else h_vc[i].x = h_vc[i].y = 0.;
	}
      }

      // compute Vlocal(G)
      thrust::host_vector<Real2> h_vl(cd->nfft3d);
      v_local(h_vl);

      // compute Vxc(R,ms) ~= alpha*(dn(R,1)+dn(R,ispin))^1/3
      int n2ft3d = 2*cd->nfft3d;
      thrust::host_vector<Real> h_xcp(2*ispin*cd->nfft3d);
      thrust::host_vector<Real> h_xce(2*ispin*cd->nfft3d);
      v_exc(ispin,n2ft3d,(Real*)&h_dn[0],(Real*) &h_xcp[0],(Real*) &h_xce[0]); 

      // apply rspace operators
      thrust::host_vector<Real> h_vall(2*cd->nfft3d);
      // compute vall(R) = FFT(Vcoulomb(G)+vlocal(G) C->R
      //   cc_SMul(0,scal2,vl,vall);
      for (int i=0; i< h_vl.size(); ++i) {
	h_vall[2*i]   = scal2 * h_vl[i].x;
	h_vall[2*i+1] = scal2 * h_vl[i].y;
      }
      //cc_Sum2(0,vc,vall);
      for (int i=0; i< h_vl.size(); ++i) {
	h_vall[2*i]   += h_vc[i].x;
	h_vall[2*i+1] += h_vc[i].y;
      }
      //mygrid->cr_fft3d(vall);
      vall = h_vall; //TODO eliminate this transfer
      Real *pt = (Real*) thrust::raw_pointer_cast(&vall[0]);
      if(cufftExecC2R(fftPlan_C2R, (cufftComplex*)pt,(cufftReal*) pt) != CUFFT_SUCCESS) {
	cerr << "D C2R FFT failed!" << endl;
	exit(1);
      }
      h_vall = vall; //TODO eliminate this transfer

      // Hpsi(G,n,ms) += FFT(Vall(R)+Vxc(R,ms))*psi(R,n,ms)) r->C
      {
	int indx1 = 0;
	int indx2 = 0;
	const Real scal1 = 1./((Real) cd->nx * (Real) cd->ny * (Real) cd->nz);
	thrust::host_vector<Real> h_tmp(n2ft3d);
	thrust::host_vector<Real> h_vpsi(n2ft3d);
	thrust::device_vector<Real> vpsi(n2ft3d);
	thrust::host_vector<Real2> h_Hpsi = Hpsi;
	for (int ms=0; ms < ispin; ++ms) {
	  //mygrid->rrr_Sum(vall,&xcp[ms*n2ft3d],tmp);
	  for(int i=0; i < n2ft3d; i++) {
	    h_tmp[i] = h_vall[i] + h_xcp[ms*n2ft3d+i];
	  }
	  //TODO: make this a multiplan fft
	  for (int ne=0; ne< psiHeader.get_ne(ms); ++ne) {
	    //mygrid->rrr_Mul(tmp,&psi_r[indx2],vpsi);
	    for(int i=0; i < n2ft3d; i++) {
	      h_vpsi[i] = h_psi_r[indx2+i] * h_tmp[i];
	    }
	    //mygrid->rc_fft3d(vpsi);
	    vpsi = h_vpsi; //TODO eliminate this transfer
	    Real *pt = (Real*) thrust::raw_pointer_cast(&vpsi[0]);
	    if(cufftExecR2C(fftPlan_R2C, (cufftReal*)pt,(cufftComplex*) pt) != CUFFT_SUCCESS) {
	      cerr << "E R2C FFT failed!" << endl;
	      exit(1);
	    }
	    h_vpsi = vpsi; //TODO eliminate this transfer
	    //mygrid->c_pack(1,vpsi);
	    //mygrid->cc_daxpy(1,(-scal1),vpsi,&Hpsi[indx1]);
	    for(int i=0; i < cd->nfft3d; i++) {
	      h_Hpsi[indx1+i].x += (-scal1)* h_vpsi[2*i];
	      h_Hpsi[indx1+i].y += (-scal1)* h_vpsi[2*i+1];
	    }
	    indx1 += cd->nfft3d;
	    indx2 += n2ft3d;
	  }
	}
	Hpsi = h_Hpsi; //TODO remove this transfer
      }
          
      //*** update ****

      // psi2(G,n,ms) = psi1(G,n,ms) + dte*Hpsi(G,n,ms)
      // steepest descents gg_SMul(dte,Hpsi,psi2); gg_Sum2(psi1,psi2);
      {
	thrust::host_vector<Real2> h_Hpsi = Hpsi;
	thrust::host_vector<Real2> h_psi1 = *psi1;
	thrust::host_vector<Real2> h_psi2(psi1->size());
	for(int i=0; i < h_Hpsi.size(); i++) {
	  // gg_SMul(dte,Hpsi,psi2); 
	  h_psi2[i].x = cd->dte * h_Hpsi[i].x;
	  h_psi2[i].y = cd->dte * h_Hpsi[i].y;
	  // gg_Sum2(psi1,psi2);
	  h_psi2[i].x += h_psi1[i].x;
	  h_psi2[i].y += h_psi1[i].y;
	}
	*psi2 = h_psi2;
      }
      // compute lagrange multipliers, Lambda(n,n,ms) = Function(dte,psi1,psi2)
      // psi2(G,n,ms) += Lambda(n,n',ms) * psi1(G,n',ms)
      ggm_lambda();

      /* loop finished */
    }

    Ke = kineticEnergy();
    cerr << "END Ke is " << Ke << endl;
    //enonlocal =  nonlocalEnergy()
  }

#include "ggmLambda.h"
#ifdef FOO
#endif


  void inline initFFTs()
  {
    //TODO: double check variables
    int dim[] = { cd->nx+2, cd->ny, cd->nz};
    // Create a batched 3D plan
    cufftPlanMany(&fftPlanMany_C2R, 3, dim, NULL, 1,0, NULL,1,0,CUFFT_C2R,nWavefunc);
    cufftPlan3d(&fftPlan_C2R, dim[0], dim[1], dim[2], CUFFT_C2R);
    cufftPlanMany(&fftPlanMany_R2C, 3, dim, NULL, 1,0, NULL,1,0,CUFFT_R2C,nWavefunc);
    cufftPlan3d(&fftPlan_R2C, dim[0], dim[1], dim[2], CUFFT_R2C);
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
      for (n=0; n < nWavefunc; ++n) {
	for (k=0; k < 2*cd->nfft3d; ++k) {
	  dn[indx0+k] += alpha*psi_r[indx1+k]*psi_r[indx1+k];
	}
	indx1 += 2*cd->nfft3d;
      }
      indx0 += 2*cd->nfft3d;
    }
  }

  // methods for v_nonlocal (slammed into the code without thought)
  void inline tcc_Mul(const int nb, Real *a, Real *b, Real *c)
  {
    int i,ii;
    
    ii = 0;
    for (i=0; i< cd->nfft3d; ++i) {
      c[ii]   = b[ii]  *a[i];
      c[ii+1] = b[ii+1]*a[i];
      ii += 2;
    }
  }
  void inline tcc_iMul(const int nb, Real *a, Real *b, Real *c)
  {
    int i,ii;
    
    ii = 0;
    for (i=0; i< cd->nfft3d; ++i) {
      c[ii]   = -b[ii+1]  * a[i];
      c[ii+1] = b[ii] * a[i];
      ii += 2;
    }
  }

  inline static void mySaxpy(const int n, const Real alpha, const Real* x, Real *y)
  {
    for(int i=0; i < n; i++) y[i] += alpha * x[i];
  }

  //TODO Look at doing this in one call on GPU
  static void Multiply_Gijl_sw1(int nn,
				const int nprj,
				const int nmax,
				const int lmax,
				int *n_prj,
				int *l_prj,
				int *m_prj,
				Real *G,
				Real *sw1,
				Real *sw2)
  {
    int a,b,na,nb;
    int nmax2 = nmax*nmax;
    int nnn = nn*nprj;
    int nna = nn;
    //Real rzero = 0.0;

    // fill sw2 with zeros
    if(sizeof(Real) != sizeof(float)) {
      cerr << "Single precision only at the moment" << endl;
      exit(1);
    }
    //cublasScopy(nnn,&rzero,0,sw2,1);
    for(int i=0; i < nnn; i++) sw2[i]=0.;

    for (b=0; b<nprj; ++b)
      for (a=0; a<nprj; ++a)
	if ((l_prj[a]==l_prj[b]) && (m_prj[a]==m_prj[b])) {
	  na = n_prj[a]-1;
	  nb = n_prj[b]-1;
	  //cublasSaxpy(nna, G[nb + na*nmax + nmax2 * l_prj[a]], &sw1[a*nn],1, &sw2[b*nn],1);
	  mySaxpy(nna,
		  G[nb + na*nmax + nmax2 * l_prj[a]],
		  &sw1[a*nn],
		  &sw2[b*nn]);
	}
  }
  
  void v_nonlocal()
  {
    int ii,ia,l,sd_function;
    int neall=psiHeader.get_ne(0)+psiHeader.get_ne(1);

    int nshift0 = cd->nfft3d;
    int nshift = 2*cd->nfft3d;
    thrust::host_vector<Real> h_sw1(neall*h_pseudoData.nprj_max);
    thrust::host_vector<Real> h_sw2(neall*h_pseudoData.nprj_max);
    // exi and prj should be real2
    thrust::host_vector<Real> h_prjtmp(h_pseudoData.nprj_max*nshift);
    thrust::host_vector<Real> h_exi(nshift);
    
    for (ii=0; ii<(h_pseudoData.myIonPtr->nion); ++ii) {
      ia = h_pseudoData.myIonPtr->katm[ii];
      if (h_pseudoData.nprj[ia] > 0) {
	/* structure factor */
	h_strfac.strfac_exi(ii, &h_exi[0]);
	
	/* generate sw1's and projectors */
	for (l=0; l<h_pseudoData.nprj[ia]; ++ l) {
	  sd_function = !(h_pseudoData.l_projector[ia][l] & 1);
	  Real *vnlprj = &h_pseudoData.vnl[ia][l*nshift0];
	  Real *prj =  &h_prjtmp[l*nshift];
	  // tcc
	  if (sd_function) {
	    tcc_Mul(1,vnlprj,&h_exi[0],prj);
	  } else {
	    tcc_iMul(1,vnlprj,&h_exi[0],prj);
	  }
	  //TODO: Remove this copy from device memory
	  thrust::host_vector<Real2> h_tmp = *psi1;
	  //cc_pack_indot(1,nn,psi,prj,&h_sw1[l*nn]);
          for (int n=0; n<neall; ++n) {
	     int i1=0, i2=1;
             Real sum=0.;
	     for(int i=0; i < cd->nfft3d; i++) {
	        sum += prj[i1] * cd->h_summer[1][i] * h_tmp[i+n*cd->nfft3d].x
	             + prj[i2] * cd->h_summer[1][i] * h_tmp[i+n*cd->nfft3d].y;
	        i1 += 2; i2 +=2;
	     }
	     h_sw1[l*neall+n] = sum;
          }
	}
	// don't worry about this: parall->Vector_SumAll(1,nn*nprj[ia],sw1);
	
	/* sw2 = Gijl*sw1 */
	//Multiply_Gijl_sw1(nn,h_pseudoData.nprj[ia],nmax[ia],lmax[ia],
		//	  n_projector[ia],l_projector[ia],m_projector[ia],
		//	  Gijl[ia],sw1,sw2);
	Multiply_Gijl_sw1(neall,h_pseudoData.nprj[ia],
			  h_pseudoData.nmax[ia],
			  h_pseudoData.lmax[ia],
			  &h_pseudoData.n_projector[ia][0],
			  &h_pseudoData.l_projector[ia][0],
			  &h_pseudoData.m_projector[ia][0],
			  &h_pseudoData.Gijl[ia][0],
			  &h_sw1[0],
			  &h_sw2[0]);

	/* do Kleinman-Bylander Multiplication */
	//int ntmp = nn* h_pseudoData.nprj[ia];
	//dscal_(&ntmp,&scal,sw2,&one);
	Real scal = 1.0/cd->omega;
	for(int i=0; i < h_sw2.size(); i++) h_sw2[i] *= scal;
	
        int ntmp = h_pseudoData.nprj[ia];
	
	thrust::device_vector<Real> sw2=h_sw2;

	thrust::device_vector<Real> prjtmp = h_prjtmp;
	float* A = (float*) thrust::raw_pointer_cast(&prjtmp[0]);
	float* B = (float*) thrust::raw_pointer_cast(&sw2[0]);
	float* C = (float*) thrust::raw_pointer_cast(&Hpsi[0]);
	static float rone  = 1.0;
	static float rmone = -1.0;
        cublasSgemm('N', 'T', nshift, neall, ntmp, rmone, A, nshift,
		    B, neall, rone, C, nshift);
      } /*if nprj>0*/
    } /*ii*/
  }
  inline void v_local(thrust::host_vector<Real2> &vout) {
    thrust::host_vector<Real2> h_exi(cd->nfft3d);
    for (int ii=0; ii<(h_pseudoData.myIonPtr->nion); ++ii) {
      int ia = h_pseudoData.myIonPtr->katm[ii];
      h_strfac.strfac_exi(ii, (Real*) &h_exi[0]);
      //tcc_MulSum2(1,vl[ia],exi,vout);
      for(int i=0; i < h_exi.size(); i++) {
	vout[i].x += h_exi[i].x * h_pseudoData.vl[ia][i];
	vout[i].y += h_exi[i].y * h_pseudoData.vl[ia][i];
      }
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
