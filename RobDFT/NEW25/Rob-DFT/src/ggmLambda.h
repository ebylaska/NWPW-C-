void fmf_Multiply(const int mb, thrust::device_vector<Real> &d_hml, Real alpha, Real beta)
{
   int ms,ms1,ms2,shift1,mshift1,ishift2;
   //int ng  = 2*npack(1); int ng0 = 2*nzero(1);
   int ng  = cd->n2ft3d; //int ng0 = 2*1;

   if (mb == -1) {
     ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0); 
     shift1  = 0;
     mshift1 = 0;
   } else {
     ms1=mb; ms2=mb+1; ishift2=0; 
     shift1  = mb*psiHeader.get_ne(0)*ng;
     mshift1 = 0;
   }
   for (ms=ms1; ms < ms2; ++ms) {
     int n = psiHeader.get_ne(ms);
     //dgemm_("N","N",&ng,&n,&n, &alpha, &psi1[shift1],&ng, &hml[mshift1],&n, &beta, &psi2[shift1],&ng);
     blasGemm('N', 'N', ng, n, n, alpha,
	      (Real*) thrust::raw_pointer_cast(&(*psi1)[shift1/2]), ng,
	      (Real*) thrust::raw_pointer_cast(&d_hml[mshift1]), n, beta,
	      (Real*) thrust::raw_pointer_cast(&(*psi2)[shift1/2]),ng);
     shift1  += psiHeader.get_ne(0)*ng;
     mshift1 += ishift2;
   }
}

inline int myIdamax(int n, Real* x)
{
  int ret=0;
  Real m=x[0];
  for(int i=1; i < n; i++) 
    if(m < x[i]) { m = x[i]; ret=i; }
  cerr << "Max is " << m << " index " << ret << endl;
  return(ret+1); // fortran notation
}
inline Real m_dmax(const int mb, Real *hml)
{
  int nn;
  if (mb==-1)
     nn = psiHeader.get_ne(0)*psiHeader.get_ne(0) + psiHeader.get_ne(1)*psiHeader.get_ne(1);
  else 
     nn = psiHeader.get_ne(mb)*psiHeader.get_ne(mb);

  //return( fabs(hml[idamax_(&nn,hml,&one)-1]) );
  return( fabs(hml[myIdamax(nn,hml)-1]) );
}

void m_scale_s22(const int mb, const Real dte, Real *s22)
{
   int j,k,ms,ms1,ms2,ishift2,indx0,indx,indxt;

   if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
   else {   ms1=mb; ms2=mb+1; ishift2=0;}

   for (ms=ms1; ms<ms2; ++ms) {
     int ne_ms = psiHeader.get_ne(ms);
     indx0 = ms*ishift2;
     for (k=0; k < ne_ms; ++k) {
       s22[indx0] = (1.0-s22[indx0])*(0.5/dte);
       indx  = indx0 + 1;
       indxt = indx0 + ne_ms;
       for (j=(k+1); j<ne_ms; ++j) {
	 s22[indx]  *= (-0.5/dte);
	 s22[indxt] *= (-0.5/dte);
	 indx  += 1;
	 indxt += ne_ms;
       }
       indx0 += (ne_ms+1);
     }
   }
}

void m_scale_s21(const int mb, const Real dte, Real *s21)
{
   int ms1,ms2,ishift2,indx0,indx,indxt;
   if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
   else {   ms1=mb; ms2=mb+1; ishift2=0;}

   for (int ms=ms1; ms<ms2; ++ms) {
     int ne_ms = psiHeader.get_ne(ms);
     indx0 = ms*ishift2;
     for (int k=0; k < ne_ms; ++k) {
       s21[indx0] = (1.0-s21[indx0])*(0.5);
       indx  = indx0 + 1;
       indxt = indx0 + ne_ms;
       for (int j=k+1; j < ne_ms; ++j) {
	 s21[indx]  *= -0.5;
	 s21[indxt] *= -0.5;
	 indx  += 1;
	 indxt += ne_ms;
       }
       indx0 += (ne_ms+1);
     }
   }
}

void mmm_Multiply(const int mb,
		  thrust::host_vector<Real> &a,
		  thrust::host_vector<Real> &b,
		  Real alpha,
		  thrust::host_vector<Real> &c,
		  Real beta)
{
  int ms1,ms2,ishift2,shift2;
  thrust::device_vector<Real> d_a = a;
  thrust::device_vector<Real> d_b = b;
  thrust::device_vector<Real> d_c = c;


  if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
  else { ms1=mb; ms2=mb+1; ishift2=0;}
  for (int ms=ms1; ms<ms2; ++ms) {
    int n = psiHeader.get_ne(ms);
    if (n>0) {
      shift2 = ms*ishift2;
      //dgemm_("N","N",&n,&n,&n, &alpha, &a[shift2], &n, &b[shift2], &n, &beta, &c[shift2], &n);
      blasGemm('N', 'N', n, n, n, alpha,
	       (Real*) thrust::raw_pointer_cast(&d_a[shift2]), n,
	       (Real*) thrust::raw_pointer_cast(&d_b[shift2]), n, beta,
	       (Real*) thrust::raw_pointer_cast(&d_c[shift2]), n);
    }
  }
  c = d_c;
}


void m_scale_s11(const int mb, const Real dte, Real *s11)
{
  int ms1,ms2,ishift2,indx0,indx,indxt;
  if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
  else {   ms1=mb; ms2=mb+1; ishift2=0;}
  
  for (int ms=ms1; ms<ms2; ++ms) {
    int ne_ms = psiHeader.get_ne(ms);
    indx0 = ms*ishift2;
    for (int k=0; k< ne_ms; ++k) {
      s11[indx0] *= -0.5*dte;
      indx  = indx0 + 1;
      indxt = indx0 + ne_ms;
      for (int j=k+1; j < ne_ms; ++j) {
	s11[indx]  *= -0.5*dte;
	s11[indxt] *= -0.5*dte;
	indx  += 1;
	indxt += ne_ms;
      }
      indx0 += (ne_ms+1);
    }
  }
}

#define ITERLMD         30
#define CONVGLMD        1e-15

int m_size(const int mb) {
  int nsize = (mb == -1)? 
    psiHeader.get_ne(0)*psiHeader.get_ne(0) + psiHeader.get_ne(1)*psiHeader.get_ne(1)
    : psiHeader.get_ne(mb)*psiHeader.get_ne(mb);
  return nsize;
}

inline void myDcopy(const int n, Real* a, Real* b)
{
  for(int i=0; i < n; i++) b[i] = a[i]; 
}
inline void blasGemm(char transa, char transb, int m , int n, int k, 
		     float alpha, const float *A, int lda, const float *B,
		     int ldb, float beta, float *C, int ldc) {
    cublasSgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
inline void blasGemm(char transa, char transb, int m , int n, int k, 
		     double alpha, const double *A, int lda, const double *B,
		     int ldb, double beta, double *C, int ldc) {
    cublasDgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    switch(cublasGetError()) {
    case CUBLAS_STATUS_NOT_INITIALIZED: cerr << "CUBLAS NOT Initialized!" << endl; exit(1);
    case CUBLAS_STATUS_INVALID_VALUE: cerr << "CUBLAS Invalid Value!" << endl; exit(1);
    case CUBLAS_STATUS_ARCH_MISMATCH: cerr << "CUBLAS Arch Mismatch!" << endl; exit(1);
    case CUBLAS_STATUS_EXECUTION_FAILED: cerr << "CUBLAS Execution Failed!" << endl; exit(1);
    } 
}

#ifdef CHECK_DGEMM
#include "RMF_Dgemm.h"
#endif
void ffm_sym_Multiply(const int mb, thrust::device_vector<Real> &d_hml)
{
  int ms,ms1,ms2,j,k,n,shift0,shift1,mshift0,mshift1;
  //int ishift2,nn;

   //int one = 1;
   //int ng  = 2*npack(1);
   //int ng0 = 2*nzero(1);
   //DEBUG*****************
   //int ng = 2 * 28874; int ng0= 2 * 1;
   int ng = cd->n2ft3d; int ng0= 2 * 1;

   if (mb==-1) {
     ms1=0; 
     // ms2=ispin;
     ms2 = psiHeader.get_ispin();
     //ishift2=ne[0]*ne[0]; 
     //nn=ne[0]*ne[0]+ne[1]*ne[1];
     //ishift2=psiHeader.get_ne(0) * psiHeader.get_ne(0);
     //nn = psiHeader.get_ne(0) * psiHeader.get_ne(0)
       //+ psiHeader.get_ne(1) * psiHeader.get_ne(1);
     shift0  = 0;
     mshift0 = 0;
   } else {
     ms1=mb; ms2=mb+1; //ishift2=0;
     //nn = ne[mb]*ne[mb];
     //nn = psiHeader.get_ne(mb) * psiHeader.get_ne(mb);
     //shift0  = mb*ne[0]*ng;
     shift0  = mb*psiHeader.get_ne(0)*ng;
     mshift0 = 0;
   }

   for (ms=ms1; ms<ms2; ++ms) {
     n       = psiHeader.get_ne(ms);
     shift1  = shift0;
     mshift1 = mshift0;
     for (k=1; k<=n; ++k) {
       //dgemm_("T","N",&k,&one,&ng, &rtwo,
	      //&psi1[shift0],&ng, &psi2[shift1],&ng, &rzero, &hml[mshift1],&k);
       //cerr << "DEBUG*** shift1 " << mshift1 << " d_hml.size() " << d_hml.size() << " k " << k << endl;
       blasGemm('T', 'N', k, 1, ng, (Real) 2.,
		(Real*) thrust::raw_pointer_cast(&(*psi1)[shift0/2]), ng,
		(Real*) thrust::raw_pointer_cast(&(*psi2)[shift1/2]), ng, (Real) 0.,
		(Real*) thrust::raw_pointer_cast(&d_hml[mshift1]), k);
#ifdef CHECK_DGEMM
       {
	 int one =1;
	 Real rzero = 0.0;
	 Real rtwo  = 2.0;
	 Real rone =  1.0;
	 Real rmone = -1.0;
	 thrust::host_vector<Real> h_hml = d_hml;
	 thrust::host_vector<Real> t_hml(h_hml.size());
	 thrust::host_vector<Real2> h_psi1 = *psi1;
	 thrust::host_vector<Real2> h_psi2 = *psi2;
	 Real* pt_psi1 = (Real*) &h_psi1[0].x;
	 Real* pt_psi2 = (Real*) &h_psi2[0].x;
	 Rmf_Dgemm("T","N",&k,&one,&ng, &rtwo,
		&pt_psi1[shift0],&ng, &pt_psi2[shift1],&ng, &rzero, &t_hml[mshift1],&k);
	 for(int i=0; i < k; i++) cerr << " hml[" << mshift1 << "+" << i << "] (" 
				       << h_hml[mshift1+i] << "," << t_hml[mshift1+i] << ")" << endl;
	 exit(1);
       }
#endif
       
       //dgemm_("T","N",&k,&one,&ng0, &rmone,
	      //&psi1[shift0],&ng, &psi2[shift1],&ng, &rone, &hml[mshift1],&k);
       blasGemm('T', 'N', k, 1, ng0, (Real) -1.,
		(Real*) thrust::raw_pointer_cast(&(*psi1)[shift0/2]), ng,
		(Real*) thrust::raw_pointer_cast(&(*psi2)[shift1/2]), ng, (Real) 1.,
		(Real*) thrust::raw_pointer_cast(&d_hml[mshift1]), k);
       shift1  += ng;
       mshift1 += n;
     }
     thrust::host_vector<Real> h_hml=d_hml;
     for (k=0; k<n; ++k)
       for (j=k+1; j<n; ++j)
	 h_hml[mshift0 + j + k*n] = h_hml[mshift0 + k + j*n];
     d_hml=h_hml;
     
     shift0  += ng*psiHeader.get_ne(0);
     mshift0 += psiHeader.get_ne(0)*psiHeader.get_ne(0);
   }
}

// stuff for calc_psi_r
struct calc_psi_r_functor {
  const Real2* psi;
  const Real* summer;
  const int nfft3d;
  const int neall;
  Real2* psi_r;
    
calc_psi_r_functor(const Real2* _psi, const Real* _summer,
		   const int _nfft3d, const int _neall, Real2* _psi_r)
: psi(_psi), summer(_summer), nfft3d(_nfft3d), neall(_neall), psi_r(_psi_r) {};
  
  __device__
  void operator()(unsigned int tid)
  {
    register Real reg_summer = summer[tid];
    
    for(int i=0; i < neall; i++) {
      register Real2 reg_psi = psi[tid + i*nfft3d];
      register Real2 reg_psi_r;
      reg_psi_r.x = reg_psi.x * reg_summer;
      reg_psi_r.y = reg_psi.y * reg_summer;
      psi_r[tid + i*nfft3d] = reg_psi_r;
    }
  }
};

void calc_psi_r( thrust::device_vector<Real2> &psi,
		 thrust::device_vector<Real> &summer,
		 thrust::device_vector<Real2> &psi_r)
{
  calc_psi_r_functor my_op(thrust::raw_pointer_cast(&psi[0]),
			   thrust::raw_pointer_cast(&summer[0]),
			   cd->nfft3d,
			   neall,
			   thrust::raw_pointer_cast(&psi_r[0]) );
  thrust::for_each(thrust::counting_iterator<unsigned int>(0),
		   thrust::counting_iterator<unsigned int>(cd->nfft3d),
		   my_op);
}

void ggm_lambda()
{
  Real dte = cd->dte;
  
  thrust::host_vector<Real> h_lambda(7*psiHeader.get_ne(0)*psiHeader.get_ne(0));//TODO: is this correct?
  thrust::host_vector<Real> h_s22(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s21(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s12(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s11(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_sa1(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_sa0(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_st1(psiHeader.get_ne(0)*psiHeader.get_ne(0));

  thrust::device_vector<Real> d_s22(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::device_vector<Real> d_s21(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::device_vector<Real> d_s11(psiHeader.get_ne(0)*psiHeader.get_ne(0));

  thrust::device_vector<Real2> d_psi_r(psi1->size());
  
  for (int ms=0; ms < psiHeader.get_ispin(); ++ms) {
    Real adiff;
    int nn = m_size(ms);
    
    calc_psi_r(*psi2, cd->d_summer[1], d_psi_r);

    // create S22 
    int offset = ms * psiHeader.get_ne(ms) * cd->nfft3d;
    blasGemm('T', 'N',  psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi2)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s22[0]), 
	     psiHeader.get_ne(ms));
    h_s22 = d_s22;

    // create S21 
    blasGemm('T', 'N', psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi1)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s21[0]), 
	     psiHeader.get_ne(ms));
    h_s21 = d_s21;

    calc_psi_r(*psi1, cd->d_summer[1], d_psi_r);

    // create S11 
    blasGemm('T', 'N', psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi1)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s11[0]), 
	     psiHeader.get_ne(ms));
    h_s11 = d_s11;
    //normCheck();
    //for(int indx=0; indx < 6*6; indx++) cout << h_s11[indx] << " ";
    //for(int indx=0; indx < 6*6; indx += 7) cout << h_s11[indx] << " ";
    //exit(1);

    m_scale_s22(ms,dte,&h_s22[0]);
    m_scale_s21(ms,dte,&h_s21[0]);
    //dcopy_(&nn,s21,&one,s12,&one);
    h_s12 = h_s21;
    m_scale_s11(ms,dte,&h_s11[0]);

    int ii   = 0;
    int done = 0;
    //dcopy_(&nn,s22,&one,sa0,&one);
    h_sa0 = h_s22;
    while ((!done) && ((ii++)<ITERLMD)) {
      //dcopy_(&nn,s22,&one,sa1,&one);
      h_sa1 = h_s22;

      mmm_Multiply(ms,h_s21,h_sa0,1.0,h_sa1,1.0);
      mmm_Multiply(ms,h_sa0,h_s12,1.0,h_sa1,1.0);
      mmm_Multiply(ms,h_s11,h_sa0,1.0,h_st1,0.0);
      mmm_Multiply(ms,h_sa0,h_st1,1.0,h_sa1,1.0);
      //dcopy_(&nn,sa1,&one,st1,&one);
      h_st1 = h_sa1;
      //daxpy_(&nn,&rmone,sa0,&one,st1,&one);
      myAxpy(nn, -1., &h_sa0[0], &h_st1[0]);
      
      adiff = m_dmax(ms,&h_st1[0]);
      if (adiff<CONVGLMD) done = 1;
      else { 
	//dcopy_(&nn,sa1,&one,sa0,&one);
	h_sa0 = h_sa1;
      }
    }
    if (!done) printf("ierr=10 adiff=%lf\n",adiff);
    //dcopy_(&nn,sa1,&one,&h_lambda[ms*psiHeader.get_ne(0)*psiHeader.get_ne(0)],&one);
    myDcopy(nn,&h_sa1[0],&h_lambda[ms*psiHeader.get_ne(0)*psiHeader.get_ne(0)]);
  }
  
  /* correction due to contraint */
  thrust::device_vector<Real> d_lambda = h_lambda;
  fmf_Multiply(-1,d_lambda,dte,1.0);

    calc_psi_r(*psi2, cd->d_summer[1], d_psi_r);

    // create S22 
    int ms =0;
    int offset = ms * psiHeader.get_ne(ms) * cd->nfft3d;
    blasGemm('T', 'N',  psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi2)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s22[0]), 
	     psiHeader.get_ne(ms));
    h_s22 = d_s22;
    cout << "S22 = " << endl;
    for (int i=0; i<5; ++i) {
       for (int j=0; j<5; ++j) cout <<  h_s22[i+6*j] << "  ";
       cout << endl;
    }
    cout << endl;
    
}

