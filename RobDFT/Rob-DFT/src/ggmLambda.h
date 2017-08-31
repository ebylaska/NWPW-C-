#ifdef NOT_USED
void fmf_Multiply(const int mb, Real *psi1, Real *hml, Real alpha, Real *psi2, Real beta)
{
   int ms,ms1,ms2,shift1,mshift1,ishift2;
   int ng  = 2*npack(1);
   int ng0 = 2*nzero(1);

   if (mb==-1) {
     ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0); 
     shift1  = 0;
     mshift1 = 0;
   } else {
     ms1=mb; ms2=mb+1; ishift2=0; 
     shift1  = mb*psiHeader.get_ne(0)*ng;
     mshift1 = 0;
   }
   for (ms=ms1; ms<ms2; ++ms) {
     int n = psiHeader.get_ne(ms);
     dgemm_("N","N",&ng,&n,&n, &alpha, &psi1[shift1],&ng, &hml[mshift1],&n, &beta, &psi2[shift1],&ng);
     shift1  += psiHeader.get_ne(0)*ng;
     mshift1 += ishift2;
   }
}
#endif

inline int myIdamax(int n, Real* x)
{
  int ret=0;
  Real m=x[0];
  for(int i=1; i < n; i++) 
    if(m < x[i]) { m = x[i]; ret=i; }
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
  thrust::device_vector<Real> d_c(c.size());

  if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
  else { ms1=mb; ms2=mb+1; ishift2=0;}
  for (int ms=ms1; ms<ms2; ++ms) {
    int n = psiHeader.get_ne(ms);
    if (n>0) {
      shift2 = ms*ishift2;
      //dgemm_("N","N",&n,&n,&n, &alpha, &a[shift2], &n, &b[shift2], &n, &beta, &c[shift2], &n);
      float* A = (float*) thrust::raw_pointer_cast(&d_a[shift2]);
      float* B = (float*) thrust::raw_pointer_cast(&d_b[shift2]);
      float* C = (float*) thrust::raw_pointer_cast(&d_c[shift2]);
      cublasSgemm('N', 'N', n, n, n, alpha, A, n, B, n, beta, C, n);
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

void ggm_lambda()
{
  Real dte = cd->dte;
  int n2ft3d = 2*cd->nfft3d;
  float *A, *B, *C;
  
  thrust::host_vector<Real> h_lambda(7*psiHeader.get_ne(0)*psiHeader.get_ne(0));//TODO: is this correct?
  thrust::host_vector<Real> h_s22(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s21(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s12(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_s11(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_sa1(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_sa0(psiHeader.get_ne(0)*psiHeader.get_ne(0));
  thrust::host_vector<Real> h_st1(psiHeader.get_ne(0)*psiHeader.get_ne(0));

  thrust::host_vector<Real2> h_psi_r(psi1->size());
  thrust::host_vector<Real2> h_psi1 = *psi1;
  thrust::host_vector<Real2> h_psi2 = *psi2;
  
  for (int ms=0; ms < psiHeader.get_ispin(); ++ms) {
    Real adiff;
    int nn = m_size(ms);

    // nWavefunc is equal to (psiHeader.get_ne(0) + psiHeader.get_ne(1)
    for(int i=0; i < nWavefunc*cd->nfft3d; i++) {
      h_psi_r[i].x = h_psi2[i].x * cd->h_summer[1][i%cd->nfft3d];
      h_psi_r[i].y = h_psi2[i].y * cd->h_summer[1][i%cd->nfft3d];
    }

    // create S22 
    thrust::device_vector<Real2> psi_r = h_psi_r;
    thrust::device_vector<Real> s22(psiHeader.get_ne(0)*psiHeader.get_ne(0));
    A = (float*) thrust::raw_pointer_cast(&psi_r[0]);
    B = (float*) thrust::raw_pointer_cast(&psi2[0]);
    C = (float*) thrust::raw_pointer_cast(&s22[0]);
    cublasSgemm('T', 'N', n2ft3d, psiHeader.get_ne(ms), psiHeader.get_ne(ms), 1., A, n2ft3d,
		B, n2ft3d, 0., C, psiHeader.get_ne(0));
    h_s22 = s22;
    
    // create S21 
    //thrust::device_vector<Real2> psi_r = h_psi_r;
    thrust::device_vector<Real> s21(psiHeader.get_ne(0)*psiHeader.get_ne(0));
    A = (float*) thrust::raw_pointer_cast(&psi_r[0]);
    B = (float*) thrust::raw_pointer_cast(&psi1[0]);
    C = (float*) thrust::raw_pointer_cast(&s21[0]);
    cublasSgemm('T', 'N', n2ft3d, psiHeader.get_ne(ms), psiHeader.get_ne(ms), 1., A, n2ft3d,
		B, n2ft3d, 0., C, psiHeader.get_ne(0));
    h_s21 = s21;
    
    for(int i=0; i < nWavefunc*cd->nfft3d; i++) {
      h_psi_r[i].x = h_psi1[i].x * cd->h_summer[1][i%cd->nfft3d];
      h_psi_r[i].y = h_psi1[i].y * cd->h_summer[1][i%cd->nfft3d];
    }

    // create S11 
    //thrust::device_vector<Real2> psi_r = h_psi_r;
    thrust::device_vector<Real> s11(psiHeader.get_ne(0)*psiHeader.get_ne(0));
    A = (float*) thrust::raw_pointer_cast(&psi_r[0]);
    B = (float*) thrust::raw_pointer_cast(&psi1[0]);
    C = (float*) thrust::raw_pointer_cast(&s11[0]);
    cublasSgemm('T', 'N', n2ft3d, psiHeader.get_ne(ms), psiHeader.get_ne(ms), 1., A, n2ft3d,
		B, n2ft3d, 0., C, psiHeader.get_ne(0));
    h_s11 = s11;

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
      mySaxpy(nn, -1., &h_sa0[0], &h_st1[0]);
      
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
  
#ifdef FOO
  /* correction due to contraint */
  //fmf_Multiply(-1,psi1,lambda,dte,psi2,1.0);
  int shift1=0;
  int mshift1=0;
  for (int ms=ms1; ms < ms2; ++ms) {
    int n = psiHeader.get_ne(ms);
    A = (float*) thrust::raw_pointer_cast(&psi1[shift1]);
    B = (float*) thrust::raw_pointer_cast(&h_lambda[mshift1]);
    C = (float*) thrust::raw_pointer_cast(&psi2[shift1]);
    cublasSgemm('N', 'N', n2ft3d, psiHeader.get_ne(ms), cd->nfft3d, 1., A, n2ft3d,
	      B, n2ft3d, 0., C, psiHeader.get_ne(ms));
    //dgemm_("N","N",n2fft3,n,n, alpha, psi1[shift1],n2fft3d, &hml[mshift1],n,
	    //beta, psi2[shift1],n2ft3d);
     shift1  += psiHeader.get_ne(0)*ng;
     mshift1 += psiHeader.get_ne(0)*psiHeader.get_ne(0);
   }
#endif
}

