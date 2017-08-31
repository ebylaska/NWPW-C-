void fmf_Multiply(const int mb, thrust::device_vector<Real> &d_hml, Real alpha, Real beta)
{
   int ms,ms1,ms2,shift1,mshift1,ishift2;
   int ng  = cd->n2ft3d; 

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
  Real m=fabs(x[0]);
  for(int i=1; i < n; i++) if(m < fabs(x[i])) { m = fabs(x[i]); ret=i; }
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


/**************************************************
 *                                                *
 *                 d_dmax                         *
 *                                                *
 **************************************************/
/*
    
*/  
  struct d_dmax_functor {
     const Real* mat1;

     d_dmax_functor(const Real* _mat1) : mat1(_mat1) {};

    __device__
    Real operator()(unsigned int tid)
    {
      register Real tmp  = mat1[tid];
      register Real sum1 = tmp*tmp;
      return(sum1);
    }
  };
inline Real d_dmax(const int mb, const thrust::device_vector<Real> &d_mat1)
{
   int nn;
   Real sum1 = 0.0;
   if (mb==-1)
      nn = psiHeader.get_ne(0)*psiHeader.get_ne(0)
         + psiHeader.get_ne(1)*psiHeader.get_ne(1);
   else
      nn = psiHeader.get_ne(mb)*psiHeader.get_ne(mb);

   d_dmax_functor my_op(thrust::raw_pointer_cast(&d_mat1[0]));
   sum1 = thrust::transform_reduce(thrust::counting_iterator<unsigned int>(0),
                                   thrust::counting_iterator<unsigned int>(nn),
                                   my_op,
                                   (Real) 0.,
                                   thrust::greater<Real>());
   return sqrt(sum1);
}





/********************************************* 
 *                                           *
 *              d_scale_s22                  *
 *                                           *
 *********************************************/
struct d_scale_s22_functor
{
     const Real  dte;
     Real*       id;
     Real*       s22;

    d_scale_s22_functor(const Real _dte, Real* _id, Real* _s22) : dte(_dte), id(_id), s22(_s22) {};
    __device__ 
    void operator()(unsigned int tid)
    {
      s22[tid] = (id[tid] - s22[tid])*(0.5/dte);
    }
};
inline void d_scale_s22(const int mb,  
                        const Real dte,
                        thrust::device_vector<Real> &s22)
{
   int msize,ishift;
   if (mb==-1)  
   {
      msize = psiHeader.get_ne(0)*psiHeader.get_ne(0)  
            + psiHeader.get_ne(1)*psiHeader.get_ne(1);
      ishift = 0;
   }
   else
   {
      msize=psiHeader.get_ne(mb)*psiHeader.get_ne(mb);
      ishift = mb*psiHeader.get_ne(0)*psiHeader.get_ne(0);
   }
  
   d_scale_s22_functor myop(dte,
                            thrust::raw_pointer_cast(&(cd->d_id[ishift])), 
                            thrust::raw_pointer_cast(&s22[0]));
   thrust::for_each(thrust::counting_iterator<unsigned int>(0),
                    thrust::counting_iterator<unsigned int>(msize),
                      myop);
}
    



/********************************************* 
 *                                           *
 *              d_scale_s21                  *
 *                                           *
 *********************************************/
struct d_scale_s21_functor
{
     Real*       id;
     Real*       s21;

    d_scale_s21_functor(Real* _id, Real* _s21) : id(_id), s21(_s21) {};
    __device__
    void operator()(unsigned int tid)
    {
      s21[tid] = (id[tid] - s21[tid])*(0.5);
    }
};
inline void d_scale_s21(const int mb,
                        const Real dte,
                        thrust::device_vector<Real> &s21)
{
   int msize,ishift;
   if (mb==-1)
   {
      msize = psiHeader.get_ne(0)*psiHeader.get_ne(0)
            + psiHeader.get_ne(1)*psiHeader.get_ne(1);
      ishift = 0;
   }
   else
   {
      msize=psiHeader.get_ne(mb)*psiHeader.get_ne(mb);
      ishift = mb*psiHeader.get_ne(0)*psiHeader.get_ne(0);
   }

   d_scale_s21_functor myop(thrust::raw_pointer_cast(&(cd->d_id[ishift])),
                            thrust::raw_pointer_cast(&s21[0]));
   thrust::for_each(thrust::counting_iterator<unsigned int>(0),
                    thrust::counting_iterator<unsigned int>(msize),
                      myop);
}



/*********************************************
 *                                           *
 *              d_scale_s11                  *
 *                                           *
 *********************************************/
struct d_scale_s11_functor
{
     const Real  dte;
     Real*       s11;

    d_scale_s11_functor(const Real _dte, Real* _s11) : dte(_dte), s11(_s11) {};
    __device__
    void operator()(unsigned int tid)
    {
      s11[tid] *= -0.5*dte;
    }
};
inline void d_scale_s11(const int mb, 
                        const Real dte,
                        thrust::device_vector<Real> &s11)
{
   int msize;
   if (mb==-1) 
      msize = psiHeader.get_ne(0)*psiHeader.get_ne(0) 
            + psiHeader.get_ne(1)*psiHeader.get_ne(1);
   else
      msize=psiHeader.get_ne(mb)*psiHeader.get_ne(mb);

   d_scale_s11_functor myop(dte,thrust::raw_pointer_cast(&s11[0]));
   thrust::for_each(thrust::counting_iterator<unsigned int>(0),
                    thrust::counting_iterator<unsigned int>(msize),
                      myop);
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



void d_mmm_Multiply(const int mb,
                    const thrust::device_vector<Real> &d_a,
                    const thrust::device_vector<Real> &d_b,
                    const Real alpha,
                    thrust::device_vector<Real>       &d_c,
                    const Real beta)
{
  int ms1,ms2,ishift2,shift2;

  if (mb==-1) {   ms1=0; ms2=psiHeader.get_ispin(); ishift2=psiHeader.get_ne(0)*psiHeader.get_ne(0);}
  else { ms1=mb; ms2=mb+1; ishift2=0;}
  for (int ms=ms1; ms<ms2; ++ms) {
    int n = psiHeader.get_ne(ms);
    if (n>0) {
      shift2 = ms*ishift2;
      blasGemm('N', 'N', n, n, n, alpha,
               (Real*) thrust::raw_pointer_cast(&d_a[shift2]), n,
               (Real*) thrust::raw_pointer_cast(&d_b[shift2]), n, beta,
               (Real*) thrust::raw_pointer_cast(&d_c[shift2]), n);
    }
  }
}



#define ITERLMD         50
#define CONVGLMD        3.0e-7

int m_size(const int mb) {
  int nsize = (mb == -1)? 
    psiHeader.get_ne(0)*psiHeader.get_ne(0) + psiHeader.get_ne(1)*psiHeader.get_ne(1)
    : psiHeader.get_ne(mb)*psiHeader.get_ne(mb);
  return nsize;
}

//TODO get rid of this
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

void ffm_sym_Multiply(const int mb, thrust::device_vector<Real> &d_hml)
{
  int ms,ms1,ms2,j,k,n,shift0,shift1,mshift0,mshift1;
   int ng = cd->n2ft3d; int ng0= 2 * 1;

   if (mb==-1) {
     ms1=0; 
     ms2 = psiHeader.get_ispin();
     shift0  = 0;
     mshift0 = 0;
   } else {
     ms1=mb; ms2=mb+1;
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
       blasGemm('T', 'N', k, 1, ng, (Real) 2.,
		(Real*) thrust::raw_pointer_cast(&(*psi1)[shift0/2]), ng,
		(Real*) thrust::raw_pointer_cast(&(*psi2)[shift1/2]), ng, (Real) 0.,
		(Real*) thrust::raw_pointer_cast(&d_hml[mshift1]), k);
       
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

/************************************************
 *                                              *
 *          ff_Multiply_Summer                  *
 *                                              *
 ************************************************/

struct ff_Multiply_Summer_functor {
    const Real2* psi;
    const Real* summer;
    const int nfft3d;
    const int neall;
    Real2* psi_r;

   ff_Multiply_Summer_functor(const Real2* _psi, const Real* _summer,
                              const int _nfft3d, const int _neall, 
                              Real2* _psi_r)
                             : psi(_psi), summer(_summer), 
                               nfft3d(_nfft3d), neall(_neall), 
                               psi_r(_psi_r) {};

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

inline void ff_Multiply_Summer(const thrust::device_vector<Real2> &psi,
                               const thrust::device_vector<Real> &summer,
                               thrust::device_vector<Real2> &psi_out)
{
   ff_Multiply_Summer_functor my_op(thrust::raw_pointer_cast(&psi[0]),
                                    thrust::raw_pointer_cast(&summer[0]),
                                    cd->nfft3d,
                                    neall,
                                    thrust::raw_pointer_cast(&psi_out[0]) );
   thrust::for_each(thrust::counting_iterator<unsigned int>(0),
                    thrust::counting_iterator<unsigned int>(cd->nfft3d),
                    my_op);
}


/************************************************
 *                                              *
 *                ffm_Multiply                  *
 *                                              *
 ************************************************/
void ffm_Multiply(const int mb, 
                  const thrust::device_vector<Real2> &d_v1,
                  const thrust::device_vector<Real2> &d_v2,
                  thrust::device_vector<Real> &d_hml)
{
   int ms1,ms2,mshift0;
   thrust::device_vector<Real2> d_tmp(d_v1.size());

   ff_Multiply_Summer(d_v1, cd->d_summer[1], d_tmp);

   if (mb==-1) 
   {
      ms1=0; 
      ms2 = psiHeader.get_ispin();
      mshift0 = psiHeader.get_ne(0)*psiHeader.get_ne(0);
   } 
   else 
   {
      ms1=mb; ms2=mb+1;
      mshift0 = 0;
   }

   for (int ms=ms1; ms<ms2; ++ms) 
   {
      int n       = psiHeader.get_ne(ms);
      int offset  = ms * psiHeader.get_ne(0) * cd->nfft3d;
      int moffset = ms*mshift0;
      blasGemm('T', 'N',  n, n, cd->n2ft3d, (Real) 1.,
	      (Real*) thrust::raw_pointer_cast(&d_tmp[offset]), cd->n2ft3d,
	      (Real*) thrust::raw_pointer_cast(&d_v2[offset]),  cd->n2ft3d, (Real) 0.,
	      (Real*) thrust::raw_pointer_cast(&d_hml[moffset]), n);
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
  
  int nn2[2];
  nn2[0] = psiHeader.get_ne(0)*psiHeader.get_ne(0);
  nn2[1] = psiHeader.get_ne(1)*psiHeader.get_ne(1);
  thrust::host_vector<Real> h_lambda(nn2[0]+nn2[1]);

  thrust::device_vector<Real> d_lambda(nn2[0]+nn2[1]);
  thrust::device_vector<Real> d_s22(nn2[0]);
  thrust::device_vector<Real> d_s21(nn2[0]);
  thrust::device_vector<Real> d_s12(nn2[0]);
  thrust::device_vector<Real> d_s11(nn2[0]);
  thrust::device_vector<Real> d_sa1(nn2[0]);
  thrust::device_vector<Real> d_sa0(nn2[0]);
  thrust::device_vector<Real> d_st1(nn2[0]);

  thrust::device_vector<Real2> d_psi_r(psi1->size());
  
  for (int ms=0; ms < psiHeader.get_ispin(); ++ms) {
    Real adiff;
    int nn = m_size(ms);
    
    // create S22 
    calc_psi_r(*psi2, cd->d_summer[1], d_psi_r);
    int offset = ms * psiHeader.get_ne(ms) * cd->nfft3d;
    blasGemm('T', 'N',  psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi2)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s22[0]), 
	     psiHeader.get_ne(ms));
    
    // create S21 
    blasGemm('T', 'N', psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi1)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s21[0]), 
	     psiHeader.get_ne(ms));
    
    // create S11 
    calc_psi_r(*psi1, cd->d_summer[1], d_psi_r);
    blasGemm('T', 'N', psiHeader.get_ne(ms), psiHeader.get_ne(ms), cd->n2ft3d, (Real) 1.,
	     (Real*) thrust::raw_pointer_cast(&d_psi_r[offset]), cd->n2ft3d,
	     (Real*) thrust::raw_pointer_cast(&(*psi1)[offset]), cd->n2ft3d, (Real) 0.,
	     (Real*) thrust::raw_pointer_cast(&d_s11[0]), 
	     psiHeader.get_ne(ms));

    d_scale_s22(ms,dte,d_s22);
    d_scale_s21(ms,dte,d_s21);
    d_s12 = d_s21;
    d_scale_s11(ms,dte,d_s11);
    
    int ii   = 0;
    int done = 0;
    d_sa0 = d_s22;
    while ((!done) && ((ii++)<ITERLMD)) {
      //dcopy_(&nn,s22,&one,sa1,&one);
      d_sa1 = d_s22;
      d_mmm_Multiply(ms,d_s21,d_sa0,1.0,d_sa1,1.0);
      d_mmm_Multiply(ms,d_sa0,d_s12,1.0,d_sa1,1.0);
      d_mmm_Multiply(ms,d_s11,d_sa0,1.0,d_st1,0.0);
      d_mmm_Multiply(ms,d_sa0,d_st1,1.0,d_sa1,1.0);
      d_st1 = d_sa1;
      d_daxpy(nn, -1.0, d_sa0,d_st1);
      
      adiff = d_dmax(ms,d_st1);
      if (adiff<CONVGLMD) 
         done = 1;
      else 
	d_sa0 = d_sa1;
    }
    if (!done) printf("ierr=10 adiff=%le\n",adiff);

    offset = ms*nn2[0];
    thrust::copy(d_sa1.begin(), d_sa1.end(),d_lambda.begin()+offset);
  }
  
  /* correction due to contraint */
  fmf_Multiply(-1,d_lambda,dte,1.0);
}

