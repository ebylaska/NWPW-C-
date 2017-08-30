
#include	<iostream>
#include	<cstdio>
#include	<stdio.h>
#include	<cmath>
#include	<cstdlib>
using namespace std;

#include	"Parallel.h"
#include	"d3db.h"
#include	"util_date.h"


int main(int argc, char *argv[])
{
   Parallel myparallel(argc,argv);

   int nfft[3],ne;
   int i,ii,maptype,ishift;
   char date[26];
   double scale,sum1,sum2;
   double cpu1,cpu2,cpu3,cpu4;
   double fftx_flop,ffty_flop,fftz_flop,fft3d_flop,cr_flop,rc_flop;
   double E[50],deltae,deltac,deltar,viral,unita[9];
   double *psi1,*psi_r;

   maptype = 2;
   ne = 100;
   nfft[0] = nfft[1] = nfft[2] = 64;

   ii = 1;
   while (ii<argc)
   {
      if (strcmp("-ne",argv[ii])==0)
      {
         ++ii; if (ii<=argc) sscanf(argv[ii], "%d",&ne);
      }
      else if (strcmp("-nfft",argv[ii])==0)
      {
         ++ii; if (ii<=argc) sscanf(argv[ii], "%d",&nfft[0]);
         ++ii; if (ii<=argc) sscanf(argv[ii], "%d",&nfft[1]);
         ++ii; if (ii<=argc) sscanf(argv[ii], "%d",&nfft[2]);
      }
      else if (strcmp("-maptype",argv[ii])==0)
      {
         ++ii; if (ii<=argc) sscanf(argv[ii], "%d",&maptype);
      }
      ++ii;
   }


   if (myparallel.is_master())
   {
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *         Parallel 3d FFT test calculation          *\n";
      cout << "          *                                                   *\n";
      cout << "          *            version #1.00   07/20/09               *\n";
      cout << "          *                                                   *\n";
      cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      cout << "          *                                                   *\n";
      cout << "          *****************************************************\n";
      cout << "          >>> job started at       " << util_date() << " <<<\n";
   }
   myparallel.init2d(1);
   d3db myd3db(&myparallel,maptype,nfft[0],nfft[1],nfft[2]);
   psi_r = myd3db.r_alloc();
   psi1 = myd3db.r_nalloc(ne);


//                 |**************************|
// *****************   summary of input data  **********************
//                 |**************************|
   fftx_flop = 2.5*nfft[0]*log(nfft[0])/log(2.0);
   ffty_flop = 5.0*nfft[1]*log(nfft[1])/log(2.0);
   fftz_flop = 5.0*nfft[2]*log(nfft[2])/log(2.0);

   fft3d_flop = fftx_flop * nfft[1]*nfft[2]
              + ffty_flop * (nfft[0]/2+1) * nfft[2]
              + fftz_flop * (nfft[0]/2+1) * nfft[1];

    cr_flop = ne*fft3d_flop;
    rc_flop = ne*fft3d_flop + ((nfft[0]+2)*nfft[1]*nfft[2]);

   if (myparallel.is_master())
   {
      cout << "\n\n\n";
      cout << "          ==============  summary of input  ==================\n";
      cout << "\n";
      cout << " number of processors used: " << myparallel.np() << "\n";
      cout << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << "\n";
      if (myd3db.maptype==1) cout << " parallel mapping         : slab"    << "\n";
      if (myd3db.maptype==2) cout << " parallel mapping         : hilbert" << "\n";


      printf("\n\n number of orbitals = %d\n",ne);
      printf("\n fft grid           = %4d x %4d x %4d (nfft3d=%d nfft3d_all=%d)\n",
              nfft[0],nfft[1],nfft[2],
              myd3db.nfft3d,
              (nfft[0]/2+1)*nfft[1]*nfft[2]);

      printf("\n\n number flops per fftx = %le\n",fftx_flop);
      printf(" number flops per ffty = %le\n",ffty_flop);
      printf(" number flops per fftz = %le\n",fftz_flop);
      printf(" number flops per 3d fft = %le\n\n",fft3d_flop);
      printf(" number flops for (complex to real) cr_ffts = %le\n",cr_flop);
      printf(" number flops for (real to complex) rc_ffts = %le\n",rc_flop);


   }

//                 |**************************|
// *****************   start iterations       **********************
//                 |**************************|
   if (myparallel.is_master()) 
      cout << "\n\n---- intializing " << ne << " orbitals ----\n ";

   scale = 1.0/((double) myd3db.nx*myd3db.ny*myd3db.nz);

   for (i=0; i<myd3db.n2ft3d; ++i)
      psi_r[i] = 2.0*rand()/((double) RAND_MAX) - 1.0;
   myd3db.r_zero_ends(psi_r);
   myd3db.rc_fft3d(psi_r);

   sum2 = 0.0;
   for (i=0; i<myd3db.n2ft3d; ++i)
       sum2 += (psi_r[i])*(psi_r[i]);
   myd3db.r_SMul(sqrt(1.0/sum2), psi_r);

   ishift = 0;
   for (ii=0; ii<ne; ++ii)
      for (i=0; i<myd3db.n2ft3d; ++i)
          psi1[ishift++] = psi_r[i];


   if (myparallel.is_master()) 
   {
      cout << "\n\n\n";
      cout << "          =============  running the ffts  =================\n";

   }
   if (myparallel.is_master()) 
   {
      cout << "\n\n---- Fourier Transforming (complex to real) " << ne << " orbitals ----\n ";
   }

   if (myparallel.is_master()) seconds(&cpu1);
   for (ii=0; ii<ne; ++ii)
   {
       myd3db.cr_fft3d(   &psi1[ii*(myd3db.n2ft3d)]);
       myd3db.r_zero_ends(&psi1[ii*(myd3db.n2ft3d)]);
   }
   if (myparallel.is_master()) seconds(&cpu2);


   if (myparallel.is_master()) 
   {
      cout << "\n\n---- Fourier Transforming (real to complex) " << ne << " orbitals ----\n ";
   }

   if (myparallel.is_master()) seconds(&cpu3);
   for (ii=0; ii<ne; ++ii)
   {
       myd3db.r_SMul(scale, &psi1[ii*(myd3db.n2ft3d)]);
       myd3db.rc_fft3d(     &psi1[ii*(myd3db.n2ft3d)]);
   }
   if (myparallel.is_master()) seconds(&cpu4);


//                  |***************************|
// ****************** report summary of results **********************
//                  |***************************|
   if (myparallel.is_master()) 
   {
      cout << "\n\n\n";
      cout << "          =============  checking results  =================\n";

   }
   for (ii=0; ii<ne; ++ii)
   {
      ishift = ii*myd3db.n2ft3d;
      sum1 = 0.0;
      for (i=0; i<myd3db.n2ft3d; ++i)
          sum1 += (psi1[i+ishift] - psi_r[i])*(psi1[i+ishift] - psi_r[i]);

      sum2 = 0.0;
      for (i=0; i<myd3db.n2ft3d; ++i)
          sum2 += (psi1[i+ishift])*(psi1[i+ishift]);
      
      sum1 = myparallel.SumAll(0,sum1);
      sum2 = myparallel.SumAll(0,sum2);

      if (sum1 > 1.0e-9)
      if (myparallel.is_master())
         cout << "checking ffts of " << ii << "  orbital: error = " << sum1 << " norm = " << sum2 << "\n";
    }


   /* deallocate memory */
   myd3db.r_dealloc(psi1);

//                 |**************************|
// *****************   report consumed time   **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      double t1 = cpu2-cpu1;
      double t2 = cpu4-cpu3;
      cout.setf(ios::scientific);
      cout << "\n";
      cout << "          =============  timing results  =================\n\n";
      cout << "time for " << ne << " cr_fft3d in seconds = " << t1 << " (" << t1/((double) ne) << " per fft)\n";
      cout << "time for " << ne << " rc_fft3d in seconds = " << t2 << " (" << t2/((double) ne) << " per fft)\n";
      cout << "\n";
      cout << "flop rate for cr_fft3d = " << cr_flop/t1 << " flops/second (" << cr_flop/t1*1.0e-9 << " Gflops)\n";
      cout << "flop rate for rc_fft3d = " << rc_flop/t2 << " flops/second (" << rc_flop/t2*1.0e-9 << " Gflops)\n\n";
      cout << "          >>> job completed at     " << util_date() << " <<<\n";

   }
}
