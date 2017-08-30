/* d3db.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/
/*
#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;
*/

#include	"Parallel.h"
#include	"d3db.h"


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

d3db::d3db(Parallel& inparall,const int inmaptype, const int nx, const int ny, const int nz)
   : Mapping3(inmaptype,inparall.np_i[1],inparall.taskid_i[1],nx,ny,nz)
{
   parall = inparall;
}


/********************************
 *                              *
 *         d3db::r_alloc        *
 *                              *
 ********************************/
double * d3db::r_alloc()
{
   double *ptr = new double[n2ft3d];
   return ptr;
}

/********************************
 *                              *
 *         d3db::r_dealloc      *
 *                              *
 ********************************/
void d3db::r_dealloc(double *ptr)
{
   delete [] ptr;
}

/********************************
 *                              *
 *         d3db::r_zero         *
 *                              *
 ********************************/
void d3db::r_zero(double *ptr)
{
   int i;
   int m = n2ft3d%7;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr[i] = 0.0;
   if (n2ft3d<7) 
      return;

   for (i=m; i<n2ft3d; i+=7)
   {
      ptr[i]   = 0.0;
      ptr[i+1] = 0.0;
      ptr[i+2] = 0.0;
      ptr[i+3] = 0.0;
      ptr[i+4] = 0.0;
      ptr[i+5] = 0.0;
      ptr[i+6] = 0.0;
   }
   return;
}

/********************************
 *                              *
 *        d3db::rr_copy         *
 *                              *
 ********************************/
void d3db::rr_copy(const double *ptr1, double *ptr2)
{
   int i;
   int m = n2ft3d%7;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = ptr1[i];
   if (n2ft3d<7) 
      return;

   for (i=m; i<n2ft3d; i+=7)
   {
      ptr2[i]   = ptr1[i];
      ptr2[i+1] = ptr1[i+1];
      ptr2[i+2] = ptr1[i+2];
      ptr2[i+3] = ptr1[i+3];
      ptr2[i+4] = ptr1[i+4];
      ptr2[i+5] = ptr1[i+5];
      ptr2[i+6] = ptr1[i+6];
   }
   return;
}

/********************************
 *                              *
 *         d3db::rr_SMul         *
 *                              *
 ********************************/
void d3db::rr_SMul(const double da, const double *ptr1, double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = da*ptr1[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   = da*ptr1[i];
      ptr2[i+1] = da*ptr1[i+1];
      ptr2[i+2] = da*ptr1[i+2];
      ptr2[i+3] = da*ptr1[i+3];
      ptr2[i+4] = da*ptr1[i+4];
   }
   return;
}




/********************************
 *                              *
 *         d3db::r_SMul        *
 *                              *
 ********************************/
void d3db::r_SMul(const double da, double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] *= da;
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   = da;
      ptr2[i+1] = da;
      ptr2[i+2] = da;
      ptr2[i+3] = da;
      ptr2[i+4] = da;
   }
   return;
}


/********************************
 *                              *
 *         d3db::t_alloc        *
 *                              *
 ********************************/
double * d3db::t_alloc()
{
   double *ptr = new double[nfft3d];
   return ptr;
}

/********************************
 *                              *
 *         d3db::t_dealloc      *
 *                              *
 ********************************/
void d3db::t_dealloc(double *ptr)
{
   delete [] ptr;
}


/********************************
 *                              *
 *         d3db::c_read         *
 *                              *
 ********************************/
void d3db::c_read(const int iunit, double *a, const int jcol)
{
   int jstart,jend,fillcolumn;
   
   if (jcol<0) 
   {
      jstart = 0;
      jend   = np-1;
      fillcolumn = 1;
   }
   else
   {
      jstart = jend = jcol;
      fillcolumn = (taskid==jcol);
   }

   /**********************
    **** slab mapping ****
	**********************/
   if (maptype==1)
   { 
   
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
      {
         for (int k=0; k<nz; ++k)
         {
            dread(iunit,tmp,(nx+2)*ny)

            index = 2*parall.ijktoindex(0,0,k);
            ii    = parall.ijktop(0,0,k);
            for (int jj=jstart jj<=jend; ++jj)
            {
               p_to = parall.convert_taskid_ij(ii,jj)
               if (p_to==MASTER)
               { 
                  index = (q-1)*(nx/2+1)*ny;
                  call dcopy((nx+2)*ny,tmp,1,a(index),1)
               }
               else
               {
                  msglen = (nx/2+1)*ny
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               }
            }
         }
      }
   
      /**** not master node ****/
      else if (fillcolumn) 
      {
         for (int k=0; k<nz; ++k)
		 {
            call D3dB_ktoqp(nb,k,q,ii)
            p_here = parall.convert_taskid_ij(ii,taskid_j)
            if (p_here==taskid)
		{ 

               msglen = (nx(nb)/2+1)*ny(nb)
               source  = MASTER

               call RCV(9+MSGDBL,tmp,mdtob(2*msglen),rcv_len,
     >                  source,rcv_proc,1)

               index = (q-1)*(nx(nb)/2+1)*ny(nb) + 1
               call dcopy((nx(nb)+2)*ny(nb),tmp,1,A(index),1)
               
            }
         }
      }
	  
   }
   
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
         
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
      {
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {

            call dread(iunit,tmp,(nx(nb)+2))

            index  = ijktoindex2(0,j,k)
            ii     = ijktop2(0,j,k)
            for (int jj=jstart; jj<=jend; ++jj)
            {
               p_to = parall.convert_taskid_ij(ii,jj)

               if (p_to.eq.MASTER) then
                  index = (q-1)*(nx(nb)/2+1) + 1
                  call dcopy((nx(nb)+2),tmp,1,A(index),1)
               else
                  msglen = (nx(nb)/2+1)
                  call SND(9+MSGDBL,tmp,mdtob(2*msglen),p_to,1)
               end if
            }
         }
      }
      
      /**** not master node ****/
      else if (fillcolumn) 
      {
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {
            
            index  = ijktoindex2(0,j,k)
            ii     = ijktop2(0,j,k)
            p_here = parall.convert_taskid_ij(ii,taskid_j)
            if (p_here==taskid) 
            {
               msglen = (nx+2)
               source  = MASTER
               call RCV(9+MSGDBL,tmp,mdtob(msglen),rcv_len,source,rcv_proc,1)
               //index = (q-1)*(nx(nb)/2+1) + 1
               call dcopy((nx+2),tmp,1,A(index),1)
            }
         }
      }

   }

}

