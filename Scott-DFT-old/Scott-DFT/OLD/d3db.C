/* d3db.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/
#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;

#include	<mpi.h>
#include	"d3db.h"

int d3db_get_np(MPI_Comm comm)
{
    int np;
    MPI_Comm_size(comm,&np);
    return np;
}
int d3db_get_taskid(MPI_Comm comm)
{
    int taskid;
    MPI_Comm_rank(comm,&taskid);
    return taskid;
}

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

d3db::d3db(const int inmaptype, MPI_Comm comm, const int nx, const int ny, const int nz)
   : Mapping3(inmaptype,d3db_get_np(comm),d3db_get_taskid(comm),nx,ny,nz)
{
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
 *         d3db::rr_copy         *
 *                              *
 ********************************/
void d3db::r_copy(const double *ptr1, double *ptr2)
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

