/* Parallel.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/

#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;

#include	"mpi.h"
#include	"Parallel.h"
//#include	"control.h"

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
Parallel::Parallel(int argc, char *argv[]) 
{

   int indx,d,dd,ii,np0,done;

   MPI::Init(argc,argv);

   dim = 1;

   comm_i[0]   = MPI::COMM_WORLD;
   npi[0]     = comm_i[0].Get_size();
   taskidi[0] = comm_i[0].Get_rank();

   comm_i[1]   = comm_i[0];
   comm_i[2]   = comm_i[0];
   npi[1]     = npi[0];
   npi[2]     = 1;
   taskidi[1] = taskidi[0];
   taskidi[2] = 0;

   procNd      = new int [npi[0]];

   npi[1] = npi[0];
   npi[2] = 1;
   for (int i=0; i<npi[0]; ++i)
      procNd[i] = i;

}

void Parallel::init2d(const int ncolumns)
{
   int ii;

   if (ncolumns>1)
   {
      dim = 2;
      npi[1] = npi[0]/ncolumns;
      npi[2] = ncolumns;

      int icount = 0;
      for (int j=0; j<npi[2]; ++j)
      for (int i=0; i<npi[1]; ++i)
      {
         if (icount==taskidi[0])
         {
            taskidi[1] = i;
            taskidi[2] = j;
         }
         procNd[i+j*npi[1]] = icount;
         icount = (icount+1)%npi[0];
      }

      int *itmp = new int[npi[0]];

      for (int i=0; i<npi[1]; ++i) itmp[i] = procNd[i+taskidi[2]*npi[1]];
      group_i[1] = MPI::COMM_WORLD.Get_group().Incl(npi[1],itmp);
      comm_i[1]  = MPI::COMM_WORLD.Create(group_i[1]);

      for (int j=0; j<npi[2]; ++j) itmp[j] = procNd[taskidi[1]+j*npi[1]];
      group_i[2] = MPI::COMM_WORLD.Get_group().Incl(npi[2],itmp);
      comm_i[2]  = MPI::COMM_WORLD.Create(group_i[2]);

      delete [] itmp;
   }

   //ii = 3+control_pfft3_qsize();
   ii = 3+8;
   reqcnt  = new int[ii];
   request = new MPI::Request*[ii];

}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Parallel::~Parallel()
{
   if (dim>1)
   {
      for (int d=1; d<=dim; ++d)
      {
         group_i[d].Free();
         comm_i[d].Free();
      }
   }

   delete [] procNd;

   delete [] reqcnt;
   delete [] request;


   MPI::Finalize();
}
/********************************
 *                              *
 *          MaxAll              *
 *                              *
 ********************************/
double Parallel::MaxAll(const int d, const double sum)
{
   double sumout;
   if (npi[d]>1) 
      comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_MAX);
   else 
      sumout = sum;
   return sumout;
}

/********************************
 *                              *
 *          SumAll              *
 *                              *
 ********************************/
double Parallel::SumAll(const int d, const double sum)
{
   double sumout;

   if (npi[d]>1) 
      comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_SUM);
   else
      sumout = sum;
   return sumout;
}

/********************************
 *                              *
 *         ISumAll              *
 *                              *
 ********************************/
int Parallel::ISumAll(const int d, const int sum)
{
   int sumout;

   if (npi[d]>1) 
      comm_i[d].Allreduce(&sum,&sumout,1,MPI_INTEGER,MPI_SUM);
   else
      sumout = sum;
   return sumout;
}


/********************************
 *                              *
 *       Vector_SumAll          *
 *                              *
 ********************************/
void Parallel::Vector_SumAll(const int d, const int n, double *sum)
{
   double *sumout;
   if (npi[d]>1)
   {
      sumout = new double [n];
      comm_i[d].Allreduce(sum,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM);
      for (int i=0; i<n; ++i) sum[i] = sumout[i];
      delete [] sumout;
   }
}

/********************************
 *                              *
 *       Vector_ISumAll         *
 *                              *
 ********************************/
void Parallel::Vector_ISumAll(const int d, const int n, int *sum)
{
   int *sumout;
   if (npi[d]>1)
   {
      sumout = new int [n];
      comm_i[d].Allreduce(sum,sumout,n,MPI_INTEGER,MPI_SUM);
      for (int i=0; i<n; ++i) sum[i] = sumout[i];
      delete [] sumout;
   }
}


/********************************
 *                              *
 *       Brdcst_Values          *
 *                              *
 ********************************/
void Parallel::Brdcst_Values(const int d, const int root, const int n, double *sum)
{
   if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_DOUBLE_PRECISION,root);
}

/********************************
 *                              *
 *       Brdcst_iValues         *
 *                              *
 ********************************/
void Parallel::Brdcst_iValues(const int d, const int root, const int n, int *sum)
{
   if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_INTEGER,root);
}

/********************************
 *                              *
 *       Brdcst_iValue          *
 *                              *
 ********************************/
void Parallel::Brdcst_iValue(const int d, const int root, int *sum)
{
   if (npi[d]>1) comm_i[d].Bcast(sum,1,MPI_INTEGER,root);
}

/********************************
 *                              *
 *       Brdcst_cValues         *
 *                              *
 ********************************/
void Parallel::Brdcst_cValues(const int d, const int root, const int n, void *sum)
{
   if (npi[d]>1) comm_i[d].Bcast(sum,n,MPI_CHAR,root);
}

/********************************
 *                              *
 *       Parallel::dsend        *
 *                              *
 ********************************/
void Parallel::dsend(const int d, const int tag, const int procto, const int n, double *sum)
{
   if (npi[d]>1) comm_i[d].Send(sum,n,MPI_DOUBLE_PRECISION,procto,tag);
}


/********************************
 *                              *
 *       Parallel::dreceive     *
 *                              *
 ********************************/
void Parallel::dreceive(const int d, const int tag, const int procfrom, const int n, double *sum)
{
   MPI::Status status;
   if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);
}



/********************************
 *                              *
 *       Parallel::isend        *
 *                              *
 ********************************/
void Parallel::isend(const int d, const int tag, const int procto, const int n, int *sum)
{
   if (npi[d]>1) comm_i[d].Send(sum,n,MPI_INTEGER,procto,tag);
}  


/********************************
 *                              *
 *       Parallel::ireceive     *
 *                              *
 ********************************/
void Parallel::ireceive(const int d, const int tag, const int procfrom, const int n, int *sum)
{
   MPI::Status status;
   if (npi[d]>1) comm_i[d].Recv(sum,n,MPI_INTEGER,procfrom,tag);
}  




/********************************
 *                              *
 *       Parallel::astart       *
 *                              *
 ********************************/
void Parallel::astart(const int d, const int sz)
{
   if (npi[d]>1)
   {
      reqcnt[d]  = 0;
      request[d] = new MPI::Request[sz];
   }
}
/********************************
 *                              *
 *       Parallel::aend         *
 *                              *
 ********************************/
void Parallel::aend(const int d)
{
   //MPI::Status status[reqcnt[d]];
   //request[d][0].Waitall(reqcnt[d],request[d],status);
   if (npi[d]>1)
   {
      request[d][0].Waitall(reqcnt[d],request[d]);
      delete [] request[d];
   }
}

/********************************
 *                              *
 *       Parallel::adreceive     *
 *                              *
 ********************************/
void Parallel::adreceive(const int d, const int tag, const int procfrom, const int n, double *sum)
{
   MPI::Status status;
   if (npi[d]>1) request[d][reqcnt[d]++] = comm_i[d].Irecv(sum,n,MPI_DOUBLE_PRECISION,procfrom,tag);
}

/********************************
 *                              *
 *       Parallel::adsend        *
 *                              *
 ********************************/
void Parallel::adsend(const int d, const int tag, const int procto, const int n, double *sum)
{
   if (npi[d]>1) request[d][reqcnt[d]++] = comm_i[d].Isend(sum,n,MPI_DOUBLE_PRECISION,procto,tag);
}

