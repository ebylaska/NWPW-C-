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

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
Parallel::Parallel(int argc, char *argv[]) 
{

   int indim, np_dim[40],i1[40];
   int indx,d,dd,ii,np0,done;
   int *itmp;

   MPI::Init(argc,argv);

   indim = 1;
   for (ii=0; ii<40; ++ii) np_dim[ii] = -1;

   ii=0; done = 0;
   while (!done)
   {
      if ( strcmp(argv[ii],"-procgeom")==0)
      {
         done = 1;
         sscanf(argv[ii+1],"%d",&indim);
         for (int jj=0; jj<indim; ++jj) sscanf(argv[ii+2+jj],"%d",&(np_dim[jj]));
      }
      ++ii;
      if (ii>=argc) done = 1;
   }

   dim = indim;

   np_i     = new int [dim+1];
   taskid_i = new int [dim+1];
   comm_i   = new MPI::Intracomm [dim+1];
   group_i  = new MPI::Group [dim+1];

   comm_i[0]   = MPI::COMM_WORLD;
   np_i[0]     = comm_i[0].Get_size();
   taskid_i[0] = comm_i[0].Get_rank();

   procNd      = new int [np_i[0]];

   for (ii=0; ii<dim; ++ii)
   {
      np_i[ii+1] = np_dim[ii];
      if (np_i[ii+1] < 1) 
         np_i[ii+1] = 1;
   }

   np0 = np_i[0];
   for (ii=dim; ii>1; --ii)
   {
      while (((np0%np_i[ii]) > 0) && (np_i[ii] > 1))
        np_i[ii] -= 1;
      np0 /= np_i[ii];
   }
   np_i[1] = np0;

   np_dim[0] = 1;
   for (d=1; d<=dim; ++d)
      np_dim[d] = np_dim[d-1]*np_i[d];

   for (d=0; d<=dim; ++d) 
      i1[d] = 0;

   for (ii=0; ii<np_i[0]; ++ii)
   {
      if (i1[0]==taskid_i[0])
      {
         for (d=1; d<=dim; ++d)
            taskid_i[d] = i1[d];
      }
      indx = i1[1];
      for (d=2; d<=dim; ++d)
         indx += i1[d]*np_dim[d-1];
      procNd[indx] = i1[0];

      i1[0] += 1;
      i1[1] += 1;
      d = 1;
      while (d<dim)
      {
         if (i1[d] >= np_i[d])
         {
            i1[d]    = 0;
            i1[d+1] += 1;
         }
         ++d;
      }
   }
   
   /* create communicators along each dimension */
   //np0     = np_i[0]; 
   //np_i[0] = 1;
   for (d=1; d<=dim; ++d)
   {
      itmp = new int [np_i[d]];
      for (ii=0; ii<np_i[d]; ++ii)
      {
         indx = ii*np_dim[d-1];
         for (dd=1; dd<=dim; ++dd)
            if (dd != d)
               indx += taskid_i[dd]*np_dim[dd-1];
         itmp[ii] = procNd[indx];
      }
      cout << "taskid=" << taskid_i[0];
      cout << " d=" << d;
      cout << " itmp=";
      for (ii=0; ii<np_i[d]; ++ii)
         cout << itmp[ii] << " ";
      cout << "\n";

      //comm_i[d] = MPI::COMM_WORLD.Create(MPI::COMM_WORLD.Get_group().Incl(np_i[d],itmp));
      group_i[d] = MPI::COMM_WORLD.Get_group().Incl(np_i[d],itmp);
      comm_i[d]  = MPI::COMM_WORLD.Create(group_i[d]);

      delete [] itmp;
   }
   //np_i[0] = np0;
}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Parallel::~Parallel()
{
   for (int d=1; d<=dim; ++d)
   {
     group_i[d].Free();
     comm_i[d].Free();
   }
   delete [] np_i;
   delete [] taskid_i;
   delete [] comm_i;
   delete [] group_i;
   delete [] procNd;

   MPI::Finalize();
}

/********************************
 *                              *
 *          SumAll              *
 *                              *
 ********************************/
double Parallel::SumAll(const int d, const double sum)
{
   double sumout;

   comm_i[d].Allreduce(&sum,&sumout,1,MPI_DOUBLE_PRECISION,MPI_SUM);
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
   sumout = new double [n];
   comm_i[d].Allreduce(sum,sumout,n,MPI_DOUBLE_PRECISION,MPI_SUM);
   for (int i=0; i<n; ++i) sum[i] = sumout[i];
   delete [] sumout;
}



/********************************
 *                              *
 *       Brdcst_Values          *
 *                              *
 ********************************/
void Parallel::Brdcst_Values(const int d, const int root, const int n, double *sum)
{
   comm_i[d].Bcast(sum,n,MPI_DOUBLE_PRECISION,root);
}

/********************************
 *                              *
 *       Brdcst_iValues         *
 *                              *
 ********************************/
void Parallel::Brdcst_iValues(const int d, const int root, const int n, int *sum)
{
   comm_i[d].Bcast(sum,n,MPI_INTEGER,root);
}
