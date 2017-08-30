/* Parallel.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/
#include        <iostream>
#include        <cmath>
#include        <cstdlib>
using namespace std;


#include	<mpi.h>
#include	"Parallel.h"

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
Parallel::Parallel(int argc, char *argv[]) 
{

   int indim, np_dim[20];
   int ii,np0,done;

   MPI::Init(argc,argv);

   indim = 1;
   for (ii=0; ii<indim; ++ii) np_dim[ii] = -1;

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

   cout << "commworld=" << MPI::COMM_WORLD << "\n";

   comm = MPI::COMM_WORLD;
   np     = comm.Get_size();
   taskid = comm.Get_rank();

   np_i     = new int [dim];
   taskid_i = new int [dim];
   comm_i   = new MPI::Intracomm [dim];

   for (ii=0; ii<dim; ++ii)
   {
      np_i[ii] = np_dim[ii];
      if (np_i[ii] < 1) 
         np_i[ii] = 1;
   }
   np0 = np;

   ii = dim;
   for (ii=dim-1; ii>0; --ii)
   {
      while (((np0%np_i[ii]) > 0) && (np_i[ii] > 1))
        np_i[ii] -= 1;
      np0 /= np_i[ii];
   }
   np_i[0] = np0;
}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
Parallel::~Parallel()
{
   delete [] np_i;
   delete [] taskid_i;
   delete [] comm_i;

   MPI::Finalize();
}
