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
Parallel::Parallel(int argc, char *argv[], 
                   const int indim, int *np_dim)
{
   int ii,np0;

   done = 0;
   while (!done)
   {
      if (strcmp("-pgeom"argv[1],"%s",
   MPI_Init(&argc,&argv);
   dim = indim;

   comm = MPI_COMM_WORLD;
   MPI_Comm_size(MPI_COMM_WORLD,&np);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

   np_i     = new int [dim];
   taskid_i = new int [dim];
   comm_i   = new MPI_Comm [dim];

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

   MPI_Finalize();
}
