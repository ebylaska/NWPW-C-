#ifndef _Parallel_H_
#define _Parallel_H_
/* Parallel.h
   Author - Eric Bylaska

	this class is used defining nd parallel geometries
*/

#include	"mpi.h"

#define	MASTER	0

class Parallel {

    int *procNd;
    MPI::Intracomm *comm_i;
    MPI::Group     *group_i;

public:
        int dim;
        int *np_i,*taskid_i;

	/* Constructors */
	Parallel(int, char **);

        /* destructor */
	~Parallel();

        int convert_taskid_ij(const int, const int);

       /* SumAll */
      double SumAll(const int, const double);
      void Vector_SumAll(const int, const int, double *);
      void Brdcst_Values(const int, const int, const int, double *);
      void Brdcst_iValues(const int, const int, const int, int *);
};

#endif
