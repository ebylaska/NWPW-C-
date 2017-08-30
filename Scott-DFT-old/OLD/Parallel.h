#ifndef _MAPPING_H_
#define _MAPPING_H_
/* Parallel.h
   Author - Eric Bylaska

	this class is used defining 3d parallel maps
*/

#include	<mpi.h>

class Parallel {

    int dim;
    int *procNd;

public:
        int np,taskid;
        int *np_i,*taskid_i;
        MPI_Comm comm,*comm_i;

	/* Constructors */
	Parallel(int, char **, const int, int *);

        /* destructor */
	~Parallel();

};

#endif
