#ifndef _MAPPING_H_
#define _MAPPING_H_
/* Parallel.h
   Author - Eric Bylaska

	this class is used defining 3d parallel maps
*/

#include	<mpi.h>

class Parallel {

    int *procNd;

public:
        int dim,np,taskid;
        int *np_i,*taskid_i;
        MPI::Intracomm comm,*comm_i;

	/* Constructors */
	Parallel(int, char **);

        /* destructor */
	~Parallel();

};

#endif
