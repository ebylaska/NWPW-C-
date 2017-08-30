#ifndef _D3dB_H_
#define _D3dB_H_
/* d3db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include	<mpi.h>
#include	"Mapping3.h"

class d3db : public Mapping3 {


public:

        /* constructor */
	d3db(const int,MPI_Comm,const int,const int,const int);

        /* destructor */
        ~d3db() {}

        /* r array operators */
        double * r_alloc();
        void     r_dealloc(double *);
        void     r_zero(double *);
        void     rr_copy(const double *,double *);
        void     rr_SMul(const double, const double *, double *);
        void	 r_SMul(const double, double *);
        double	 r_dsum(const double *);
        void  	 r_zero_ends(double *);
        void  	 r_abs(double *);
        void  	 rrr_Sum(const double *, const double *, double *);
        void  	 rrr_Minus(const double *, const double *, double *);
        void  	 rrr_Divide(const double *, const double *, double *);
        void  	 rr_daxpy(const double, const double *, double *);
        //void  	 r_read(const int, const int, double *);

        /* t array operators */
        double * t_alloc();
        void     t_dealloc(double *);
};

#endif
