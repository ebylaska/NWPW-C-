#ifndef _PROJECTION_H_
#define _PROJECTION_H_
/* Projection.h -
   Author - Eric Bylaska

    This class is used for finding the sigma and pi projections
between two atoms.
*/
#include	"Vector3.h"

class Projection {
public:

       /* Constructors */
	inline	Projection();
	inline	Projection(Vector3, Vector3);

	inline	Projection operator =(Projection);

 	double  r;
	Vector3	sigma,
		pi1,
		pi2;
};

#include "Projection.C"

#endif
