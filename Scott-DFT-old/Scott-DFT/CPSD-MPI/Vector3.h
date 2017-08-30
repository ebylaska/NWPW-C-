#ifndef _VECTOR3_H_
#define _VECTOR3_H_
/* Vector3.h
   Author - Eric Bylaska

	this class is used for storing and manipulating
	xyz objects.
*/
using namespace std;

#include	<iostream.h>

class Vector3 {

	double	xyz[3];

public:

	/* Constructors */
	inline Vector3()
        {  xyz[0]=0.0;
 	   xyz[1]=0.0;
	   xyz[2]=0.0;
        }
	inline Vector3(const double *v) 
        {  xyz[0]=v[0];
 	   xyz[1]=v[1];
	   xyz[2]=v[2];
        }
	inline Vector3(const double x, const double y, const double z)
        {  xyz[0]=x;
 	   xyz[1]=y;
	   xyz[2]=z;
        }

        /* copy constructor */
        inline Vector3(const Vector3& source)
        {
           xyz[0] = source.xyz[0];
           xyz[1] = source.xyz[1];
           xyz[2] = source.xyz[2];
        }
           


	inline double  *array() {return xyz;}
	inline double	operator [](const int i) {return xyz[i];}
        inline double	x() { return xyz[0];}
	inline double	y() { return xyz[1];}
	inline double	z() { return xyz[2];}
        double   r();

	inline double   l2norm() { return( xyz[0]*xyz[0] +
			   	           xyz[1]*xyz[1] +
					   xyz[2]*xyz[2]);   }
	Vector3 operator =(Vector3);
	int	 operator ==(Vector3);
	int	 operator !=(Vector3);

        /* dot product */
	double	 operator ,(Vector3);  

        /* Cross product */
	friend	Vector3 operator *(Vector3, Vector3);

        /* Scalar products */
	friend	Vector3 operator *(const double,  Vector3);
  	friend	Vector3 operator *(Vector3, const double);
  	friend	Vector3 operator *(const int,     Vector3);
  	friend	Vector3 operator *(Vector3, const int);

        /* addition operators */
  	friend	Vector3 operator +(Vector3,       Vector3);
  	friend	Vector3 operator +(const double,  Vector3);
  	friend	Vector3 operator +(Vector3,       const double);
  	friend	Vector3 operator +(const int,     Vector3);
  	friend	Vector3 operator +(Vector3,       const int);

        /* substraction operators */
  	friend	Vector3 operator -(Vector3,       Vector3);
  	friend	Vector3 operator -(const double,  Vector3);
  	friend	Vector3 operator -(Vector3,       const double);
  	friend	Vector3 operator -(const int,     Vector3);
  	friend	Vector3 operator -(Vector3,       const int);

	/* io operations */
  	friend	ostream& operator <<(ostream&, Vector3);
  	friend	istream& operator >>(istream&, Vector3&);
};

#include	"Vector3.C"

#endif
