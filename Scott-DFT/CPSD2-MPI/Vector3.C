/* Vector3.C
   Author - Eric Bylaska

	this class is used for storing and manipulating
	xyz objects.
*/
#include	<iostream.h>
#include	<math.h>
#include	"Vector3.h"

/********************************
 *				*
 *	       r		*
 *				*
 ********************************/

/* returns the length of the vector 

*/
inline double	Vector3::r()
{
   double tmp;

   tmp = xyz[0]*xyz[0]
       + xyz[1]*xyz[1]
       + xyz[2]*xyz[2];
   tmp = sqrt(tmp);

   return tmp;
} /*r*/


/********************************
 *				*
 *         dot product		*
 *				*
 ********************************/

inline double	 Vector3::operator ,(Vector3 source)
{
   double tmp;

   tmp = xyz[0]*(source.x())
       + xyz[1]*(source.y())
       + xyz[2]*(source.z());

   return tmp;
} /* dot product */


/********************************
 *				*
 *      Equality operator       *
 *				*
 ********************************/
inline Vector3 Vector3::operator =(Vector3 source)
{
   if (this != &source)
   {
      xyz[0] = source.x();
      xyz[1] = source.y();
      xyz[2] = source.z();
   }
   return *this;
} /* = */


/********************************
 *				*
 *      boolean operators 	*
 *				*
 ********************************/
inline int	Vector3::operator==(Vector3 source)
{
   int bl;

   bl   = ((xyz[0]==source.x()) 
        && (xyz[1]==source.y())
        && (xyz[2]==source.z()) );

   return bl;
} /* == */

inline int	Vector3::operator !=(Vector3 source)
{
   int bl;
   
   bl   = ((xyz[0]==source.x()) 
        && (xyz[1]==source.y())
        && (xyz[2]==source.z()) );
   bl = !bl;

   return bl;
} /* != */


/********************************
 *				*
 *        Cross Product		*
 *				*
 ********************************/
inline Vector3 operator *(Vector3 v1, Vector3 v2)
{
   double x,y,z;

   x = v1.y()*v2.z() - v1.z()*v2.y();
   y = v1.z()*v2.x() - v1.x()*v2.z();
   z = v1.x()*v2.y() - v1.y()*v2.x();


   return Vector3(x,y,z);
}

/********************************
 *				*
 *	  Scalar Products	*
 *				*
 ********************************/

inline Vector3 operator *(const double a,   Vector3 v1)
{
   double x,y,z;

   x = a*v1.x();
   y = a*v1.y();
   z = a*v1.z();

   return Vector3(x,y,z);
}

inline Vector3 operator *(Vector3 v1, const double a)
{
   double x,y,z;

   x = a*v1.x();
   y = a*v1.y();
   z = a*v1.z();

   return Vector3(x,y,z);
}

inline Vector3 operator *(const int a,  Vector3 v1)
{
   double x,y,z;

   x = a*v1.x();
   y = a*v1.y();
   z = a*v1.z();

   return Vector3(x,y,z);
}

inline Vector3 operator *(Vector3 v1, const int a)
{
   double x,y,z;

   x = a*v1.x();
   y = a*v1.y();
   z = a*v1.z();

   return Vector3(x,y,z);
} /* scalar products */




/********************************
 *				*
 *     Addition operators 	*
 *				*
 ********************************/
inline Vector3 operator +(Vector3 v1, Vector3 v2)
{
   double x,y,z;

   x =  v1.x() + v2.x(); 
   y =  v1.y() + v2.y(); 
   z =  v1.z() + v2.z(); 

   return Vector3(x,y,z);
}

inline Vector3 operator +(const double a,   Vector3 v1)
{
   double x,y,z;

   x = a + v1.x();
   y = a + v1.y();
   z = a + v1.z();

   return Vector3(x,y,z);
}
inline Vector3 operator +(Vector3 v1, const double a)
{
   double x,y,z;

   x = v1.x() + a;
   y = v1.y() + a;
   z = v1.z() + a;

   return Vector3(x,y,z);
}

inline Vector3 operator +(const int a,   Vector3 v1)
{
   double x,y,z;

   x = v1.x() + a;
   y = v1.y() + a;
   z = v1.z() + a;

   return Vector3(x,y,z);
}

inline Vector3 operator +(Vector3 v1, const int a)
{
   double x,y,z;

   x = v1.x() + a;
   y = v1.y() + a;
   z = v1.z() + a;

   return Vector3(x,y,z);
}



/********************************
 *				*
 *      subtraction operators   *
 *				*
 ********************************/

inline Vector3 operator -(Vector3 v1, Vector3 v2)
{
   double x,y,z;

   x = v1.x() - v2.x();
   y = v1.y() - v2.y();
   z = v1.z() - v2.z();

   return Vector3(x,y,z);
}

inline Vector3 operator -(const double a,   Vector3 v1)
{
   double x,y,z;

   x = a - v1.x();
   y = a - v1.y();
   z = a - v1.z();

   return Vector3(x,y,z);
}

inline Vector3 operator -(Vector3 v1, const double a)
{
   double x,y,z;

   x = v1.x() - a;
   y = v1.y() - a;
   z = v1.z() - a;

   return Vector3(x,y,z);
}

inline Vector3 operator -(const int a,   Vector3 v1)
{
   double x,y,z;

   x = a - v1.x();
   y = a - v1.y();
   z = a - v1.z();

   return Vector3(x,y,z);
}

inline Vector3 operator -(Vector3 v1, const int a)
{
   double x,y,z;

   x = v1.x() - a;
   y = v1.y() - a;
   z = v1.z() - a;

   return Vector3(x,y,z);
}


/********************************
 *				*
 *	io operations		*
 *				*
 ********************************/

inline ostream& operator <<(ostream& s, Vector3 v1)
{
    s << v1.x() << " "
      << v1.y() << " "
      << v1.z() << " ";

   return s;
}

inline istream& operator >>(istream& s, Vector3& v1)
{
   s >> v1.xyz[0] >> v1.xyz[1] >> v1.xyz[2];

   return s;
}

