/* Projection.C -
   Author - Eric Bylaska

    This class is used for finding the sigma and pi projections
between two atoms.
*/

#include	<math.h>
#include	"Projection.h"

#define zero	1.0e-6

/********************************
 *				*
 *	   Constructors		*
 *				*
 ********************************/
inline Projection::Projection()
{
   r = 1.0;
   sigma = Vector3(1.0,0.0,0.0);
   pi1   = Vector3(0.0,1.0,0.0);
   pi2   = Vector3(0.0,0.0,1.0);
}

inline Projection::Projection(Vector3 r1, Vector3 r2)
{
   double x,y,z;
  
   sigma = r2 - r1;
   x = sigma.x();
   y = sigma.y();
   z = sigma.z();
   r   = sqrt(x*x + y*y + z*z);

   /* define sigma */
   sigma = sigma * (1.0/r);

   /* define pi1 */
   if (fabs(x) < zero)
      pi1 = Vector3(1.0,0.0,0.0);
   else if (fabs(y) < zero)
      pi1 = Vector3(0.0,1.0,0.0);
   else if (fabs(z) < zero)
      pi1 = Vector3(0.0,0.0,1.0);
   else
      pi1 = Vector3(0.0,sqrt((z*z)/(y*y+z*z)), -y*sqrt(1.0/(y*y+z*z)) );

   /* define pi2 = sigma X pi1 */
   pi2 = sigma*pi1;

   /* define pi1 = pi2 x sigma */
   pi1 = pi2*sigma;


} /* Constructor */


/********************************
 *				*
 *	operator =		*
 *				*
 ********************************/

inline Projection Projection::operator =(Projection source)
{
   if (this != &source)
   { 
      r = source.r;
      sigma = source.sigma;
      pi1   = source.pi1;
      pi2   = source.pi2;
   }

   return *this;
}
