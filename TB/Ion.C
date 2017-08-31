/* Ion.C - 
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>
#include	"Convert_Name.h"
#include	"Ion.h"


/* Constructors */

inline Ion::Ion() : Vector3()
{
   strcpy(name,"X");
   charge = NameToCharge(name);
   mass   = NameToMass(name);
}

inline Ion::Ion(const char* atom, const Vector3 v1) : Vector3(v1)
{
   strcpy(name,atom);
   charge = NameToCharge(name);
   mass   = NameToMass(name);
}
  
inline Ion::Ion(const char* atom, 
         const double x, const double y, const double z) : Vector3(x,y,z)
{
   strcpy(name,atom);
   charge = NameToCharge(name);
   mass   = NameToMass(name);
}

inline Ion::Ion(const char* atom, const double* v1) : Vector3(v1)
{
   strcpy(name,atom);
   charge = NameToCharge(name);
   mass   = NameToMass(name);
}

inline Ion::Ion(Ion& i1) : Vector3(i1)
{
   strcpy(name,i1.name);
   charge  = i1.charge;
   mass    = i1.mass;

}


inline Ion& Ion::operator =(Vector3 v1)
{
    double *v2;

    v2 = this->array();
    v2[0] = v1[0];
    v2[1] = v1[1];
    v2[2] = v1[2];

   return *this;
}

inline int	Ion::Charge()
{
   return charge;
}

inline double	Ion::Mass()
{
   return mass;
}

inline char*	Ion::Name()
{
   return name;
}


           
/* io operations */
inline ostream& operator <<(ostream& s, Ion i1)
{
   Vector3 v1 = i1;

   
   s << i1.Name() << " ";
   s.setf(ios::fixed,ios::floatfield);
   s << v1;
   s << i1.Charge() << "\n";

   return s;

}

inline istream& operator >>(istream& s, Ion& i1)
{
   char    name[5];
   Vector3 v1;
   double  not_used;

   s >> ws >> name;
   s >> v1; 
   s >> not_used; /* charge from xyz file not used */

   i1 = v1;
   strcpy(i1.name,name);
   i1.charge = NameToCharge(name);
   i1.mass   = NameToMass(name);

   return s;
}

