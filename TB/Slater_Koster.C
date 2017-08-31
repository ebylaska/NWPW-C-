/* Slater_Koster.C -
*/

using namespace std;
#include	<iostream.h>
#include	<math.h>
#include	<string.h>
#include	"Slater_Koster.h"
#include	"Convert_Name.h"


/********************************
 *				*
 *         Constructors 	*
 *				*
 ********************************/

inline Slater_Koster::Slater_Koster(const Slater_Koster& source)
{
   int i;
   strcpy(atom1,source.atom1);
   strcpy(atom2,source.atom2);
   z1 = source.z1;
   z2 = source.z2;
   for (i=0; i<5; ++i)
   {
      f0[i] = source.f0[i];
      r0[i] = source.r0[i];
      rc[i] = source.rc[i];
      na[i] = source.na[i];
      nb[i] = source.nb[i];
      nc[i] = source.nc[i];
   }
   for (i=0; i<4; ++i)
     taa[i] = source.taa[i];
}





/********************************
 *				*
 *	     operator =		*
 *				*
 ********************************/

inline Slater_Koster& Slater_Koster::operator=(const Slater_Koster& source)
{

   if (this != &source)
   {
      int i;
      strcpy(atom1,source.atom1);
      strcpy(atom2,source.atom2);
      z1 = source.z1;
      z2 = source.z2;
      for (i=0; i<5; ++i)
      {
         f0[i] = source.f0[i];
         r0[i] = source.r0[i];
         rc[i] = source.rc[i];
         na[i] = source.na[i];
         nb[i] = source.nb[i];
         nc[i] = source.nc[i];
      }
      for (i=0; i<4; ++i)
        taa[i] = source.taa[i];
   }

   return *this;
}

/********************************
 *				*
 *	        F		*
 *				*
 ********************************/
/* returns 
     F = Slater-Koster functional(r)

*/

inline double Slater_Koster::F(const int i, const double r)
{
   double tmp;
   tmp = f0[i]*pow(r0[i]/r,na[i])
          *exp(-nb[i]*pow(r/rc[i],nc[i]) + nb[i]*pow(r0[i]/rc[i],nc[i]));

    return tmp;
}

/********************************
 *				*
 *	     dFdr		*
 *				*
 ********************************/
/* returns 
     dFdr = d/dr(Slater-Koster functional(r))

*/


inline double Slater_Koster::dFdr(const int i, const double r)
{
   double tmp;

   tmp = (-(na[i]/r) - (nb[i]*nc[i]/r)*pow(r/rc[i],nc[i]));
   tmp = tmp*
         f0[i]*pow(r0[i]/r,na[i])
         *exp(-nb[i]*pow(r/rc[i],nc[i]) + nb[i]*pow(r0[i]/rc[i],nc[i]));
   
   return tmp;
}


/********************************
 *				*
 *           Ecore		*
 *				*
 ********************************/

inline double Slater_Koster::Ecore(const double r)
{
   double tmp;
   tmp = this->F(0,r);
   return tmp;
}

/********************************
 *				*
 *	    dEcoredr		*
 *				*
 ********************************/

inline double Slater_Koster::dEcoredr(const double r)
{
   double tmp;
   tmp = this->dFdr(0,r);
   return tmp;
}


/********************************
 *				*
 *	interaction_kind	*
 *				*
 ********************************/

/* 
    Interaction labels
     0: S-S 
     1: S-P
     2: P-S
     3: P-P
*/

inline int	interaction_kind(const int alpha, const int beta)
{
   int ii;

   if ((alpha == 1) && (beta == 1))  ii = 0;
   if ((alpha == 1) && (beta >  1))  ii = 1;
   if ((alpha >  1) && (beta == 1))  ii = 2;
   if ((alpha >  1) && (beta >  1))  ii = 3;

   return ii;
} /* interaction_kind */

       
    

/********************************
 * 				*
 *            t_ab		*
 *				*
 ********************************/

inline double Slater_Koster::t_ab(const int alpha, const int beta, 
                           const Projection& p)
{
   double tmp;
   int    ik;
   double r12 = p.r;
   Vector3 sigma = p.sigma;
   Vector3 pi1   = p.pi1;
   Vector3 pi2   = p.pi2;

   ik  = interaction_kind(alpha,beta);

   /* S-S interaction */
   if (ik == 0)  
      tmp = (this->F(1,r12)); 
 
    /* S-P interaction */
    if (ik == 1)
       tmp = (sigma[beta-2])*(this->F(2,r12));

    /* P-S interaction */
    if (ik == 2)
       tmp = -(sigma[alpha-2])*(this->F(2,r12));
 
    /* P-P interaction */
    if (ik == 3)
       tmp = (sigma[alpha-2])*(sigma[beta-2])*(this->F(3,r12))
           + ((pi1[alpha-2])*(pi1[beta-2])
           +  (pi2[alpha-2])*(pi2[beta-2]))
             *(this->F(4,r12));

   return tmp;
} /* t_ab */
 
/********************************
 * 				*
 *         dt_abdr		*
 *				*
 ********************************/

inline double Slater_Koster::dt_abdr(const int alpha, const int beta, 
                              const Projection& p)
{
   double tmp;
   int    ik;
   double r12 = p.r;
   Vector3 sigma = p.sigma;
   Vector3 pi1   = p.pi1;
   Vector3 pi2   = p.pi2;


   ik  = interaction_kind(alpha,beta);

   /* S-S interaction */
   if (ik == 0)  
      tmp = (this->dFdr(1,r12)); 
 
    /* S-P interaction */
    if (ik == 1)
       tmp = (sigma[beta-2])*(this->dFdr(2,r12));
       
    /* P-S interaction */
    if (ik == 2)
       tmp = -(sigma[alpha-2])*(this->dFdr(2,r12));
 
    /* P-P interaction */
    if (ik == 3)
       tmp = (sigma[alpha-2])*(sigma[beta-2])*(this->dFdr(3,r12))
           + ((pi1[alpha-2])*(pi1[beta-2])
           +  (pi2[alpha-2])*(pi2[beta-2]))
             *(this->dFdr(4,r12));

   return tmp;

} /* dt_abdr */
 


/********************************
 *				*
 *	    io operations	*
 *				*
 ********************************/

inline ostream& operator << (ostream& s, Slater_Koster& source)
{
   
   int i;
   s  << source.Atom1() << " " << source.Atom2() << "\n";
   for (i=0; i<4; ++i)
      s << source.taa[i] << "\n";

   for (i=0; i<5; ++i)
   {
      s  << source.F0(i) << " ";
      s  << source.R0(i) << " ";
      s  << source.Rc(i) << " ";
      s  << source.Na(i) << " ";
      s  << source.Nb(i) << " ";
      s  << source.Nc(i) << "\n";
   }

   return s;
}
inline istream& operator >> (istream& s, Slater_Koster& source)
{

    int i;
    s >> source.atom1;
    s >> source.atom2;

    source.z1 = NameToCharge(source.atom1);
    source.z2 = NameToCharge(source.atom2);

    for (i=0; i<4; ++i)
       s >> source.taa[i];

    for (i=0; i<5; ++i)
    {
       s >> source.f0[i];
       s >> source.r0[i];
       s >> source.rc[i];
       s >> source.na[i];
       s >> source.nb[i];
       s >> source.nc[i];
    }
    if (source.z1 > source.z2)
    {
       int itmp  = source.z1;
       source.z1 = source.z2;
       source.z2 = itmp;

       char	stmp[5];
       strcpy(stmp,source.atom1);
       strcpy(source.atom1,source.atom2);
       strcpy(source.atom2,stmp);

     }


    return s;

}

