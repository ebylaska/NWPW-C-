#ifndef	_Slater_Koster_H_
#define	_Slater_Koster_H_
/* Slater_Koster.h -
*/

using namespace std;
#include	<iostream.h>
#include	<string.h>
#include	"Projection.h"

/* interaction index
        Data structure stores:

		0:core
       		1:SS
       		2:S-sigma
       		3:sigma-sigma
       		4:pi-pi 

        ++++++++++++++++++++++++++++++++
	t_ab(alpha, beta, projection) assumes

	 	alpha,beta
		1: s
		2: px
		3: py
		4: pz
*/

class Slater_Koster {

        int    z1;
        int    z2;
	char   atom1[5];
	char   atom2[5];
        
	double	taa[4];
        int	nf;     /* degree of polynomial expansions */
        double cf[10];  /* coefficients of polynomial expansion */
	double f0[5],r0[5],rc[5],
               na[5],nb[5],nc[5];


        
public:
        /* Constructors */
        inline Slater_Koster() 
        {
           int i;
	   strcpy(atom1,"X");
	   strcpy(atom2,"X"); 
           z1 = 0;   z2 = 0;
           for (i=0; i<5; ++i)
           {
             f0[i]=0.0;
             r0[i]=0.0;
             rc[i]=0.0;
             na[i]=0.0;
             nb[i]=0.0;
             nc[i]=0.0;
           }
   	   for (i=0; i<4; ++i)
             taa[i] = 0.0;
          nf = 0;
          for (i=0; i<10; ++i)
            cf[i] = 0.0;
        }
	inline Slater_Koster(const Slater_Koster&);	
        inline Slater_Koster& operator=(const Slater_Koster&);

        /* Slater-Koster function objects */
	inline double F(const int, const double);
        inline double dFdr(const int, const double);

	/* Tight Binding objects */
        inline double Ecore(const double);
        inline double dEcoredr(const double);
        inline double t_ab(const int,    const int, const Projection&);
        inline double dt_abdr(const int, const int, const Projection&);
        inline double t_aa(const int alpha) { return taa[alpha-1]; }

        /* Mostly for checking purposes */
        inline double F0(const int i) { return f0[i]; }
        inline double R0(const int i) { return r0[i]; }
        inline double Rc(const int i) { return rc[i]; }
        inline double Na(const int i) { return na[i]; }
        inline double Nb(const int i) { return nb[i]; }
        inline double Nc(const int i) { return nc[i]; }
        inline int    Z1() {return z1;}
        inline int    Z2() {return z2;}
        inline char*  Atom1() {return atom1;}
        inline char*  Atom2() {return atom2;}

       

        /* stdio operations */
        inline friend ostream& operator << (ostream&, Slater_Koster&);
        inline friend istream& operator >> (istream&, Slater_Koster&);


};
#include	"Slater_Koster.C"

#endif
