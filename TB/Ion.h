#ifndef	_ION_H_
#define _ION_H_

using namespace std;

#include	<iostream.h>
#include	"Vector3.h"


#define	inline	

class	Ion : public Vector3 {
	char		name[5];
	int   		charge;
	double		mass;

public:
        /* Constructors */
        inline Ion();
        inline Ion(const char*, const Vector3);
        inline Ion(const char*, const double, const double, const double);
        inline Ion(const char*, const double*);
        inline Ion(Ion&);

        inline int Charge();
	inline double Mass();
        inline char*  Name();
        inline Ion& operator =(Vector3);

           
        /* io operations */
	inline friend ostream& operator <<(ostream&, Ion);
	inline friend istream& operator >>(istream&, Ion&);
};

#undef inline
#include	"Ion.C"

#endif
