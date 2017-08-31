#ifndef	_Slater_Koster_List_H_
#define _Slater_Koster_List_H_

#include	<iostream.h>
#include	"Slater_Koster.h"

class	Slater_Koster_List {

        int		size;
	Slater_Koster*	array;

public:
        /* Constructors */
        inline Slater_Koster_List();
	inline Slater_Koster_List(const int);

        /* Destructor */
        inline ~Slater_Koster_List();
           

        inline Slater_Koster* operator()(const int);
        inline Slater_Koster_List& operator=(Slater_Koster_List&);
        inline int Size();

        /* stdio operations */
        inline friend ostream& operator << (ostream&, Slater_Koster_List&);
        inline friend istream& operator >> (istream&, Slater_Koster_List&);
            

};
#include	"Slater_Koster_List.C"

#endif
