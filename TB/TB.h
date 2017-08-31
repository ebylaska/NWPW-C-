#ifndef	_TB_H_
#define _TB_H_
/* TB.h -
   Author - Eric Bylaska

*/

#include	"Ion.h"
#include	"TB_Eigen.h"

class	TB : public TB_Eigen {
	int		neigs;
	
	double		coreE,
			totalE;
        int		nions;
	Ion*		ions;
	Ion*		fions;

	inline void	Setup();
public:
     	/* Constructor */
	inline TB(const int, const int, const double, 
	   const int, Ion*, Ion*);

	/* Destructor */
//	~TB();

	/* reset the ions */
	inline void Reset();

        
	inline double	Energy()      {return totalE;}
	inline double	Core_Energy() {return coreE;}

	inline void Force();
};

#include	"TB.C"

#endif
