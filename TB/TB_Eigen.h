#ifndef _TB_EIGEN_H_
#define _TB_EIGEN_H_
/* TB_Eigen.h -
   Author - Eric Bylaska

	This class is used for setting up and solving
a a tigt-binding.

*/
#include	<iostream.h>


class	TB_Eigen {

	int	N;
	int 	ispin,nelc;
        int	spin;
	int	*fill;
	double	*H,
		*psi,
		*eigens;
	double	Evalence,
		U,
 		Uterm;

public:
	/* Constructors */
	inline TB_Eigen();
	inline TB_Eigen(const int, const int, const int, const double);

	/* Destructor */
	inline ~TB_Eigen();

	/* routines for setting up the matrix and viewing it */
	inline void	Set_Matrix(const int, const int, const double);
	inline double	operator () (const int, const int);

        inline int Number_Eigenvalues() {return N;}
	inline double	operator () (const int);
  	inline int   	Fill(const int);
        inline int	Spin() {return spin;}

  	inline double	U_Energy();
  	inline double	Valence_Energy();

	/* routines for Solving matrix */
	inline void	Diagonalize();

	/* equality operator */
	inline TB_Eigen& operator =(TB_Eigen&);

	/* io operations */
        inline friend void	operator <<(ostream&, TB_Eigen&);
};

#include "TB_Eigen.C"
#endif
