#define EIGSOLVER	eigen_

extern "C" {
void eigen_(const int*, const int*,
	   const double*, const double*, const double*, 
           const int*);
}

