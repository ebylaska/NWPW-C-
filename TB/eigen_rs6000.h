#define EIGSOLVER	eigen

extern "C" {
void eigen(const int*, const int*,
	   const double*, const double*, const double*, 
           const int*);
}

