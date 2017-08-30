/* evector.h - 4/12/94 */

#include	<iostream>
#include	<cmath>
#include	<cstdlib>
using namespace std;

template<class T> class evector{
	T*  v;
        int sz;
public:
        /* Constructors */
  	evector() { v =((T*) 0);
                   sz=0; }
  	evector(const int);

  	~evector() { if (v) delete [] v; }

        int	  size();
        T&         operator [] (const int);
        evector<T>& operator =  (const evector<T>&);
        evector<T>& operator +  (const evector<T>&, const evector<T>&);
        ostream&   harry(ostream&);
};

#define MAX(x,y)	((x>y) ? x : y)
#define error(X) cout << X

template<class T>
   evector<T>::evector(const int ss)
   {
       v = new T[ss];
       sz = ss;
    }
template<class T>
   int evector<T>::size()
   {
      return sz;
   }


template<class T> 
   T& evector<T>::operator[](const int i)
   {
       if (i<0 || sz <=i) error("evector: range error");
          return v[i];
   }

template<class T> 
   evector<T>& evector<T>::operator=(const evector<T>& v1)
   {

      if (this != &v1)
      {
         int ss = MAX(sz,(v1.sz));
         if (sz!=0)
            delete [] v;
         v  = new  T[ss];
         sz = ss;
         for (int i=0; i<ss; ++i)
            v[i] = v1.v[i];
      }
      return *this;
   }

template<class T> 
   evector<T>& evector<T>::operator + (const evector<T>& v1, const evector<T> v2)
   {

      if (this != &v1)
      {
         int ss = MAX(sz,(v1.sz));
         if (sz!=0)
            delete [] v;
         v  = new  T[ss];
         sz = ss;
         for (int i=0; i<ss; ++i)
            v[i] = v1.v[i] + v2.v[i];
      }
      return *this;
   }


   

template<class T>
  ostream& evector<T>::harry(ostream& s)
  {
       int i;
       int size = sz;

       s << "(";
       if (size > 0)
       {
          for (i=0; i<(size-1); ++i)
             s << v[i] << ",";
          s << v[size-1];
       }
       s << ")";

       return s;
   }


