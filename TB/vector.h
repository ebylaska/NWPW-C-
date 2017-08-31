/* vector.h - 4/12/94 */


template<class T> class vector{
	T*  v;
        int sz;
public:
        /* Constructors */
  	vector() { v =((T*) 0);
                   sz=0; }
  	vector(const int);

  	~vector() { if (v) delete [] v; }

        int	  size();
        T&        operator [] (const int);
        vector<T>& operator =  (const vector<T>&);
//        friend ostream&  operator << (ostream&, vector<T>&);
};

#define MAX(x,y)	((x>y) ? x : y)
#define error(X) cout << X

template<class T>
   vector<T>::vector(const int ss)
   {
       v = new T[ss];
       sz = ss;
    }
template<class T>
   vector<T>::size()
   {
      return sz;
   }


template<class T> 
   T& vector<T>::operator[](const int i)
   {
       if (i<0 || sz <=i) error("vector: range error");
          return v[i];
   }

template<class T> 
   vector<T>& vector<T>::operator=(const vector<T>& v1)
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
  ostream& operator << (ostream& s, vector<T>& v1)
  {
       int i;
       int size = v1.size();

       s << "(";
       if (size > 0)
       {
          for (i=0; i<(size-1); ++i)
             s << v1[i] << ",";
          s << v1[size-1];
       }
       s << ")";

       return s;
   }
