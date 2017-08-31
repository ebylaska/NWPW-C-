
#include	"Slater_Koster_List.h"


        /* Constructors */
inline Slater_Koster_List::Slater_Koster_List()
{ 
   size  = 0;
}

inline Slater_Koster_List::Slater_Koster_List(const int sz)
{
   size = sz;
  array = new Slater_Koster[size];
}

        /* Destructor */
inline Slater_Koster_List::~Slater_Koster_List()
{
   delete [] array;
}
           

inline Slater_Koster* Slater_Koster_List::operator()(const int i)
{
    return &array[i];
}

inline Slater_Koster_List& Slater_Koster_List::operator=(Slater_Koster_List& source)
{
   if (this != &source)
   {
      if (size > 0)
         delete [] array;
      size  = source.size;
      array = new Slater_Koster[size];
      for (int i=0; i<size; ++i)
         array[i] = source.array[i];
   }
   return *this;
}

inline int Slater_Koster_List::Size()
{
   return size;
}

        /* stdio operations */
inline ostream& operator << (ostream& s, Slater_Koster_List& source)
{   
   s << source.Size() << "\n";
   for (int i=0; i<source.Size(); ++i)
   {
      s << *source(i);
      s << "\n";
   }

    return s;
 }

inline istream& operator >> (istream& s, Slater_Koster_List& source)
{
     if (source.size > 0)
        delete [] source.array;

     s >> source.size;
      
     source.array = new Slater_Koster[source.size];

     for (int i=0; i<source.Size(); ++i)
       s >> *source(i);

   return s;
 }
            
        



