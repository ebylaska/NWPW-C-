#ifndef PSIHEADER_H
#define PSIHEADER_H

#include <stdint.h>
template <typename Real>
class PsiHeader {
  int nBytes;
  int version;
  int nfft[3];
  Real unita[9];
  int ispin;
  int ne[2];
  int occupation;

  inline int _bRead(void* buf, int* pt, int n)
    {
      for(int i=0; i < n; i++) { pt[i] = (int) ((uint64_t*)buf)[i]; }
      return(n*sizeof(uint64_t)); 
    }
  inline int _bRead(void* buf, Real* pt, int n)
    {
      for(int i=0; i < n; i++) { pt[i] = (Real) ((double*)buf)[i]; }
      return(n*sizeof(double)); 
    }
  char* read(char* buf)
    {
      buf += _bRead(buf, &version,1);
      buf += _bRead(buf, nfft,3);
      buf += _bRead(buf, unita,9);
      buf += _bRead(buf, &ispin,1);
      buf += _bRead(buf, ne,2);
      buf += _bRead(buf, &occupation,1);
      return(buf);
    }

 public:
  PsiHeader() {nBytes=0;};
  PsiHeader(char* buf)
    {
      char *pt = read(buf);
      nBytes = pt-buf;
    }
  inline int get_nBytes() {return(nBytes); }
  inline int get_ispin() {return(ispin); }
  inline int get_ne(int index) {return(ne[index]); }
  inline Real get_unita(int index) {return(unita[index]); }
  inline int get_nfft(int index) {return(nfft[index]); }
  inline Real* get_unita_ptr() {return(unita); }
  void reportHeader()
  {
    cerr << "nBytes " << nBytes << endl;
    cerr << "ispin " << ispin << endl;
    cerr << "ne " ; for(int i=0; i < 2; i++) cerr << " " << ne[i]; cerr << endl;
    cerr << "unita " ; for(int i=0; i < 9; i++) cerr << " " << unita[i]; cerr << endl;
    cerr << "nfft " ; for(int i=0; i < 3; i++) cerr << " " << nfft[i]; cerr << endl;
    cerr << "occupation " << occupation << endl;
  }
};

#endif
