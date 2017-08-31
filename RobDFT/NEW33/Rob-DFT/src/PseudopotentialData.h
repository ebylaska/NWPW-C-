#ifndef PSEUDOPOTENTIAL_DATA_H
#define PSEUDOPOTENTIAL_DATA_H

#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
using namespace std;

#include <thrust/host_vector.h>
#include <vector>
#include <string>

template <typename Real>
class PseudoPotentialData {
 public:
  int nn;
  int npsp;
  Ion *myIonPtr;
  int nprj_max;
  int nfft[3];
  Real unita[9];
  char atom[2];
  int version;
  thrust::host_vector<int> psp_type,nprj;
  thrust::host_vector<Real> zv,amass,rcore,semicore;
  thrust::host_vector<int> lmax,lmmax, locp,nmax;
  vector< thrust::host_vector<Real> > rc;
  vector< thrust::host_vector<int> > n_projector;
  vector< thrust::host_vector<int> > l_projector;
  vector< thrust::host_vector<int> > m_projector;
  vector< thrust::host_vector<Real> > h_Gijl;
  vector< thrust::host_vector<Real> > vl;

  vector< thrust::host_vector<Real> > h_vnl;
  vector< thrust::device_vector<Real> > d_vnl;

  // device versions
  vector< thrust::device_vector<Real> > d_Gijl;
  
  vector<string> comment;

  void errorExit(char* s) { cerr << s << endl; exit(1);}

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
  inline int _bRead(void* buf, char* pt, int n)
    {
      for(int i=0; i < n; i++) { pt[i] = (char) ((char*)buf)[i]; }
      return(n*sizeof(char)); 
    }

  PseudoPotentialData() {};
  PseudoPotentialData(char* rootFilename, Ion* myIon, int nfft3d)
    {
      myIonPtr = myIon;
      npsp=myIon->nkatm;
      nprj_max = 0;
      psp_type = thrust::host_vector<int>(npsp);
      amass = thrust::host_vector<Real>(npsp);
      zv = thrust::host_vector<Real>(npsp);
      lmax = thrust::host_vector<int>(npsp);
      lmmax = thrust::host_vector<int>(npsp);
      locp = thrust::host_vector<int>(npsp);
      nmax = thrust::host_vector<int>(npsp);
      nprj = thrust::host_vector<int>(npsp);
      rcore = thrust::host_vector<Real>(npsp);
      semicore = thrust::host_vector<Real>(npsp);

      for(int ia=0; ia < npsp; ia++) {
	int fd;
	struct stat statbuf;
	cerr << myIon->atom(ia) << endl; 
	string filename = rootFilename + string(".") + string(myIon->atom(ia)) + ".vpp";
	cerr << "PseudoPotential filename " << filename << endl;
	if( (fd=open(filename.c_str(), O_RDONLY)) < 0) errorExit((char*)"cannot open file");
	if(fstat(fd,&statbuf) < 0) errorExit((char*)"cannot stat file");
	char *mmapBuf = (char *) mmap(0,statbuf.st_size,PROT_READ,MAP_SHARED,fd,0);
	char *buf = mmapBuf;

	buf += 80; // step past comment

	buf += _bRead(buf,&psp_type[ia],1);
	buf += _bRead(buf,&version,1);
	buf += _bRead(buf,nfft,3);
	buf += _bRead(buf,unita,9);
	buf += _bRead(buf,atom,2);
	cerr << atom[0] << atom[1] << " atom " << endl;
	buf += _bRead(buf,&amass[ia],1);
	cerr << amass[ia] << " amass " << endl;
	buf += _bRead(buf,&zv[ia],1);
	buf += _bRead(buf,&lmax[ia],1);
	buf += _bRead(buf,&locp[ia],1);
	buf += _bRead(buf,&nmax[ia],1);

	lmmax[ia] =((lmax[ia])+1)*((lmax[ia])+1) - (2*(locp[ia])+1);
	
	rc.push_back(thrust::host_vector<Real>(lmax[ia]+1));
	buf += _bRead(buf,&rc[ia][0], lmax[ia]+1); 
	buf += _bRead(buf,&nprj[ia], 1); 
	cerr << nprj[ia] << " nprj " << endl;
	if (nprj[ia]>nprj_max) nprj_max = nprj[ia];

	if (nprj[ia] > 0)  {
	  n_projector.push_back(thrust::host_vector<int>(nprj[ia]));
	  buf += _bRead(buf,&n_projector[ia][0], nprj[ia]); 
	  l_projector.push_back(thrust::host_vector<int>(nprj[ia]));
	  buf += _bRead(buf,&l_projector[ia][0], nprj[ia]); 
	  m_projector.push_back(thrust::host_vector<int>(nprj[ia]));
	  buf += _bRead(buf,&m_projector[ia][0], nprj[ia]); 
	}
	nn = (nmax[ia])*(nmax[ia])*(lmax[ia]+1);
	h_Gijl.push_back(thrust::host_vector<Real>(nn));
	buf += _bRead(buf,&h_Gijl[ia][0], nn); 
	
	buf += _bRead(buf,&rcore[ia],1);
	if (rcore[ia] > 0.0) semicore[ia] = 1.;
	else semicore[ia] = 0.;

	if (nprj[ia] > 0) {
	  vl.push_back(thrust::host_vector<Real>(nfft3d));
	  buf += _bRead(buf,&vl[ia][0], nfft3d); 
          cerr << " nfft3d=" << nfft3d << endl;
          cerr << " vl =" << vl[ia][0] << " " << vl[ia][1] << endl;
	  
	  h_vnl.push_back(thrust::host_vector<Real>(nprj[ia]*nfft3d));
	  buf += _bRead(buf,&h_vnl[ia][0], nprj[ia]*nfft3d); 
          cerr << " vnl1 =" << h_vnl[ia][0] << " " << h_vnl[ia][1] << endl;
          cerr << " vnl2 =" << h_vnl[ia][nfft3d] << " " << h_vnl[ia][nfft3d+1] << endl;
          cerr << " vnl3 =" << h_vnl[ia][2*nfft3d] << " " << h_vnl[ia][2*nfft3d+1] << endl;
          cerr << " vnl4 =" << h_vnl[ia][3*nfft3d] << " " << h_vnl[ia][3*nfft3d+1] << endl;
	}
	if (semicore[ia] > 0) {
	  cerr << "Semicore not implemented!" << endl;
	  exit(1);
	}

	munmap(mmapBuf,statbuf.st_size);
	close(fd);
	// allocate and transfer to GPU
	for(int i=0; i < h_Gijl.size(); i++) d_Gijl.push_back(thrust::device_vector<Real>(h_Gijl[i]));
	for(int i=0; i < h_vnl.size(); i++) d_vnl.push_back(thrust::device_vector<Real>(h_vnl[i]));
      }
    }
  ~PseudoPotentialData()
    {
    }
  void reportHeader()
  {
    cerr << "npsp " << npsp << endl;
  }
};

#endif
