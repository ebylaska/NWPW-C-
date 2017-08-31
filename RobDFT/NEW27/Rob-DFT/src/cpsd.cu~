#include <iostream>
using namespace std;
#include "boxArray.h"

int main(int argc, char* argv[])
{
  if(argc < 2) {
    cerr << "Use: rootFilename\n" << endl;
    return(1);
  }

  //boxArray<float, float2> ba(argv[1]);
  boxArray<double, double2> ba(argv[1]);
  ba.inner_loop();
  
  return(0);
}

