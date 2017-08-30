/* d1db.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/

/*
#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>
using namespace std;
*/

#include	"Parallel.h"
#include	"Mapping1.h"
#include	"d1db.h"


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

d1db::d1db(Parallel *inparall,const int inmaptype, const int inispin, int *inne)
   : Mapping1(inmaptype,inparall->np_j(),inparall->taskid_j(),inispin,inne)
{
   parall = inparall;
}

