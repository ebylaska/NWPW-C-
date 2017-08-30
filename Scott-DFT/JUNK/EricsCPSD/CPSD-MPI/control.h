#ifndef _CONTROL_H_
#define _CONTROL_H_
/* control.h
   Author - Eric Bylaska

*/

#include	"rtdb.h"
extern void control_read(RTDB&);

extern int control_mapping();
extern int control_balance();
extern int control_np_orbital();
extern int control_ngrid(const int);
extern double control_unita(const int, const int);
extern double control_ecut();
extern double control_wcut();

#endif
