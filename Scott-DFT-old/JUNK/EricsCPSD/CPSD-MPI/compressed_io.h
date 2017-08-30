#ifndef _COMPRESSED_IO_H_
#define _COMPRESSED_IO_H_
/* compressed_io.h -
   Author - Eric Bylaska

*/

extern void cwrite(const int, char *, const int);
extern void cread(const int, char *, const int);
extern void iwrite(const int, const int *, const int);
extern void iread(const int, int *, const int);
extern void dwrite(const int, const double *, const int);
extern void dread(const int, double *, const int);
extern void openfile(const int, char *, char *, const int);
extern void closefile(const int);

#endif