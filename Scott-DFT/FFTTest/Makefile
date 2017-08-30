INCLUDES = -I../include

include ../arch.mac

OBJECTS	= ffttest.o d3db.o Mapping3.o hilbert.o Parallel.o compressed_io.o util_date.o util.o 


LIBDIR1 = ../lib
LDFLAGS += -L$(LIBDIR1)

#LIBDIR2 = /Users/bylaska/lib
#LIBDIR2 = /home/bylaska/lib
#LDFLAGS += -L$(LIBDIR2)

#LIBRARY = ../lib/libfftpack.a ../BLAS/blas_LINUX.a 
#LIBRARY = ../lib/libfftpack.a ../lib/lapack_LINUX.a ../lib/blas_LINUX.a  -lgfortran
#LIBRARY = -lfftw3 -llapack -lcblas -latlas
LIBRARY = -lfftpack -llapack -lblas

ffttest.x:	$(OBJECTS)
	$(C++) -o ffttest.x $(C++FLAGS) $(OBJECTS) $(LDFLAGS) $(LIBRARY)

ffttest.o: 	ffttest.C Parallel.h util_date.h d3db.h 

d3db.o:	d3db.C d3db.h Mapping3.h util.h

Mapping3.o:	Mapping3.C Mapping3.h

hilbert.o:	hilbert.C hilbert.h olist.h

Parallel.o:	Parallel.C Parallel.h

util.o:	util.C util.h

control.o:	control.C control.h

util_date.o:	util_date.C util_date.h

compressed_io.o:	compressed_io.C compressed_io.h


clean:
	rm -f *.o ffttest.x

run:
	mpirun -np 4 ./ffttest.x 

run2:
	mpirun -np 2 ./ffttest.x

run8:
	mpirun -np 8 ./ffttest.x
