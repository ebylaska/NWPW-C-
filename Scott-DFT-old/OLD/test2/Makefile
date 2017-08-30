
OBJECTS	= d3db.o Mapping3.o hilbert.o Parallel.o compressed_io.o

Options = -O3

CC  = mpicc
CPP = mpicxx



test.x:	$(OBJECTS)
	$(CPP) $(Options) test.C $(OBJECTS) -o test.x


d3db.o:	d3db.C d3db.h
	$(CPP) -c $(Options) d3db.C

Mapping3.o:	Mapping3.C Mapping3.h
	$(CPP) -c $(Options) Mapping3.C

hilbert.o:	hilbert.C hilbert.h
	$(CPP) -c $(Options) hilbert.C

Parallel.o:	Parallel.C Parallel.h
	$(CPP) -c $(Options) Parallel.C

compressed_io.o:	compressed_io.c compressed_io.h
	$(CPP) -c $(Options) compressed_io.c


clean:
	rm -f *.o test.x
