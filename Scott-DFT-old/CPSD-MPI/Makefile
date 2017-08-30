include ../arch.abe

OBJECTS	= cpsd.o Pseudopotential.o Ion.o PGrid.o Balance.o lattice.o d3db.o Mapping3.o hilbert.o Parallel.o util_date.o control.o rtdb.o rtdb_seq.o hdbm.o compressed_io.o


LIBDIR = ../lib
LDFLAGS += -L$(LIBDIR)
LIBRARY = -lfftpack -llapack -lblas

cpsd.x:	$(OBJECTS)
	$(C++) -o cpsd.x $(C++FLAGS) $(OBJECTS) $(LIBRARY) $(LDFLAGS)

CPSD.o: 	PGrid.C PGrid.h Parallel.h rtdb.h control.h util_date.h d3db.h lattice.h

Pseudopotential.o: Pseudopotential.h Pseudopotential.C Ion.h Ion.C PGrid.C PGrid.h

Ion.o:		Ion.C Ion.h

PGrid.o:	PGrid.C PGrid.h

Balance.o:	Balance.C Balance.h

lattice.o:	lattice.C lattice.h

d3db.o:	d3db.C d3db.h

Mapping3.o:	Mapping3.C Mapping3.h

hilbert.o:	hilbert.C hilbert.h

Parallel.o:	Parallel.C Parallel.h

rtdb.o:	rtdb.C 

control.o:	control.C 

util_date.o:	util_date.c util_date.h

compressed_io.o:	compressed_io.c compressed_io.h

rtdb_seq.o:	rtdb_seq.c 

hdbm.o:	hdbm.c hdbm.h
#	$(C++) $(C++FLAGS) -c $<

clean:
	rm -f *.o cpsd.x

run:
	/usr/local/bin/mpiexec -n 4 ./cpsd.x

run2:
	/usr/local/bin/mpiexec -n 2 ./cpsd.x

run8:
	/usr/local/bin/mpiexec -n 8 ./cpsd.x
