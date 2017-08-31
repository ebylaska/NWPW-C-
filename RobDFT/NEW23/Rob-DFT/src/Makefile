INCLUDES = -I../inc

include ../arch.gpu

OBJECTS	= cpsd.o  control.o Ion.o rtdb_seq.o hdbm.o rtdb.o

LIBDIR1 = ../lib 
LDFLAGS += -L$(LIBDIR1) -L/usr/lib64 -lcublas -lcufft
LIBRARY = -lm

cpsd.x:	$(OBJECTS) cpsd.o
	$(C++) -o cpsd.x $(C++FLAGS) $(OBJECTS) $(LDFLAGS) $(LIBRARY)

cpsd.o: cpsd.cu boxArray.h psiHeader.h CommonData.h PseudopotentialData.h Strfac.h ggmLambda.h
	$(C++) $(C++FLAGS) -c cpsd.cu
rtdb.o:	rtdb.cpp rtdb.h
	$(C++) $(C++FLAGS) -c rtdb.cpp 
control.o: control.cpp control.h
	$(C++) $(C++FLAGS) -c control.cpp 
Ion.o: 	control.cpp Ion.h
	$(C++) $(C++FLAGS) -c Ion.cpp 
rdtb_seq.o: rtdb_seq.c rtdb_seq.h
	$(CC) $(CFLAGS) -c rtdb_seq.c
hdbm.o: hdbm.c hdbm.h
	$(CC) $(CFLAGS) -c hdbm.c

clean:
	rm -f *.o cpsd.x
