
OBJECTS	= tbsd.o Interactions.o \
	xyz_io.o \
	extra.o \
	error.o seconds.o \
        message.o convert_name.o

OBJECTS_MD = tbmd.o Interactions.o \
	xyz_io.o \
	extra.o \
	error.o seconds.o \
        message.o convert_name.o


LIBDIR = /home/bylaska/lib
LIBS = -L$(LIBDIR) -lblas -lc -lm
Options = -O3
LIB_i486   = 
#LIB_sun4   = -lF77
#IB_sun4   = 
#LIB_rs6000 =
LIBARCH    = $(LIB_$(ARCH))

CC  = gcc
CPP = g++
FC  = g77



tbsd.x:	$(ARCH) $(OBJECTS)
	$(CPP) $(Options) $(OBJECTS) $(LIBS) $(LIBARCH) -o tbsd.x

tbmd.x:	$(ARCH) $(OBJECTS_MD)
	$(CPP) $(Options) $(OBJECTS_MD) $(LIBS) $(LIBARCH) -o tbmd.x

$(ARCH):
	cp eigen_$(ARCH).h eigen.h

tbsd.o:	tbsd.C
	$(CPP) -c $(Options) tbsd.C

tbmd.o:	tbmd.C
	$(CPP) -c $(Options) tbmd.C

Slater_Koster.o:	Slater_Koster.h Slater_Koster.C
	$(CPP) -c $(Options) Slater_Koster.C

Slater_Koster_List.o:	Slater_Koster_List.h Slater_Koster_List.C
	$(CPP) -c $(Options) Slater_Koster_List.C

Interactions.o:	Interactions.h Interactions.C
	$(CPP) -c $(Options) Interactions.C

xyz_io.o:	xyz_io.h xyz_io.C
	$(CPP) -c $(Options) xyz_io.C

convert_name.o:	convert_name.h convert_name.c
	$(CC) -c $(Options) convert_name.c

extra.o:	extra.f
	$(FC) -c $(Options) extra.f

TB_Eigen.o:	TB_Eigen.h TB_Eigen.C
	$(CPP) -c $(Options) TB_Eigen.C

TB.o:	TB.h TB.C
	$(CPP) -c $(Options) TB.C

error.o:	error.h error.C
	$(CPP) -c $(Options) error.C

seconds.o:	seconds.h seconds.C
	$(CPP) -c $(Options) seconds.C

message.o:	message.h message.c
	$(CC) -c $(Options) message.c




clean:
	rm -f *.o tbsd.x tbmd.x core
