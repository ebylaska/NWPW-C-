include Make.inc

mySRCdir = $(SRCdir)/src/threads
dep = $(INCSdir)/atlas_misc.h $(INCAdir)/atlas_pthreads.h \
      $(INCSdir)/atlas_threads.h

obj = ATL_thread_start.o ATL_thread_exit.o ATL_thread_join.o \
      ATL_log2tlaunch.o

all : lib
lib : lib.grd l3thr l3pt
slib3 : lib.grd $(sextthr)
	cd blas/level3 ; $(MAKE) slib
	cd $(GMMdir) ; $(MAKE) stlib
somp:
	cd blas/level3/omp ; $(MAKE) slib
spt:
	cd $(BLDdir)/src/pthreads/blas/level3/ ; $(MAKE) slib
	cd $(BLDdir)/src/pthreads/misc/ ; $(MAKE) lib
clib3 : lib.grd $(cextthr)
	cd blas/level3 ; $(MAKE) clib
	cd $(GMMdir) ; $(MAKE) ctlib
comp:
	cd blas/level3/omp ; $(MAKE) clib
cpt:
	cd $(BLDdir)/src/pthreads/blas/level3/ ; $(MAKE) clib
	cd $(BLDdir)/src/pthreads/misc/ ; $(MAKE) lib
dlib3 : lib.grd $(dextthr)
	cd blas/level3 ; $(MAKE) dlib
	cd $(GMMdir) ; $(MAKE) dtlib
domp:
	cd blas/level3/omp ; $(MAKE) dlib
dpt:
	cd $(BLDdir)/src/pthreads/blas/level3/ ; $(MAKE) dlib
	cd $(BLDdir)/src/pthreads/misc/ ; $(MAKE) lib
zlib3 : lib.grd $(zextthr)
	cd blas/level3 ; $(MAKE) zlib
	cd $(GMMdir) ; $(MAKE) ztlib
zomp:
	cd blas/level3/omp ; $(MAKE) zlib
zpt:
	cd $(BLDdir)/src/pthreads/blas/level3/ ; $(MAKE) zlib
	cd $(BLDdir)/src/pthreads/misc/ ; $(MAKE) lib

lib.grd : $(obj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(obj)
	$(RANLIB) $(ATLASlib)
	touch lib.grd
l3thr: dlib3 zlib3 clib3 slib3
l3pt:
	cd $(BLDdir)/src/pthreads/blas/level3/ ; $(MAKE) lib
	cd $(BLDdir)/src/pthreads/misc/ ; $(MAKE) lib

ATL_log2tlaunch.o : $(dep) $(mySRCdir)/ATL_log2tlaunch.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/ATL_log2tlaunch.c
ATL_thread_join.o : $(dep) $(mySRCdir)/ATL_thread_join.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/ATL_thread_join.c
ATL_thread_exit.o : $(dep) $(mySRCdir)/ATL_thread_exit.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/ATL_thread_exit.c
ATL_thread_start.o : $(dep) $(mySRCdir)/ATL_thread_start.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/ATL_thread_start.c
