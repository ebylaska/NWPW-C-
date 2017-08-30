include Make.inc
mySRCdir = $(SRCdir)/src/threads/blas/level3/omp

all : lib
lib : ilib dlib slib zlib clib
double    : dlib
single    : slib
complex   : clib
complex16 : zlib
clean : iclean dclean sclean cclean zclean

deps = $(INCSdir)/atlas_misc.h $(INCAdir)/atlas_type.h \
       $(INCSdir)/atlas_lvl3.h $(INCSdir)/atlas_tlvl3.h \
       $(INCAdir)/atlas_pthreads.h $(INCSdir)/atlas_threads.h

iobj = ATL_Xtompsyrk.o

iclean:
	rm -f $(iobj) ilib.grd
ilib : ilib.grd
ilib.grd : $(iobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(iobj)
	$(RANLIB) $(ATLASlib)
	touch ilib.grd

ATL_Xtompsyrk.o : $(deps) $(mySRCdir)/ATL_Xtompsyrk.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/ATL_Xtompsyrk.c



sINCdep = $(deps) $(INCAdir)/smm.h

sobj = ATL_stomptrmm.o ATL_stomptrsm.o ATL_stompsymm.o ATL_stompsyrk.o \
       ATL_stompsyr2k.o ATL_stompgemm.o

slib : slib.grd ilib 
slib.grd : $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(sobj)
	$(RANLIB) $(ATLASlib)
	touch slib.grd
sclean:
	rm -f $(sobj) slib.grd

ATL_stompgemm.o : $(deps) $(mySRCdir)/ATL_tompgemm.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tompgemm.c
ATL_stompsyr2k.o : $(deps) $(mySRCdir)/ATL_tompsyr2k.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tompsyr2k.c
ATL_stompsyrk.o : $(deps) $(mySRCdir)/ATL_tompsyrk.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tompsyrk.c
ATL_stompsymm.o : $(deps) $(mySRCdir)/ATL_tompsymm.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tompsymm.c
ATL_stomptrsm.o : $(deps) $(mySRCdir)/ATL_tomptrsm.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tomptrsm.c
ATL_stomptrmm.o : $(deps) $(mySRCdir)/ATL_tomptrmm.c
	$(ICC) $(ICCFLAGS) -DSREAL -o $@ -c $(mySRCdir)/ATL_tomptrmm.c
dINCdep = $(deps) $(INCAdir)/dmm.h

dobj = ATL_dtomptrmm.o ATL_dtomptrsm.o ATL_dtompsymm.o ATL_dtompsyrk.o \
       ATL_dtompsyr2k.o ATL_dtompgemm.o

dlib : dlib.grd ilib 
dlib.grd : $(dobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(dobj)
	$(RANLIB) $(ATLASlib)
	touch dlib.grd
dclean:
	rm -f $(dobj) dlib.grd

ATL_dtompgemm.o : $(deps) $(mySRCdir)/ATL_tompgemm.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tompgemm.c
ATL_dtompsyr2k.o : $(deps) $(mySRCdir)/ATL_tompsyr2k.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tompsyr2k.c
ATL_dtompsyrk.o : $(deps) $(mySRCdir)/ATL_tompsyrk.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tompsyrk.c
ATL_dtompsymm.o : $(deps) $(mySRCdir)/ATL_tompsymm.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tompsymm.c
ATL_dtomptrsm.o : $(deps) $(mySRCdir)/ATL_tomptrsm.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tomptrsm.c
ATL_dtomptrmm.o : $(deps) $(mySRCdir)/ATL_tomptrmm.c
	$(ICC) $(ICCFLAGS) -DDREAL -o $@ -c $(mySRCdir)/ATL_tomptrmm.c
cINCdep = $(deps) $(INCAdir)/cmm.h

cobj = ATL_ctomptrmm.o ATL_ctomptrsm.o ATL_ctompsymm.o ATL_ctompsyrk.o \
       ATL_ctompsyr2k.o ATL_ctompgemm.o ATL_ctompherk.o ATL_ctompher2k.o \
       ATL_ctomphemm.o

clib : clib.grd ilib 
clib.grd : $(cobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(cobj)
	$(RANLIB) $(ATLASlib)
	touch clib.grd
cclean:
	rm -f $(cobj) clib.grd

ATL_ctomphemm.o : $(deps) $(mySRCdir)/ATL_tomphemm.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tomphemm.c
ATL_ctompher2k.o : $(deps) $(mySRCdir)/ATL_tompher2k.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompher2k.c
ATL_ctompherk.o : $(deps) $(mySRCdir)/ATL_tompherk.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompherk.c
ATL_ctompgemm.o : $(deps) $(mySRCdir)/ATL_tompgemm.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompgemm.c
ATL_ctompsyr2k.o : $(deps) $(mySRCdir)/ATL_tompsyr2k.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompsyr2k.c
ATL_ctompsyrk.o : $(deps) $(mySRCdir)/ATL_tompsyrk.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompsyrk.c
ATL_ctompsymm.o : $(deps) $(mySRCdir)/ATL_tompsymm.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tompsymm.c
ATL_ctomptrsm.o : $(deps) $(mySRCdir)/ATL_tomptrsm.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tomptrsm.c
ATL_ctomptrmm.o : $(deps) $(mySRCdir)/ATL_tomptrmm.c
	$(ICC) $(ICCFLAGS) -DSCPLX -o $@ -c $(mySRCdir)/ATL_tomptrmm.c
zINCdep = $(deps) $(INCAdir)/zmm.h

zobj = ATL_ztomptrmm.o ATL_ztomptrsm.o ATL_ztompsymm.o ATL_ztompsyrk.o \
       ATL_ztompsyr2k.o ATL_ztompgemm.o ATL_ztompherk.o ATL_ztompher2k.o \
       ATL_ztomphemm.o

zlib : zlib.grd ilib 
zlib.grd : $(zobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(zobj)
	$(RANLIB) $(ATLASlib)
	touch zlib.grd
zclean:
	rm -f $(zobj) zlib.grd

ATL_ztomphemm.o : $(deps) $(mySRCdir)/ATL_tomphemm.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tomphemm.c
ATL_ztompher2k.o : $(deps) $(mySRCdir)/ATL_tompher2k.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompher2k.c
ATL_ztompherk.o : $(deps) $(mySRCdir)/ATL_tompherk.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompherk.c
ATL_ztompgemm.o : $(deps) $(mySRCdir)/ATL_tompgemm.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompgemm.c
ATL_ztompsyr2k.o : $(deps) $(mySRCdir)/ATL_tompsyr2k.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompsyr2k.c
ATL_ztompsyrk.o : $(deps) $(mySRCdir)/ATL_tompsyrk.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompsyrk.c
ATL_ztompsymm.o : $(deps) $(mySRCdir)/ATL_tompsymm.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tompsymm.c
ATL_ztomptrsm.o : $(deps) $(mySRCdir)/ATL_tomptrsm.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tomptrsm.c
ATL_ztomptrmm.o : $(deps) $(mySRCdir)/ATL_tomptrmm.c
	$(ICC) $(ICCFLAGS) -DDCPLX -o $@ -c $(mySRCdir)/ATL_tomptrmm.c
