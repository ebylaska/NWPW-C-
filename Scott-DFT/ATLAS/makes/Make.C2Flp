include Make.inc
mySRCdir = $(SRCdir)/interfaces/lapack/C2F/src

sobj = ATL_C2Fsgels.o ATL_C2Fsormlq.o ATL_C2Fsormql.o ATL_C2Fsormqr.o \
       ATL_C2Fsormrq.o
dobj = ATL_C2Fdgels.o ATL_C2Fdormlq.o ATL_C2Fdormql.o ATL_C2Fdormqr.o \
       ATL_C2Fdormrq.o
cobj = ATL_C2Fcgels.o ATL_C2Fcunmlq.o ATL_C2Fcunmql.o ATL_C2Fcunmqr.o \
       ATL_C2Fcunmrq.o
zobj = ATL_C2Fzgels.o ATL_C2Fzunmlq.o ATL_C2Fzunmql.o ATL_C2Fzunmqr.o \
       ATL_C2Fzunmrq.o

all : lib
lib : dlib slib zlib clib
clean : dclean sclean zclean cclean iclean

#
# NOTE: only build C2F interface if F77 LAPACK exists; otherwise, we will
#       not be able to build dynamic libs due to unsatisfied calls to F77LAPACK
#
sclean :
	rm -f $(sobj) slib.grd
slib : 
	if [ -d "$(FLAdir)" ]; then \
           $(MAKE) slib.grd ; \
        fi
slib.grd : $(sobj)
	if [ -d "$(FLAdir)" ]; then \
           $(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(sobj) ; \
           $(RANLIB) $(LAPACKlib) ; \
        fi
	touch slib.grd
dclean :
	rm -f $(dobj) dlib.grd
dlib : 
	if [ -d "$(FLAdir)" ]; then \
           $(MAKE) dlib.grd ; \
        fi
dlib.grd : $(dobj)
	if [ -d "$(FLAdir)" ]; then \
           $(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(dobj) ; \
           $(RANLIB) $(LAPACKlib) ; \
        fi
	touch dlib.grd
cclean :
	rm -f $(cobj) clib.grd
clib : 
	if [ -d "$(FLAdir)" ]; then \
           $(MAKE) clib.grd ; \
        fi
clib.grd : $(cobj)
	if [ -d "$(FLAdir)" ]; then \
           $(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(cobj) ; \
           $(RANLIB) $(LAPACKlib) ; \
        fi
	touch clib.grd
zclean :
	rm -f $(zobj) zlib.grd
zlib : 
	if [ -d "$(FLAdir)" ]; then \
           $(MAKE) zlib.grd ; \
        fi
zlib.grd : $(zobj)
	if [ -d "$(FLAdir)" ]; then \
           $(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(zobj) ; \
           $(RANLIB) $(LAPACKlib) ; \
        fi
	touch zlib.grd

ATL_C2Fcunmlq.o : $(mySRCdir)/ATL_C2Funmlq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/ATL_C2Funmlq.c
ATL_C2Fcunmrq.o : $(mySRCdir)/ATL_C2Funmrq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/ATL_C2Funmrq.c
ATL_C2Fcunmql.o : $(mySRCdir)/ATL_C2Funmql.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/ATL_C2Funmql.c
ATL_C2Fcunmqr.o : $(mySRCdir)/ATL_C2Funmqr.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/ATL_C2Funmqr.c
ATL_C2Fzunmlq.o : $(mySRCdir)/ATL_C2Funmlq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/ATL_C2Funmlq.c
ATL_C2Fzunmrq.o : $(mySRCdir)/ATL_C2Funmrq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/ATL_C2Funmrq.c
ATL_C2Fzunmql.o : $(mySRCdir)/ATL_C2Funmql.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/ATL_C2Funmql.c
ATL_C2Fzunmqr.o : $(mySRCdir)/ATL_C2Funmqr.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/ATL_C2Funmqr.c
ATL_C2Fsormlq.o : $(mySRCdir)/ATL_C2Formlq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/ATL_C2Formlq.c
ATL_C2Fsormrq.o : $(mySRCdir)/ATL_C2Formrq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/ATL_C2Formrq.c
ATL_C2Fsormql.o : $(mySRCdir)/ATL_C2Formql.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/ATL_C2Formql.c
ATL_C2Fsormqr.o : $(mySRCdir)/ATL_C2Formqr.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/ATL_C2Formqr.c
ATL_C2Fdormlq.o : $(mySRCdir)/ATL_C2Formlq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/ATL_C2Formlq.c
ATL_C2Fdormrq.o : $(mySRCdir)/ATL_C2Formrq.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/ATL_C2Formrq.c
ATL_C2Fdormql.o : $(mySRCdir)/ATL_C2Formql.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/ATL_C2Formql.c
ATL_C2Fdormqr.o : $(mySRCdir)/ATL_C2Formqr.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/ATL_C2Formqr.c

ATL_C2Fsgels.o : $(mySRCdir)/ATL_C2Fgels.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/ATL_C2Fgels.c
ATL_C2Fdgels.o : $(mySRCdir)/ATL_C2Fgels.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/ATL_C2Fgels.c
ATL_C2Fcgels.o : $(mySRCdir)/ATL_C2Fgels.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/ATL_C2Fgels.c
ATL_C2Fzgels.o : $(mySRCdir)/ATL_C2Fgels.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/ATL_C2Fgels.c

