include Make.inc
mySRCdir = $(SRCdir)/interfaces/lapack/F77/src

ilaobj = icmax1.o ilaver.o iparmq.o izmax1.o lsamen.o

sobj = ATL_f77wrap_sgelqf.o ATL_f77wrap_sgeqlf.o ATL_f77wrap_sgeqrf.o \
       ATL_f77wrap_sgerqf.o ATL_f77wrap_sgesv.o ATL_f77wrap_sgetnb.o \
       ATL_f77wrap_sgetrf.o ATL_f77wrap_sgetri.o ATL_f77wrap_sgetrs.o \
       ATL_f77wrap_slarfb.o ATL_f77wrap_slarft.o ATL_f77wrap_slauum.o \
       ATL_f77wrap_sposv.o ATL_f77wrap_spotrf.o ATL_f77wrap_spotrs.o \
       ATL_f77wrap_strtri.o ieeeck.o sgelqf.o sgeqlf.o sgeqrf.o sgerqf.o \
       sgesv.o sgetrf.o sgetri.o sgetrs.o slamc3.o slamch.o slarfb.o slarft.o \
       slauum.o sposv.o spotrf.o spotri.o spotrs.o strtri.o
dobj = ATL_f77wrap_dgelqf.o ATL_f77wrap_dgeqlf.o ATL_f77wrap_dgeqrf.o \
       ATL_f77wrap_dgerqf.o ATL_f77wrap_dgesv.o ATL_f77wrap_dgetnb.o \
       ATL_f77wrap_dgetrf.o ATL_f77wrap_dgetri.o ATL_f77wrap_dgetrs.o \
       ATL_f77wrap_dlarfb.o ATL_f77wrap_dlarft.o ATL_f77wrap_dlauum.o \
       ATL_f77wrap_dposv.o ATL_f77wrap_dpotrf.o ATL_f77wrap_dpotrs.o \
       ATL_f77wrap_dtrtri.o dgelqf.o dgeqlf.o dgeqrf.o dgerqf.o dgesv.o \
       dgetrf.o dgetri.o dgetrs.o dlamc3.o dlamch.o dlarfb.o dlarft.o \
       dlauum.o dposv.o dpotrf.o dpotri.o dpotrs.o dtrtri.o ieeeck.o
cobj = ATL_f77wrap_cgelqf.o ATL_f77wrap_cgeqlf.o ATL_f77wrap_cgeqrf.o \
       ATL_f77wrap_cgerqf.o ATL_f77wrap_cgesv.o ATL_f77wrap_cgetnb.o \
       ATL_f77wrap_cgetrf.o ATL_f77wrap_cgetri.o ATL_f77wrap_cgetrs.o \
       ATL_f77wrap_clarfb.o ATL_f77wrap_clarft.o ATL_f77wrap_clauum.o \
       ATL_f77wrap_cposv.o ATL_f77wrap_cpotrf.o ATL_f77wrap_cpotrs.o \
       ATL_f77wrap_ctrtri.o cgelqf.o cgeqlf.o cgeqrf.o cgerqf.o cgesv.o \
       cgetrf.o cgetri.o cgetrs.o clarfb.o clarft.o clauum.o cposv.o cpotrf.o \
       cpotri.o cpotrs.o ctrtri.o ieeeck.o
zobj = ATL_f77wrap_zgelqf.o ATL_f77wrap_zgeqlf.o ATL_f77wrap_zgeqrf.o \
       ATL_f77wrap_zgerqf.o ATL_f77wrap_zgesv.o ATL_f77wrap_zgetnb.o \
       ATL_f77wrap_zgetrf.o ATL_f77wrap_zgetri.o ATL_f77wrap_zgetrs.o \
       ATL_f77wrap_zlarfb.o ATL_f77wrap_zlarft.o ATL_f77wrap_zlauum.o \
       ATL_f77wrap_zposv.o ATL_f77wrap_zpotrf.o ATL_f77wrap_zpotrs.o \
       ATL_f77wrap_ztrtri.o ieeeck.o zgelqf.o zgeqlf.o zgeqrf.o zgerqf.o \
       zgesv.o zgetrf.o zgetri.o zgetrs.o zlarfb.o zlarft.o zlauum.o zposv.o \
       zpotrf.o zpotri.o zpotrs.o ztrtri.o

all : lib
lib : reflapack.grd dlib slib zlib clib
clean : dclean sclean zclean cclean iclean
iclean :
	rm -f $(ilaobj) ieeeck.o

sclean :
	rm -f $(sobj) slib.grd reflapack.grd
slib : slib.grd 
slib.grd : $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(sobj)
	$(RANLIB) $(LAPACKlib)
	touch slib.grd
dclean :
	rm -f $(dobj) dlib.grd reflapack.grd
dlib : dlib.grd 
dlib.grd : $(dobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(dobj)
	$(RANLIB) $(LAPACKlib)
	touch dlib.grd
cclean :
	rm -f $(cobj) clib.grd reflapack.grd
clib : clib.grd 
clib.grd : $(cobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(cobj)
	$(RANLIB) $(LAPACKlib)
	touch clib.grd
zclean :
	rm -f $(zobj) zlib.grd reflapack.grd
zlib : zlib.grd 
zlib.grd : $(zobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(zobj)
	$(RANLIB) $(LAPACKlib)
	touch zlib.grd

spotri.o : $(mySRCdir)/spotri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/spotri.f
sgetri.o : $(mySRCdir)/sgetri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgetri.f
ATL_f77wrap_sgetri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
slauum.o : $(mySRCdir)/slauum.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/slauum.f
ATL_f77wrap_slauum.o : $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
strtri.o : $(mySRCdir)/strtri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/strtri.f
ATL_f77wrap_strtri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
sgetrs.o : $(mySRCdir)/sgetrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgetrs.f
ATL_f77wrap_sgetrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
sgetrf.o : $(mySRCdir)/sgetrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgetrf.f
ATL_f77wrap_sgetrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
sgesv.o : $(mySRCdir)/sgesv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgesv.f
ATL_f77wrap_sgesv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
spotrs.o : $(mySRCdir)/spotrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/spotrs.f
ATL_f77wrap_spotrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
spotrf.o : $(mySRCdir)/spotrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/spotrf.f
ATL_f77wrap_spotrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
sposv.o : $(mySRCdir)/sposv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sposv.f
ATL_f77wrap_sposv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
slarft.o : $(mySRCdir)/slarft.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/slarft.f
ATL_f77wrap_slarft.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
slarfb.o : $(mySRCdir)/slarfb.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/slarfb.f
ATL_f77wrap_slarfb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
sgelqf.o : $(mySRCdir)/sgelqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgelqf.f
ATL_f77wrap_sgelqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
sgerqf.o : $(mySRCdir)/sgerqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgerqf.f
ATL_f77wrap_sgerqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
sgeqlf.o : $(mySRCdir)/sgeqlf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgeqlf.f
ATL_f77wrap_sgeqlf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
sgeqrf.o : $(mySRCdir)/sgeqrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/sgeqrf.f
ATL_f77wrap_sgeqrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
ATL_f77wrap_sgetnb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c

dpotri.o : $(mySRCdir)/dpotri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dpotri.f
dgetri.o : $(mySRCdir)/dgetri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgetri.f
ATL_f77wrap_dgetri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
dlauum.o : $(mySRCdir)/dlauum.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dlauum.f
ATL_f77wrap_dlauum.o : $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
dtrtri.o : $(mySRCdir)/dtrtri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dtrtri.f
ATL_f77wrap_dtrtri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
dgetrs.o : $(mySRCdir)/dgetrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgetrs.f
ATL_f77wrap_dgetrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
dgetrf.o : $(mySRCdir)/dgetrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgetrf.f
ATL_f77wrap_dgetrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
dgesv.o : $(mySRCdir)/dgesv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgesv.f
ATL_f77wrap_dgesv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
dpotrs.o : $(mySRCdir)/dpotrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dpotrs.f
ATL_f77wrap_dpotrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
dpotrf.o : $(mySRCdir)/dpotrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dpotrf.f
ATL_f77wrap_dpotrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
dposv.o : $(mySRCdir)/dposv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dposv.f
ATL_f77wrap_dposv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
dlarft.o : $(mySRCdir)/dlarft.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dlarft.f
ATL_f77wrap_dlarft.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
dlarfb.o : $(mySRCdir)/dlarfb.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dlarfb.f
ATL_f77wrap_dlarfb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
dgelqf.o : $(mySRCdir)/dgelqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgelqf.f
ATL_f77wrap_dgelqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
dgerqf.o : $(mySRCdir)/dgerqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgerqf.f
ATL_f77wrap_dgerqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
dgeqlf.o : $(mySRCdir)/dgeqlf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgeqlf.f
ATL_f77wrap_dgeqlf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
dgeqrf.o : $(mySRCdir)/dgeqrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/dgeqrf.f
ATL_f77wrap_dgeqrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
ATL_f77wrap_dgetnb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c

cpotri.o : $(mySRCdir)/cpotri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cpotri.f
cgetri.o : $(mySRCdir)/cgetri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgetri.f
ATL_f77wrap_cgetri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
clauum.o : $(mySRCdir)/clauum.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/clauum.f
ATL_f77wrap_clauum.o : $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
ctrtri.o : $(mySRCdir)/ctrtri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/ctrtri.f
ATL_f77wrap_ctrtri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
cgetrs.o : $(mySRCdir)/cgetrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgetrs.f
ATL_f77wrap_cgetrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
cgetrf.o : $(mySRCdir)/cgetrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgetrf.f
ATL_f77wrap_cgetrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
cgesv.o : $(mySRCdir)/cgesv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgesv.f
ATL_f77wrap_cgesv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
cpotrs.o : $(mySRCdir)/cpotrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cpotrs.f
ATL_f77wrap_cpotrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
cpotrf.o : $(mySRCdir)/cpotrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cpotrf.f
ATL_f77wrap_cpotrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
cposv.o : $(mySRCdir)/cposv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cposv.f
ATL_f77wrap_cposv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
clarft.o : $(mySRCdir)/clarft.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/clarft.f
ATL_f77wrap_clarft.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
clarfb.o : $(mySRCdir)/clarfb.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/clarfb.f
ATL_f77wrap_clarfb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
cgelqf.o : $(mySRCdir)/cgelqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgelqf.f
ATL_f77wrap_cgelqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
cgerqf.o : $(mySRCdir)/cgerqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgerqf.f
ATL_f77wrap_cgerqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
cgeqlf.o : $(mySRCdir)/cgeqlf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgeqlf.f
ATL_f77wrap_cgeqlf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
cgeqrf.o : $(mySRCdir)/cgeqrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/cgeqrf.f
ATL_f77wrap_cgeqrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
ATL_f77wrap_cgetnb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c

zpotri.o : $(mySRCdir)/zpotri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zpotri.f
zgetri.o : $(mySRCdir)/zgetri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgetri.f
ATL_f77wrap_zgetri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getri.c
zlauum.o : $(mySRCdir)/zlauum.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zlauum.f
ATL_f77wrap_zlauum.o : $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_lauum.c
ztrtri.o : $(mySRCdir)/ztrtri.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/ztrtri.f
ATL_f77wrap_ztrtri.o : $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_trtri.c
zgetrs.o : $(mySRCdir)/zgetrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgetrs.f
ATL_f77wrap_zgetrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getrs.c
zgetrf.o : $(mySRCdir)/zgetrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgetrf.f
ATL_f77wrap_zgetrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getrf.c
zgesv.o : $(mySRCdir)/zgesv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgesv.f
ATL_f77wrap_zgesv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gesv.c
zpotrs.o : $(mySRCdir)/zpotrs.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zpotrs.f
ATL_f77wrap_zpotrs.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_potrs.c
zpotrf.o : $(mySRCdir)/zpotrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zpotrf.f
ATL_f77wrap_zpotrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_potrf.c
zposv.o : $(mySRCdir)/zposv.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zposv.f
ATL_f77wrap_zposv.o : $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_posv.c
zlarft.o : $(mySRCdir)/zlarft.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zlarft.f
ATL_f77wrap_zlarft.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_larft.c
zlarfb.o : $(mySRCdir)/zlarfb.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zlarfb.f
ATL_f77wrap_zlarfb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_larfb.c
zgelqf.o : $(mySRCdir)/zgelqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgelqf.f
ATL_f77wrap_zgelqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gelqf.c
zgerqf.o : $(mySRCdir)/zgerqf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgerqf.f
ATL_f77wrap_zgerqf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_gerqf.c
zgeqlf.o : $(mySRCdir)/zgeqlf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgeqlf.f
ATL_f77wrap_zgeqlf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_geqlf.c
zgeqrf.o : $(mySRCdir)/zgeqrf.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/zgeqrf.f
ATL_f77wrap_zgeqrf.o : $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_geqrf.c
ATL_f77wrap_zgetnb.o : $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDCPLX $(mySRCdir)/f77wrap/ATL_f77wrap_getnb.c

slamch.o : $(mySRCdir)/lamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DSREAL $(mySRCdir)/lamch.c
slamc3.o : $(mySRCdir)/slamc3.f
	$(F77) $(F77FLAGS) $(F77NOOPT) -o $@ -c $(mySRCdir)/slamc3.f
dlamch.o : $(mySRCdir)/lamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c -DDREAL $(mySRCdir)/lamch.c
dlamc3.o : $(mySRCdir)/dlamc3.f
	$(F77) $(F77FLAGS) $(F77NOOPT) -o $@ -c $(mySRCdir)/dlamc3.f


ieeeck.o : $(mySRCdir)/ieeeck.f
	$(F77) -o $@ -c $(F77FLAGS) $(mySRCdir)/ieeeck.f

reflapack.grd : 
	if [ -d "$(FLAdir)" ]; then \
           $(MAKE) reflabuild ; \
        fi
	touch reflapack.grd
reflabuild : $(FLAdir)/lapack_$(ARCH).a
	mkdir flapack 
	cd flapack ; $(ARCHIVER) x $(FLAdir)/lapack_$(ARCH).a
	cd flapack ; rm -f ?posv.o ?potrf.o ?potrs.o ?gesv.o ?getrf.o ?getrs.o \
           ?potri.o ?lauum.o ?getri.o ?trtri.o ieeeck.o ilaenv.o [s,d]lamch.o \
           ?geqrf.o ?gerqf.o ?geqlf.o ?gelqf.o ?larfb.o ?larft.o
	cd flapack ; $(ARCHIVER) $(ARFLAGS) $(LAPACKlib) *.o
	rm -rf flapack
	$(RANLIB) $(LAPACKlib)
makeincex = $(FLAdir)/make.inc.example
$(FLAdir)/make.inc : $(makeincex)
	rm -f $(FLAdir)/make.inc
	cat $(makeincex) | sed -e 's@PLAT *=.*@PLAT = _$(ARCH)@' \
         -e 's@^FORTRAN *=.*@FORTRAN = $(F77)@' \
         -e 's@^OPTS *=.*@OPTS = $(F77FLAGS)@' \
         -e 's@^NOOPT *=.*@NOOPT = $(F77FLAGS) $(F77NOOPT)@' \
         -e 's@^LOADER *=.*@LOADER = $(F77)@' \
         -e 's@^LOADOPTS *=.*@LOADOPTS = $(F77FLAGS)@' \
         -e 's@^TIMER *=.*@TIMER = NONE@' \
         > $(FLAdir)/make.inc

$(FLAdir)/lapack_$(ARCH).a : $(FLAdir)/make.inc
	cd $(FLAdir) ; $(PMAKE) lib
