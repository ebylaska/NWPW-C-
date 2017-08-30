include Make.inc
mySRCdir = $(SRCdir)/interfaces/lapack/C/src

sobj = clapack_sgelqf.o clapack_sgeqlf.o clapack_sgeqrf.o clapack_sgerqf.o \
       clapack_sgesv.o clapack_sgetrf.o clapack_sgetri.o clapack_sgetrs.o \
       clapack_slamch.o clapack_slauum.o clapack_sposv.o clapack_spotrf.o \
       clapack_spotri.o clapack_spotrs.o clapack_strtri.o

dobj = clapack_dgelqf.o clapack_dgeqlf.o clapack_dgeqrf.o clapack_dgerqf.o \
       clapack_dgesv.o clapack_dgetrf.o clapack_dgetri.o clapack_dgetrs.o \
       clapack_dlamch.o clapack_dlauum.o clapack_dposv.o clapack_dpotrf.o \
       clapack_dpotri.o clapack_dpotrs.o clapack_dtrtri.o

cobj = clapack_cgelqf.o clapack_cgeqlf.o clapack_cgeqrf.o clapack_cgerqf.o \
       clapack_cgesv.o clapack_cgetrf.o clapack_cgetri.o clapack_cgetrs.o \
       clapack_clauum.o clapack_cposv.o clapack_cpotrf.o clapack_cpotri.o \
       clapack_cpotrs.o clapack_ctrtri.o

zobj = clapack_zgelqf.o clapack_zgeqlf.o clapack_zgeqrf.o clapack_zgerqf.o \
       clapack_zgesv.o clapack_zgetrf.o clapack_zgetri.o clapack_zgetrs.o \
       clapack_zlauum.o clapack_zposv.o clapack_zpotrf.o clapack_zpotri.o \
       clapack_zpotrs.o clapack_ztrtri.o

all : lib
lib : dlib slib zlib clib
clean : dclean sclean zclean cclean

sclean :
	rm -f $(sobj) slib.grd
slib :
	$(PMAKE) slib.grd
slib.grd : $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(sobj)
	$(RANLIB) $(LAPACKlib)
	touch slib.grd

dclean :
	rm -f $(dobj) dlib.grd
dlib :
	$(PMAKE) dlib.grd
dlib.grd : $(dobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(dobj)
	$(RANLIB) $(LAPACKlib)
	touch dlib.grd

cclean :
	rm -f $(cobj) clib.grd
clib :
	$(PMAKE) clib.grd
clib.grd : $(cobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(cobj)
	$(RANLIB) $(LAPACKlib)
	touch clib.grd

zclean :
	rm -f $(zobj) zlib.grd
zlib :
	$(PMAKE) zlib.grd
zlib.grd : $(zobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(zobj)
	$(RANLIB) $(LAPACKlib)
	touch zlib.grd


clapack_sgetri.o : $(mySRCdir)/clapack_sgetri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgetri.c
clapack_spotri.o : $(mySRCdir)/clapack_spotri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_spotri.c
clapack_slauum.o : $(mySRCdir)/clapack_slauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_slauum.c
clapack_strtri.o : $(mySRCdir)/clapack_strtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_strtri.c
clapack_sgetrs.o : $(mySRCdir)/clapack_sgetrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgetrs.c
clapack_sgetrf.o : $(mySRCdir)/clapack_sgetrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgetrf.c
clapack_sgesv.o : $(mySRCdir)/clapack_sgesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgesv.c
clapack_spotrs.o : $(mySRCdir)/clapack_spotrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_spotrs.c
clapack_spotrf.o : $(mySRCdir)/clapack_spotrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_spotrf.c
clapack_sposv.o : $(mySRCdir)/clapack_sposv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sposv.c
clapack_sgelqf.o : $(mySRCdir)/clapack_sgelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgelqf.c
clapack_sgerqf.o : $(mySRCdir)/clapack_sgerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgerqf.c
clapack_sgeqlf.o : $(mySRCdir)/clapack_sgeqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgeqlf.c
clapack_sgeqrf.o : $(mySRCdir)/clapack_sgeqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_sgeqrf.c
clapack_slamch.o : $(mySRCdir)/clapack_slamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_slamch.c

clapack_dgetri.o : $(mySRCdir)/clapack_dgetri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgetri.c
clapack_dpotri.o : $(mySRCdir)/clapack_dpotri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dpotri.c
clapack_dlauum.o : $(mySRCdir)/clapack_dlauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dlauum.c
clapack_dtrtri.o : $(mySRCdir)/clapack_dtrtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dtrtri.c
clapack_dgetrs.o : $(mySRCdir)/clapack_dgetrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgetrs.c
clapack_dgetrf.o : $(mySRCdir)/clapack_dgetrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgetrf.c
clapack_dgesv.o : $(mySRCdir)/clapack_dgesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgesv.c
clapack_dpotrs.o : $(mySRCdir)/clapack_dpotrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dpotrs.c
clapack_dpotrf.o : $(mySRCdir)/clapack_dpotrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dpotrf.c
clapack_dposv.o : $(mySRCdir)/clapack_dposv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dposv.c
clapack_dgelqf.o : $(mySRCdir)/clapack_dgelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgelqf.c
clapack_dgerqf.o : $(mySRCdir)/clapack_dgerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgerqf.c
clapack_dgeqlf.o : $(mySRCdir)/clapack_dgeqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgeqlf.c
clapack_dgeqrf.o : $(mySRCdir)/clapack_dgeqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dgeqrf.c
clapack_dlamch.o : $(mySRCdir)/clapack_dlamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_dlamch.c

clapack_cgetri.o : $(mySRCdir)/clapack_cgetri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgetri.c
clapack_cpotri.o : $(mySRCdir)/clapack_cpotri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cpotri.c
clapack_clauum.o : $(mySRCdir)/clapack_clauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_clauum.c
clapack_ctrtri.o : $(mySRCdir)/clapack_ctrtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_ctrtri.c
clapack_cgetrs.o : $(mySRCdir)/clapack_cgetrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgetrs.c
clapack_cgetrf.o : $(mySRCdir)/clapack_cgetrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgetrf.c
clapack_cgesv.o : $(mySRCdir)/clapack_cgesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgesv.c
clapack_cpotrs.o : $(mySRCdir)/clapack_cpotrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cpotrs.c
clapack_cpotrf.o : $(mySRCdir)/clapack_cpotrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cpotrf.c
clapack_cposv.o : $(mySRCdir)/clapack_cposv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cposv.c
clapack_cgelqf.o : $(mySRCdir)/clapack_cgelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgelqf.c
clapack_cgerqf.o : $(mySRCdir)/clapack_cgerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgerqf.c
clapack_cgeqlf.o : $(mySRCdir)/clapack_cgeqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgeqlf.c
clapack_cgeqrf.o : $(mySRCdir)/clapack_cgeqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_cgeqrf.c
clapack_clamch.o : $(mySRCdir)/clapack_clamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_clamch.c

clapack_zgetri.o : $(mySRCdir)/clapack_zgetri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgetri.c
clapack_zpotri.o : $(mySRCdir)/clapack_zpotri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zpotri.c
clapack_zlauum.o : $(mySRCdir)/clapack_zlauum.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zlauum.c
clapack_ztrtri.o : $(mySRCdir)/clapack_ztrtri.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_ztrtri.c
clapack_zgetrs.o : $(mySRCdir)/clapack_zgetrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgetrs.c
clapack_zgetrf.o : $(mySRCdir)/clapack_zgetrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgetrf.c
clapack_zgesv.o : $(mySRCdir)/clapack_zgesv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgesv.c
clapack_zpotrs.o : $(mySRCdir)/clapack_zpotrs.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zpotrs.c
clapack_zpotrf.o : $(mySRCdir)/clapack_zpotrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zpotrf.c
clapack_zposv.o : $(mySRCdir)/clapack_zposv.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zposv.c
clapack_zgelqf.o : $(mySRCdir)/clapack_zgelqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgelqf.c
clapack_zgerqf.o : $(mySRCdir)/clapack_zgerqf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgerqf.c
clapack_zgeqlf.o : $(mySRCdir)/clapack_zgeqlf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgeqlf.c
clapack_zgeqrf.o : $(mySRCdir)/clapack_zgeqrf.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zgeqrf.c
clapack_zlamch.o : $(mySRCdir)/clapack_zlamch.c
	$(ICC) $(ICCFLAGS) -o $@ -c $(mySRCdir)/clapack_zlamch.c

