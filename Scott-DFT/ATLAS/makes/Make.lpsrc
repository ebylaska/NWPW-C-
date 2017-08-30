include Make.inc
mySRCdir = $(SRCdir)/src/lapack

all : lib
lib : dlib slib zlib clib laenv
clean : sclean dclean cclean zclean iclean
iclean :
	rm -f $(iobj)

#
# We build ATL_ilaenv into ATLAS rather than LAPACK, because it is dependent
# on the BLAS, and more importantly, the thread count, and objects put in
# lapack should be independent of both the BLAS and the nthreads.  Therefore,
# the C/F77 interface for ilaenv can be found in the blas/ subdir rather than
# the lapack (this way the serial/threaded version is determined by what
# lib you link in, just as with the BLAS)
#
ptlib : tenvlib.grd
tenvobj = ATL_itlaenv.o
tlaenv: tenvlib.grd
tenvlib.grd : $(tenvobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(tenvobj)
	$(RANLIB) $(ATLASlib)
	touch tenvlib.grd
envobj = ATL_ilaenv.o
laenv: envlib.grd
envlib.grd : $(envobj)
	$(ARCHIVER) $(ARFLAGS) $(ATLASlib) $(envobj)
	$(RANLIB) $(ATLASlib)
	touch envlib.grd
	
sobj = ATL_sgelq2.o ATL_sgelqf.o ATL_sgelqr.o ATL_sgeql2.o ATL_sgeqlf.o \
       ATL_sgeqlr.o ATL_sgeqr2.o ATL_sgeqrf.o ATL_sgeqrr.o ATL_sgerq2.o \
       ATL_sgerqf.o ATL_sgerqr.o ATL_sgetrf.o ATL_sgetrfC.o ATL_sgetrfR.o \
       ATL_sgetri.o ATL_sgetriC.o ATL_sgetriR.o ATL_sgetrs.o ATL_slamch.o \
       ATL_slapy2.o ATL_slarf.o ATL_slarfb.o ATL_slarfg.o ATL_slarfp.o \
       ATL_slarft.o ATL_slaswp.o ATL_slauum.o ATL_slauumCL.o ATL_slauumCU.o \
       ATL_slauumRL.o ATL_slauumRU.o ATL_spotrf.o ATL_spotrfL.o ATL_spotrfU.o \
       ATL_spotrs.o ATL_strtri.o ATL_strtriCL.o ATL_strtriCU.o ATL_strtriRL.o \
       ATL_strtriRU.o

incdep = $(INCSdir)/atlas_lapack.h
slib : 
	$(PMAKE) slib.grd
slib.grd : $(sobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(sobj)
	$(RANLIB) $(LAPACKlib)
	touch slib.grd

   
ATL_slamch.o : $(mySRCdir)/ATL_lamch.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slamch.o -DSREAL \
              $(mySRCdir)/ATL_lamch.c
ATL_slapy2.o : $(mySRCdir)/ATL_lapy2.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slapy2.o -DSREAL \
              $(mySRCdir)/ATL_lapy2.c
ATL_slarfp.o : $(mySRCdir)/ATL_larfp.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slarfp.o -DSREAL \
              $(mySRCdir)/ATL_larfp.c
ATL_slarft.o : $(mySRCdir)/ATL_larft.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slarft.o -DSREAL \
              $(mySRCdir)/ATL_larft.c
ATL_slarfg.o : $(mySRCdir)/ATL_larfg.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slarfg.o -DSREAL \
              $(mySRCdir)/ATL_larfg.c
ATL_slarf.o : $(mySRCdir)/ATL_larf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slarf.o -DSREAL \
              $(mySRCdir)/ATL_larf.c
ATL_slarfb.o : $(mySRCdir)/ATL_larfb.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slarfb.o -DSREAL \
              $(mySRCdir)/ATL_larfb.c
ATL_sgelqr.o : $(mySRCdir)/ATL_gelqr.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgelqr.o -DSREAL \
              $(mySRCdir)/ATL_gelqr.c
ATL_sgelqf.o : $(mySRCdir)/ATL_gelqf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgelqf.o -DSREAL \
              $(mySRCdir)/ATL_gelqf.c
ATL_sgelq2.o : $(mySRCdir)/ATL_gelq2.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgelq2.o -DSREAL \
              $(mySRCdir)/ATL_gelq2.c
ATL_sgeqlr.o : $(mySRCdir)/ATL_geqlr.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeqlr.o -DSREAL \
              $(mySRCdir)/ATL_geqlr.c
ATL_sgeqlf.o : $(mySRCdir)/ATL_geqlf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeqlf.o -DSREAL \
              $(mySRCdir)/ATL_geqlf.c
ATL_sgeql2.o : $(mySRCdir)/ATL_geql2.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeql2.o -DSREAL \
              $(mySRCdir)/ATL_geql2.c
ATL_sgerqr.o : $(mySRCdir)/ATL_gerqr.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgerqr.o -DSREAL \
              $(mySRCdir)/ATL_gerqr.c
ATL_sgerqf.o : $(mySRCdir)/ATL_gerqf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgerqf.o -DSREAL \
              $(mySRCdir)/ATL_gerqf.c
ATL_sgerq2.o : $(mySRCdir)/ATL_gerq2.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgerq2.o -DSREAL \
              $(mySRCdir)/ATL_gerq2.c
ATL_sgeqrr.o : $(mySRCdir)/ATL_geqrr.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeqrr.o -DSREAL \
              $(mySRCdir)/ATL_geqrr.c
ATL_sgeqrf.o : $(mySRCdir)/ATL_geqrf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeqrf.o -DSREAL \
              $(mySRCdir)/ATL_geqrf.c
ATL_sgeqr2.o : $(mySRCdir)/ATL_geqr2.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgeqr2.o -DSREAL \
              $(mySRCdir)/ATL_geqr2.c
ATL_spotrfL.o : $(mySRCdir)/ATL_potrfL.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_spotrfL.o -DSREAL \
              $(mySRCdir)/ATL_potrfL.c
ATL_spotrfU.o : $(mySRCdir)/ATL_potrfU.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_spotrfU.o -DSREAL \
              $(mySRCdir)/ATL_potrfU.c
ATL_spotrs.o : $(mySRCdir)/ATL_potrs.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_spotrs.o -DSREAL \
              $(mySRCdir)/ATL_potrs.c
ATL_spotrf.o : $(mySRCdir)/ATL_potrf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_spotrf.o -DSREAL \
              $(mySRCdir)/ATL_potrf.c
ATL_sgetrfR.o : $(mySRCdir)/ATL_getrfR.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetrfR.o -DSREAL \
              $(mySRCdir)/ATL_getrfR.c
ATL_sgetrfC.o : $(mySRCdir)/ATL_getrfC.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetrfC.o -DSREAL \
              $(mySRCdir)/ATL_getrfC.c
ATL_sgetrs.o : $(mySRCdir)/ATL_getrs.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetrs.o -DSREAL \
              $(mySRCdir)/ATL_getrs.c
ATL_sgetrf.o : $(mySRCdir)/ATL_getrf.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetrf.o -DSREAL \
              $(mySRCdir)/ATL_getrf.c
ATL_slaswp.o : $(mySRCdir)/ATL_laswp.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slaswp.o -DSREAL \
              $(mySRCdir)/ATL_laswp.c
ATL_sgetriC.o : $(mySRCdir)/ATL_getriC.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetriC.o -DSREAL \
              $(mySRCdir)/ATL_getriC.c
ATL_sgetriR.o : $(mySRCdir)/ATL_getriR.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetriR.o -DSREAL \
              $(mySRCdir)/ATL_getriR.c
ATL_sgetri.o : $(mySRCdir)/ATL_getri.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_sgetri.o -DSREAL \
              $(mySRCdir)/ATL_getri.c
ATL_slauum.o : $(mySRCdir)/ATL_lauum.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slauum.o -DSREAL \
              $(mySRCdir)/ATL_lauum.c
ATL_strtriRL.o : $(mySRCdir)/ATL_trtriRL.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_strtriRL.o -DSREAL \
              $(mySRCdir)/ATL_trtriRL.c
ATL_strtriRU.o : $(mySRCdir)/ATL_trtriRU.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_strtriRU.o -DSREAL \
              $(mySRCdir)/ATL_trtriRU.c
ATL_strtriCL.o : $(mySRCdir)/ATL_trtriCL.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_strtriCL.o -DSREAL \
              $(mySRCdir)/ATL_trtriCL.c
ATL_strtriCU.o : $(mySRCdir)/ATL_trtriCU.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_strtriCU.o -DSREAL \
              $(mySRCdir)/ATL_trtriCU.c
ATL_strtri.o : $(mySRCdir)/ATL_trtri.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_strtri.o -DSREAL \
              $(mySRCdir)/ATL_trtri.c
ATL_slauumCL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slauumCL.o -DSREAL -DColMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_slauumRL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slauumRL.o -DSREAL -DRowMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_slauumCU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slauumCU.o -DSREAL -DColMajor_ \
              $(mySRCdir)/ATL_lauumU.c
ATL_slauumRU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(SKC) $(sKCFLAGS) -c -o ATL_slauumRU.o -DSREAL -DRowMajor_ \
              $(mySRCdir)/ATL_lauumU.c

sclean :
	rm -f $(sobj)
dobj = ATL_dgelq2.o ATL_dgelqf.o ATL_dgelqr.o ATL_dgeql2.o ATL_dgeqlf.o \
       ATL_dgeqlr.o ATL_dgeqr2.o ATL_dgeqrf.o ATL_dgeqrr.o ATL_dgerq2.o \
       ATL_dgerqf.o ATL_dgerqr.o ATL_dgetrf.o ATL_dgetrfC.o ATL_dgetrfR.o \
       ATL_dgetri.o ATL_dgetriC.o ATL_dgetriR.o ATL_dgetrs.o ATL_dlamch.o \
       ATL_dlapy2.o ATL_dlarf.o ATL_dlarfb.o ATL_dlarfg.o ATL_dlarfp.o \
       ATL_dlarft.o ATL_dlaswp.o ATL_dlauum.o ATL_dlauumCL.o ATL_dlauumCU.o \
       ATL_dlauumRL.o ATL_dlauumRU.o ATL_dpotrf.o ATL_dpotrfL.o ATL_dpotrfU.o \
       ATL_dpotrs.o ATL_dtrtri.o ATL_dtrtriCL.o ATL_dtrtriCU.o ATL_dtrtriRL.o \
       ATL_dtrtriRU.o

incdep = $(INCSdir)/atlas_lapack.h
dlib : 
	$(PMAKE) dlib.grd
dlib.grd : $(dobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(dobj)
	$(RANLIB) $(LAPACKlib)
	touch dlib.grd

   
ATL_dlamch.o : $(mySRCdir)/ATL_lamch.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlamch.o -DDREAL \
              $(mySRCdir)/ATL_lamch.c
ATL_dlapy2.o : $(mySRCdir)/ATL_lapy2.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlapy2.o -DDREAL \
              $(mySRCdir)/ATL_lapy2.c
ATL_dlarfp.o : $(mySRCdir)/ATL_larfp.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlarfp.o -DDREAL \
              $(mySRCdir)/ATL_larfp.c
ATL_dlarft.o : $(mySRCdir)/ATL_larft.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlarft.o -DDREAL \
              $(mySRCdir)/ATL_larft.c
ATL_dlarfg.o : $(mySRCdir)/ATL_larfg.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlarfg.o -DDREAL \
              $(mySRCdir)/ATL_larfg.c
ATL_dlarf.o : $(mySRCdir)/ATL_larf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlarf.o -DDREAL \
              $(mySRCdir)/ATL_larf.c
ATL_dlarfb.o : $(mySRCdir)/ATL_larfb.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlarfb.o -DDREAL \
              $(mySRCdir)/ATL_larfb.c
ATL_dgelqr.o : $(mySRCdir)/ATL_gelqr.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgelqr.o -DDREAL \
              $(mySRCdir)/ATL_gelqr.c
ATL_dgelqf.o : $(mySRCdir)/ATL_gelqf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgelqf.o -DDREAL \
              $(mySRCdir)/ATL_gelqf.c
ATL_dgelq2.o : $(mySRCdir)/ATL_gelq2.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgelq2.o -DDREAL \
              $(mySRCdir)/ATL_gelq2.c
ATL_dgeqlr.o : $(mySRCdir)/ATL_geqlr.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeqlr.o -DDREAL \
              $(mySRCdir)/ATL_geqlr.c
ATL_dgeqlf.o : $(mySRCdir)/ATL_geqlf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeqlf.o -DDREAL \
              $(mySRCdir)/ATL_geqlf.c
ATL_dgeql2.o : $(mySRCdir)/ATL_geql2.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeql2.o -DDREAL \
              $(mySRCdir)/ATL_geql2.c
ATL_dgerqr.o : $(mySRCdir)/ATL_gerqr.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgerqr.o -DDREAL \
              $(mySRCdir)/ATL_gerqr.c
ATL_dgerqf.o : $(mySRCdir)/ATL_gerqf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgerqf.o -DDREAL \
              $(mySRCdir)/ATL_gerqf.c
ATL_dgerq2.o : $(mySRCdir)/ATL_gerq2.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgerq2.o -DDREAL \
              $(mySRCdir)/ATL_gerq2.c
ATL_dgeqrr.o : $(mySRCdir)/ATL_geqrr.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeqrr.o -DDREAL \
              $(mySRCdir)/ATL_geqrr.c
ATL_dgeqrf.o : $(mySRCdir)/ATL_geqrf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeqrf.o -DDREAL \
              $(mySRCdir)/ATL_geqrf.c
ATL_dgeqr2.o : $(mySRCdir)/ATL_geqr2.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgeqr2.o -DDREAL \
              $(mySRCdir)/ATL_geqr2.c
ATL_dpotrfL.o : $(mySRCdir)/ATL_potrfL.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dpotrfL.o -DDREAL \
              $(mySRCdir)/ATL_potrfL.c
ATL_dpotrfU.o : $(mySRCdir)/ATL_potrfU.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dpotrfU.o -DDREAL \
              $(mySRCdir)/ATL_potrfU.c
ATL_dpotrs.o : $(mySRCdir)/ATL_potrs.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dpotrs.o -DDREAL \
              $(mySRCdir)/ATL_potrs.c
ATL_dpotrf.o : $(mySRCdir)/ATL_potrf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dpotrf.o -DDREAL \
              $(mySRCdir)/ATL_potrf.c
ATL_dgetrfR.o : $(mySRCdir)/ATL_getrfR.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetrfR.o -DDREAL \
              $(mySRCdir)/ATL_getrfR.c
ATL_dgetrfC.o : $(mySRCdir)/ATL_getrfC.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetrfC.o -DDREAL \
              $(mySRCdir)/ATL_getrfC.c
ATL_dgetrs.o : $(mySRCdir)/ATL_getrs.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetrs.o -DDREAL \
              $(mySRCdir)/ATL_getrs.c
ATL_dgetrf.o : $(mySRCdir)/ATL_getrf.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetrf.o -DDREAL \
              $(mySRCdir)/ATL_getrf.c
ATL_dlaswp.o : $(mySRCdir)/ATL_laswp.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlaswp.o -DDREAL \
              $(mySRCdir)/ATL_laswp.c
ATL_dgetriC.o : $(mySRCdir)/ATL_getriC.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetriC.o -DDREAL \
              $(mySRCdir)/ATL_getriC.c
ATL_dgetriR.o : $(mySRCdir)/ATL_getriR.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetriR.o -DDREAL \
              $(mySRCdir)/ATL_getriR.c
ATL_dgetri.o : $(mySRCdir)/ATL_getri.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dgetri.o -DDREAL \
              $(mySRCdir)/ATL_getri.c
ATL_dlauum.o : $(mySRCdir)/ATL_lauum.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlauum.o -DDREAL \
              $(mySRCdir)/ATL_lauum.c
ATL_dtrtriRL.o : $(mySRCdir)/ATL_trtriRL.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dtrtriRL.o -DDREAL \
              $(mySRCdir)/ATL_trtriRL.c
ATL_dtrtriRU.o : $(mySRCdir)/ATL_trtriRU.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dtrtriRU.o -DDREAL \
              $(mySRCdir)/ATL_trtriRU.c
ATL_dtrtriCL.o : $(mySRCdir)/ATL_trtriCL.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dtrtriCL.o -DDREAL \
              $(mySRCdir)/ATL_trtriCL.c
ATL_dtrtriCU.o : $(mySRCdir)/ATL_trtriCU.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dtrtriCU.o -DDREAL \
              $(mySRCdir)/ATL_trtriCU.c
ATL_dtrtri.o : $(mySRCdir)/ATL_trtri.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dtrtri.o -DDREAL \
              $(mySRCdir)/ATL_trtri.c
ATL_dlauumCL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlauumCL.o -DDREAL -DColMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_dlauumRL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlauumRL.o -DDREAL -DRowMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_dlauumCU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlauumCU.o -DDREAL -DColMajor_ \
              $(mySRCdir)/ATL_lauumU.c
ATL_dlauumRU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(DKC) $(dKCFLAGS) -c -o ATL_dlauumRU.o -DDREAL -DRowMajor_ \
              $(mySRCdir)/ATL_lauumU.c

dclean :
	rm -f $(dobj)
cobj = ATL_cgelq2.o ATL_cgelqf.o ATL_cgelqr.o ATL_cgeql2.o ATL_cgeqlf.o \
       ATL_cgeqlr.o ATL_cgeqr2.o ATL_cgeqrf.o ATL_cgeqrr.o ATL_cgerq2.o \
       ATL_cgerqf.o ATL_cgerqr.o ATL_cgetrf.o ATL_cgetrfC.o ATL_cgetrfR.o \
       ATL_cgetri.o ATL_cgetriC.o ATL_cgetriR.o ATL_cgetrs.o ATL_clacgv.o \
       ATL_cladiv.o ATL_clapy2.o ATL_clapy3.o ATL_clarf.o ATL_clarfb.o \
       ATL_clarfg.o ATL_clarfp.o ATL_clarft.o ATL_claswp.o ATL_clauum.o \
       ATL_clauumCL.o ATL_clauumCU.o ATL_clauumRL.o ATL_clauumRU.o \
       ATL_cpotrf.o ATL_cpotrfL.o ATL_cpotrfRL.o ATL_cpotrfRU.o ATL_cpotrfU.o \
       ATL_cpotrs.o ATL_ctrtri.o ATL_ctrtriCL.o ATL_ctrtriCU.o ATL_ctrtriRL.o \
       ATL_ctrtriRU.o

incdep = $(INCSdir)/atlas_lapack.h
clib : 
	$(PMAKE) clib.grd
clib.grd : $(cobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(cobj)
	$(RANLIB) $(LAPACKlib)
	touch clib.grd

   
ATL_clacgv.o : $(mySRCdir)/ATL_lacgv.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clacgv.o -DSCPLX \
              $(mySRCdir)/ATL_lacgv.c
ATL_clapy3.o : $(mySRCdir)/ATL_lapy3.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clapy3.o -DSCPLX \
              $(mySRCdir)/ATL_lapy3.c
ATL_cladiv.o : $(mySRCdir)/ATL_ladiv.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cladiv.o -DSCPLX \
              $(mySRCdir)/ATL_ladiv.c
ATL_cpotrfRU.o : $(mySRCdir)/ATL_potrfRU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrfRU.o -DSCPLX \
              $(mySRCdir)/ATL_potrfRU.c
ATL_cpotrfRL.o : $(mySRCdir)/ATL_potrfRL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrfRL.o -DSCPLX \
              $(mySRCdir)/ATL_potrfRL.c
ATL_clapy2.o : $(mySRCdir)/ATL_lapy2.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clapy2.o -DSCPLX \
              $(mySRCdir)/ATL_lapy2.c
ATL_clarfp.o : $(mySRCdir)/ATL_larfp.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clarfp.o -DSCPLX \
              $(mySRCdir)/ATL_larfp.c
ATL_clarft.o : $(mySRCdir)/ATL_larft.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clarft.o -DSCPLX \
              $(mySRCdir)/ATL_larft.c
ATL_clarfg.o : $(mySRCdir)/ATL_larfg.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clarfg.o -DSCPLX \
              $(mySRCdir)/ATL_larfg.c
ATL_clarf.o : $(mySRCdir)/ATL_larf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clarf.o -DSCPLX \
              $(mySRCdir)/ATL_larf.c
ATL_clarfb.o : $(mySRCdir)/ATL_larfb.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clarfb.o -DSCPLX \
              $(mySRCdir)/ATL_larfb.c
ATL_cgelqr.o : $(mySRCdir)/ATL_gelqr.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgelqr.o -DSCPLX \
              $(mySRCdir)/ATL_gelqr.c
ATL_cgelqf.o : $(mySRCdir)/ATL_gelqf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgelqf.o -DSCPLX \
              $(mySRCdir)/ATL_gelqf.c
ATL_cgelq2.o : $(mySRCdir)/ATL_gelq2.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgelq2.o -DSCPLX \
              $(mySRCdir)/ATL_gelq2.c
ATL_cgeqlr.o : $(mySRCdir)/ATL_geqlr.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeqlr.o -DSCPLX \
              $(mySRCdir)/ATL_geqlr.c
ATL_cgeqlf.o : $(mySRCdir)/ATL_geqlf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeqlf.o -DSCPLX \
              $(mySRCdir)/ATL_geqlf.c
ATL_cgeql2.o : $(mySRCdir)/ATL_geql2.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeql2.o -DSCPLX \
              $(mySRCdir)/ATL_geql2.c
ATL_cgerqr.o : $(mySRCdir)/ATL_gerqr.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgerqr.o -DSCPLX \
              $(mySRCdir)/ATL_gerqr.c
ATL_cgerqf.o : $(mySRCdir)/ATL_gerqf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgerqf.o -DSCPLX \
              $(mySRCdir)/ATL_gerqf.c
ATL_cgerq2.o : $(mySRCdir)/ATL_gerq2.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgerq2.o -DSCPLX \
              $(mySRCdir)/ATL_gerq2.c
ATL_cgeqrr.o : $(mySRCdir)/ATL_geqrr.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeqrr.o -DSCPLX \
              $(mySRCdir)/ATL_geqrr.c
ATL_cgeqrf.o : $(mySRCdir)/ATL_geqrf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeqrf.o -DSCPLX \
              $(mySRCdir)/ATL_geqrf.c
ATL_cgeqr2.o : $(mySRCdir)/ATL_geqr2.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgeqr2.o -DSCPLX \
              $(mySRCdir)/ATL_geqr2.c
ATL_cpotrfL.o : $(mySRCdir)/ATL_potrfL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrfL.o -DSCPLX \
              $(mySRCdir)/ATL_potrfL.c
ATL_cpotrfU.o : $(mySRCdir)/ATL_potrfU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrfU.o -DSCPLX \
              $(mySRCdir)/ATL_potrfU.c
ATL_cpotrs.o : $(mySRCdir)/ATL_potrs.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrs.o -DSCPLX \
              $(mySRCdir)/ATL_potrs.c
ATL_cpotrf.o : $(mySRCdir)/ATL_potrf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cpotrf.o -DSCPLX \
              $(mySRCdir)/ATL_potrf.c
ATL_cgetrfR.o : $(mySRCdir)/ATL_getrfR.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetrfR.o -DSCPLX \
              $(mySRCdir)/ATL_getrfR.c
ATL_cgetrfC.o : $(mySRCdir)/ATL_getrfC.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetrfC.o -DSCPLX \
              $(mySRCdir)/ATL_getrfC.c
ATL_cgetrs.o : $(mySRCdir)/ATL_getrs.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetrs.o -DSCPLX \
              $(mySRCdir)/ATL_getrs.c
ATL_cgetrf.o : $(mySRCdir)/ATL_getrf.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetrf.o -DSCPLX \
              $(mySRCdir)/ATL_getrf.c
ATL_claswp.o : $(mySRCdir)/ATL_laswp.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_claswp.o -DSCPLX \
              $(mySRCdir)/ATL_laswp.c
ATL_cgetriC.o : $(mySRCdir)/ATL_getriC.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetriC.o -DSCPLX \
              $(mySRCdir)/ATL_getriC.c
ATL_cgetriR.o : $(mySRCdir)/ATL_getriR.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetriR.o -DSCPLX \
              $(mySRCdir)/ATL_getriR.c
ATL_cgetri.o : $(mySRCdir)/ATL_getri.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_cgetri.o -DSCPLX \
              $(mySRCdir)/ATL_getri.c
ATL_clauum.o : $(mySRCdir)/ATL_lauum.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clauum.o -DSCPLX \
              $(mySRCdir)/ATL_lauum.c
ATL_ctrtriRL.o : $(mySRCdir)/ATL_trtriRL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_ctrtriRL.o -DSCPLX \
              $(mySRCdir)/ATL_trtriRL.c
ATL_ctrtriRU.o : $(mySRCdir)/ATL_trtriRU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_ctrtriRU.o -DSCPLX \
              $(mySRCdir)/ATL_trtriRU.c
ATL_ctrtriCL.o : $(mySRCdir)/ATL_trtriCL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_ctrtriCL.o -DSCPLX \
              $(mySRCdir)/ATL_trtriCL.c
ATL_ctrtriCU.o : $(mySRCdir)/ATL_trtriCU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_ctrtriCU.o -DSCPLX \
              $(mySRCdir)/ATL_trtriCU.c
ATL_ctrtri.o : $(mySRCdir)/ATL_trtri.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_ctrtri.o -DSCPLX \
              $(mySRCdir)/ATL_trtri.c
ATL_clauumCL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clauumCL.o -DSCPLX -DColMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_clauumRL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clauumRL.o -DSCPLX -DRowMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_clauumCU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clauumCU.o -DSCPLX -DColMajor_ \
              $(mySRCdir)/ATL_lauumU.c
ATL_clauumRU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(CKC) $(cKCFLAGS) -c -o ATL_clauumRU.o -DSCPLX -DRowMajor_ \
              $(mySRCdir)/ATL_lauumU.c

cclean :
	rm -f $(cobj)
zobj = ATL_zgelq2.o ATL_zgelqf.o ATL_zgelqr.o ATL_zgeql2.o ATL_zgeqlf.o \
       ATL_zgeqlr.o ATL_zgeqr2.o ATL_zgeqrf.o ATL_zgeqrr.o ATL_zgerq2.o \
       ATL_zgerqf.o ATL_zgerqr.o ATL_zgetrf.o ATL_zgetrfC.o ATL_zgetrfR.o \
       ATL_zgetri.o ATL_zgetriC.o ATL_zgetriR.o ATL_zgetrs.o ATL_zlacgv.o \
       ATL_zladiv.o ATL_zlapy2.o ATL_zlapy3.o ATL_zlarf.o ATL_zlarfb.o \
       ATL_zlarfg.o ATL_zlarfp.o ATL_zlarft.o ATL_zlaswp.o ATL_zlauum.o \
       ATL_zlauumCL.o ATL_zlauumCU.o ATL_zlauumRL.o ATL_zlauumRU.o \
       ATL_zpotrf.o ATL_zpotrfL.o ATL_zpotrfRL.o ATL_zpotrfRU.o ATL_zpotrfU.o \
       ATL_zpotrs.o ATL_ztrtri.o ATL_ztrtriCL.o ATL_ztrtriCU.o ATL_ztrtriRL.o \
       ATL_ztrtriRU.o

incdep = $(INCSdir)/atlas_lapack.h
zlib : 
	$(PMAKE) zlib.grd
zlib.grd : $(zobj)
	$(ARCHIVER) $(ARFLAGS) $(LAPACKlib) $(zobj)
	$(RANLIB) $(LAPACKlib)
	touch zlib.grd

   
ATL_zlacgv.o : $(mySRCdir)/ATL_lacgv.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlacgv.o -DDCPLX \
              $(mySRCdir)/ATL_lacgv.c
ATL_zlapy3.o : $(mySRCdir)/ATL_lapy3.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlapy3.o -DDCPLX \
              $(mySRCdir)/ATL_lapy3.c
ATL_zladiv.o : $(mySRCdir)/ATL_ladiv.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zladiv.o -DDCPLX \
              $(mySRCdir)/ATL_ladiv.c
ATL_zpotrfRU.o : $(mySRCdir)/ATL_potrfRU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrfRU.o -DDCPLX \
              $(mySRCdir)/ATL_potrfRU.c
ATL_zpotrfRL.o : $(mySRCdir)/ATL_potrfRL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrfRL.o -DDCPLX \
              $(mySRCdir)/ATL_potrfRL.c
ATL_zlapy2.o : $(mySRCdir)/ATL_lapy2.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlapy2.o -DDCPLX \
              $(mySRCdir)/ATL_lapy2.c
ATL_zlarfp.o : $(mySRCdir)/ATL_larfp.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlarfp.o -DDCPLX \
              $(mySRCdir)/ATL_larfp.c
ATL_zlarft.o : $(mySRCdir)/ATL_larft.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlarft.o -DDCPLX \
              $(mySRCdir)/ATL_larft.c
ATL_zlarfg.o : $(mySRCdir)/ATL_larfg.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlarfg.o -DDCPLX \
              $(mySRCdir)/ATL_larfg.c
ATL_zlarf.o : $(mySRCdir)/ATL_larf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlarf.o -DDCPLX \
              $(mySRCdir)/ATL_larf.c
ATL_zlarfb.o : $(mySRCdir)/ATL_larfb.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlarfb.o -DDCPLX \
              $(mySRCdir)/ATL_larfb.c
ATL_zgelqr.o : $(mySRCdir)/ATL_gelqr.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgelqr.o -DDCPLX \
              $(mySRCdir)/ATL_gelqr.c
ATL_zgelqf.o : $(mySRCdir)/ATL_gelqf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgelqf.o -DDCPLX \
              $(mySRCdir)/ATL_gelqf.c
ATL_zgelq2.o : $(mySRCdir)/ATL_gelq2.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgelq2.o -DDCPLX \
              $(mySRCdir)/ATL_gelq2.c
ATL_zgeqlr.o : $(mySRCdir)/ATL_geqlr.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeqlr.o -DDCPLX \
              $(mySRCdir)/ATL_geqlr.c
ATL_zgeqlf.o : $(mySRCdir)/ATL_geqlf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeqlf.o -DDCPLX \
              $(mySRCdir)/ATL_geqlf.c
ATL_zgeql2.o : $(mySRCdir)/ATL_geql2.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeql2.o -DDCPLX \
              $(mySRCdir)/ATL_geql2.c
ATL_zgerqr.o : $(mySRCdir)/ATL_gerqr.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgerqr.o -DDCPLX \
              $(mySRCdir)/ATL_gerqr.c
ATL_zgerqf.o : $(mySRCdir)/ATL_gerqf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgerqf.o -DDCPLX \
              $(mySRCdir)/ATL_gerqf.c
ATL_zgerq2.o : $(mySRCdir)/ATL_gerq2.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgerq2.o -DDCPLX \
              $(mySRCdir)/ATL_gerq2.c
ATL_zgeqrr.o : $(mySRCdir)/ATL_geqrr.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeqrr.o -DDCPLX \
              $(mySRCdir)/ATL_geqrr.c
ATL_zgeqrf.o : $(mySRCdir)/ATL_geqrf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeqrf.o -DDCPLX \
              $(mySRCdir)/ATL_geqrf.c
ATL_zgeqr2.o : $(mySRCdir)/ATL_geqr2.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgeqr2.o -DDCPLX \
              $(mySRCdir)/ATL_geqr2.c
ATL_zpotrfL.o : $(mySRCdir)/ATL_potrfL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrfL.o -DDCPLX \
              $(mySRCdir)/ATL_potrfL.c
ATL_zpotrfU.o : $(mySRCdir)/ATL_potrfU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrfU.o -DDCPLX \
              $(mySRCdir)/ATL_potrfU.c
ATL_zpotrs.o : $(mySRCdir)/ATL_potrs.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrs.o -DDCPLX \
              $(mySRCdir)/ATL_potrs.c
ATL_zpotrf.o : $(mySRCdir)/ATL_potrf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zpotrf.o -DDCPLX \
              $(mySRCdir)/ATL_potrf.c
ATL_zgetrfR.o : $(mySRCdir)/ATL_getrfR.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetrfR.o -DDCPLX \
              $(mySRCdir)/ATL_getrfR.c
ATL_zgetrfC.o : $(mySRCdir)/ATL_getrfC.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetrfC.o -DDCPLX \
              $(mySRCdir)/ATL_getrfC.c
ATL_zgetrs.o : $(mySRCdir)/ATL_getrs.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetrs.o -DDCPLX \
              $(mySRCdir)/ATL_getrs.c
ATL_zgetrf.o : $(mySRCdir)/ATL_getrf.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetrf.o -DDCPLX \
              $(mySRCdir)/ATL_getrf.c
ATL_zlaswp.o : $(mySRCdir)/ATL_laswp.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlaswp.o -DDCPLX \
              $(mySRCdir)/ATL_laswp.c
ATL_zgetriC.o : $(mySRCdir)/ATL_getriC.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetriC.o -DDCPLX \
              $(mySRCdir)/ATL_getriC.c
ATL_zgetriR.o : $(mySRCdir)/ATL_getriR.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetriR.o -DDCPLX \
              $(mySRCdir)/ATL_getriR.c
ATL_zgetri.o : $(mySRCdir)/ATL_getri.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zgetri.o -DDCPLX \
              $(mySRCdir)/ATL_getri.c
ATL_zlauum.o : $(mySRCdir)/ATL_lauum.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlauum.o -DDCPLX \
              $(mySRCdir)/ATL_lauum.c
ATL_ztrtriRL.o : $(mySRCdir)/ATL_trtriRL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_ztrtriRL.o -DDCPLX \
              $(mySRCdir)/ATL_trtriRL.c
ATL_ztrtriRU.o : $(mySRCdir)/ATL_trtriRU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_ztrtriRU.o -DDCPLX \
              $(mySRCdir)/ATL_trtriRU.c
ATL_ztrtriCL.o : $(mySRCdir)/ATL_trtriCL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_ztrtriCL.o -DDCPLX \
              $(mySRCdir)/ATL_trtriCL.c
ATL_ztrtriCU.o : $(mySRCdir)/ATL_trtriCU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_ztrtriCU.o -DDCPLX \
              $(mySRCdir)/ATL_trtriCU.c
ATL_ztrtri.o : $(mySRCdir)/ATL_trtri.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_ztrtri.o -DDCPLX \
              $(mySRCdir)/ATL_trtri.c
ATL_zlauumCL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlauumCL.o -DDCPLX -DColMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_zlauumRL.o : $(mySRCdir)/ATL_lauumL.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlauumRL.o -DDCPLX -DRowMajor_ \
              $(mySRCdir)/ATL_lauumL.c
ATL_zlauumCU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlauumCU.o -DDCPLX -DColMajor_ \
              $(mySRCdir)/ATL_lauumU.c
ATL_zlauumRU.o : $(mySRCdir)/ATL_lauumU.c $(incdep)
	$(ZKC) $(zKCFLAGS) -c -o ATL_zlauumRU.o -DDCPLX -DRowMajor_ \
              $(mySRCdir)/ATL_lauumU.c

zclean :
	rm -f $(zobj)

tenvdep = $(INCSdir)/atlas_lapack.h $(INCAdir)/atlas_stGetNB_gerqf.h \
          $(INCAdir)/atlas_stGetNB_geqrf.h $(INCAdir)/atlas_stGetNB_geqlf.h \
          $(INCAdir)/atlas_stGetNB_gelqf.h $(INCAdir)/atlas_dtGetNB_gerqf.h \
          $(INCAdir)/atlas_dtGetNB_geqrf.h $(INCAdir)/atlas_dtGetNB_geqlf.h \
          $(INCAdir)/atlas_dtGetNB_gelqf.h $(INCAdir)/atlas_ctGetNB_gerqf.h \
          $(INCAdir)/atlas_ctGetNB_geqrf.h $(INCAdir)/atlas_ctGetNB_geqlf.h \
          $(INCAdir)/atlas_ctGetNB_gelqf.h $(INCAdir)/atlas_ztGetNB_gerqf.h \
          $(INCAdir)/atlas_ztGetNB_geqrf.h $(INCAdir)/atlas_ztGetNB_geqlf.h \
          $(INCAdir)/atlas_ztGetNB_gelqf.h
ATL_itlaenv.o : $(tenvdep) $(mySRCdir)/ATL_ilaenv.c
	$(DKC) $(DKCFLAGS) $(CDEFS) -c -o ATL_itlaenv.o -DATL_USEPTHREADS \
               $(mySRCdir)/ATL_ilaenv.c
envdep = $(INCSdir)/atlas_lapack.h $(INCAdir)/atlas_sGetNB_gerqf.h \
         $(INCAdir)/atlas_sGetNB_geqrf.h $(INCAdir)/atlas_sGetNB_geqlf.h \
         $(INCAdir)/atlas_sGetNB_gelqf.h $(INCAdir)/atlas_dGetNB_gerqf.h \
         $(INCAdir)/atlas_dGetNB_geqrf.h $(INCAdir)/atlas_dGetNB_geqlf.h \
         $(INCAdir)/atlas_dGetNB_gelqf.h $(INCAdir)/atlas_cGetNB_gerqf.h \
         $(INCAdir)/atlas_cGetNB_geqrf.h $(INCAdir)/atlas_cGetNB_geqlf.h \
         $(INCAdir)/atlas_cGetNB_gelqf.h $(INCAdir)/atlas_zGetNB_gerqf.h \
         $(INCAdir)/atlas_zGetNB_geqrf.h $(INCAdir)/atlas_zGetNB_geqlf.h \
         $(INCAdir)/atlas_zGetNB_gelqf.h
ATL_ilaenv.o : $(envdep) $(mySRCdir)/ATL_ilaenv.c
	$(DKC) $(DKCFLAGS) $(CDEFS) -c -o ATL_ilaenv.o  \
               $(mySRCdir)/ATL_ilaenv.c


