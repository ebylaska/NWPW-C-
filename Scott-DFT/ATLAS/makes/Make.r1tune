include Make.inc
mySRCdir = $(SRCdir)/tune/blas/ger
r1SRCdir = $(SRCdir)/src/blas/ger

sR1CC = $(SKC)
sR1CFLAGS = $(SKCFLAGS)
cR1CC = $(SKC)
cR1CFLAGS = $(SKCFLAGS)
dR1CC = $(DKC)
dR1CFLAGS = $(DKCFLAGS)
zR1CC = $(DKC)
zR1CFLAGS = $(DKCFLAGS)
cas=0
xu=16
yu=2
l1mul=75
iflag=0
Mt=297
Nt=177
ldat=$(Mt)
M=1000
N=1000
lda=$(M)
incy=3
flushKB = 4192
align = -Fx 16 -Fy 16 -Fa 16
tflags =
ssegenargs=

xr1search : xemit_r1h r1search.o l1lib tstlib
	$(XCC) $(XCCFLAGS) -o xr1search r1search.o
r1search.o : $(mySRCdir)/r1search.c scases.dsc dcases.dsc ccases.dsc zcases.dsc
	$(XCC) -c $(XCCFLAGS) -o r1search.o $(mySRCdir)/r1search.c

xr1hgen : r1hgen.o
	$(XCC) $(XCCFLAGS) -o $@ r1hgen.o
r1hgen.o : $(INCAdir)/atlas_type.h $(mySRCdir)/r1hgen.c $(INCSdir)/atlas_r1parse.h
	$(XCC) -c $(XCCFLAGS) -o $@ $(mySRCdir)/r1hgen.c
xr1ksearch : r1ksearch.o
	$(XCC) $(XCCFLAGS) -o $@ r1ksearch.o
r1ksearch.o : $(INCAdir)/atlas_type.h $(mySRCdir)/r1ksearch.c $(INCSdir)/atlas_r1parse.h $(INCSdir)/atlas_r1testtime.h
	$(XCC) -c $(XCCFLAGS) -o $@ $(mySRCdir)/r1ksearch.c
xr1gen_sse : r1gen_sse.o
	$(XCC) $(XCCFLAGS) -o $@ r1gen_sse.o
r1gen_sse.o : $(INCAdir)/atlas_type.h $(mySRCdir)/r1gen_sse.c  
	$(XCC) -c $(XCCFLAGS) -o $@ $(mySRCdir)/r1gen_sse.c
xemit_r1h : emit_r1h.o
	$(XCC) $(XCCFLAGS) -o $@ emit_r1h.o
emit_r1h.o : $(INCAdir)/atlas_type.h $(mySRCdir)/emit_r1h.c  
	$(XCC) -c $(XCCFLAGS) -o $@ $(mySRCdir)/emit_r1h.c
$(INCAdir)/atlas_type.h :
	cd $(SYSdir) ; $(MAKE) RunTyp
l1lib : sl1lib dl1lib cl1lib zl1lib 
tstlib : ststlib dtstlib ctstlib ztstlib

force_build :

res/L1CacheSize : $(SYSdir)/res/L1CacheSize
	cd $(SYSdir) ; $(MAKE) res/L1CacheSize
	ln -s $(SYSdir)/res/L1CacheSize res/L1CacheSize

ATL_ger.c :
	rm -f ATL_ger.c
	echo "#define ATL_TUNING" > ATL_ger.c
	cat $(r1SRCdir)/ATL_ger.c >> ATL_ger.c

ststobj = ATL_sger1.o ATL_sgerk_Mlt16.o ATL_sgerk_axpy.o
sobj = ATL_sger.o ATL_sger1.o \
       ATL_sgerk_Mlt16.o ATL_sgerk_axpy.o

#
# This target is for directly timing a call to the kernel, w/o cache flushing
#
sr1ktime : sr1ktime.o ATL_sgerk_Mlt16.o ATL_sgerk_axpy.o
	rm -f ATL_sger1.c
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_sger1.c
	$(sR1CC) -c $(CDEFS) $(sR1CFLAGS) -I$(R1Tdir) -o ATL_sger1.o \
           -DSREAL -DATL_sger1_a1_x1_yX=ATL_UGERK ATL_sger1.c
	$(CLINKER) $(CLINKFLAGS) -o xsr1ktime sr1ktime.o \
                   ATL_sger1.o \
                   ATL_sgerk_Mlt16.o ATL_sgerk_axpy.o \
                   $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xsr1ktime -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
#
# This target calls the GER driver (which blocks), with cache flushing
#
sr1time : sr1time.o xr1hgen
	rm -f $(sobj) ATL_sger1.c atlas_sr1.h r1h.bak
	if [ -e "$(INCAdir)/atlas_sr1.h" ]; then \
           mv $(INCAdir)/atlas_sr1.h r1h.bak ; \
        fi
	./xr1hgen -p s -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_sger1.c
	$(MAKE) $(sobj)
	$(CLINKER) $(CLINKFLAGS) -o xsr1time sr1time.o \
                   $(sobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xsr1time -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_sr1.h ; \
        fi

sr1ktest : sr1ktest.o xr1hgen
	rm -f $(ststobj) ATL_sger1.c 
	if [ -e "$(INCAdir)/atlas_sr1.h" ]; then \
           mv $(INCAdir)/atlas_sr1.h r1h.bak ; \
        fi
	./xr1hgen -p s -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_sger1.c
	$(MAKE) $(ststobj)
	$(CLINKER) $(CLINKFLAGS) -o xsr1ktest sr1ktest.o \
                   $(ststobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xsr1ktest -m 1 $(Mt) -n 1 $(Nt) -y 1 $(incy) \
                  -l 1 $(ldat) $(align)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_sr1.h ; \
        fi


sr1time.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL $(mySRCdir)/r1ktime.c
sr1ktime.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL -DTIME_KERNEL $(mySRCdir)/r1ktime.c
sr1ktest.o : $(mySRCdir)/r1ktest.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL $(mySRCdir)/r1ktest.c
sr1test.o : $(mySRCdir)/r1test.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL $(mySRCdir)/r1test.c
sger1tune.o : $(mySRCdir)/ger1tune.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL $(mySRCdir)/ger1tune.c
ATL_sger.o : ATL_ger.c atlas_sr1.h
	$(ICC) -c $(ICCFLAGS) -o $@ -DSREAL -I$(R1Tdir) ATL_ger.c
ATL_sgerk_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DSREAL $(r1SRCdir)/ATL_gerk_axpy.c
ATL_sgerk_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DSREAL $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_sger1.o : ATL_sger1.c
	$(sR1CC) -c $(CDEFS) $(sR1CFLAGS) -I$(R1Tdir) -o $@ -DSREAL ATL_sger1.c
ATL_sger1c.o : ATL_sger1.c
	$(sR1CC) -c $(CDEFS) $(sR1CFLAGS) -I$(R1Tdir) -o $@ -DSREAL -DConj_ \
              -DATL_UGERK=ATL_GERK_Conj ATL_sger1.c

res/sR1SUMM : xr1ksearch force_build
	./xr1ksearch -p s
	$(MAKE) sinstall
sr1outd:
	mkdir sr1outd
sinstall : xccobj sr1outd xr1hgen force_build
	rm -f sr1outd/*
	if [ -s "$(GR1dir)/Make_sr1" ]; then \
           cd  $(GR1dir) ; $(MAKE) -f Make_sr1 killall ; \
           rm -f $(GR1dir)/Make_sr1 ; \
        fi
	./xr1hgen -p s -F res/sR1SUMM -d sr1outd
	cp -f sr1outd/Make_sr1 $(GR1dir)/.
	cp sr1outd/*.c $(GR1dir)/.
	cp -f sr1outd/*.h $(INCAdir)/.
	if [ -s "xccobj.exe" ] ; then \
	   cp -f xccobj.exe $(GR1dir)/. ; \
	elif [ -s "xccobj" ]; then \
           cp -f xccobj $(GR1dir)/. ; \
        fi
	touch $(INCSdir)/atlas_r1.h
	rm -f sr1outd/* ; rmdir sr1outd

CASES/ATL_sssegen.c : xr1gen_sse force_build
	rm -f CASES/ATL_sssegen.c
	./xr1gen_sse $(ssegenargs) > CASES/ATL_sssegen.c

sl1lib :
	cd $(BINdir) ; $(MAKE) sl1lib
ststlib :
	cd $(BINdir) ; $(MAKE) ststlib
dtstobj = ATL_dger1.o ATL_dgerk_Mlt16.o ATL_dgerk_axpy.o
dobj = ATL_dger.o ATL_dger1.o \
       ATL_dgerk_Mlt16.o ATL_dgerk_axpy.o

#
# This target is for directly timing a call to the kernel, w/o cache flushing
#
dr1ktime : dr1ktime.o ATL_dgerk_Mlt16.o ATL_dgerk_axpy.o
	rm -f ATL_dger1.c
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_dger1.c
	$(dR1CC) -c $(CDEFS) $(dR1CFLAGS) -I$(R1Tdir) -o ATL_dger1.o \
           -DDREAL -DATL_dger1_a1_x1_yX=ATL_UGERK ATL_dger1.c
	$(CLINKER) $(CLINKFLAGS) -o xdr1ktime dr1ktime.o \
                   ATL_dger1.o \
                   ATL_dgerk_Mlt16.o ATL_dgerk_axpy.o \
                   $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xdr1ktime -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
#
# This target calls the GER driver (which blocks), with cache flushing
#
dr1time : dr1time.o xr1hgen
	rm -f $(dobj) ATL_dger1.c atlas_dr1.h r1h.bak
	if [ -e "$(INCAdir)/atlas_dr1.h" ]; then \
           mv $(INCAdir)/atlas_dr1.h r1h.bak ; \
        fi
	./xr1hgen -p d -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_dger1.c
	$(MAKE) $(dobj)
	$(CLINKER) $(CLINKFLAGS) -o xdr1time dr1time.o \
                   $(dobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xdr1time -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_dr1.h ; \
        fi

dr1ktest : dr1ktest.o xr1hgen
	rm -f $(dtstobj) ATL_dger1.c 
	if [ -e "$(INCAdir)/atlas_dr1.h" ]; then \
           mv $(INCAdir)/atlas_dr1.h r1h.bak ; \
        fi
	./xr1hgen -p d -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_dger1.c
	$(MAKE) $(dtstobj)
	$(CLINKER) $(CLINKFLAGS) -o xdr1ktest dr1ktest.o \
                   $(dtstobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xdr1ktest -m 1 $(Mt) -n 1 $(Nt) -y 1 $(incy) \
                  -l 1 $(ldat) $(align)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_dr1.h ; \
        fi


dr1time.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL $(mySRCdir)/r1ktime.c
dr1ktime.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL -DTIME_KERNEL $(mySRCdir)/r1ktime.c
dr1ktest.o : $(mySRCdir)/r1ktest.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL $(mySRCdir)/r1ktest.c
dr1test.o : $(mySRCdir)/r1test.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL $(mySRCdir)/r1test.c
dger1tune.o : $(mySRCdir)/ger1tune.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL $(mySRCdir)/ger1tune.c
ATL_dger.o : ATL_ger.c atlas_dr1.h
	$(ICC) -c $(ICCFLAGS) -o $@ -DDREAL -I$(R1Tdir) ATL_ger.c
ATL_dgerk_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DDREAL $(r1SRCdir)/ATL_gerk_axpy.c
ATL_dgerk_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DDREAL $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_dger1.o : ATL_dger1.c
	$(dR1CC) -c $(CDEFS) $(dR1CFLAGS) -I$(R1Tdir) -o $@ -DDREAL ATL_dger1.c
ATL_dger1c.o : ATL_dger1.c
	$(dR1CC) -c $(CDEFS) $(dR1CFLAGS) -I$(R1Tdir) -o $@ -DDREAL -DConj_ \
              -DATL_UGERK=ATL_GERK_Conj ATL_dger1.c

res/dR1SUMM : xr1ksearch force_build
	./xr1ksearch -p d
	$(MAKE) dinstall
dr1outd:
	mkdir dr1outd
dinstall : xccobj dr1outd xr1hgen force_build
	rm -f dr1outd/*
	if [ -s "$(GR1dir)/Make_dr1" ]; then \
           cd  $(GR1dir) ; $(MAKE) -f Make_dr1 killall ; \
           rm -f $(GR1dir)/Make_dr1 ; \
        fi
	./xr1hgen -p d -F res/dR1SUMM -d dr1outd
	cp -f dr1outd/Make_dr1 $(GR1dir)/.
	cp dr1outd/*.c $(GR1dir)/.
	cp -f dr1outd/*.h $(INCAdir)/.
	if [ -s "xccobj.exe" ] ; then \
	   cp -f xccobj.exe $(GR1dir)/. ; \
	elif [ -s "xccobj" ]; then \
           cp -f xccobj $(GR1dir)/. ; \
        fi
	touch $(INCSdir)/atlas_r1.h
	rm -f dr1outd/* ; rmdir dr1outd

CASES/ATL_dssegen.c : xr1gen_sse force_build
	rm -f CASES/ATL_dssegen.c
	./xr1gen_sse $(ssegenargs) > CASES/ATL_dssegen.c

dl1lib :
	cd $(BINdir) ; $(MAKE) dl1lib
dtstlib :
	cd $(BINdir) ; $(MAKE) dtstlib
ctstobj = ATL_cger1.o  ATL_cger1c.o  \
          ATL_cgerk_Mlt16.o ATL_cgerk_axpy.o \
          ATL_cgerck_Mlt16.o ATL_cgerck_axpy.o 
cobj = ATL_cger.o ATL_cger1.o \
       ATL_cgerk_Mlt16.o ATL_cgerk_axpy.o \
       ATL_cgerck_Mlt16.o ATL_cgerck_axpy.o 

#
# This target is for directly timing a call to the kernel, w/o cache flushing
#
cr1ktime : cr1ktime.o ATL_cgerk_Mlt16.o ATL_cgerk_axpy.o
	rm -f ATL_cger1.c
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_cger1.c
	$(cR1CC) -c $(CDEFS) $(cR1CFLAGS) -I$(R1Tdir) -o ATL_cger1.o \
           -DSCPLX -DATL_cger1_a1_x1_yX=ATL_UGERK ATL_cger1.c
	$(CLINKER) $(CLINKFLAGS) -o xcr1ktime cr1ktime.o \
                   ATL_cger1.o \
                   ATL_cgerk_Mlt16.o ATL_cgerk_axpy.o \
                   $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xcr1ktime -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
#
# This target calls the GER driver (which blocks), with cache flushing
#
cr1time : cr1time.o xr1hgen
	rm -f $(cobj) ATL_cger1.c atlas_cr1.h r1h.bak
	if [ -e "$(INCAdir)/atlas_cr1.h" ]; then \
           mv $(INCAdir)/atlas_cr1.h r1h.bak ; \
        fi
	./xr1hgen -p c -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_cger1.c
	$(MAKE) $(cobj)
	$(CLINKER) $(CLINKFLAGS) -o xcr1time cr1time.o \
                   $(cobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xcr1time -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_cr1.h ; \
        fi

cr1ktest : cr1ktest.o xr1hgen
	rm -f $(ctstobj) ATL_cger1.c 
	if [ -e "$(INCAdir)/atlas_cr1.h" ]; then \
           mv $(INCAdir)/atlas_cr1.h r1h.bak ; \
        fi
	./xr1hgen -p c -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_cger1.c
	$(MAKE) $(ctstobj)
	$(CLINKER) $(CLINKFLAGS) -o xcr1ktest cr1ktest.o \
                   $(ctstobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xcr1ktest -m 1 $(Mt) -n 1 $(Nt) -y 1 $(incy) \
                  -l 1 $(ldat) $(align)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_cr1.h ; \
        fi


cr1time.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX $(mySRCdir)/r1ktime.c
cr1ktime.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX -DTIME_KERNEL $(mySRCdir)/r1ktime.c
cr1ktest.o : $(mySRCdir)/r1ktest.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX $(mySRCdir)/r1ktest.c
cr1test.o : $(mySRCdir)/r1test.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX $(mySRCdir)/r1test.c
cger1tune.o : $(mySRCdir)/ger1tune.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX $(mySRCdir)/ger1tune.c
ATL_cger.o : ATL_ger.c atlas_cr1.h
	$(ICC) -c $(ICCFLAGS) -o $@ -DSCPLX -I$(R1Tdir) ATL_ger.c
ATL_cgerk_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DSCPLX $(r1SRCdir)/ATL_gerk_axpy.c
ATL_cgerck_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -DConj_ -I$(R1Tdir) \
                       -o $@ -DSCPLX $(r1SRCdir)/ATL_gerk_axpy.c
ATL_cgerk_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DSCPLX $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_cgerck_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(SKC) -c $(SKCFLAGS) $(CDEFS) -DConj_ -I$(R1Tdir) \
                       -o $@ -DSCPLX $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_cger1.o : ATL_cger1.c
	$(cR1CC) -c $(CDEFS) $(cR1CFLAGS) -I$(R1Tdir) -o $@ -DSCPLX ATL_cger1.c
ATL_cger1c.o : ATL_cger1.c
	$(cR1CC) -c $(CDEFS) $(cR1CFLAGS) -I$(R1Tdir) -o $@ -DSCPLX -DConj_ \
              -DATL_UGERK=ATL_GERK_Conj ATL_cger1.c

res/cR1SUMM : xr1ksearch force_build
	./xr1ksearch -p c
	$(MAKE) cinstall
cr1outd:
	mkdir cr1outd
cinstall : xccobj cr1outd xr1hgen force_build
	rm -f cr1outd/*
	if [ -s "$(GR1dir)/Make_cr1" ]; then \
           cd  $(GR1dir) ; $(MAKE) -f Make_cr1 killall ; \
           rm -f $(GR1dir)/Make_cr1 ; \
        fi
	./xr1hgen -p c -F res/cR1SUMM -d cr1outd
	cp -f cr1outd/Make_cr1 $(GR1dir)/.
	cp cr1outd/*.c $(GR1dir)/.
	cp -f cr1outd/*.h $(INCAdir)/.
	if [ -s "xccobj.exe" ] ; then \
	   cp -f xccobj.exe $(GR1dir)/. ; \
	elif [ -s "xccobj" ]; then \
           cp -f xccobj $(GR1dir)/. ; \
        fi
	touch $(INCSdir)/atlas_r1.h
	rm -f cr1outd/* ; rmdir cr1outd

CASES/ATL_cssegen.c : xr1gen_sse force_build
	rm -f CASES/ATL_cssegen.c
	./xr1gen_sse $(ssegenargs) > CASES/ATL_cssegen.c

cl1lib :
	cd $(BINdir) ; $(MAKE) cl1lib
ctstlib :
	cd $(BINdir) ; $(MAKE) ctstlib
ztstobj = ATL_zger1.o  ATL_zger1c.o  \
          ATL_zgerk_Mlt16.o ATL_zgerk_axpy.o \
          ATL_zgerck_Mlt16.o ATL_zgerck_axpy.o 
zobj = ATL_zger.o ATL_zger1.o \
       ATL_zgerk_Mlt16.o ATL_zgerk_axpy.o \
       ATL_zgerck_Mlt16.o ATL_zgerck_axpy.o 

#
# This target is for directly timing a call to the kernel, w/o cache flushing
#
zr1ktime : zr1ktime.o ATL_zgerk_Mlt16.o ATL_zgerk_axpy.o
	rm -f ATL_zger1.c
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_zger1.c
	$(zR1CC) -c $(CDEFS) $(zR1CFLAGS) -I$(R1Tdir) -o ATL_zger1.o \
           -DDCPLX -DATL_zger1_a1_x1_yX=ATL_UGERK ATL_zger1.c
	$(CLINKER) $(CLINKFLAGS) -o xzr1ktime zr1ktime.o \
                   ATL_zger1.o \
                   ATL_zgerk_Mlt16.o ATL_zgerk_axpy.o \
                   $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xzr1ktime -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
#
# This target calls the GER driver (which blocks), with cache flushing
#
zr1time : zr1time.o xr1hgen
	rm -f $(zobj) ATL_zger1.c atlas_zr1.h r1h.bak
	if [ -e "$(INCAdir)/atlas_zr1.h" ]; then \
           mv $(INCAdir)/atlas_zr1.h r1h.bak ; \
        fi
	./xr1hgen -p z -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_zger1.c
	$(MAKE) $(zobj)
	$(CLINKER) $(CLINKFLAGS) -o xzr1time zr1time.o \
                   $(zobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xzr1time -m $(M) -n $(N) -l $(lda) \
                  -C $(flushKB) $(align) $(tflags)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_zr1.h ; \
        fi

zr1ktest : zr1ktest.o xr1hgen
	rm -f $(ztstobj) ATL_zger1.c 
	if [ -e "$(INCAdir)/atlas_zr1.h" ]; then \
           mv $(INCAdir)/atlas_zr1.h r1h.bak ; \
        fi
	./xr1hgen -p z -l $(l1mul) -x $(xu) -y $(yu) -f $(iflag)
	cp -f $(mySRCdir)/CASES/$(r1rout) ATL_zger1.c
	$(MAKE) $(ztstobj)
	$(CLINKER) $(CLINKFLAGS) -o xzr1ktest zr1ktest.o \
                   $(ztstobj) $(TESTlib) $(ATLASlib)
	$(ATLRUN) $(R1Tdir) xzr1ktest -m 1 $(Mt) -n 1 $(Nt) -y 1 $(incy) \
                  -l 1 $(ldat) $(align)
	if [ -e "r1h.bak" ]; then \
           mv r1h.bak $(INCAdir)/atlas_zr1.h ; \
        fi


zr1time.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX $(mySRCdir)/r1ktime.c
zr1ktime.o : $(mySRCdir)/r1ktime.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX -DTIME_KERNEL $(mySRCdir)/r1ktime.c
zr1ktest.o : $(mySRCdir)/r1ktest.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX $(mySRCdir)/r1ktest.c
zr1test.o : $(mySRCdir)/r1test.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX $(mySRCdir)/r1test.c
zger1tune.o : $(mySRCdir)/ger1tune.c
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX $(mySRCdir)/ger1tune.c
ATL_zger.o : ATL_ger.c atlas_zr1.h
	$(ICC) -c $(ICCFLAGS) -o $@ -DDCPLX -I$(R1Tdir) ATL_ger.c
ATL_zgerk_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DDCPLX $(r1SRCdir)/ATL_gerk_axpy.c
ATL_zgerck_axpy.o : $(r1SRCdir)/ATL_gerk_axpy.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -DConj_ -I$(R1Tdir) \
                       -o $@ -DDCPLX $(r1SRCdir)/ATL_gerk_axpy.c
ATL_zgerk_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -I$(R1Tdir) \
                       -o $@ -DDCPLX $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_zgerck_Mlt16.o : $(r1SRCdir)/ATL_gerk_Mlt16.c
	$(DKC) -c $(DKCFLAGS) $(CDEFS) -DConj_ -I$(R1Tdir) \
                       -o $@ -DDCPLX $(r1SRCdir)/ATL_gerk_Mlt16.c
ATL_zger1.o : ATL_zger1.c
	$(zR1CC) -c $(CDEFS) $(zR1CFLAGS) -I$(R1Tdir) -o $@ -DDCPLX ATL_zger1.c
ATL_zger1c.o : ATL_zger1.c
	$(zR1CC) -c $(CDEFS) $(zR1CFLAGS) -I$(R1Tdir) -o $@ -DDCPLX -DConj_ \
              -DATL_UGERK=ATL_GERK_Conj ATL_zger1.c

res/zR1SUMM : xr1ksearch force_build
	./xr1ksearch -p z
	$(MAKE) zinstall
zr1outd:
	mkdir zr1outd
zinstall : xccobj zr1outd xr1hgen force_build
	rm -f zr1outd/*
	if [ -s "$(GR1dir)/Make_zr1" ]; then \
           cd  $(GR1dir) ; $(MAKE) -f Make_zr1 killall ; \
           rm -f $(GR1dir)/Make_zr1 ; \
        fi
	./xr1hgen -p z -F res/zR1SUMM -d zr1outd
	cp -f zr1outd/Make_zr1 $(GR1dir)/.
	cp zr1outd/*.c $(GR1dir)/.
	cp -f zr1outd/*.h $(INCAdir)/.
	if [ -s "xccobj.exe" ] ; then \
	   cp -f xccobj.exe $(GR1dir)/. ; \
	elif [ -s "xccobj" ]; then \
           cp -f xccobj $(GR1dir)/. ; \
        fi
	touch $(INCSdir)/atlas_r1.h
	rm -f zr1outd/* ; rmdir zr1outd

CASES/ATL_zssegen.c : xr1gen_sse force_build
	rm -f CASES/ATL_zssegen.c
	./xr1gen_sse $(ssegenargs) > CASES/ATL_zssegen.c

zl1lib :
	cd $(BINdir) ; $(MAKE) zl1lib
ztstlib :
	cd $(BINdir) ; $(MAKE) ztstlib

$(BINdir)/xsubstr :
	cd $(BINdir) ; $(MAKE) xsubstr
$(BINdir)/xccobj :
	cd $(BINdir) ; $(MAKE) xccobj
xccobj : $(BINdir)/xccobj
	-@ cp -f $(BINdir)/xccobj .
	-@ cp -f $(BINdir)/xccobj.exe .
