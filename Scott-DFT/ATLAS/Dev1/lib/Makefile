LIBTOOL=libtool
DESTDIR=/home/bylaska/Codes/Scott-DFT
INCINSTdir=$(DESTDIR)/include
LIBINSTdir=$(DESTDIR)/lib
include Make.inc
mySRCdir = $(SRCdir)/lib
#
# override with libatlas.so only when atlas is built to one lib
#
DYNlibs = liblapack.so libf77blas.so libcblas.so libatlas.so 
PTDYNlibs = liblapack.so libptf77blas.so libptcblas.so libatlas.so 
CDYNlibs = liblapack.so libcblas.so libatlas.so 
CPTDYNlibs = liblapack.so libptcblas.so libatlas.so 

VER=3.9.23
tmpd = RCW_tMp
tarnam = atlas$(VER)_$(ARCH)
tar : tarfile
tarfile : $(tarnam).tgz
$(tarnam).tgz :
	mkdir $(ARCH)
	cd $(ARCH) ; mkdir include lib
	cp $(TOPdir)/doc/LibReadme.txt $(ARCH)/README
	cp $(TOPdir)/Make.$(ARCH) $(ARCH)/.
	cp $(BINdir)/INSTALL_LOG/SUMMARY.LOG $(ARCH)/.
	cp $(INCSdir)/cblas.h $(ARCH)/include/.
	cp $(INCSdir)/clapack.h $(ARCH)/include/.
	cp $(LIBdir)/libatlas.a $(ARCH)/lib/.
	cp $(LIBdir)/libf77blas.a $(ARCH)/lib/.
	cp $(LIBdir)/libcblas.a $(ARCH)/lib/.
	cp $(LIBdir)/liblapack.a $(ARCH)/lib/.
	- cp $(LIBdir)/libptcblas.a $(ARCH)/lib/.
	- cp $(LIBdir)/libptf77blas.a $(ARCH)/lib/.
	$(TAR) cf $(tarnam).tar $(ARCH)
	rm -rf $(ARCH)
	$(GZIP) --best $(tarnam).tar
	mv $(tarnam).tar.gz $(tarnam).tgz

# ===================================================================
# The following commands are to build dynamic/shared objects on Linux
# ===================================================================
ptshared :
	- rm -f libatlas.so liblapack.so
	$(MAKE) libatlas.so libptf77blas.so libf77blas.so \
                libptcblas.so libcblas.so liblapack.so
shared :
	- rm -f libatlas.so liblapack.so
	$(MAKE) libatlas.so libf77blas.so libcblas.so liblapack.so
cptshared :
	- rm -f libatlas.so libclapack.so
	$(MAKE) libatlas.so libclapack.so libptcblas.so libcblas.so
cshared :
	- rm -f libatlas.so libclapack.so
	$(MAKE) libatlas.so libclapack.so libcblas.so

libatlas.so : libatlas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libatlas.so -o libatlas.so \
           -rpath-link $(LIBINSTdir) \
           --whole-archive libatlas.a --no-whole-archive -lc $(LIBS)
liblapack.so : liblapack.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/liblapack.so \
            -o liblapack.so  -rpath-link $(LIBINSTdir) --whole-archive \
           liblapack.a --no-whole-archive $(F77SYSLIB)
libclapack.so : libclapack.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libclapack.so \
           -o libclapack.so -rpath-link $(LIBINSTdir) \
           --whole-archive libclapack.a --no-whole-archive -lc $(LIBS)
libptf77blas.so : libptf77blas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libptf77blas.so \
           -o libptf77blas.so -rpath-link $(LIBINSTdir) \
           --whole-archive libptf77blas.a --no-whole-archive $(F77SYSLIB)
libf77blas.so : libf77blas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libf77blas.so \
           -o libf77blas.so -rpath-link $(LIBINSTdir) \
           --whole-archive libf77blas.a --no-whole-archive $(F77SYSLIB)
libptcblas.so : libptcblas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libptcblas.so \
           -o libptcblas.so -rpath-link $(LIBINSTdir) \
          --whole-archive libptcblas.a --no-whole-archive -lc $(LIBS)
libcblas.so : libcblas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libcblas.so \
           -o libcblas.so -rpath-link $(LIBINSTdir) \
          --whole-archive libcblas.a --no-whole-archive -lc $(LIBS)

#
# Builds one shared lib from all ATLAS files
#
fat_ptshared :                              # threaded target
	ld $(LDFLAGS) -shared -soname libatlas.so -o libatlas.so \
        --whole-archive liblapack.a libptf77blas.a libptcblas.a libatlas.a \
           --no-whole-archive -lc $(F77SYSLIB) $(LIBS)
fat_shared :                                # unthreaded target
	ld $(LDFLAGS) -shared -soname libatlas.so -o libatlas.so \
        --whole-archive liblapack.a libf77blas.a libcblas.a libatlas.a \
           --no-whole-archive -lc $(F77SYSLIB) $(LIBS)
#
# Builds shared lib, not include fortran codes from LAPACK
#
fat_cptshared : libclapack.a                # threaded target
	ld $(LDFLAGS) -shared -soname libcatlas.so -o libcatlas.so \
        --whole-archive libclapack.a libptcblas.a libatlas.a \
        --no-whole-archive -lc $(LIBS)
fat_cshared : libclapack.a                  # unthreaded target
	ld $(LDFLAGS) -shared -soname libcatlas.so -o libcatlas.so \
        --whole-archive libclapack.a libcblas.a libatlas.a \
        --no-whole-archive -lc $(LIBS)

libclapack.a : liblapack.a
	rm -rf clapack libclapack.a
	mkdir clapack
	cd clapack ; ar x ../liblapack.a
	rm -f clapack/*f77wrap* clapack/*C2F*
	ar r libclapack.a clapack/ATL_* clapack/clapack_*
	rm -rf clapack

#  ============================================
#  The following commands build DLLs on Windows
#  ============================================
dlls:   libatlas.dll libf77blas.dll libcblas.dll liblapack.dll
ptdlls: libatlas.dll libptf77blas.dll libptcblas.dll liblapack.dll
cdlls: libatlas.dll libcblas.dll libclapack.dll
ptcdlls: libatlas.dll libptcblas.dll libclapack.dll

libatlas.dll : libatlas.a
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libatlas.dll \
           -o libatlas.dll -rpath-link $(LIBINSTdir) \
           --whole-archive libatlas.a --no-whole-archive -lc $(LIBS)
libcblas.dll : libcblas.a libatlas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libcblas.dll \
           -o libcblas.dll -rpath-link $(LIBINSTdir) \
          --whole-archive libcblas.a --no-whole-archive libatlas.dll -lc $(LIBS)
libf77blas.dll : libf77blas.a libatlas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libf77blas.dll \
           -o libf77blas.dll -rpath-link $(LIBINSTdir) --whole-archive \
           libf77blas.a --no-whole-archive libatlas.dll $(F77SYSLIB) -lc $(LIBS)
libptcblas.dll : libptcblas.a libatlas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libptcblas.dll \
           -o libptcblas.dll -rpath-link $(LIBINSTdir) --whole-archive \
           libptcblas.a --no-whole-archive libatlas.dll -lc $(LIBS)
libptf77blas.dll : libptf77blas.a libatlas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libptf77blas.dll \
           -o libptf77blas.dll -rpath-link $(LIBINSTdir) --whole-archive \
           libptf77blas.a --no-whole-archive libatlas.a $(F77SYSLIB) -lc $(LIBS)
liblapack.dll : liblapack.a libatlas.dll libf77blas.dll libcblas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/liblapack.dll \
            -o liblapack.dll  -rpath-link $(LIBINSTdir) --whole-archive \
           liblapack.a --no-whole-archive libf77blas.dll libcblas.dll \
           libatlas.dll $(F77SYSLIB) -lc $(LIBS)
libclapack.dll : libclapack.a libatlas.dll libcblas.dll
	ld $(LDFLAGS) -shared -soname $(LIBINSTdir)/libclapack.dll \
           -o libclapack.dll -rpath-link $(LIBINSTdir) --whole-archive \
           libclapack.a --no-whole-archive libcblas.dll libatlas.dll -lc $(LIBS)

#  =======================================================================
#  The following commands are to build dynamib libraries on OS X (in BETA)
#  =======================================================================
dylib :
	$(MAKE) libatlas.dylib libf77blas.dylib libcblas.dylib liblapack.dylib
ptdylib :
	$(MAKE) libatlas.dylib libptf77blas.dylib libptcblas.dylib \
                liblapack.dylib
cdylib :
	$(MAKE) libatlas.dylib libcblas.dylib libclapack.dylib
ptcdylib :
	$(MAKE) libatlas.dylib libptcblas.dylib libclapack.dylib

libatlas.dylib : libatlas.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../libatlas.a
	rm -f $(tmpd)ATL_[z,c,d,s]ref*.o
	cd $(tmpd) ; $(LIBTOOL) -dynamic -o ../libatlas.dylib \
        -install_name $(LIBINSTdir)/libatlas.dylib -current_version $(VER) \
        -compatibility_version $(VER) *.o $(LIBS)
	rm -rf $(tmpd)
libptcblas.dylib : libatlas.dylib libptcblas.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../libptcblas.a
	cd $(tmpd) ; $(LIBTOOL) -dynamic -o ../libptcblas.dylib \
           -install_name $(LIBINSTdir)/libptcblas.dylib \
           -compatibility_version $(VER) -current_version $(VER) \
           *.o ../libatlas.dylib $(LIBS)
	rm -rf $(tmpd)
libcblas.dylib : libatlas.dylib libcblas.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../libcblas.a
	cd $(tmpd) ; $(LIBTOOL) -dynamic -o ../libcblas.dylib \
           -install_name $(LIBINSTdir)/libcblas.dylib \
           -compatibility_version $(VER) -current_version $(VER) \
           *.o ../libatlas.dylib $(LIBS)
	rm -rf $(tmpd)
libptf77blas.dylib : libatlas.dylib libptf77blas.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../libptf77blas.a
	cd $(tmpd) ; $(LIBTOOL) -dynamic -o ../libptf77blas.dylib \
           -install_name $(LIBINSTdir)/libptf77blas.dylib \
           -compatibility_version $(VER) -current_version $(VER) \
           *.o ../libatlas.dylib $(F77SYSLIB) $(LIBS)
	rm -rf $(tmpd)
libf77blas.dylib : libatlas.dylib libf77blas.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../libf77blas.a
	cd $(tmpd) ; libtool -dynamic -o ../libf77blas.dylib \
           -install_name $(LIBINSTdir)/libf77blas.dylib \
           -compatibility_version 3.9.23 -current_version 3.9.23 \
           *.o ../libatlas.dylib $(F77SYSLIB) $(LIBS)
	rm -rf $(tmpd)
liblapack.dylib : libf77blas.dylib libcblas.dylib libatlas.dylib liblapack.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../liblapack.a
	cd $(tmpd) ; libtool -dynamic -o ../liblapack.dylib \
           -install_name $(LIBINSTdir)/liblapack.dylib \
           -compatibility_version $(VER) -current_version $(VER) \
           *.o ../libf77blas.dylib ../libcblas.dylib ../libatlas.dylib \
           $(F77SYSLIB) $(LIBS)
	rm -rf $(tmpd)
libclapack.dylib : libcblas.dylib libatlas.dylib liblapack.a
	rm -rf $(tmpd) ; mkdir $(tmpd)
	cd $(tmpd) ; ar x ../liblapack.a
	rm -f $(tmpd)/*C2F $(tmpd)/*f77wrap*
	cd $(tmpd) ; libtool -dynamic -o ../libclapack.dylib \
           -install_name $(LIBINSTdir)/libclapack.dylib \
           -compatibility_version $(VER) -current_version $(VER) \
           *.o ../libcblas.dylib ../libatlas.dylib $(LIBS)
	rm -rf $(tmpd)
xtst_lp: $(DYNlibs)
	$(ICC) $(CDEFS) -o $@ $(mySRCdir)/qr.c $(DYNlibs) -Wl,--rpath ./
xtst : $(DYNlibs)
	$(ICC) $(CDEFS) -o $@ $(mySRCdir)/test_dynlink.c $(DYNlibs) \
           -Wl,--rpath ./

xtry_lp:
	$(ICC) $(CDEFS) -o $@ $(mySRCdir)/qr.c libatlas.so -Wl,--rpath ./
xtry :
	$(ICC) $(CDEFS) -o $@ $(mySRCdir)/test_dynlink.c libatlas.so \
           -Wl,--rpath ./
xtry_c :
	$(ICC) $(CDEFS) -o $@ $(mySRCdir)/test_dynlink.c libcatlas.so \
           -Wl,--rpath ./
