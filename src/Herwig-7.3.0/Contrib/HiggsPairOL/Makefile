# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files

OPENLOOPSPREFIX=/Users/user/phd/herwig_gh/./opt/OpenLoops-2.1.2
THEPEGINCLUDE=-I/opt/homebrew/Cellar/thepeg/2.3.0_1//include
HERWIGINCLUDE=-I/Users/user/phd/herwig_gh/./include/
HERWIGINSTALL=/Users/user/phd/herwig_gh/.
GSLINCLUDE=-I/opt/homebrew/Cellar/gsl/2.8//include 
FASTJETINCLUDE=-I/opt/homebrew/Cellar/fastjet/3.4.3/include
LDFLAGS= 
SHARED_FLAG=-bundle -Wl,-undefined,dynamic_lookup 
FASTJETLIB = `echo $(FASTJETINCLUDE) | sed "s/-I//" | sed "s/ //"`/../lib/
LDFLAGS= -L$(FASTJETLIB) -lfastjet 
INCLUDE = $(THEPEGINCLUDE) $(HERWIGINCLUDE) $(GSLINCLUDE) $(FASTJETINCLUDE)
OPENLOOPSLIB=$(OPENLOOPSPREFIX)/lib
OPENLOOPSPROCLIB=$(OPENLOOPSPREFIX)/proclib
OL = -L$(OPENLOOPSLIB) -lopenloops -L$(OPENLOOBSPROCLIB) -lopenloops
#
# C++ flags
# 
CXX=/usr/bin/clang++ -std=c++14
CXXFLAGS=-O2 -DBOOST_UBLAS_NDEBUG
FC=/opt/homebrew/bin/gfortran
FCFLAGS =
FCLIBS= -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15/../../.. -lemutls_w -lheapt_w -lgfortran -lquadmath

default : HiggsPair.so AlpGenHandlerOL.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c $< -o $@  

%.o : %.f
	$(FC) -fPIC $(CPPFLAGS) $(INCLUDE) $(FCFLAGS) -c $< -o $@  

HiggsPair.so : HiggsPair.o MEHiggsPairOL.o MEHiggsPairJet.o scalint.o 
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) $(SHARED_FLAG) $(OL) $(LDFLAGS)\
	$(FCLIBS) $^ -o $@

AlpGenHandlerOL.so : AlpGenHandlerOL.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) $(SHARED_FLAG) $(LDFLAGS) \
	-Wl,-undefined -Wl,dynamic_lookup \
	AlpGenHandlerOL.o -o AlpGenHandlerOL.so

clean:
	rm -f HiggsPair.so MEHiggsPairJet.o MEHiggsPairOL.o scalint.o AlpGenHandlerOL.o
	
install: 
	cp *.so $(HERWIGINSTALL)/lib/Herwig
