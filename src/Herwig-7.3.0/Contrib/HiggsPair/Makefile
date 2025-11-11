# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files
THEPEGINCLUDE=-I/opt/homebrew/Cellar/thepeg/2.3.0_1//include
HERWIGINCLUDE=-I/Users/user/phd/herwig_gh/./include/
HERWIGINSTALL=/Users/user/phd/herwig_gh/.
GSLINCLUDE=-I/opt/homebrew/Cellar/gsl/2.8//include 
LDFLAGS= 
SHARED_FLAG=-bundle -Wl,-undefined,dynamic_lookup 
INCLUDE = $(THEPEGINCLUDE) $(HERWIGINCLUDE) $(GSLINCLUDE)
#
# C++ flags
# 
CXX=/usr/bin/clang++ -std=c++14
CXXFLAGS=-O2 -DBOOST_UBLAS_NDEBUG
FC=/opt/homebrew/bin/gfortran
FCFLAGS =
FCLIBS= -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15/../../.. -lemutls_w -lheapt_w -lgfortran -lquadmath

default : MEHiggsPair.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c $< -o $@  

%.o : %.f
	$(FC) -fPIC $(CPPFLAGS) $(INCLUDE) $(FCFLAGS) -c $< -o $@  

MEHiggsPair.so : MEHiggsPair.o scalint.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) $(SHARED_FLAG) $(LDFLAGS) \
	$(FCLIBS) $^ -o $@

clean:
	rm -f MEHiggsPair.so MEHiggsPair.o scalint.o

install:
	cp *.so $(HERWIGINSTALL)/lib/Herwig
        
        
        
        