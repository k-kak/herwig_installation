# -*- Makefile -*- (for emacs)

#
# This Makefile is intended for compiling Herwig++ plugins
# You can find plugins here: INSERT URL
#
# This Makefile received very little testing, 
# any bug reports are very welcome!
#

# location of include files
THEPEGINCLUDE=-I/opt/homebrew/Cellar/thepeg/2.3.0_1//include
GSLINCLUDE=-I/opt/homebrew/Cellar/gsl/2.8//include
HERWIGINCLUDE=-I/Users/user/phd/herwig_gh/./include/
INCLUDE = $(THEPEGINCLUDE) $(GSLINCLUDE) $(HERWIGINCLUDE)
#
# C++ flags
# 
CXX=/usr/bin/clang++ -std=c++14 
CXXFLAGS=-O2 -DBOOST_UBLAS_NDEBUG 
LDFLAGS= 
SHARED_FLAG=-bundle -Wl,-undefined,dynamic_lookup 

ALLCCFILES=$(shell echo *.cc)

default : HwShowerVeto.so

%.o : %.cc %.h
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) -c -shared $< -o $@  

HwShowerVeto.so: NonBShowerVeto.o
	$(CXX) -fPIC $(CPPFLAGS) $(INCLUDE) $(CXXFLAGS) \
		NonBShowerVeto.o $(SHARED_FLAG) $(LDFLAGS)  -o HwShowerVeto.so

clean:
	rm -f $(ALLCCFILES:.cc=.o) HwShowerVeto.so
