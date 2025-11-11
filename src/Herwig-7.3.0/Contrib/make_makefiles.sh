#!/bin/bash
# script to put the HERWIG and THEPEG paths into the simple makefiles,
# for everything else the user is on their own!
# loop over all the files

OPENLOOPSPREFIX="/Users/user/phd/herwig_gh/./opt/OpenLoops-2.1.2"
THEPEGINCLUDE="-I/opt/homebrew/Cellar/thepeg/2.3.0_1//include"
GSLINCLUDE="-I/opt/homebrew/Cellar/gsl/2.8//include"
FASTJETINCLUDE="-I/opt/homebrew/Cellar/fastjet/3.4.3/include"
HERWIGINCLUDE=-I/Users/user/phd/herwig_gh/./include/
HERWIGINSTALL="/Users/user/phd/herwig_gh/."
FC="/opt/homebrew/bin/gfortran"
FCLIBS=" -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15 -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc -L/opt/homebrew/Cellar/gcc/15.2.0/bin/../lib/gcc/current/gcc/aarch64-apple-darwin24/15/../../.. -lemutls_w -lheapt_w -lgfortran -lquadmath"
CXX="/usr/bin/clang++ -std=c++14"
CXXFLAGS="-O2 -DBOOST_UBLAS_NDEBUG"
LDFLAGS="" 
SHARED_FLAG="-bundle -Wl,-undefined,dynamic_lookup" 

for i in *
do
# if a directory
  if [ -d $i ]; then
# check input files exists
      file=$i/Makefile.in
      if [ -e $file ]; then
	  file2=`echo $file | sed s!\.in!!`
	  echo 'Making ' $file2
	  sed "s!THEPEGINCLUDE *\=!THEPEGINCLUDE=$THEPEGINCLUDE!" < $file | \
	  sed "s!OPENLOOPSPREFIX *\=!OPENLOOPSPREFIX=$OPENLOOPSPREFIX!" | \
	  sed "s!FC *\=!FC=$FC!" | \
	  sed "s!FCLIBS *\=!FCLIBS=$FCLIBS!" | \
          sed "s!CXX *\=!CXX=$CXX!" | \
          sed "s!SHARED_FLAG *\=!SHARED_FLAG=$SHARED_FLAG!" | \
          sed "s!LDFLAGS *\=!LDFLAGS=$LDFLAGS!" | \
          sed "s!CXXFLAGS *\=!CXXFLAGS=$CXXFLAGS!" | \
	  sed "s!HERWIGINCLUDE *\=!HERWIGINCLUDE=$HERWIGINCLUDE!" | \
          sed "s!HERWIGINSTALL *\=!HERWIGINSTALL=$HERWIGINSTALL!" | \
	  sed "s!GSLINCLUDE *\=!GSLINCLUDE=$GSLINCLUDE!" | \
	  sed "s!FASTJETINCLUDE *\=!FASTJETINCLUDE=$FASTJETINCLUDE!" > $file2
      fi
  fi
done
