# Fighting Herwig Bootstrap

This repository contains a compiled and functional version of Herwig7 bootstrapped installation, aiming to document steps that were needed to make the installation succeed on Mac OS 15.7.2, Apple M4 chip device. 

Authors: Mohamed Aly (mohamed.aly@cern.ch), Karim Kandeel (karim.kandeel@cern.ch)

## Installation directory

Create a new directory to house herwig installation on your device. We refer to this as <my_herwig> from hereonafter.

```
mkdir /path/to/<my_herwig>
cd /path/to/<my_herwig>
```

To visualise directory structure, we recomment you also install `tree`:
```
brew install tree
```

## Python

The compilation was carried out by using a python virtual environment set by `uv`:

0. If you don't have `uv`, install with

```
curl -LsSf https://astral.sh/uv/install.sh | sh
```

1. Initialise a `uv` project:
```
uv init
```
2. Set the python version to `3.9` (this is chosen arbitrarily given experience with installing other HEP event generators):
```
uv python install 3.9
```

This installed `3.9.6` for me, but the process worked also with `3.9.24`.

3. Create a virtual environment:

```
uv venv
```

3. Check if you have a `venv` created from `uv` in `.venv/` and activate it

```bash
tree ./venv
```
returns:

```bash
.venv/
├── bin
└── lib
    └── python3.9
        └── site-packages
```
Then run:

```bash
source .venv/bin/activate
```

to enter the environment. This environment is only valid in the directory where you created it.

4. Check python version and binaries:

```bash
which python && python --version
```
returns:

```bash
/path/to/<my_herwig>/.venv/bin/python
Python 3.9.24
```

5. Install some packages we will need:

```bash
uv add six cython lib 
```

6. Sync installations
```bash
uv sync
```

7. Check what we have
```bash
tree .venv/
```

which gives:

```bash
.venv/
├── bin
└── lib
    └── python3.9
        └── site-packages
            ├── Cython
            │   ├── Build
            │   │   └── Tests
            │   ├── Compiler
            │   │   └── Tests
            │   ├── Debugger
            │   │   └── Tests
            │   ├── Distutils
            │   ├── Includes
            │   │   ├── cpython
            │   │   ├── libc
            │   │   ├── libcpp
            │   │   ├── numpy
            │   │   └── posix
            │   ├── Plex
            │   ├── Runtime
            │   ├── Tempita
            │   ├── Tests
            │   └── Utility
            ├── cython-3.2.0.dist-info
            ├── lib
            ├── lib-4.0.0.dist-info
            ├── pyximport
            └── six-1.17.0.dist-info

29 directories
```

**Now our python environment is set and we can proceed with installation!**

## The Herwig Bootstrap script

```bash
wget https://herwig.hepforge.org/downloads/herwig-bootstrap
```

will install the bootstrapping script in your current directory. This is a `python` script which organises the installation.

It is meant to be ran as an executable, so we need to set permission:

```
chmod +x herwig-bootstrap
```

The Herwig docs suggest this should run with:

```bash

./herwig-bootstrap -j4 /path/to/<my_herwig>
```

**but this didn't go well, with errors everywhere.**

### Errors we encountered:
- In compiling `fastjet` library, we encountered compilation erros regarding non-existent class attrbibutes in a print statement:

    ```bash
    ./ProtoJet.hpp:198:51: error: no member named '_Et' in 'ProtoJet<Item>'
    198 |   os<<"y phi Et = ("<<_y<<", "<<_phi<<", "<<this->_Et<<")"<<std::endl;
        |   
    ```

    I was able to circumvent this by commenting the problematic line in:
    ```
    src/fastjet-3.4.0/plugins/D0RunIICone/ProtoJet.hpp
    ```

    which in-theory is harmless (it is just a print out), but shouldn't happen.

- Our `cmake` version was `4.1.2` the latest from brew, obtained with:
    ```
    brew install cmake
    ```
    This did not make `HEPMC3` installation happy:
    
    ```bash
    cmake /Users/user/phd/herwig/./src/HepMC3-3.2.5 -DCMAKE_INSTALL_PREFIX=/Users/user/phd/herwig/. -DHEPMC3_ENABLE_ROOTIO:BOOL=OFF -DHEPMC3_ENABLE_PYTHON:BOOL=OFF
    CMake Error at CMakeLists.txt:1 (cmake_minimum_required):
    Compatibility with CMake < 3.5 has been removed from CMake.

    Update the VERSION argument <min> value.  Or, use the <min>...<max> syntax
    to tell CMake that the project requires at least <min> but has been updated
    to work with policies introduced by <max> or earlier.

    Or, add -DCMAKE_POLICY_VERSION_MINIMUM=3.5 to try configuring anyway.
    ```

    so we took the adivce and tried to configure anyway:

    ```bash
    export CMAKE_POLICY_VERSION_MINIMUM=3.5
    ```

    This, again, move us forward.

- `lhapdf` installtion then failed, because it didn't propagate correctly our python environment, and couldn't find the `lib` library:
    ```bash
    checking for Python include path... -I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/include/python3.9
    checking for Python library path... -L/Applications/Xcode.app/Contents/Developer/Library/Frameworks/Python3.framework/Versions/3.9/lib -lpython3.9
    checking for Python site-packages path... /Users/user/phd/herwig/.venv/lib/python3.9/site-packages
    checking for Python platform specific site-packages path... 
    checking python extra libraries... -ldl -lSystem  -framework CoreFoundation 
    checking python extra linking flags... -Wl,-stack_size,1000000  -framework CoreFoundation Python3.framework/Versions/3.9/Python3
    checking consistency of all components of python development environment... no
    configure: error: in `/Users/user/phd/herwig/src/LHAPDF-6.5.3':
    configure: error: 
    Could not link test program to Python. Maybe the main Python library has been
    installed in some non-standard library path. If so, pass it to configure,
    via the LIBS environment variable.
    Example: ./configure LIBS="-L/usr/non-standard-path/python/lib"
    ============================================================================
    ERROR!
    You probably have to install the development version of the Python package
    for your distribution.  The exact name of this package varies among them.
    ============================================================================

    See `config.log' for more details
    ```

    This also should be fixed. It occured if we try to compile `lhapdf` inside `src/` ourselves too, so the issue is in the `lhapdf` `MakeFile` and how `Python` is configured for it.

    The solution here was to install `lhapdf` from `brew`:

    ```
    brew install lhapdf
    ```

    **Note: if you don't have `davidchall/hep` tapped with brew, you need to first run:
    ```
    brew tap davidchall/hep
    ```

    Then, to tell the bootstrap script to use this `lhapdf` installtion, we ran:

    ```bash
    ./herwig-bootstrap -j 8 --lhapdf-location=/opt/homebrew/opt/lhapdf/ .
    ```
- We immediatly ran into the same error with configuring `python` for installting both `rivet` and `yoda`. So we installed both with brew and continued:
    ```bash
    brew install davidchall/hep/yoda # just install yoda didn't work
    brew install --formula rivet # just install rivet didn't work
    ./herwig-bootstrap -j 8 --with-lhapdf=/opt/homebrew/opt/lhapdf/ --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ .
    ```

- `boost` complained about linking erros, so we also installed it with `brew`. At this point, we decided it is best to install also `gsl` with brew:
    ```
    brew install boost gsl
    ./herwig-bootstrap -j 8 --with-lhapdf=/opt/homebrew/opt/lhapdf/ --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ --with-boost=/opt/homebrew/Cellar/boost/1.89.0/ --with-gsl=/opt/homebrew/Cellar/gsl/2.8/ .
    ```

- `THEPEG` the couldn't compile wit the default `fortran` flags, raising the error:
    ```fortran
    Making all in Looptools
    PPF77    D/libHwLooptoolsXFC_la-D0funcC.lo
    PPF77    D/libHwLooptoolsXFC_la-D0func.lo
    PPF77    util/libHwLooptoolsXFC_la-ffcxs3.lo
    PPF77    util/libHwLooptoolsXFC_la-ffcxs4.lo
    PPF77    util/libHwLooptoolsXFC_la-ffbndc.lo
    PPF77    util/libHwLooptoolsXFC_la-ffdcxs.lo
    PPF77    E/libHwLooptoolsCFC_la-Eget.lo
    PPF77    E/libHwLooptoolsCFC_la-E0func.lo
    D/D0func.F:486:35:

    486 |         parameter (nz2 = -2147483648)   ! O'20000000000'
        |                                            1
    Error: Integer too big for its kind at (1). This check can be disabled with the option '-fno-range-check'
    D/D0func.F:491:27:

    491 |         parameter (nz2p1234 = nz2 + p1234)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:495:27:

    495 |         parameter (nz2p1243 = nz2 + p1243)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:499:27:

    499 |         parameter (nz2p2134 = nz2 + p2134)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:503:27:

    503 |         parameter (nz2p2143 = nz2 + p2143)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:507:27:

    507 |         parameter (nz2p3214 = nz2 + p3214)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:511:27:

    511 |         parameter (nz2p4213 = nz2 + p4213)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0func.F:525:29:

    525 |      &    nz3p1234, nz2p1234, nz3p2143, nz2p1234,
        |                             1
    Error: Symbol 'nz2p1234' must be a PARAMETER in DATA statement at (1)
    D/D0funcC.F:428:35:

    428 |         parameter (nz2 = -2147483648)   ! O'20000000000'
        |                                            1
    Error: Integer too big for its kind at (1). This check can be disabled with the option '-fno-range-check'
    D/D0funcC.F:433:27:

    433 |         parameter (nz2p1234 = nz2 + p1234)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:437:27:

    437 |         parameter (nz2p1243 = nz2 + p1243)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:441:27:

    441 |         parameter (nz2p2134 = nz2 + p2134)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:445:27:

    445 |         parameter (nz2p2143 = nz2 + p2143)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:449:27:

    449 |         parameter (nz2p3214 = nz2 + p3214)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:453:27:

    453 |         parameter (nz2p4213 = nz2 + p4213)
        |                                  1
    Error: Parameter 'nz2' at (1) has not been declared or is a variable, which does not reduce to a constant expression
    D/D0funcC.F:467:29:

    467 |      &    nz3p1234, nz2p1234, nz3p2143, nz2p1234,
        |                             1
    Error: Symbol 'nz2p1234' must be a PARAMETER in DATA statement at (1)
    D/D0funcC.F:951:39:

    711 |      &    cLi2omrat2(q3, s, q4, m4) - pi6 +
        |                              2         
    ......
    951 |      &             cLi2omrat2(q4, s, -1D0, -1D0)
        |                                       1
    Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8)).
    D/D0funcC.F:951:45:

    711 |      &    cLi2omrat2(q3, s, q4, m4) - pi6 +
        |                                  2           
    ......
    951 |      &             cLi2omrat2(q4, s, -1D0, -1D0)
        |                                             1
    Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8)).
    D/D0func.F:1212:36:

    858 |      &        2*(Li2omrat2(q2, s, q4, t) -
        |                                    2       
    ......
    1212 |           dilogs = Li2omrat2(q3, t, x43(4), x43(2)) +
        |                                           1
    Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
    D/D0func.F:1212:44:

    858 |      &        2*(Li2omrat2(q2, s, q4, t) -
        |                                        2           
    ......
    1212 |           dilogs = Li2omrat2(q3, t, x43(4), x43(2)) +
        |                                                   1
    Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
    make[1]: *** [D/libHwLooptoolsXFC_la-D0funcC.lo] Error 1
    make[1]: *** Waiting for unfinished jobs....
    make[1]: *** [D/libHwLooptoolsXFC_la-D0func.lo] Error 1
    make: *** [all-recursive] Error 1
    ```

    In particular, the following starts the error stack:

    ```fortran
    Error: Integer too big for its kind at (1). This check can be disabled with the option '-fno-range-check'
    ```

    which we were able to fix by passing the flag during the fortran compilation in the bootstrap script that we are running.
    To do this, we modify the function:
    ```python
    def compile(config_flags=[]) :
    if os.access('Makefile',os.R_OK):
        print('Makefile exists, skipping configure')
    else:
        if not os.access('configure', os.X_OK):
            os.chmod('configure', 0o755)
        flags = ["./configure","--prefix="+base_dir]
        flags += config_flags
        check_call(flags)
    check_call(["make","-s","-j%s" % opts.ncore, "FFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check", "FCFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check"])

    # TODO ?? skip here if some file already exists in /lib or /bin ??
    check_call(["make","-s","install"])
    os.chdir(src_dir)
    ```

    and in particular the line:

    ```python
    check_call(["make","-s","-j%s" % opts.ncore]) 
    ```
    to:
    ```python
    check_call(["make","-s","-j%s" % opts.ncore, "FFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check", "FCFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check"])
    ```

    then, clean the build:
    ```bash
    cd src/Herwig-7.3.0/
    make clean
    ```
    and running again:

    ```
    ./herwig-bootstrap -j 8 --with-lhapdf=/opt/homebrew/opt/lhapdf/ --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ .
    ```
    and this finally finishes the compile:
    ```bash
    ################ /  / ^^/ ##########################
    ############### /--/   / ###########################
    ############## /  /   / ############################

    Herwig 7 bootstrap was successful.

    $ source bin/activate

        activates all required environment variables.

    $ deactivate

        returns to the original environment variables.
    ```

## TL;DR

We found the following steps to be most stable to get everything running without too many hacks (while still hacky):

```bash
brew install thepeg lhapdf rivet yoda boost gsl fastjet
mkdir bin # needed for qgraf error, see below
ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib /opt/homebrew/opt/gsl/lib/libgsl.27.dylib # needed for thepeg-gsl dylib incompatibility, see below
./herwig-bootstrap -j8 --with-lhapdf=/opt/homebrew/opt/lhapdf/ --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ --with-boost=/opt/homebrew/Cellar/boost/1.89.0/ --with-gsl=/opt/homebrew/Cellar/gsl/2.8/ --with-thepeg=/opt/homebrew/Cellar/thepeg/2.3.0_1/ --with-hepmc=/opt/homebrew/Cellar/hepmc3/3.2.7_1/ --with-fastjet /opt/homebrew/Cellar/fastjet/3.4.3/ .
````

The above erros when installing `gosam`, which we fixed by manual compilation of `gosam`:

```bash
cd src/gosam-contrib-2.0 
make clean
make -j8
cd ../../
./herwig-bootstrap -j8 --with-lhapdf=/opt/homebrew/opt/lhapdf/ --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ --with-boost=/opt/homebrew/Cellar/boost/1.89.0/ --with-gsl=/opt/homebrew/Cellar/gsl/2.8/ --with-thepeg=/opt/homebrew/Cellar/thepeg/2.3.0_1/ --with-hepmc=/opt/homebrew/Cellar/hepmc3/3.2.7_1/ --with-fastjet /opt/homebrew/Cellar/fastjet/3.4.3/ .
```

which worked finally:

```bash
################ /  / ^^/ ##########################
############### /--/   / ###########################
############## /  /   / ############################

Herwig 7 bootstrap was successful.

$ source bin/activate

    activates all required environment variables.

$ deactivate

    returns to the original environment variables.
```

### The errors in TL;DR

- `qgraf` without manually creating `bin` directory:
    ```bash
    Extracted qgraf-3.4.2 exists already
    /Users/user/phd/herwig_gh/src
    /Users/user/phd/herwig_gh/src/qgraf-3.4.2
    /opt/homebrew/bin/gfortran qgraf-3.4.2.f -o qgraf -O2
    Traceback (most recent call last):
    File "/Users/user/.local/share/uv/python/cpython-3.9.24-macos-aarch64-none/lib/python3.9/shutil.py", line 825, in move
        os.rename(src, real_dst)
    FileNotFoundError: [Errno 2] No such file or directory: 'qgraf' -> '/Users/user/phd/herwig_gh/./bin/qgraf'
    ```

- `dylib` error was encountered when using the `brew` version (latest at time, see `versions` section) of `thepeg`:
    ```
     /opt/homebrew/bin/gmkdir -p '/Users/user/phd/herwig_gh/./share/Herwig'
    /usr/bin/install -c -m 644 Makefile-UserModules '/Users/user/phd/herwig_gh/./share/Herwig'
    Creating repository
    dyld[71554]: Library not loaded: /opt/homebrew/opt/gsl/lib/libgsl.27.dylib
    Referenced from: <CDEF1F71-8686-3192-B087-C339CEA19970> /opt/homebrew/Cellar/thepeg/2.3.0_1/lib/ThePEG/libThePEG.30.dylib
    Reason: tried: '/Users/user/phd/herwig_gh/src/Herwig-7.3.0/API/.libs/libgsl.27.dylib' (no such file), '/opt/homebrew/opt/gsl/lib/libgsl.27.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/homebrew/opt/gsl/lib/libgsl.27.dylib' (no such file), '/opt/homebrew/opt/gsl/lib/libgsl.27.dylib' (no such file), '/Users/user/phd/herwig_gh/src/Herwig-7.3.0/API/.libs/libgsl.27.dylib' (no such file), '/opt/homebrew/Cellar/gsl/2.8/lib/libgsl.27.dylib' (no such file), '/System/Volumes/Preboot/Cryptexes/OS/opt/homebrew/Cellar/gsl/2.8/lib/libgsl.27.dylib' (no such file), '/opt/homebrew/Cellar/gsl/2.8/lib/libgsl.27.dylib' (no such file)
    make[5]: *** [install-data-hook] Abort trap: 6
    make[4]: *** [install-data-am] Error 2
    make[3]: *** [install-am] Error 2
    make[2]: *** [install-recursive] Error 1
    make[1]: *** [install] Error 2
    make: *** [install-recursive] Error 1
    ```

    which occurs because `thepeg` looks for an older version of `dylib`  in the `gsl` installtion than the one coming from the `brew` version of `gsl` (again, latest at the time).

    which was magically fixed by a symoblic link to fool the compiler:
    ```
    ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib /opt/homebrew/opt/gsl/lib/libgsl.27.dylib 
    ```

# Versions
This section documents all the versions of tools we installed with `brew`:

