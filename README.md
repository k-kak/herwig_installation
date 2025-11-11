# Herwig7 Installation Guide for macOS (Apple Silicon)

This repository documents a successful Herwig7 bootstrap installation on macOS 15.7.2 (Apple M4), providing working solutions to common compilation issues.

**Authors:** Mohamed Aly (mohamed.aly@cern.ch), Karim Kandeel (karim.kandeel@cern.ch)

## Table of Contents

- [Quick Start (Recommended)](#quick-start-recommended)
- [Prerequisites](#prerequisites)
- [Detailed Installation Steps](#detailed-installation-steps)
  - [1. Installation Directory Setup](#1-installation-directory-setup)
  - [2. Python Environment Setup](#2-python-environment-setup)
  - [3. Installing Dependencies via Homebrew](#3-installing-dependencies-via-homebrew)
  - [4. Running the Bootstrap Script](#4-running-the-bootstrap-script)
  - [5. Manual GoSAM Compilation (if needed)](#5-manual-gosam-compilation-if-needed)
- [Common Errors and Solutions](#common-errors-and-solutions)
  - [FastJet Compilation Error](#fastjet-compilation-error)
  - [CMake Version Compatibility](#cmake-version-compatibility)
  - [LHAPDF Python Configuration](#lhapdf-python-configuration)
  - [RIVET and YODA Python Issues](#rivet-and-yoda-python-issues)
  - [Boost Linking Errors](#boost-linking-errors)
  - [ThePEG Fortran Compilation](#thepeg-fortran-compilation)
  - [qgraf Missing bin Directory](#qgraf-missing-bin-directory)
  - [GSL dylib Version Mismatch](#gsl-dylib-version-mismatch)
- [Installed Versions](#installed-versions)

---

## Quick Start (Recommended)

**This is the fastest way to get Herwig7 running.** After encountering numerous compilation issues, we found this workflow to be the most stable.

### Complete Installation Commands

```bash
# 1. Tap HEP packages repository
brew tap davidchall/hep

# 2. Install all dependencies via Homebrew
brew install thepeg lhapdf rivet yoda boost gsl fastjet

# 3. Create and setup installation directory
mkdir -p ~/herwig_installation  # or your preferred path
cd ~/herwig_installation

# 4. Download bootstrap script
wget https://herwig.hepforge.org/downloads/herwig-bootstrap
chmod +x herwig-bootstrap

# 5. Create bin directory (prevents qgraf error)
mkdir bin

# 6. Fix GSL dylib version mismatch
ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib /opt/homebrew/opt/gsl/lib/libgsl.27.dylib

# 7. Setup Python environment with uv
curl -LsSf https://astral.sh/uv/install.sh | sh  # if uv not installed
uv init
uv python install 3.9
uv venv
source .venv/bin/activate
uv add six cython lib
uv sync

# 8. Run bootstrap (adjust version numbers to match your brew installations)
./herwig-bootstrap -j8 \
  --with-lhapdf=/opt/homebrew/opt/lhapdf/ \
  --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ \
  --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ \
  --with-boost=/opt/homebrew/Cellar/boost/1.89.0/ \
  --with-gsl=/opt/homebrew/Cellar/gsl/2.8/ \
  --with-thepeg=/opt/homebrew/Cellar/thepeg/2.3.0_1/ \
  --with-hepmc=/opt/homebrew/Cellar/hepmc3/3.2.7_1/ \
  --with-fastjet=/opt/homebrew/Cellar/fastjet/3.4.3/ .
```

**Important:** Version numbers in the bootstrap command (e.g., `yoda/1.9.10`) must match your installed versions. Check with:
```bash
brew info <package>
```

### If GoSAM Fails

If the bootstrap fails during GoSAM compilation:

```bash
cd src/gosam-contrib-2.0
make clean
make -j8
cd ../../

# Re-run the bootstrap command from step 8 above
```

### Success Indicator

You should see:

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

---

## Prerequisites

Before starting the installation, ensure you have the following:

### System Requirements
- **OS:** macOS 15.7.2 or similar (tested on Apple M4)
- **Architecture:** Apple Silicon (ARM64)
- **Xcode Command Line Tools:** Install with `xcode-select --install`
- **Disk Space:** ~5-10 GB for source files and compiled binaries

### Package Manager
- **Homebrew:** Install from [brew.sh](https://brew.sh) if not already installed
- **HEP tap:** Required for physics packages: `brew tap davidchall/hep`

### Essential Tools
- **wget:** `brew install wget` (for downloading bootstrap script)
- **tree (optional):** `brew install tree` (helpful for visualizing directory structure)

### Python Environment
- **uv:** Fast Python package/environment manager
  - Install: `curl -LsSf https://astral.sh/uv/install.sh | sh`
  - Why uv? Provides more reliable environment management than conda/virtualenv for this workflow
- **Python 3.9:** Managed via uv (tested with 3.9.6 and 3.9.24)
  - Why 3.9? Chosen based on compatibility with other HEP event generators

### Why Install Dependencies via Homebrew?

Building dependencies from source encounters numerous compilation errors on Apple Silicon:
- **Python configuration issues:** LHAPDF, RIVET, and YODA fail to detect virtual environments
- **Fortran compatibility:** ThePEG requires specific compiler flags
- **Linking errors:** Boost and other libraries have complex dependency chains

Installing pre-built binaries via Homebrew avoids these issues entirely.

---

## Detailed Installation Steps

This section provides a step-by-step breakdown with explanations.

### 1. Installation Directory Setup

Create a dedicated directory for your Herwig installation:

```bash
mkdir -p ~/herwig_installation  # Example path
cd ~/herwig_installation
```

**Note:** Throughout this guide, we refer to this as `<my_herwig>`. Replace with your actual path.

**Optional - visualize directory structure:**
```bash
brew install tree  # if not already installed
tree ./
```

### 2. Python Environment Setup

We use `uv` to manage a Python 3.9 virtual environment.

**Install uv (if not already installed):**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Create virtual environment:**
```bash
cd <my_herwig>  # Your installation directory

# Initialize uv project
uv init

# Install Python 3.9
uv python install 3.9

# Create virtual environment
uv venv

# Activate the environment
source .venv/bin/activate
```

**Verify Python setup:**
```bash
which python && python --version
# Expected output:
# <my_herwig>/.venv/bin/python
# Python 3.9.24 (or 3.9.6)
```

**Install required Python packages:**
```bash
uv add six cython lib
uv sync
```

**Expected directory structure:**
```bash
.venv/
├── bin
└── lib
    └── python3.9
        └── site-packages
            ├── Cython/
            ├── lib/
            ├── pyximport/
            └── six-*.dist-info/
```

### 3. Installing Dependencies via Homebrew

Install all required HEP dependencies using Homebrew:

```bash
# Tap the HEP repository (if not already done)
brew tap davidchall/hep

# Install all dependencies
brew install thepeg lhapdf rivet yoda boost gsl fastjet
```

**Note:** HepMC3 will be installed automatically as a dependency.

**Verify installation:**
```bash
brew list --versions thepeg yoda rivet gsl boost fastjet lhapdf
```

### 4. Running the Bootstrap Script

**Download and prepare the bootstrap script:**
```bash
cd <my_herwig>
wget https://herwig.hepforge.org/downloads/herwig-bootstrap
chmod +x herwig-bootstrap
```

**Prepare environment (critical steps):**

1. **Create bin directory** (prevents qgraf installation error):
```bash
mkdir bin
```

2. **Fix GSL dylib compatibility** (ThePEG expects libgsl.27, but Homebrew provides libgsl.28):
```bash
ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib \
      /opt/homebrew/opt/gsl/lib/libgsl.27.dylib
```

**Run the bootstrap script:**
```bash
./herwig-bootstrap -j8 \
  --with-lhapdf=/opt/homebrew/opt/lhapdf/ \
  --with-yoda=/opt/homebrew/Cellar/yoda/1.9.10/ \
  --with-rivet=/opt/homebrew/Cellar/rivet/3.1.10/ \
  --with-boost=/opt/homebrew/Cellar/boost/1.89.0/ \
  --with-gsl=/opt/homebrew/Cellar/gsl/2.8/ \
  --with-thepeg=/opt/homebrew/Cellar/thepeg/2.3.0_1/ \
  --with-hepmc=/opt/homebrew/Cellar/hepmc3/3.2.7_1/ \
  --with-fastjet=/opt/homebrew/Cellar/fastjet/3.4.3/ .
```

**Important notes:**
- The `-j8` flag uses 8 parallel jobs (adjust based on your CPU cores)
- Version numbers in paths must match your installations (see [Installed Versions](#installed-versions))
- Check installed versions with: `brew info <package_name>`
- The final `.` specifies the current directory as installation target

**Finding correct paths:**
```bash
# For packages in /opt/homebrew/opt/ (symlinks):
brew --prefix lhapdf

# For packages in Cellar (actual installations):
brew info yoda | grep Cellar
```

### 5. Manual GoSAM Compilation (if needed)

If bootstrap fails during GoSAM compilation:

```bash
cd src/gosam-contrib-2.0
make clean
make -j8
cd ../../
```

Then re-run the bootstrap command from step 4.

**Why this happens:** The GoSAM build system sometimes fails with parallel compilation. Manual compilation resolves race conditions.

---

## Common Errors and Solutions

This section documents errors encountered when building dependencies from source. **If you follow the Quick Start, you should not encounter these errors.**

### FastJet Compilation Error

**Error:**
```
./ProtoJet.hpp:198:51: error: no member named '_Et' in 'ProtoJet<Item>'
198 |   os<<"y phi Et = ("<<_y<<", "<<_phi<<", "<<this->_Et<<")"<<std::endl;
```

**Location:** `src/fastjet-3.4.0/plugins/D0RunIICone/ProtoJet.hpp:198`

**Workaround:** Comment out the problematic print statement (harmless, only affects debug output).

**Recommended solution:** Install FastJet via Homebrew: `brew install fastjet`

---

### CMake Version Compatibility

**Error:**
```
CMake Error at CMakeLists.txt:1 (cmake_minimum_required):
Compatibility with CMake < 3.5 has been removed from CMake.
```

**Cause:** Homebrew's CMake (4.1.2+) removed support for old version requirements in HepMC3.

**Workaround:**
```bash
export CMAKE_POLICY_VERSION_MINIMUM=3.5
```

**Recommended solution:** Use Homebrew HepMC3: installed automatically with other dependencies.

---

### LHAPDF Python Configuration

**Error:**
```
configure: error:
Could not link test program to Python. Maybe the main Python library has been
installed in some non-standard library path.
```

**Full error details:**
```
checking for Python include path... -I/Library/Developer/CommandLineTools/...
checking python extra libraries... -ldl -lSystem  -framework CoreFoundation
checking consistency of all components of python development environment... no
```

**Cause:** LHAPDF's configure script fails to detect the uv virtual environment's Python library.

**Recommended solution:**
```bash
brew install lhapdf
./herwig-bootstrap --with-lhapdf=/opt/homebrew/opt/lhapdf/ .
```

---

### RIVET and YODA Python Issues

**Error:** Same Python configuration issues as LHAPDF.

**Recommended solution:**
```bash
brew tap davidchall/hep
brew install davidchall/hep/yoda  # explicit tap specification needed
brew install --formula rivet      # --formula flag needed
```

Then use `--with-yoda` and `--with-rivet` flags in bootstrap command.

---

### Boost Linking Errors

**Error:** Boost library linking failures during compilation.

**Recommended solution:**
```bash
brew install boost gsl
```

Then use `--with-boost` and `--with-gsl` flags in bootstrap command.

---

### ThePEG Fortran Compilation

**Error:**
```fortran
D/D0func.F:486:35:
486 |         parameter (nz2 = -2147483648)   ! O'20000000000'
    |                                            1
Error: Integer too big for its kind at (1). This check can be disabled with
the option '-fno-range-check'
```

**Full error context:**
```fortran
D/D0func.F:491:27:
491 |         parameter (nz2p1234 = nz2 + p1234)
    |                                  1
Error: Parameter 'nz2' at (1) has not been declared or is a variable, which
does not reduce to a constant expression
```

**Cause:** Modern gfortran compilers on Apple Silicon are stricter about integer sizes in legacy Fortran code.

**Workaround:** Modify `herwig-bootstrap` script's `compile()` function.

Find this line:
```python
check_call(["make","-s","-j%s" % opts.ncore])
```

Replace with:
```python
check_call(["make", "-s", "-j%s" % opts.ncore,
           "FFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check",
           "FCFLAGS=-g -std=legacy -ffixed-line-length-none -fno-range-check"])
```

Then clean and rebuild:
```bash
cd src/Herwig-7.3.0/
make clean
cd ../../
./herwig-bootstrap ...  # re-run with same flags
```

**Recommended solution:** Use Homebrew ThePEG: `brew install thepeg`

---

### qgraf Missing bin Directory

**Error:**
```
Traceback (most recent call last):
  File ".../shutil.py", line 825, in move
    os.rename(src, real_dst)
FileNotFoundError: [Errno 2] No such file or directory: 'qgraf' -> '.../bin/qgraf'
```

**Cause:** Bootstrap script tries to move compiled qgraf binary to `bin/` before the directory exists.

**Solution:** Create the directory before running bootstrap:
```bash
mkdir bin
```

---

### GSL dylib Version Mismatch

**Error:**
```
dyld[71554]: Library not loaded: /opt/homebrew/opt/gsl/lib/libgsl.27.dylib
  Referenced from: <...> /opt/homebrew/Cellar/thepeg/2.3.0_1/lib/ThePEG/libThePEG.30.dylib
  Reason: tried: '/opt/homebrew/opt/gsl/lib/libgsl.27.dylib' (no such file)
```

**Full error context:**
```
make[5]: *** [install-data-hook] Abort trap: 6
make[4]: *** [install-data-am] Error 2
make[3]: *** [install-am] Error 2
make[2]: *** [install-recursive] Error 1
make[1]: *** [install] Error 2
make: *** [install-recursive] Error 1
```

**Cause:** ThePEG (from Homebrew) was compiled against libgsl.27.dylib, but current GSL provides libgsl.28.dylib.

**Solution:** Create a symbolic link:
```bash
ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib \
      /opt/homebrew/opt/gsl/lib/libgsl.27.dylib
```

---

## Installed Versions

This section documents the Homebrew package versions used in our successful installation:

| Package  | Version    | Installation Notes                                    |
|----------|------------|-------------------------------------------------------|
| boost    | 1.89.0     | Standard Homebrew package                             |
| fastjet  | 3.4.3      | Standard Homebrew package                             |
| gsl      | 2.8        | Requires symlink (see below)                          |
| rivet    | 3.1.10     | Install with `brew install --formula rivet`           |
| thepeg   | 2.3.0_1    | Via `brew tap davidchall/hep`, requires GSL symlink   |
| yoda     | 1.9.10     | Install with `brew install davidchall/hep/yoda`       |
| lhapdf   | 6.5.4      | Via `brew install lhapdf`                             |
| hepmc3   | 3.2.7_1    | Automatically installed as dependency                 |
| Python   | 3.9.24     | Via uv (3.9.6 also tested and works)                  |

**Check your installed versions:**
```bash
brew list --versions thepeg yoda rivet gsl boost fastjet lhapdf
```

**Important:** The version numbers in your bootstrap command paths (e.g., `/opt/homebrew/Cellar/yoda/1.9.10/`) must match your actual installed versions. Use `brew info <package>` to verify paths.

### GSL Symlink Requirement

ThePEG 2.3.0_1 was compiled against libgsl.27.dylib, but GSL 2.8 provides libgsl.28.dylib. Create a symlink to resolve this:

```bash
ln -s /opt/homebrew/opt/gsl/lib/libgsl.28.dylib \
      /opt/homebrew/opt/gsl/lib/libgsl.27.dylib
```

---

## Additional Notes

### Understanding Bootstrap Flags

- `-j8`: Use 8 parallel compilation jobs (adjust based on CPU cores)
- `--with-<package>=<path>`: Use pre-installed package instead of building from source
- Final `.`: Install to current directory

### Environment Activation

After successful installation, activate the Herwig environment:

```bash
source bin/activate
```

To deactivate:

```bash
deactivate
```

### Troubleshooting

If you encounter errors not listed above:

1. Check Homebrew package versions match the bootstrap command paths
2. Ensure Python virtual environment is activated
3. Verify `bin/` directory exists
4. Confirm GSL symlink is in place
5. Review the detailed error messages in the terminal

### Contributing

If you encounter additional errors or find better solutions, please contact the authors or contribute to this documentation.

---

## References

- **Herwig Official Documentation:** https://herwig.hepforge.org
- **Bootstrap Installation Guide:** https://herwig.hepforge.org/tutorials/installation/
- **Homebrew:** https://brew.sh
- **uv Python Manager:** https://astral.sh/uv
