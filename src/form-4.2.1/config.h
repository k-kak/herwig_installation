/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <boost/unordered_map.hpp> header file. */
/* #undef HAVE_BOOST_UNORDERED_MAP_HPP */

/* Define to 1 if you have the <boost/unordered_set.hpp> header file. */
/* #undef HAVE_BOOST_UNORDERED_SET_HPP */

/* Define to 1 if you have __builtin_popcount function. */
#define HAVE_BUILTIN_POPCOUNT 1

/* Define to 1 if you have clock_gettime. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have ftime. */
/* #undef HAVE_FTIME */

/* Define to 1 if you have gettimeofday. */
/* #undef HAVE_GETTIMEOFDAY */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <limits.h> header file. */
#define HAVE_LIMITS_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have __popcnt function. */
/* #undef HAVE_POPCNT */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/file.h> header file. */
#define HAVE_SYS_FILE_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <tr1/unordered_map> header file. */
/* #undef HAVE_TR1_UNORDERED_MAP */

/* Define to 1 if you have the <tr1/unordered_set> header file. */
/* #undef HAVE_TR1_UNORDERED_SET */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <unordered_map> header file. */
#define HAVE_UNORDERED_MAP 1

/* Define to 1 if you have the <unordered_set> header file. */
#define HAVE_UNORDERED_SET 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Compiling for ILP32 data model */
/* #undef ILP32 */

/* Compiling for a Linux system. */
/* #undef LINUX */

/* Compiling for LLP64 data model */
/* #undef LLP64 */

/* Compiling for LP64 data model */
#define LP64 /**/

/* Name of package */
#define PACKAGE "form"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/vermaseren/form/issues"

/* Define to the full name of this package. */
#define PACKAGE_NAME "FORM"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "FORM 4.2.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "form"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "4.2.1"

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `off_t', as computed by sizeof. */
#define SIZEOF_OFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Compiling for UNIX system */
#define UNIX /**/

/* Version number of package */
#define VERSION "4.2.1"

/* Compiling for WINDOWS system */
/* #undef WINDOWS */

/* Define to use GMP for long integer arithmetic. */
/* #undef WITHGMP */

/* Define to use POSIX thread clock. */
#ifdef WITHPTHREADS
#define WITHPOSIXCLOCK /**/
#endif

/* Define to use zlib for compression. */
#define WITHZLIB /**/

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
