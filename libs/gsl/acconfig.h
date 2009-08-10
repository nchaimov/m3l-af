@BOTTOM@
/* Define if you have inline */
#undef HAVE_INLINE

/* Define if you need to hide the static definitions of inline functions */
#undef HIDE_INLINE_STATIC

/* Define if you have the ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_CLOCKS_PER_SEC

/* Defined if configure has guessed a missing ansi CLOCKS_PER_SEC clock rate */
#undef HAVE_GUESSED_CLOCKS_PER_SEC

/* Use configure's best guess for CLOCKS_PER_SEC if it is unknown */
#ifndef HAVE_CLOCKS_PER_SEC
#define CLOCKS_PER_SEC HAVE_GUESSED_CLOCKS_PER_SEC
#endif

/* Defined if you have ansi EXIT_SUCCESS and EXIT_FAILURE in stdlib.h */
#undef HAVE_EXIT_SUCCESS_AND_FAILURE

/* Use 0 and 1 for EXIT_SUCCESS and EXIT_FAILURE if we don't have them */
#ifndef HAVE_EXIT_SUCCESS_AND_FAILURE
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

/* Define one of these if you have a known IEEE arithmetic interface */
#undef HAVE_SPARCLINUX_IEEE_INTERFACE
#undef HAVE_M68KLINUX_IEEE_INTERFACE
#undef HAVE_PPCLINUX_IEEE_INTERFACE
#undef HAVE_X86LINUX_IEEE_INTERFACE
#undef HAVE_SUNOS4_IEEE_INTERFACE
#undef HAVE_SOLARIS_IEEE_INTERFACE
#undef HAVE_HPUX11_IEEE_INTERFACE
#undef HAVE_HPUX_IEEE_INTERFACE
#undef HAVE_TRU64_IEEE_INTERFACE
#undef HAVE_IRIX_IEEE_INTERFACE
#undef HAVE_AIX_IEEE_INTERFACE
#undef HAVE_FREEBSD_IEEE_INTERFACE
#undef HAVE_OS2EMX_IEEE_INTERFACE
#undef HAVE_NETBSD_IEEE_INTERFACE
#undef HAVE_OPENBSD_IEEE_INTERFACE
#undef HAVE_DARWIN_IEEE_INTERFACE

/* Define this if we need to include /usr/include/float.h explicitly
   in order to get FP_RND_RN and related macros.  This is known to be
   a problem on some Compaq Tru64 unix systems when compiled with GCC. */
#undef FIND_FP_RND_IN_USR_INCLUDE_FLOAT_H

/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#undef HAVE_EXTENDED_PRECISION_REGISTERS

#ifdef HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

/* Define this if printf can handle %Lf for long double */
#undef HAVE_PRINTF_LONGDOUBLE

/* Define this is IEEE comparisons work correctly (e.g. NaN != NaN) */
#undef HAVE_IEEE_COMPARISONS

/* Define this is IEEE denormalized numbers are available */
#undef HAVE_IEEE_DENORMALS

/* Substitute gsl functions for missing system functions */

#ifndef HAVE_HYPOT
#define hypot gsl_hypot
#endif

#ifndef HAVE_LOG1P
#define log1p gsl_log1p
#endif

#ifndef HAVE_EXPM1
#define expm1 gsl_expm1
#endif

#ifndef HAVE_ACOSH
#define acosh gsl_acosh
#endif

#ifndef HAVE_ASINH
#define asinh gsl_asinh
#endif

#ifndef HAVE_ATANH
#define atanh gsl_atanh
#endif

#ifndef HAVE_ISINF
#define isinf gsl_isinf
#endif

#ifndef HAVE_ISNAN
#define isnan gsl_isnan
#endif

#ifndef HAVE_FINITE
#ifdef HAVE_ISFINITE
#define finite isfinite
#else
#define finite gsl_finite
#endif
#endif

#ifdef __GNUC__
#define DISCARD_POINTER(p) do { ; } while(p ? 0 : 0);
#else
#define DISCARD_POINTER(p) /* ignoring discarded pointer */
#endif

#ifndef RANGE_CHECK_ON
#define RANGE_CHECK_OFF  /* turn off range checking by default */
#endif
