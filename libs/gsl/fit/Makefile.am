noinst_LTLIBRARIES = libgslfit.la 

pkginclude_HEADERS = gsl_fit.h

INCLUDES= -I$(top_builddir)

libgslfit_la_SOURCES = linear.c

check_PROGRAMS = test #demo

TESTS = test

test_SOURCES = test.c
test_LDADD = libgslfit.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

#demo_SOURCES = demo.c
#demo_LDADD = libgslfit.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la ../sys/libgslsys.la

