noinst_LTLIBRARIES = libgslcombination.la 

pkginclude_HEADERS = gsl_combination.h

INCLUDES= -I$(top_builddir) -I$(top_srcdir)

libgslcombination_la_SOURCES = init.c file.c combination.c 

noinst_HEADERS = 

TESTS = test

check_PROGRAMS = test

test_SOURCES = test.c

test_LDADD = libgslcombination.la ../vector/libgslvector.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

#noinst_PROGRAMS = demo
#demo_SOURCES = demo.c
#demo_LDADD = libgslcombination.la ../vector/libgslvector.la ../err/libgslerr.la ../test/libgsltest.la ../utils/libutils.la

#CLEANFILES = test.txt test.dat


