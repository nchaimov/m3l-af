noinst_LTLIBRARIES = libgslrandist.la

pkginclude_HEADERS= gsl_randist.h

INCLUDES= -I$(top_builddir)

libgslrandist_la_SOURCES = bernoulli.c beta.c bigauss.c binomial.c cauchy.c chisq.c discrete.c erlang.c exponential.c exppow.c fdist.c flat.c gamma.c gauss.c gausstail.c geometric.c gumbel.c hyperg.c laplace.c levy.c logarithmic.c logistic.c lognormal.c nbinomial.c pareto.c pascal.c poisson.c rayleigh.c shuffle.c sphere.c tdist.c weibull.c landau.c

TESTS = test

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la


