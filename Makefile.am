#AM_CFLAGS=-O3 -fomit-frame-pointer -funroll-loops -ftree-vectorize -Wall -fopenmp
#AM_CFLAGS=-O3 -funroll-loops -Wall -g -pg -fopenmp
AM_CFLAGS=-O3 -fomit-frame-pointer -funroll-loops -Wall -ftree-vectorize -ffast-math -fopenmp

AM_LDFLAGS=-lm -lgsl


# Uncomment (i.e. remove the `#' character at the begining of) 
# the two lines below if you want to use MPI.
# Comment the two lines below if you don't want to use MPI.

# CC=mpicc
# DEFS=-DUNIX -D$(PROG) -DDEBUG -DMPI
#SCFLAGS=-finstrument-functions -pg

# Comment the line below if you want to use MPI.
# Uncomment the line below if you don't want to use MPI.

DEFS=-DUNIX -D$(PROG) -DDEBUG

bin_PROGRAMS = phyml
EXTRA_PROGRAMS = m4 mg mc
PROG = PHYML


phyml_SOURCES = main.c utilities.c optimiz.c annealing.c emig.c lk.c bionj.c models.c free.c\
options.c simu.c eigen.c pars.c alrt.c interface.c cl.c spr.c  draw.c mpi_boot.c numeric.c
phyml_LDADD = main.o utilities.o optimiz.o annealing.o lk.o emig.o bionj.o models.o free.o\
options.o simu.o eigen.o pars.o alrt.o interface.o cl.o spr.o draw.o mpi_boot.o numeric.o


SCOMPILE = $(CC) $(SDEFS) $(DEFS) $(INCLUDES) $(AM_CFLAGS) $(SCFLAGS)
SCCLD = $(CC) $(LIBS)
SLINK = $(SCCLD) $(AM_CFLAGS) $(SLDFLAGS) $(LDFLAGS) $(INCLUDES) $(AM_LDFLAGS)


intro:	
	@echo ""
	@echo ""
	@echo ":: Building PhyML version $(VERSION) ::"
	@echo ""
	@echo ""

putversion:
	sed -e 's/#define VERSION.*/#define VERSION \"$(VERSION)\"/' utilities.h > new_utilities.h
	mv new_utilities.h utilities.h

phyml$(EXEEXT): intro putversion $(phyml_LDADD) $(phyml_SOURCES)
	$(SLINK) $(phyml_LDADD) -o $@
	@echo ""
	@echo "Done."

main.o: main.c
	$(SCOMPILE) -c main.c

utilities.o: utilities.c utilities.h
	$(SCOMPILE) -c utilities.c

optimiz.o: optimiz.c optimiz.h
	$(SCOMPILE) -c optimiz.c

annealing.o: annealing.c annealing.h
	$(SCOMPILE) -c annealing.c

emig.o: emig.c emig.h
	$(SCOMPILE) -c emig.c

lk.o: lk.c lk.h
	$(SCOMPILE) -c lk.c

bionj.o: bionj.c bionj.h
	$(SCOMPILE) -c bionj.c

models.o: models.c models.h
	$(SCOMPILE) -c models.c

free.o: free.c free.h
	$(SCOMPILE) -c free.c

options.o: options.c options.h
	$(SCOMPILE) -c options.c

simu.o: simu.c simu.h
	$(SCOMPILE) -c simu.c

eigen.o: eigen.c eigen.h
	$(SCOMPILE) -c eigen.c

pars.o: pars.c pars.h
	$(SCOMPILE) -c pars.c

alrt.o: alrt.c alrt.h
	$(SCOMPILE) -c alrt.c

interface.o: interface.c interface.h
	$(SCOMPILE) -c interface.c

cl.o: cl.c cl.h
	$(SCOMPILE) -c cl.c

spr.o: spr.c spr.h
	$(SCOMPILE) -c spr.c

# m4.o: m4.c m4.h
# 	$(SCOMPILE) -c m4.c

mg.o: mg.c mg.h
	$(SCOMPILE) -c mg.c

draw.o: draw.c draw.h
	$(SCOMPILE) -c draw.c

# rates.o: rates.c rates.h
# 	$(SCOMPILE) -c rates.c

mpi_boot.o: mpi_boot.c mpi_boot.h
	$(SCOMPILE) -c mpi_boot.c

numeric.o: numeric.c numeric.h
	$(SCOMPILE) -c numeric.c
