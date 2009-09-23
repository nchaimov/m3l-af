#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "options.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef COMPRESS_SUBALIGNMENTS
#include "compress.h"
#endif

//#define RUN_TESTS // This toggles on/off the entire unit test suite.

void Write_File(char *filename, char *s, int append);
void Run_tests(int argc, char **argv);
int Test_main(int argc, char **argv, int testid);

void Test_alignment_read(option *io, seq **data, int testid);
