##################################################################
#
# PhyML-M3l Performance Visualizer
#
# Victor Hanson-Smith
# victorhs@cs.uoregon.edu
# August 2009 
#
# DESCRIPTION:
# This script parses the output from PhyML-M3l, captures lines with performance data,
# and then plots the data as PDF files. If the PhyML execution was configured with
# low verbosity, then the output will not contain performance data.  In this case,
# this script will exit gracefully and not create PDF plots.
# 
#
# INPUT:
# A text file with output from PhyML-M3l.  This file should contain lines with
# the following format:
#
# plot1 <x> <y>
#
# . . .where <x> is an iteration number and <y> is the -lnL of the current tree at that iteration.
#
# And:
#
# plot2 <x> <y>
#
# . . .where <x> is the clock time in milliseconds and <y> is the number of iterations thus completed.
#
# OUTPUT:
# Two PDFs will be written to the current working directory:
# (1) plot1.pdf shows iterations versus -lnL
# (2) plot2.pdf shows milliseconds versus iterations
#
#####################################################################