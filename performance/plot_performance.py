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
# USAGE:
# python plot_performance.py <filepath> 
#
# ... where <filepath> is the path to a text file capturing output from PhyML-M3l.  
# This file should contain some lines with the following format:
# 
# plot1 <x> <y>
#
# ...where <x> is an iteration number and <y> is the -lnL of the current tree at that iteration.
#
# And:
#
# plot2 <x> <y>
#
# ...where <x> is the clock time in milliseconds and <y> is the number of iterations thus completed.
#
# OUTPUT:
# Two CRAN scripts (i.e. 'R' scripts) will be written to the current working directory.  
# When these scripts are invoked, they will generate PDF plots.
# (1) plot1.cran plots iterations vs. -lnL
# (2) plot2.cran plots milliseconds vs. -lnL
#
#####################################################################

import os
import sys
import re

try:
    filename = sys.argv[1]
except IndexError:
    print "Yikes! I think you forgot to specify the filepath containing PhyML output."
    print "\nUsage: python plot_performance.py <filepath>"

try:
    fin = open(filename, "r")
except IOError:
    print "Yikes! The file named ", filename, " does not exist."
    exit()

plot1_xy = {} # where key = iterations as integers, value = -lnL of the current tree at that iteration
plot2_xy = {} # where key = milliseconds as doubles, value = current iteration at that time

#
# Parse the PhyML output
#
for line in fin.readlines():
    if line.startswith("plot1"):
        line = line.strip()
        tokens = line.split()
        plot1_xy[ float(tokens[1]) ] = float(tokens[2])
    elif line.startswith("plot2"):
        line = line.strip()
        tokens = line.split()
        plot2_xy[ float(tokens[1]) ] = float(tokens[2])

fin.close()

#
# Sanity check
#
early_quit = False
if plot1_xy.__len__() < 1:
    print "I did not find any data for plot1 (iterations vs. -lnL)"
    early_quit = True
if plot2_xy.__len__() < 1:
    print "I did not find any data for plot2 (time vs. iterations)"
    early_quit = True
if early_quit:
    exit(1)

#
# Write CRAN scripts and plot PDFs.
#
def plot_in_r(points, output_filename_seed, title, xlab, ylab):
    x_sorted = points.keys()
    x_sorted.sort()
    
    string = "x<-c("
    for x in x_sorted:
        string += x.__str__() + ","   
    string = re.sub(",$", "", string)
    string += ");\n"
    
    string += "y<-c("
    for x in x_sorted:
        string += points[x].__str__() + ","
    string = re.sub(",$", "", string)
    string += ");\n"
    string += "plot(x, y, type='l', main='" + title + "', xlab='" + xlab + "', ylab='" + ylab + "');\n"
    
    fout_cran = open(output_filename_seed + ".cran", "w")
    fout_cran.write("pdf('" + output_filename_seed + ".pdf', width=10, height=5);\n")
    fout_cran.write(string)
    fout_cran.write("dev.off();\n")
    fout_cran.close()
    
    os.system("r --no-save < " + output_filename_seed + ".cran")

plot_in_r(plot1_xy, "plot1", "-lnL(iterations)", "iterations", "-lnL")
plot_in_r(plot2_xy, "plot2", "iterations(time)", "seconds", "iterations")