#!/usr/bin/env python

import sys
import os
import re
from getopt import getopt, GetoptError

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

usage = """
This script will take one or more datafiles and build a histogram of
one or more columns of data from each file.  The options described
below allow control of the binning, input and output.  Items sort into
bins if they are greater then the bin's lower bound and less than or
equal to the bin's upper bound.

histogram [options] datafile1 [datafile2 datafile3...]

Options:

    --col n
      Which column of the input data should we use?  You can also pass
      a list of column numbers and all of those columns will be porcessed
      in the same manner i.e. (1 3 7 10) - output will include histograms
      of columns 1, 3, 7 and 10 all with the same binning.

    --bin_auto [n]
      Automatically makes n bins divided evenly over the data range.  The
      value of n defaults to 10 and if no bin options are given, then this
      option is the default.

    --bin_range 'min max step'
      Builds bins begining at min, continuing to max by step.

    --bin_range_log 'min max step base'
      Builds bins begining at min, continuing to max by step where the base
      is raised to the power of step.

    --bin_file filename
      Build bins based on the data given in a file.  Each value in the first
      column of the file becomes a bin.

   --find filename
     Use the linux command 'find' to search for all files of a given name below
     the current working directory and then add these files to any passed on the
     command line.

  --accuracy
     Determines the accuracy at which to compare floating point numbers.  Integers
     are converted to floating point.  Default is 6 digits after the decimal place.

  --graph
     Print a ASCII graph of the histogram.
"""

bin_auto = 10

bin_range = None
bin_range_log = None
bin_file = None
col_string = None
find = None
accuracy = 6
graph = False

try:
    opts, args = getopt(sys.argv[1:], '', ['bin_auto=', 'bin_range=', 'bin_range_log=', 'bin_file=', 'col=', 'find=', 'accuracy=', 'graph'])
except GetoptError:
    print(usage)
    sys.exit(2)

for opt, arg in opts:
    if opt == '--bin_auto':
        bin_auto = int(arg)
    elif opt == '--bin_range':
        bin_range = arg
    elif opt == '--bin_range_log':
        bin_range_log = arg
    elif opt == '--bin_file':
        bin_file = arg
    elif opt == '--col':
        col_string = arg
    elif opt == '--find':
        find = arg
    elif opt == '--accuracy':
        accuracy = int(arg)
    elif opt == '--graph':
        graph = True

cols = []
if col_string:
    cols = [int(x) - 1 for x in re.split(r'[\s,;-_]+', col_string)]
if not cols:
    cols = [-1]

files = args
if find:
    find_files = os.popen(f'find ./ -name {find}').read().strip().split('\n')
    files.extend(find_files)
if not files:
    files = ['-']

def parse_data(files):
    data = []
    data_min = float('inf')
    data_max = float('-inf')
    for file in files:
        if file == '-':
            f = sys.stdin
        else:
            f = open(file, 'r')
        for line in f:
            fields = [float(x) for x in line.strip().split()]
            for col in cols:
                value = fields[col]
                data.append(value)
                data_min = min(data_min, value)
                data_max = max(data_max, value)
        if file != '-':
            f.close()
    return data, data_min, data_max

def build_bin(data_min, data_max, bin_auto, bin_range, bin_file):
    if bin_range:
        min_val, max_val, step = [float(x) for x in bin_range.split()]
        bins = [min_val + i*step for i in range(int((max_val - min_val) / step) + 1)]
    elif bin_range_log:
        min_val, max_val, step, base = [float(x) for x in bin_range_log.split()]
        bins = [base**(min_val + i*step) for i in range(int((max_val - min_val) / step) + 1)]
    elif bin_file:
        bins = [float(x) for x in open(bin_file, 'r').read().strip().split()]
    else:
        bins = [data_min + i*(data_max - data_min)/bin_auto for i in range(bin_auto + 1)]
    return bins

def build_histogram(bins, data):
    histogram = [0] * len(bins)
    for value in data:
        for i, bin_val in enumerate(bins):
            if i == 0 and value <= bin_val:
                histogram[i] += 1
                break
            elif i > 0 and bins[i-1] < value <= bin_val:
                histogram[i] += 1
                break
    return histogram

def print_graph(histogram, bin_range):
    if bin_range:
        min_val, max_val, step = [float(x) for x in bin_range.split()]
        bin_labels = [f'{min_val + i*step:.{accuracy}f}' for i in range(len(histogram))]
    else:
        bin_labels = [f'{x:.{accuracy}f}' for x in bins]

    max_count = max(histogram)
    for i, count in enumerate(histogram):
        print(f'{bin_labels[i]:>12} | {"*" * int(count/max_count*50)}')

def print_histogram(histogram):
    for i, count in enumerate(histogram):
        print(f'Bin {i+1}: {count}')

data, data_min, data_max = parse_data(files)
bins = build_bin(data_min, data_max, bin_auto, bin_range, bin_file)
histogram = build_histogram(bins, data)

if graph:
    print_graph(histogram, bin_range)
else:
    print_histogram(histogram)

sys.exit(0)
