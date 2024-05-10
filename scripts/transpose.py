#!/usr/bin/env python

import sys
import argparse
import pandas as pd

usage_msg = """
This script transposes a CSV or TSV file.

Positional Arguments:
file            The path to the file or STDIN if omitted.

Keyword Arguments:
--in_format, -i Specify the input file format ('csv' or 'tsv').
--out_format, -o Specify the output file format ('csv' or 'tsv').
--header, -d    Number of header rows (default 1).
--keep, -k      Keep header rows as columns in the output (not yet implimented).
--index, -n     Add an incrementing index to the output (not yet implimented).
--pad, -p       Pad the index to a fixed width with zeroes (not yet implimented).

Synopsis:
transpose.py --in_format csv --out_format tsv --header 1 --keep --index --pad 4 file.csv
"""

parser = argparse.ArgumentParser(description='Transpose a CSV/TSV file', usage=usage_msg, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file', nargs='?', type=str, default=None, help='Input file or STDIN if not specified')
parser.add_argument('--in_format', '-i', choices=['csv', 'tsv'], required=False, help='Input file format')
parser.add_argument('--out_format', '-o', choices=['csv', 'tsv'], required=False, help='Output file format')
parser.add_argument('--header', '-d', type=int, default=1, help='Number of header rows')
parser.add_argument('--keep', '-k', action='store_true', help='Keep header row in output (not yet implemented)')
parser.add_argument('--index', '-n', action='store_true', help='Add an index column to output (not yet implemented)')
parser.add_argument('--pad', '-p', type=int, help='Padding length for index column (not yet implemented)')
args = parser.parse_args()

# Guess input format from file extension unless specified with --in_format

if args.file:
    if args.in_format is None:
        if args.file.endswith('.csv'):
            args.in_format = 'csv'
        elif args.file.endswith('.tsv'):
            args.in_format = 'tsv'
        else:
            sys.exit('Error: Could not determine input format from file extension')
else:
    if args.in_format is None:
        sys.exit('FATAL : unknown_file_format : Must specify --in_format when reading from STDIN')

# Handle input from file or STDIN
if args.file:
    fh = open(args.file, 'r')
else:
    fh = sys.stdin

# Determine delimiter based on format
delimiter = ',' if args.in_format == 'csv' else '\t'

# Load the DataFrame
header_rows = None if args.header == 0 else args.header-1
df = pd.read_csv(fh, delimiter=delimiter, header=header_rows)

# Transpose the DataFrame
df = df.T

# Handle --keep which should keep the header rows as columns when transposing
if args.keep:
    df.reset_index(drop=False, inplace=True)
    df.rename(columns={df.columns[0]: 'header'}, inplace=True)
    df.set_index('header', inplace=True)

# Handle --keep
if args.keep:
    exit('FATAL : keep_not_yet_implemented')

# Handle --index
if args.index:
    exit('FATAL : index_not_yet_implemented')

# Handle --pad
if args.pad:
    exit('FATAL : pad_not_yet_implemented')

# Output formatting
out_delimiter = ',' if args.out_format == 'csv' else '\t'
df.to_csv(sys.stdout, index=False, header=bool(args.keep), sep=out_delimiter)

# Close the file handle if one was opened
if fh is not sys.stdin:
    fh.close()
