#!/usr/bin/env python

"""
Split a FASTA file based on the occurrence of markers.


The output is FASTA-like, depending on the replacement defined in the library.
If the replacement is identical to the marker, the output is in FASTA format.

Per marker, two files are created:
- markername.txt         : all sequences that have the marker as substring
- markername_counted.txt : unique sequences (counts are placed in the header)

Format of the library file:
name marker replacement
"""

import argparse
import sys

from Bio import SeqIO

from . import version


def add_to_library(library, line):
    """
    Add a line from the library file to the library.

    A library entry contains:
    - name of the marker
    - sequence of the marker
    - replacement sequence
    - open handle to the .txt file
    - open handle to the _counted.txt file
    - number of times the marker has been seen
    - sequences containing the marker, its header and its count

    :arg list library: Library (see above).
    :arg str line: Line of the library file.
    """
    library.append(line.split())
    library[-1].append(open(library[-1][0] + '.txt', 'w'))
    library[-1].append(open(library[-1][0] + '_counted.txt', 'w'))
    library[-1].append(0)
    library[-1].append({})


def write_output(entry, description, sequence):
    """
    Write a sequence and its header to its .txt file.

    While writing, also increase the marker counter and count the number of
    times this particular sequence has been found.

    :arg list entry: Library entry.
    :arg str description: Sequence header.
    :arg str sequence:  The sequence.
    """
    entry[3].write('>{}\n'.format(description))
    entry[3].write('{}\n'.format(sequence))
    entry[5] += 1
    if entry[6].has_key(sequence):
        entry[6][sequence][1] += 1
    else:
        entry[6][sequence] = [description, 1]


def split_fasta(input_handle, library_handle, output_handle):
    """
    Split a FASTA file based on the occurrence of markers.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream library_handle: Open readable handle to the library file.
    :arg stream output_handle: Open writable handle for summary data.
    """
    # Read the library file.
    library = []
    line = library_handle.readline()
    while line:
        add_to_library(library, line)
        line = library_handle.readline()
    add_to_library(library, 'Unrecognised _INVALID_ _INVALID_')

    # Split the fasta file.
    total = 0
    recognised = 0
    for record in SeqIO.parse(input_handle, 'fasta'):
        hit = False
        for entry in library:
            if entry[1] in record.seq:
                flanks = str(record.seq).split(entry[1])
                write_output(entry, record.description,
                    flanks[0] + entry[2] + flanks[1])
                recognised += 1
                hit = True
            if entry[1] == '_INVALID_' and not hit:
                write_output(entry, record.description, str(record.seq))
        total += 1

    # Write summary information and write the _counted.txt files.
    output_handle.write('Total: {}\n'.format(total))
    for entry in library:
        output_handle.write('{}: {}\n'.format(entry[0], entry[5]))
        unique = 0
        maximum = 0
        for i in sorted(entry[6], key = lambda item: -entry[6][item][1]):
            entry[4].write('>COUNT {:04d} {}\n'.format(entry[6][i][1],
                entry[6][i][0]))
            entry[4].write('{}\n'.format(i))
            if entry[6][i][1] > maximum:
                maximum = entry[6][i][1]
            unique += 1
        output_handle.write('  unique: {}\n'.format(unique))
        output_handle.write('  max: {}\n'.format(maximum))
        entry[3].close()
        entry[4].close()
    output_handle.write('Unrecognised: {}\n'.format(total - recognised))


def main():
    """
    Main entry point.
    """
    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])

    parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='input file')
    parser.add_argument('library_handle', metavar='LIBRARY',
        type=argparse.FileType('r'), help='file containing markers')
    parser.add_argument('-o', dest='output_handle', metavar='OUTPUT',
        type=argparse.FileType('w'), default=sys.stdout,
        help='output file (default=stdout)')
    parser.add_argument('-v', action='version', version=version(parser.prog))

    args = parser.parse_args()

    split_fasta(args.input_handle, args.library_handle, args.output_handle)


if __name__ == '__main__':
    main()
