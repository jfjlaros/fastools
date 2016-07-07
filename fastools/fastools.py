#!/usr/bin/python
import argparse
import itertools
import re
import Levenshtein
import urllib2
import random
import sys

from collections import defaultdict

from Bio import Seq, SeqIO, Entrez, pairwise2, Restriction
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from . import doc_split, version, usage


def guess_file_format(handle):
    """
    Guess the file type of an NGS data file.

    :arg file handle: Open readable handle to an NGS data file.

    :return str: Either 'fasta' or 'fastq'.
    """
    try:
        extension = getattr(handle, 'name').split('.')[-1]
    except AttributeError:
        pass
    else:
        if extension in ('fastq', 'fq'):
            return 'fastq'
        elif extension in ('fasta', 'fa'):
            return 'fasta'

    try:
        position = handle.tell()
        handle.seek(0)
    except IOError:
        sys.stderr.write('Cannot deterine file type in stream, assuming FASTA')
        return 'fasta'

    token = handle.read(1)
    handle.seek(position)

    if token == '>':
        return 'fasta'
    return 'fastq'


def bed_read(handle):
    """
    Parse a four-column BED file (chrom start end replacement).

    :arg file handle: Open readable handle to a BED file.

    :returns dict: A list of edits (ranges and replacements) per chromosome.
    """
    records = defaultdict(list)

    for line in handle.readlines():
        record = line.split()
        records[record[0]].append([int(record[1]), int(record[2]), record[3]])
    for reference in records:
        records[reference].sort(reverse=True)
    return records


def sanitise(input_handle, output_handle):
    """
    Convert a FASTA/FASTQ file to a standard FASTA/FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        SeqIO.write(record, output_handle, file_format)


def fa2fq(input_handle, output_handle, quality):
    """
    Convert a FASTA file to a FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a FASTQ file.
    :arg int quality: Quality score.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        record.letter_annotations = {'phred_quality':
            [quality] * len(record.seq)}
        SeqIO.write(record, output_handle, 'fastq')


def fq2fa(input_handle, output_handle):
    """
    Convert a FASTQ file to a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    try:
        for record in SeqIO.parse(input_handle, 'fastq'):
            SeqIO.write(record, output_handle, 'fasta')
    except ValueError, error:
        print 'Error: %s' % error
        emptyRecord = SeqRecord(Seq.Seq(''), '', '', '')
        SeqIO.write(emptyRecord, output_handle, 'fasta')


def add(input_handle, output_handle, sequence, quality):
    """
    Add a sequence to the 5' end of each read in a FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTQ file.
    :arg str sequence: Sequence to be added to the 5' end of the read.
    :arg int quality: Quality score.
    """
    addition_q = [quality] * len(sequence)

    for record in SeqIO.parse(input_handle, 'fastq'):
        qual = record.letter_annotations['phred_quality']
        record.letter_annotations = {}
        record.seq = sequence + record.seq
        record.letter_annotations = {'phred_quality': addition_q + qual}
        SeqIO.write(record, output_handle, 'fastq')


def aln(input_handles):
    """
    Calculate the Levenshtein distance between two FASTA files.

    :arg stream input_handles: Open readable handles to FASTA files.

    :return list(tuple(str, str, int)): List of distances.
    """
    distances = []

    for i in SeqIO.parse(input_handles[0], 'fasta'):
        for j in SeqIO.parse(input_handles[1], 'fasta'):
            distances.append((i.name, j.name,
                Levenshtein.distance(str(i.seq), str(j.seq))))

    return distances


def maln(input_handle):
    """
    Calculate the Hamming distance between all sequences in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :return list(tuple(str, str, int)): List of distances.
    """
    distances = []

    data = {x.name: str(x.seq) for x in SeqIO.parse(input_handle, 'fasta')}
    for j in data:
        print j,
    print 
    for i in data:
        print i,
        for j in data:
            print Levenshtein.hamming(data[i], data[j]),
        print

    return distances


def length(input_handle):
    """
    Report the lengths of all FASTA records in a file.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :return list(int): List of lengths.
    """
    lengths = []

    for record in SeqIO.parse(input_handle, 'fasta'):
        lengths.append(len(record.seq))

    return lengths


def restrict(input_handle, enzymes):
    """
    Fragment a genome with restriction enzymes.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg list(str) enzymes: List of restiction enzymes.

    :return list(int): List of fragment sizes.
    """
    restriction_batch = Restriction.RestrictionBatch(enzymes)
    lengths = []

    for record in SeqIO.parse(input_handle, 'fasta'):
        positions = sorted(set([0, len(record.seq)] + 
            sum(restriction_batch.search(record.seq).values(), [])))

        for i in range(len(positions) - 1):
            lengths.append(positions[i + 1] - positions[i])

    return lengths


def collapse(word, max_stretch):
    """
    Collapse stretches of single letters in a word that exceed a certain
    length.

    :arg str word: Non empty input string.
    :arg int max_stretch: Maximum stretch of single letters, must be larger
        than 1.

    :return tuple(str, int): The collapsed word and the number of collapsed
        stretches.
    """
    stretch = 0
    collapsed_word = word[0]
    number_of_collapses = 0

    for i in range(1, len(word)):
        if word[i - 1] == word[i]:
            stretch += 1
        else:
            stretch = 0
        if stretch < max_stretch:
            collapsed_word += word[i]
        if stretch == max_stretch:
            number_of_collapses += 1

    return collapsed_word, number_of_collapses


def collapse_fasta(input_handle, output_handle, stretch):
    """
    Remove all mononucleotide stretches from a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writeable handle to a FASTA file.
    :arg int stretch: Maximum stretch of single letters, must be larger than 1.

    :return int: Number of collapsed stretches.
    """
    total_collapses = 0

    for record in SeqIO.parse(input_handle, 'fasta'):
        sequence, collapses = collapse(record.seq, stretch)
        record.seq = Seq.Seq(sequence)
        SeqIO.write(record, output_handle, 'fasta')
        total_collapses += collapses

    return total_collapses


def s2i(input_handle, output_handle):
    """
    Convert sanger FASTQ to illumina FASTQ.

    :arg stream input_handle: Open readable handle to a FASTQ file.
    :arg stream output_handle: Open writeable handle to a FASTQ file.
    """
    return SeqIO.convert(input_handle, 'fastq', output_handle,
        'fastq-illumina')


def count_tags(input_handle, sequence, mismatches):
    """
    Count tags in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg str sequence: The sequence that needs to be counted.
    :arg int mismatches: The number of mismatches allowed.

    :return int: Number of occurrences.
    """
    count = 0

    for record in SeqIO.parse(input_handle, 'fasta'):
        alignment = pairwise2.align.localms(str(record.seq),
            sequence, 1, -1, 0, -1)

        if alignment and len(sequence) - alignment[0][2] <= mismatches:
            count += 1

    return count


def select(input_handle, output_handle, first, last):
    """
    Select a substring from every read.
    Positions are one-based and inclusive.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    :arg int first: First base of the selection.
    :arg int last: Last base of the selection.
    """
    file_format = guess_file_format(input_handle)
    real_first = first - 1

    for record in SeqIO.parse(input_handle, file_format):
        SeqIO.write([record[real_first:last]], output_handle, file_format)


def rselect(input_handle, output_handle, name, first, last):
    """
    Select a substring from every read.
    Positions are one-based and inclusive.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    :arg str name: Accession number.
    :arg int first: First base of the selection.
    :arg int last: Last base of the selection.
    """
    file_format = guess_file_format(input_handle)
    real_first = first - 1

    for record in SeqIO.parse(input_handle, file_format):
        full_acc_no = record.name

        if '|' in record.name:
            full_acc_no = record.name.split('|')[3]

        accno = full_acc_no.split('.')[0]

        if accno == name:
            SeqIO.write([record[real_first:last]], output_handle, file_format)


def fa2gb(input_handle, output_handle, name):
    """
    Convert a FASTA file to a GenBank file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a GenBank file.
    :arg str name: A GenBank accession number.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        record.seq.alphabet = IUPAC.unambiguous_dna
        record.id = name
        record.name = name
        SeqIO.write(record, output_handle, 'genbank')


def gb2fa(input_handle, output_handle):
    """
    Convert a GenBank file to a FASTA file.

    :arg stream input_handle: Open readable handle to a GenBank file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'genbank'):
        SeqIO.write(record, output_handle, 'fasta')


def mangle(input_handle, output_handle):
    """
    Calculate the complement (not reverse-complement) of a FASTA sequence.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        seq = ''

        for i in record.seq:
            if i in ['A', 'a']:
                seq += 'T'
            if i in ['C', 'c']:
                seq += 'G'
            if i in ['G', 'g']:
                seq += 'C'
            if i in ['T', 't']:
                seq += 'A'

        new_record = SeqRecord(Seq.Seq(seq), record.id + 'C', '',
            'Complement (not reverse-complement) of the non-N part of %s' %
            record.id)
        SeqIO.write(new_record, output_handle, 'fasta')


def generate_dna(length, output_handle, name, description):
    """
    Generate a DNA sequence in FASTA format.

    :arg int length: Length of the DNA sequence.
    :arg stream output_handle: Open writable handle to a FASTA file.
    :arg str name: Name of the DNA sequence.
    :arg str description: Description of the DNA sequence.
    """
    dna = ['A', 'C', 'G', 'T']
    seq = ''

    for i in range(length):
        seq += dna[random.randint(0, 3)]

    record = SeqRecord(Seq.Seq(seq), name, '', description)
    SeqIO.write(record, output_handle, 'fasta')


def get_reference(name, email, output_handle, start=0, stop=0, orientation=0):
    """
    Retrieve a reference sequence and find the location of a specific gene.

    :arg str name: An accession number.
    :arg str email: An email address.
    :arg stream output_handle: An open writable handle.
    :arg int start: Start of the area of interest.
    :arg int stop: End of the area of interest.
    :arg int orientation: Orientation (1=forward, 2=reverse).
    """
    Entrez.email = email

    try:
        if start:
            handle = Entrez.efetch(db='nuccore', rettype='fasta', id=name,
                seq_start=start, seq_stop=stop, strand=orientation)
        else:
            handle = Entrez.efetch(db='nuccore', rettype='fasta', id=name)
    except urllib2.HTTPError:
        sys.stderr.write('Error: could not retrieve %s\n' % name)
        return

    output_handle.write(handle.read())


def cat(input_handle):
    """
    Return the sequence content of a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        print record.seq


def descr(input_handle):
    """
    Return the description of all records in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        print record.description


def splitseq(input_handle, output_handles, sequence):    
    """
    split a FASTA/FASTQ file based on containing part of the sequence 
 
    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handles: List of open writable handles to FASTA/FASTQ
        files.
    :arg str sequence: filter reads containing this sequence.
    """

    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format) :
        if sequence in record.seq:
            SeqIO.write(record,output_handles[0], file_format)
        else:
            SeqIO.write(record,output_handles[1], file_format)



def length_split(input_handle, output_handles, length):
    """
    Split a FASTA/FASTQ file on length.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handles: List of open writable handles to FASTA/FASTQ
        files.
    :arg int length: Length threshold.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        if len(record.seq) >= length:
            SeqIO.write([record], output_handles[0], file_format)
        else:
            SeqIO.write([record], output_handles[1], file_format)


def reverse(input_handle, output_handle):
    """
    Make the reverse complement a FASTA/FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        reverse_record = record.reverse_complement()
        reverse_record.id = record.id
        reverse_record.description = record.description
        SeqIO.write([reverse_record], output_handle, file_format)


def merge(input_handles, output_handle, fill):
    """
    Merge two FASTA files.

    :arg stream input_handles: List of open readable handle to FASTA/FASTQ
        files.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    :arg int fill: Amount of 'N's to be added between the reads.
    """
    for records in itertools.izip(SeqIO.parse(input_handles[0], 'fasta'),
            SeqIO.parse(input_handles[1], 'fasta')):
        record = SeqRecord(Seq.Seq(
            str(records[0].seq) + 'N' * fill + str(records[1].seq)),
            records[0].name, records[0].id, records[0].description)
        SeqIO.write([record], output_handle, 'fasta')


def find_motif(record, motif):
    """
    Find a certain sequence in a FASTA record.

    :arg SeqRecord record: Seq object which will be searched.
    :arg str motif: The sequence to be found.

    :returns generator(tuple(int, int)): tuple of start and end of matches in
        record.
    """
    regex = re.compile(motif.strip(), re.IGNORECASE)

    for match in regex.finditer(str(record.seq)):
        yield (int(match.start()), int(match.end()))


def fa_motif2bed(input_handle, output_handle, motif):
    """
    Find a given sequence in a FASTA file and write the results to a Bed file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a BED file.
    :arg str motif: The sequence to be found.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        for m in find_motif(record, motif):
            output_handle.write('\t'.join(map(str, [record.id, m[0], m[1]])) +
                '\n')


def edit(input_handle, bed_handle, output_handle):
    """
    Replace regions in a reference sequence. The BED file is structured as
    follows: "chrom start end replacement".

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream bed_handle: Open readable_handle to a BED file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    bed_records = bed_read(bed_handle)

    for record in SeqIO.parse(input_handle, 'fasta'):
        for bed_record in bed_records[record.name]:
            record.seq = (record.seq[:bed_record[0] - 1] + bed_record[2] +
                record.seq[bed_record[1]:])
        SeqIO.write([record], output_handle, 'fasta')


def main():
    """
    Main entry point.
    """
    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument('input_handle', metavar='INPUT',
        type=argparse.FileType('r'), help='input file')

    input2_parser = argparse.ArgumentParser(add_help=False)
    input2_parser.add_argument('input_handles', metavar='INPUT',
        type=argparse.FileType('r'), nargs=2, help='input files')

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument('output_handle', metavar='OUTPUT',
        type=argparse.FileType('w'), help='output file')

    output2_parser = argparse.ArgumentParser(add_help=False)
    output2_parser.add_argument('output_handles', metavar='OUTPUT',
        type=argparse.FileType('w'), nargs=2, help='output files')

    file_parser = argparse.ArgumentParser(add_help=False,
        parents=[input_parser, output_parser])

    qual_parser = argparse.ArgumentParser(add_help=False)
    qual_parser.add_argument('-q', dest='quality', type=int, default=40,
        help='quality score (%(type)s default=%(default)s)')

    seq_parser = argparse.ArgumentParser(add_help=False)
    seq_parser.add_argument('sequence', metavar='SEQ', type=str,
        help='a sequence (%(type)s)')

    range_parser = argparse.ArgumentParser(add_help=False)
    range_parser.add_argument('first', metavar='FIRST', type=int,
        help='first base of the selection (%(type)s)')
    range_parser.add_argument('last', metavar='LAST', type=int,
        help='last base of the selection (%(type)s)')

    name_parser = argparse.ArgumentParser(add_help=False)
    name_parser.add_argument('name', metavar='ACCNO', type=str,
        help='accession number')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action='version', version=version(parser.prog))
    subparsers = parser.add_subparsers(dest='subcommand')

    parser_sanitise = subparsers.add_parser('sanitise',
        parents=[file_parser], description=doc_split(sanitise))
    parser_sanitise.set_defaults(func=sanitise)

    parser_fa2fq = subparsers.add_parser('fa2fq',
        parents=[file_parser, qual_parser], description=doc_split(fa2fq))
    parser_fa2fq.set_defaults(func=fa2fq)

    parser_fq2fa = subparsers.add_parser('fq2fa', parents=[file_parser],
        description=doc_split(fq2fa))
    parser_fq2fa.set_defaults(func=fq2fa)

    parser_add = subparsers.add_parser('add', parents=[file_parser, seq_parser,
        qual_parser], description=doc_split(add))
    parser_add.set_defaults(func=add)

    parser_aln = subparsers.add_parser('aln', parents=[input2_parser],
        description=doc_split(aln))
    parser_aln.set_defaults(func=aln)

    parser_maln = subparsers.add_parser('maln', parents=[input_parser],
        description=doc_split(maln))
    parser_maln.set_defaults(func=maln)

    parser_len = subparsers.add_parser('len', parents=[input_parser],
        description=doc_split(length))
    parser_len.set_defaults(func=length)

    parser_restrict = subparsers.add_parser('restrict', parents=[input_parser],
        description=doc_split(restrict))
    parser_restrict.add_argument('-r', dest='enzymes', type=str, nargs='+',
        default=['EcoRI', 'MseI'], help='restriction enzymes')
    parser_restrict.set_defaults(func=restrict)

    parser_collapse = subparsers.add_parser('collapse', parents=[file_parser],
        description=doc_split(collapse_fasta))
    parser_collapse.add_argument('-s', '--stretch', dest='max_stretch',
        default=3, type=int,
        help='Length of the stretch (%(type)s default: %(default)s)')
    parser_collapse.set_defaults(func=collapse_fasta)

    parser_s2i = subparsers.add_parser('s2i', parents=[file_parser],
        description=doc_split(s2i))
    parser_s2i.set_defaults(func=s2i)

    parser_tagcount = subparsers.add_parser('tagcount', parents=[input_parser,
        seq_parser], description=doc_split(count_tags))
    parser_tagcount.add_argument('-m', dest='mismatches', type=int, default=2,
        help='amount of mismatches allowed (%(type)s default=%(default)s)')
    parser_tagcount.set_defaults(func=count_tags)

    parser_select = subparsers.add_parser('select', parents=[file_parser,
        range_parser], description=doc_split(select))
    parser_select.set_defaults(func=select)

    parser_rselect = subparsers.add_parser('rselect', parents=[file_parser,
        name_parser, range_parser], description=doc_split(rselect))
    parser_rselect.set_defaults(func=rselect)

    parser_fa2gb = subparsers.add_parser('fa2gb', parents=[file_parser,
        name_parser], description=doc_split(fa2gb))
    parser_fa2gb.set_defaults(func=fa2gb)

    parser_gb2fa = subparsers.add_parser('gb2fa', parents=[file_parser],
        description=doc_split(gb2fa))
    parser_gb2fa.set_defaults(func=gb2fa)

    parser_mangle = subparsers.add_parser('mangle', parents=[file_parser],
        description=doc_split(mangle))
    parser_mangle.set_defaults(func=mangle)

    parser_gen = subparsers.add_parser('gen', parents=[output_parser,
        name_parser], description=doc_split(generate_dna))
    parser_gen.add_argument('length', metavar='LENGTH', type=int,
        help='length of the DNA sequence')
    parser_gen.add_argument('description', metavar='DESCR', type=str,
        help='descriptino of the DNA sequence')
    parser_gen.set_defaults(func=generate_dna)

    parser_get = subparsers.add_parser('get', parents=[output_parser,
        name_parser], description=doc_split(get_reference))
    parser_get.add_argument('email', metavar='EMAIL', type=str,
        help='email address')
    parser_get.add_argument('-s', dest='start', type=int,
        help='start of the area of interest')
    parser_get.add_argument('-p', dest='stop', type=int,
        help='end of the area of interest')
    parser_get.add_argument('-o', dest='orientation', type=int,
        help='orientation (1=forward, 2=reverse)')
    parser_get.set_defaults(func=get_reference)

    parser_cat = subparsers.add_parser('cat', parents=[input_parser],
        description=doc_split(cat))
    parser_cat.set_defaults(func=cat)

    parser_descr = subparsers.add_parser('descr', parents=[input_parser],
        description=doc_split(descr))
    parser_descr.set_defaults(func=descr)

    parser_splitseq = subparsers.add_parser('splitseq', parents=[input_parser,
        output2_parser, seq_parser], description=doc_split(splitseq))
    parser_splitseq.set_defaults(func=splitseq)

    parser_lenfilt = subparsers.add_parser('lenfilt', parents=[input_parser,
        output2_parser], description=doc_split(length_split))
    parser_lenfilt.add_argument('-l', dest='length', type=int, default=25, 
        help='length threshold (%(type)s default: %(default)s)')
    parser_lenfilt.set_defaults(func=length_split)

    parser_reverse = subparsers.add_parser('reverse', parents=[file_parser],
        description=doc_split(reverse))
    parser_reverse.set_defaults(func=reverse)

    parser_merge = subparsers.add_parser('merge', parents=[input2_parser,
        output_parser], description=doc_split(merge))
    parser_merge.add_argument('-f', dest='fill', type=int, default=0,
        help="Add 'N's between the reads (%(type)s default: %(default)s)")
    parser_merge.set_defaults(func=merge)

    parser_fa_mot2bed = subparsers.add_parser('famotif2bed',
        parents=[file_parser], description=doc_split(fa_motif2bed))
    parser_fa_mot2bed.add_argument('motif', metavar='MOTIF', type=str,
        help='The sequence to be found')
    parser_fa_mot2bed.set_defaults(func=fa_motif2bed)

    parser_edit = subparsers.add_parser('edit',
        parents=[input_parser, output_parser], description=doc_split(edit))
    parser_edit.add_argument('bed_handle', metavar='BED',
        type=argparse.FileType('r'), help='BED file')
    parser_edit.set_defaults(func=edit)

    try:
        args = parser.parse_args()
    except IOError, error:
        parser.error(error)

    if args.subcommand == 'aln':
        for i in aln(args.input_handles):
            print '%s %s %i' % i

    elif args.subcommand == 'len':
        print ' '.join(map(lambda x: str(x), length(args.input_handle)))

    elif args.subcommand == 'restrict':
        print ' '.join(map(lambda x: str(x), restrict(args.input_handle,
            args.enzymes)))

    elif args.subcommand == 'collapse':
        print 'Collapsed %i stretches longer than %i.' % (collapse_fasta(
            args.input_handle, args.output_handle, args.max_stretch),
            args.max_stretch)

    elif args.subcommand == 's2i':
        print 'converted %i records' % s2i(args.input_handle,
            args.output_handle)

    elif args.subcommand == 'tagcount':
        print count_tags(args.input_handle, args.sequence, args.mismatches)

    else:
        try:
            args.func(**{k: v for k, v in vars(args).items()
                if k not in ('func', 'subcommand')})
        except ValueError, error:
            parser.error(error)


if __name__ == "__main__":
    main()
