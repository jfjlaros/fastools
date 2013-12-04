#!/usr/bin/python

import sys
import random
import urllib2
import argparse
import itertools
import Levenshtein
from Bio import Seq, SeqIO, Entrez, pairwise2, Restriction
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from . import docSplit, version, usage

def guessFileType(handle):
    """
    Guess the file type of an NGS data file.

    We assume that the stream is rewinded before use, after use, the input
    stream will be rewinded.

    @arg handle: Open readable handle to an NGS data file.
    @type handle: stream

    @returns: Either "fasta" or "fastq".
    @rtype: str
    """
    if handle.isatty():
        sys.stderr.write("Cannot deterine file type in stream, assuming FASTA")
        return "fasta"
    #if

    token = handle.read(1)
    handle.seek(0)

    if token == '>':
        return "fasta"
    return "fastq"
#guessFileType

def sanitise(inputHandle, outputHandle):
    """
    Convert a FASTA/FASTQ file to a standard FASTA/FASTQ file.

    @arg inputHandle: Open readable handle to a FASTA/FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA/FASTQ file.
    @type outputHandle: stream
    """
    fileFormat = guessFileType(inputHandle)

    for record in SeqIO.parse(inputHandle, fileFormat):
        SeqIO.write(record, outputHandle, fileFormat)
#sanitise

def fa2fq(inputHandle, outputHandle, quality):
    """
    Convert a FASTA file to a FASTQ file.

    @arg inputHandle: Open readable handle to a FASTA file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTQ file.
    @type outputHandle: stream
    @arg quality: Quality score.
    @type quality: int
    """
    for record in SeqIO.parse(inputHandle, "fasta"):
        record.letter_annotations = {"phred_quality":
            [quality] * len(record.seq)}
        SeqIO.write(record, outputHandle, "fastq")
    #for
#fa2fq

def fq2fa(inputHandle, outputHandle):
    """
    Convert a FASTQ file to a FASTA file.

    @arg inputHandle: Open readable handle to a FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA file.
    @type outputHandle: stream
    """
    try:
        for record in SeqIO.parse(inputHandle, "fastq"):
            SeqIO.write(record, outputHandle, "fasta")
    except ValueError, error:
        print "Error: %s" % error
        emptyRecord = SeqRecord(Seq.Seq(""), "", "", "")
        SeqIO.write(emptyRecord, outputHandle, "fasta")
    #except
#fq2fa

def add(inputHandle, outputHandle, sequence, quality):
    """
    Add a sequence to the 5' end of each read in a FASTQ file.

    @arg inputHandle: Open readable handle to a FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTQ file.
    @type outputHandle: stream
    @arg sequence: Sequence to be added to the 5' end of the read.
    @type sequence: str
    @arg quality: Quality score.
    @type quality: int
    """
    addition_q = [quality] * len(sequence)

    for record in SeqIO.parse(inputHandle, "fastq"):
        qual = record.letter_annotations["phred_quality"]
        record.letter_annotations = {}
        record.seq = sequence + record.seq
        record.letter_annotations = {"phred_quality": addition_q + qual}
        SeqIO.write(record, outputHandle, "fastq")
    #for
#add

def aln(fastaHandle1, fastaHandle2):
    """
    Calculate the Levenshtein distance between two FASTA files.

    @arg fastaHandle1: Open readable handle to a FASTA file.
    @type fastaHandle1: stream
    @arg fastaHandle2: Open readable handle to a FASTA file.
    @type fastaHandle2: stream

    @returns: List of distances.
    @rtype: list[tuple(str, str, int)]
    """
    distances = []

    for i in SeqIO.parse(fastaHandle1, "fasta"):
        for j in SeqIO.parse(fastaHandle2, "fasta"):
            distances.append((i.name, j.name,
                Levenshtein.distance(str(i.seq), str(j.seq))))

    return distances
#aln

def length(handle):
    """
    Report the lengths of all FASTA records in a file.

    @arg handle: Open readable handle to a FASTA file.
    @type handle: stream

    @returns: List of lengths.
    @rtype: list[int]
    """
    lengths = []

    for record in SeqIO.parse(handle, "fasta"):
        lengths.append(len(record.seq))

    return lengths
#length

def restrict(handle, enzymes):
    """
    Fragment a genome with restriction enzymes.

    @arg handle: Open readable handle to a FASTA file.
    @type handle: stream
    @arg enzymes: List of restiction enzymes.
    @type enzymes: list(str)

    @returns: List of fragment sizes.
    @rtype: list[int]
    """
    restrictionBatch = Restriction.RestrictionBatch(enzymes)
    lengths = []

    for record in SeqIO.parse(handle, "fasta"):
        positions = sorted(set([0, len(record.seq)] + 
            sum(restrictionBatch.search(record.seq).values(), [])))

        for i in range(len(positions) - 1):
            lengths.append(positions[i + 1] - positions[i])
    #for

    return lengths
#restrict

def collapse(word, maxStretch):
    """
    Collapse stretches of single letters in a word that exceed a certain
    length.

    @arg word: Non empty input string.
    @type word: str
    @arg maxStretch: Maximum stretch of single letters, must be larger than 1.
    @type maxStretch: int

    @returns: The collapsed word and the number of collapsed stretches.
    @rtype: tuple(srt, int)
    """
    stretch = 0
    collapsedWord = word[0]
    numberOfCollapses = 0

    for i in range(1, len(word)):
        if word[i - 1] == word[i]:
            stretch += 1
        else:
            stretch = 0
        if stretch < maxStretch:
            collapsedWord += word[i]
        if stretch == maxStretch:
            numberOfCollapses += 1
    #for

    return collapsedWord, numberOfCollapses
#collapse

def collapseFasta(inputHandle, outputHandle, stretch):
    """
    Remove all mononucleotide stretches from a FASTA file.

    @arg inputHandle: Open readable handle to a FASTA file.
    @type inputHandle: stream
    @arg outputHandle: Open writeable handle to a FASTA file.
    @type outputHandle: stream
    @arg stretch: Maximum stretch of single letters, must be larger than 1.
    @type stretch: int

    @returns: Number of collapsed stretches.
    @rtype: int
    """
    totalCollapses = 0

    for record in SeqIO.parse(inputHandle, "fasta"):
        sequence, collapses = collapse(record.seq, stretch)
        record.seq = Seq.Seq(sequence)
        SeqIO.write(record, outputHandle, "fasta")
        totalCollapses += collapses
    #for

    return totalCollapses
#collapseFasta

def s2i(inputHandle, outputHandle):
    """
    Convert sanger FASTQ to illumina FASTQ.

    @arg inputHandle: Open readable handle to a FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writeable handle to a FASTQ file.
    @type outputHandle: stream
    """
    return SeqIO.convert(inputHandle, "fastq", outputHandle, "fastq-illumina")
#s2i

def countTags(inputHandle, tag, mismatches):
    """
    Count tags in a FASTA file.

    @arg inputHandle: Open readable handle to a FASTA file.
    @type inputHandle: stream
    @arg tag: The tag that needs to be counted.
    @type tag: str
    @arg mismatches: The number of mismatches allowed.
    @type mismatches: int

    @returns: Number of occurrences.
    @rtype: int
    """
    count = 0

    for record in SeqIO.parse(inputHandle, "fasta"):
        alignment = pairwise2.align.localms(str(record.seq),
            tag, 1, -1, 0, -1)

        if alignment and len(tag) - alignment[0][2] <= mismatches:
            count += 1
    #for

    return count
#countTags

def select(inputHandle, outputHandle, first, last):
    """
    Select a substring from every read.
    Positions are one-based and inclusive.

    @arg inputHandle: Open readable handle to a FASTA/FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA/FASTQ file.
    @type outputHandle: stream
    @arg first: First base of the selection.
    @type first: int
    @arg last: Last base of the selection.
    @type last: int
    """
    fileFormat = guessFileType(inputHandle)
    realFirst = first - 1

    for record in SeqIO.parse(inputHandle, fileFormat):
        SeqIO.write([record[realFirst:last]], outputHandle, fileFormat)
#select

def rselect(inputHandle, outputHandle, name, first, last):
    """
    Select a substring from every read.
    Positions are one-based and inclusive.

    @arg inputHandle: Open readable handle to a FASTA/FASTQ file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA/FASTQ file.
    @type outputHandle: stream
    @arg name: Accession number.
    @type name: str
    @arg first: First base of the selection.
    @type first: int
    @arg last: Last base of the selection.
    @type last: int
    """
    fileFormat = guessFileType(inputHandle)
    realFirst = first - 1

    for record in SeqIO.parse(inputHandle, fileFormat):
        fullAccNo = record.name

        if '|' in record.name:
            fullAccNo = record.name.split('|')[3]

        accno = fullAccNo.split('.')[0]

        if accno == name:
            SeqIO.write([record[realFirst:last]], outputHandle, fileFormat)
    #for
#rselect

def fa2gb(inputHandle, outputHandle, accno):
    """
    Convert a FASTA file to a GenBank file.

    @arg inputHandle: Open readable handle to a FASTA file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a GenBank file.
    @type outputHandle: stream
    @arg accno: A GenBank accession number.
    @type accno: str
    """
    for record in SeqIO.parse(inputHandle, "fasta"):
        record.seq.alphabet = IUPAC.unambiguous_dna
        record.id = accno
        record.name = accno
        SeqIO.write(record, outputHandle, "genbank")
    #for
#fa2gb

def gb2fa(inputHandle, outputHandle):
    """
    Convert a GenBank file to a FASTA file.

    @arg inputHandle: Open readable handle to a GenBank file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA file.
    @type outputHandle: stream
    """
    for record in SeqIO.parse(inputHandle, "genbank"):
        SeqIO.write(record, outputHandle, "fasta")
#gb2fa

def mangle(inputHandle, outputHandle):
    """
    Calculate the complement (not reverse-complement) of a FASTA sequence.

    @arg inputHandle: Open readable handle to a FASTA file.
    @type inputHandle: stream
    @arg outputHandle: Open writable handle to a FASTA file.
    @type outputHandle: stream
    """
    for record in SeqIO.parse(inputHandle, "fasta"):
        seq = ""

        for i in record.seq:
            if i in ['A', 'a']:
                seq += 'T'
            if i in ['C', 'c']:
                seq += 'G'
            if i in ['G', 'g']:
                seq += 'C'
            if i in ['T', 't']:
                seq += 'A'
        #for

        newRecord = SeqRecord(Seq.Seq(seq), record.id + 'C', "",
            "Complement (not reverse-complement) of the non-N part of %s" %
            record.id)
        SeqIO.write(newRecord, outputHandle, "fasta")
    #for
#mangle

def generateDNA(size, handle, name, description):
    """
    Generate a DNA sequence in FASTA format.

    @arg size: Length of the DNA sequence.
    @arg size: int
    @arg handle: Open writable handle to a FASTA file.
    @type handle: stream
    @arg name: Name of the DNA sequence.
    @arg name: str
    @arg description: Description of the DNA sequence.
    @arg description: str
    """
    DNA = ['A', 'C', 'G', 'T']
    seq = ""

    for i in range(size):
        seq += DNA[random.randint(0, 3)]

    record = SeqRecord(Seq.Seq(seq), name, "", description)
    SeqIO.write(record, handle, "fasta")
#generateDNA

def getReference(acc, email, outputHandle, start=0, stop=0, orientation=0):
    """
    Retrieve a reference sequence and find the location of a specific gene.

    @arg acc: An accession number.
    @type acc: str
    @arg email: An email address.
    @type email: str
    @arg outputHandle: An open writable handle.
    @type outputHandle: stream
    @arg start: Start of the area of interest.
    @type start: int
    @arg stop: End of the area of interest.
    @type stop: int
    @arg orientation: Orientation (1=forward, 2=reverse).
    @type orientation: int
    """
    Entrez.email = email

    try:
        if start:
            handle = Entrez.efetch(db="nuccore", rettype="fasta", id=acc,
                seq_start=start, seq_stop=stop, strand=orientation)
        else:
            handle = Entrez.efetch(db="nuccore", rettype="fasta", id=acc)
    except urllib2.HTTPError:
        sys.stderr.write("Error: could not retrieve %s\n" % acc)
        return
    #except

    outputHandle.write(handle.read())
#getReference

def cat(handle):
    """
    Return the sequence content of a FASTA file.

    @arg handle: Open readable handle to a FASTA file.
    @type handle: stream
    """
    for record in SeqIO.parse(handle, "fasta"):
        print record.seq
#cat

def descr(handle):
    """
    Return the description of all records in a FASTA file.

    @arg handle: Open readable handle to a FASTA file.
    @type handle: stream
    """
    for record in SeqIO.parse(handle, "fasta"):
        print record.description
#descr

def lengthSplit(inputHandle, outputHandles, length):
    """
    Split a fasta/fastq file on length.

    @arg inputHandle: Open readable handle to a fasta/fastq file.
    @type inputHandle: stream
    @arg outputHandles: List of open writable handles to fasta/fastq files.
    @type outputHandles: list[stream]
    @arg length: Length threshold.
    @type length: int
    """
    fileType = guessFileType(inputHandle)

    for record in SeqIO.parse(inputHandle, fileType):
        if len(record.seq) >= length:
            SeqIO.write([record], outputHandles[0], fileType)
        else:
            SeqIO.write([record], outputHandles[1], fileType)
#lengthSplit

def reverse(input_handle, output_handle):
    """
    Make the reverse complement a fasta/fastq file.

    @arg input_handle: Open readable handle to a fasta/fastq file.
    @type input_handle: stream
    @arg output_handle: Open writable handle to a fasta/fastq file.
    @type output_handle: stream
    """
    file_type = guessFileType(input_handle)

    for record in SeqIO.parse(input_handle, file_type):
        reverse_record = record.reverse_complement()
        reverse_record.id = record.id
        reverse_record.description = record.description
        SeqIO.write([reverse_record], output_handle, file_type)
    #for
#reverse

def merge(input_handles, output_handle, fill):
    """
    Merge two fasta files.

    @arg input_handle: List of open readable handle to fasta/fastq files.
    @type input_handle: list(stream)
    @arg output_handle: Open writable handle to a fasta/fastq file.
    @type output_handle: stream
    @arg fill: Amount of 'N's to be added between the reads.
    @type fill: int
    """
    for records in itertools.izip(SeqIO.parse(input_handles[0], "fasta"),
            SeqIO.parse(input_handles[1], "fasta")):
        record = SeqRecord(Seq.Seq(
            str(records[0].seq) + 'N' * fill + str(records[1].seq)),
            records[0].name, records[0].id, records[0].description)
        SeqIO.write([record], output_handle, "fasta")
    #for
#merge

def main():
    """
    Main entry point.
    """
    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("INPUT", type=argparse.FileType('r'),
        help="input file")

    input2_parser = argparse.ArgumentParser(add_help=False)
    input2_parser.add_argument("INPUT", type=argparse.FileType('r'), nargs=2,
        help="input files")

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("OUTPUT", type=argparse.FileType('w'),
        help="output file")

    output2_parser = argparse.ArgumentParser(add_help=False)
    output2_parser.add_argument("OUTPUT", type=argparse.FileType('w'), nargs=2,
        help="output files")

    file_parser = argparse.ArgumentParser(add_help=False,
        parents=[input_parser, output_parser])

    qual_parser = argparse.ArgumentParser(add_help=False)
    qual_parser.add_argument("-q", dest="quality", type=int, default=40,
        help="quality score (%(type)s default=%(default)s)")

    seq_parser = argparse.ArgumentParser(add_help=False)
    seq_parser.add_argument("SEQ", type=str, help="a sequence")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers(dest="subcommand")

    parser_sanitise = subparsers.add_parser("sanitise",
        parents=[file_parser], description=docSplit(sanitise))

    parser_fa2fq = subparsers.add_parser("fa2fq",
        parents=[file_parser, qual_parser], description=docSplit(fa2fq))

    parser_fq2fa = subparsers.add_parser("fq2fa", parents=[file_parser],
        description=docSplit(fq2fa))

    parser_add = subparsers.add_parser("add", parents=[file_parser, seq_parser,
        qual_parser], description=docSplit(add))

    parser_aln = subparsers.add_parser("aln", parents=[input2_parser],
        description=docSplit(aln))

    parser_len = subparsers.add_parser("len", parents=[input_parser],
        description=docSplit(length))

    parser_restrict = subparsers.add_parser("restrict", parents=[input_parser],
        description=docSplit(restrict))
    parser_restrict.add_argument("-r", dest="enzyme", type=str, nargs="+",
        default=["EcoRI", "MseI"], help="restriction enzymes")

    parser_collapse = subparsers.add_parser("collapse", parents=[file_parser],
        description=docSplit(collapseFasta))
    parser_collapse.add_argument('-s', '--stretch', dest='stretch', default=3,
        type=int, help='Length of the stretch (%(type)s default: %(default)s)')

    parser_s2i = subparsers.add_parser("s2i", parents=[file_parser],
        description=docSplit(s2i))

    parser_tagcount = subparsers.add_parser("tagcount", parents=[input_parser,
        seq_parser], description=docSplit(countTags))
    parser_tagcount.add_argument("-m", dest="mismatches", type=int, default=2,
        help="amount of mismatches allowed (%(type)s default=%(default)s)")

    parser_select = subparsers.add_parser("select", parents=[file_parser],
        description=docSplit(select))
    parser_select.add_argument("FIRST", type=int,
        help="first base of the selection (%(type)s)")
    parser_select.add_argument("LAST", type=int,
        help="last base of the selection (%(type)s)")

    parser_rselect = subparsers.add_parser("rselect", parents=[file_parser],
        description=docSplit(rselect))
    parser_rselect.add_argument("NAME", type=str,
        help="accession number")
    parser_rselect.add_argument("FIRST", type=int,
        help="first base of the selection (%(type)s)")
    parser_rselect.add_argument("LAST", type=int,
        help="last base of the selection (%(type)s)")

    parser_fa2gb = subparsers.add_parser("fa2gb", parents=[file_parser],
        description=docSplit(fa2gb))
    parser_fa2gb.add_argument("ACCNO", type=str,
        help="GenBank accession number")

    parser_gb2fa = subparsers.add_parser("gb2fa", parents=[file_parser],
        description=docSplit(gb2fa))

    parser_mangle = subparsers.add_parser("mangle", parents=[file_parser],
        description=docSplit(mangle))

    parser_gen = subparsers.add_parser("gen", parents=[output_parser],
        description=docSplit(generateDNA))
    parser_gen.add_argument("LENGTH", type=int,
        help="length of the DNA sequence")
    parser_gen.add_argument("NAME", type=str,
        help="name of the DNA sequence")
    parser_gen.add_argument("DESCR", type=str,
        help="descriptino of the DNA sequence")

    parser_get = subparsers.add_parser("get", parents=[output_parser],
        description=docSplit(getReference))
    parser_get.add_argument("ACC", type=str,
        help="accession number")
    parser_get.add_argument("EMAIL", type=str,
        help="email address")
    parser_get.add_argument("-s", dest="start", type=int,
        help="start of the area of interest")
    parser_get.add_argument("-p", dest="stop", type=int,
        help="end of the area of interest")
    parser_get.add_argument("-o", dest="orientation", type=int,
        help="orientation (1=forward, 2=reverse)")

    parser_cat = subparsers.add_parser("cat", parents=[input_parser],
        description=docSplit(cat))

    parser_cat = subparsers.add_parser("descr", parents=[input_parser],
        description=docSplit(descr))

    parser_lenfilt = subparsers.add_parser("lenfilt", parents=[input_parser,
        output2_parser], description=docSplit(lengthSplit))
    parser_lenfilt.add_argument("-l", dest="length", type=int, default=25, 
        help="length threshold (%(type)s default: %(default)s)")

    parser_reverse = subparsers.add_parser("reverse", parents=[file_parser],
        description=docSplit(reverse))

    parser_merge = subparsers.add_parser("merge", parents=[input2_parser,
        output_parser], description=docSplit(merge))
    parser_merge.add_argument("-f", dest="fill", type=int, default=0,
        help="Add 'N's between the reads (%(type)s default: %(default)s)")

    args = parser.parse_args()

    if args.subcommand == "sanitise":
        sanitise(args.INPUT, args.OUTPUT)

    if args.subcommand == "fa2fq":
        fa2fq(args.INPUT, args.OUTPUT, args.quality)

    if args.subcommand == "fq2fa":
        fq2fa(args.INPUT, args.OUTPUT)

    if args.subcommand == "add":
        add(args.INPUT, args.OUTPUT, args.SEQ, args.quality)

    if args.subcommand == "aln":
        for i in aln(*args.INPUT):
            print "%s %s %i" % i

    if args.subcommand == "len":
        print ' '.join(map(lambda x: str(x), length(args.INPUT)))

    if args.subcommand == "restrict":
        print ' '.join(map(lambda x: str(x), restrict(args.INPUT,
            args.enzyme)))

    if args.subcommand == "collapse":
        collapses = collapseFasta(args.INPUT, args.OUTPUT, args.stretch)

        print "Collapsed %i stretches longer than %i." % (collapses,
            args.stretch)
    #if

    if args.subcommand == "s2i":
        print "converted %i records" % s2i(args.INPUT, args.OUTPUT)

    if args.subcommand == "tagcount":
        print countTags(args.INPUT, args.SEQ, args.mismatches)

    if args.subcommand == "select":
        select(args.INPUT, args.OUTPUT, args.FIRST, args.LAST)

    if args.subcommand == "rselect":
        rselect(args.INPUT, args.OUTPUT, args.NAME, args.FIRST, args.LAST)

    if args.subcommand == "fa2gb":
        fa2gb(args.INPUT, args.OUTPUT, args.ACCNO)

    if args.subcommand == "gb2fa":
        gb2fa(args.INPUT, args.OUTPUT)

    if args.subcommand == "mangle":
        mangle(args.INPUT, args.OUTPUT)

    if args.subcommand == "gen":
        generateDNA(args.LENGTH, args.OUTPUT, args.NAME, args.DESCR)

    if args.subcommand == "get":
        getReference(args.ACC, args.EMAIL, args.OUTPUT, args.start, args.stop,
            args.orientation)

    if args.subcommand == "cat":
        cat(args.INPUT)

    if args.subcommand == "descr":
        descr(args.INPUT)

    if args.subcommand == "lenfilt":
        lengthSplit(args.INPUT, args.OUTPUT, args.length)

    if args.subcommand == "reverse":
        reverse(args.INPUT, args.OUTPUT)

    if args.subcommand == "merge":
        merge(args.INPUT, args.OUTPUT, args.fill)
#main

if __name__ == "__main__":
    main()
