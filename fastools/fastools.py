from sys import stderr, stdout

from csv import Sniffer, reader as csv_reader
from random import randint
from urllib.error import HTTPError

from Bio import Seq, SeqIO, Entrez, pairwise2, Restriction
from Bio.SeqRecord import SeqRecord
from Levenshtein import distance, hamming

from .utils import (
    guess_file_format, _collapse, _edits_read, _find_motif, _write_seq)


def add(input_handle, output_handle, sequence, quality):
    """Add a sequence to the 5' end of each read in a FASTQ file.

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
    """Calculate the Levenshtein distance between two FASTA files.

    :arg stream input_handles: Open readable handles to FASTA files.

    :return list(tuple(str, str, int)): List of distances.
    """
    distances = []

    for i in SeqIO.parse(input_handles[0], 'fasta'):
        for j in SeqIO.parse(input_handles[1], 'fasta'):
            distances.append(
                (i.name, j.name, distance(str(i.seq), str(j.seq))))

    return distances


def cat(input_handle):
    """Return the sequence content of a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        yield str(record.seq)


def collapse(input_handle, output_handle, stretch):
    """Remove all mononucleotide stretches from a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writeable handle to a FASTA file.
    :arg int stretch: Maximum stretch of single letters, must be larger than 1.

    :return int: Number of collapsed stretches.
    """
    total_collapses = 0

    for record in SeqIO.parse(input_handle, 'fasta'):
        sequence, collapses = _collapse(record.seq, stretch)
        record.seq = Seq.Seq(sequence)
        SeqIO.write(record, output_handle, 'fasta')
        total_collapses += collapses

    return total_collapses


def csv2fa2(input_handle, output_handles, skip_header=False):
    """Convert a CSV file to two FASTA files.

    :arg stream input_handle: Open readable handle to a CSV file.
    :arg list(stream) output_handles: List of open writable handles to FASTA
        files.
    :arg bool skip_header: Ignore the first line of the CSV file.
    """
    dialect = Sniffer().sniff(input_handle.read(1024))
    input_handle.seek(0)

    reader = csv_reader(input_handle, dialect)
    if skip_header:
        next(reader)
    for record in reader:
        _write_seq(output_handles[0], record[1], record[0])
        _write_seq(output_handles[1], record[2], record[0])


def descr(input_handle):
    """Return the description of all records in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        yield record.description


def dna2rna(input_handle, output_handle):
    """Convert the FASTA/FASTQ content from DNA to RNA.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        record.seq = record.seq.transcribe()
        SeqIO.write(record, output_handle, file_format)


def edit(input_handle, edits_handle, output_handle):
    """Replace regions in a reference sequence. The header of the edits file
    must have the following strucure:
    >name chrom:start_end

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream edits_handle: Open readable_handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    edits_records = _edits_read(edits_handle)

    for record in SeqIO.parse(input_handle, 'fasta'):
        for bed_record in edits_records[record.name]:
            record.seq = (
                record.seq[:bed_record[0] - 1] + bed_record[2] +
                record.seq[bed_record[1]:])
        SeqIO.write([record], output_handle, 'fasta')


def fa2fq(input_handle, output_handle, quality):
    """Convert a FASTA file to a FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a FASTQ file.
    :arg int quality: Quality score.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        record.letter_annotations = {'phred_quality':
            [quality] * len(record.seq)}
        SeqIO.write(record, output_handle, 'fastq')


def fa2gb(input_handle, output_handle, name):
    """Convert a FASTA file to a GenBank file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a GenBank file.
    :arg str name: A GenBank accession number.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        record.annotations["molecule_type"] = "DNA"
        record.id = name
        record.name = name
        SeqIO.write(record, output_handle, 'genbank')


def famotif2bed(input_handle, output_handle, motif):
    """Find a given sequence in a FASTA file and write the results to a Bed
    file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg stream output_handle: Open writable handle to a BED file.
    :arg str motif: The sequence to be found.
    """
    for record in SeqIO.parse(input_handle, 'fasta'):
        for m in _find_motif(record, motif):
            output_handle.write(
                '\t'.join(map(str, [record.id, m[0], m[1]])) + '\n')


def fq2fa(input_handle, output_handle):
    """Convert a FASTQ file to a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    try:
        for record in SeqIO.parse(input_handle, 'fastq'):
            SeqIO.write(record, output_handle, 'fasta')
    except ValueError as error:
        stderr.write('Error: {}\n'.format(error))
        emptyRecord = SeqRecord(Seq.Seq(''), '', '', '')
        SeqIO.write(emptyRecord, output_handle, 'fasta')


def gb2fa(input_handle, output_handle):
    """Convert a GenBank file to a FASTA file.

    :arg stream input_handle: Open readable handle to a GenBank file.
    :arg stream output_handle: Open writable handle to a FASTA file.
    """
    for record in SeqIO.parse(input_handle, 'genbank'):
        SeqIO.write(record, output_handle, 'fasta')


def gen(length, output_handle, name, description):
    """Generate a DNA sequence in FASTA format.

    :arg int length: Length of the DNA sequence.
    :arg stream output_handle: Open writable handle to a FASTA file.
    :arg str name: Name of the DNA sequence.
    :arg str description: Description of the DNA sequence.
    """
    dna = ['A', 'C', 'G', 'T']
    seq = ''

    for i in range(length):
        seq += dna[randint(0, 3)]

    record = SeqRecord(Seq.Seq(seq), name, '', description)
    SeqIO.write(record, output_handle, 'fasta')


def get(name, email, output_handle, start=0, stop=0, orientation=0):
    """Retrieve a reference sequence and find the location of a specific gene.

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
            handle = Entrez.efetch(
                db='nuccore', rettype='fasta', id=name, seq_start=start,
                seq_stop=stop, strand=orientation)
        else:
            handle = Entrez.efetch(db='nuccore', rettype='fasta', id=name)
    except HTTPError: # URLError if NCBI is down.
        stderr.write('Error: could not retrieve {}\n'.format(name))
        return

    output_handle.write(handle.read())


def length(input_handle):
    """Report the lengths of all FASTA records in a file.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :return list(int): List of lengths.
    """
    lengths = []

    for record in SeqIO.parse(input_handle, 'fasta'):
        lengths.append(len(record.seq))

    return lengths


def lenfilt(input_handle, output_handles, length):
    """Split a FASTA/FASTQ file on length.

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


def list_enzymes():
    """Return a list of supported restiction enzymes."""
    return Restriction.Restriction_Dictionary.rest_dict.keys()


def maln(input_handle):
    """Calculate the Hamming distance between all sequences in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :return list(tuple(str, str, int)): List of distances.
    """
    distances = []

    data = {x.name: str(x.seq) for x in SeqIO.parse(input_handle, 'fasta')}
    for j in data:
        stdout.write(j)
    stdout.write('\n')
    for i in data:
        stdout.write(i)
        for j in data:
            stdout.write(hamming(data[i], data[j]))
        stdout.write('\n')

    return distances


def mangle(input_handle, output_handle):
    """Calculate the complement (not reverse-complement) of a FASTA sequence.

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

        new_record = SeqRecord(
            Seq.Seq(seq), record.id + 'C', '',
            'Complement (not reverse-complement) of the non-N part of' +
            '{}'.format(record.id))
        SeqIO.write(new_record, output_handle, 'fasta')


def merge(input_handles, output_handle, fill):
    """Merge two FASTA files.

    :arg stream input_handles: List of open readable handle to FASTA/FASTQ
        files.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    :arg int fill: Amount of 'N's to be added between the reads.
    """
    for records in zip(
            SeqIO.parse(input_handles[0], 'fasta'),
            SeqIO.parse(input_handles[1], 'fasta')):
        record = SeqRecord(
            Seq.Seq(str(records[0].seq) + 'N' * fill + str(records[1].seq)),
            records[0].name, records[0].id, records[0].description)
        SeqIO.write([record], output_handle, 'fasta')


def raw2fa(input_handle, output_handle, name, description):
    """Make a FASTA file from a raw sequence.

    :arg stream input_handle: Open readable handle to a raw sequence file.
    :arg stream output_handle: Open writeable handle to a FASTA file.
    :arg str name: Name of the DNA sequence.
    :arg str description: Description of the DNA sequence.
    """
    record = SeqRecord(
        Seq.Seq(input_handle.readline().strip('\n')), name, '', description)
    SeqIO.write(record, output_handle, 'fasta')


def restrict(input_handle, enzymes):
    """Fragment a genome with restriction enzymes.

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


def reverse(input_handle, output_handle):
    """Make the reverse complement a FASTA/FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        reverse_record = record.reverse_complement()
        reverse_record.id = record.id
        reverse_record.description = record.description
        SeqIO.write([reverse_record], output_handle, file_format)


def rna2dna(input_handle, output_handle):
    """Convert the FASTA/FASTQ content from RNA to DNA.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        record.seq = record.seq.back_transcribe()
        SeqIO.write(record, output_handle, file_format)


def rselect(input_handle, output_handle, name, first, last):
    """Select a substring from every read.
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


def s2i(input_handle, output_handle):
    """Convert sanger FASTQ to illumina FASTQ.

    :arg stream input_handle: Open readable handle to a FASTQ file.
    :arg stream output_handle: Open writeable handle to a FASTQ file.
    """
    return SeqIO.convert(
        input_handle, 'fastq', output_handle, 'fastq-illumina')


def sanitise(input_handle, output_handle):
    """Convert a FASTA/FASTQ file to a standard FASTA/FASTQ file.

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handle: Open writable handle to a FASTA/FASTQ file.
    """
    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        SeqIO.write(record, output_handle, file_format)


def select(input_handle, output_handle, first, last):
    """Select a substring from every read.
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


def splitseq(input_handle, output_handles, sequence):
    """Split a FASTA/FASTQ file based on containing part of the sequence

    :arg stream input_handle: Open readable handle to a FASTA/FASTQ file.
    :arg stream output_handles: List of open writable handles to FASTA/FASTQ
        files.
    :arg str sequence: filter reads containing this sequence.
    """

    file_format = guess_file_format(input_handle)

    for record in SeqIO.parse(input_handle, file_format):
        if sequence in record.seq:
            SeqIO.write(record,output_handles[0], file_format)
        else:
            SeqIO.write(record,output_handles[1], file_format)


def tagcount(input_handle, sequence, mismatches):
    """Count tags in a FASTA file.

    :arg stream input_handle: Open readable handle to a FASTA file.
    :arg str sequence: The sequence that needs to be counted.
    :arg int mismatches: The number of mismatches allowed.

    :return int: Number of occurrences.
    """
    count = 0

    for record in SeqIO.parse(input_handle, 'fasta'):
        alignment = pairwise2.align.localms(
            str(record.seq), sequence, 1, -1, -1, -1)

        if alignment and len(sequence) - alignment[0][2] <= mismatches:
            count += 1

    return count
