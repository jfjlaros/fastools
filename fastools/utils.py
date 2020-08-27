from collections import defaultdict
from re import compile as re_compile, split as re_split, IGNORECASE

from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord


def _collapse(word, max_stretch):
    """Collapse stretches of single letters in a word that exceed a certain
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


def _edits_read(handle):
    """Parse a FASTA file that contains edits.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :returns dict: A list of edits (ranges and replacements) per chromosome.
    """
    records = defaultdict(list)

    for record in SeqIO.parse(handle, 'fasta'):
         chrom, start, end = re_split(':|_', record.description.split()[-1])
         records[chrom].append([int(start), int(end), record.seq])
    for reference in records:
        records[reference].sort(reverse=True)
    return records


def _find_motif(record, motif):
    """Find a certain sequence in a FASTA record.

    :arg SeqRecord record: Seq object which will be searched.
    :arg str motif: The sequence to be found.

    :returns generator(tuple(int, int)): tuple of start and end of matches in
        record.
    """
    regex = re_compile(motif.strip(), IGNORECASE)

    for match in regex.finditer(str(record.seq)):
        yield (int(match.start()), int(match.end()))


def _write_seq(handle, seq, name, file_format='fasta'):
    record = SeqRecord(Seq.Seq(seq), name, '', '')
    SeqIO.write(record, handle, file_format)


def guess_file_format(handle):
    """Guess the file type of an NGS data file.

    :arg file handle: Open readable handle to an NGS data file.

    :return str: Either 'fasta' or 'fastq'.
    """
    if not hasattr(handle, 'name') or handle.name != '<stdin>':
        token = handle.read(1)
        handle.seek(0)
    else:
        token = handle.peek(1)

    if token == '>':
        return 'fasta'
    return 'fastq'


def guess_header_format(handle):
    """Guess the header format.

    :arg stream handle: Open readable handle to an NGS data file.

    :return str: Either 'normal', 'x' or 'unknown'.
    """
    if not hasattr(handle, 'name') or handle.name != '<stdin>':
        line = handle.readline().strip('\n')
        handle.seek(0)
    else:
        line = handle.peek(1024).split('\n')[0]

    if line.count('#') == 1 and line.split('#')[1].count('/') == 1:
        return 'normal'
    if line.count(' ') == 1 and line.split(' ')[1].count(':') == 3:
        return 'x'
    return 'unknown'
