import re
from collections import defaultdict
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from . import exceptions

class SeqExtractor(object):

    def __init__(self, input_handle, output_handle, location, bases_start, bases_end ):

        self.input_handle = input_handle
        self.output_handle = output_handle
        self.location = location
        self.bases_start = bases_start
        self.bases_end  = bases_end

        self._value_pair_checker()

    def _value_pair_checker(self):

        if self.location == 'start':

            if self.bases_start is None:
                raise exceptions.InvalidPairValueError("Ivalid Value-Pair:"+self.location+"-"+str(self.bases_start))

        if self.location == 'end':

            if self.bases_end is None:
                raise exceptions.InvalidPairValueError("Ivalid Value-Pair:"+self.location+"-"+ str(self.bases_end))

        if self.location == 'both_ends':
            if self.bases_start is None or self.bases_end is None:
                raise exceptions.InvalidPairValueError("Ivalid Value-Pair:"+self.location+"-"+ str(self.bases_start) +":" + str(self.bases_end))


def guess_file_format(handle):
    """
    Guess the file type of an NGS data file.

    :arg file handle: Open readable handle to an NGS data file.

    :return str: Either 'fasta' or 'fastq'.
    """
    if handle.name != '<stdin>':
        token = handle.read(1)
        handle.seek(0)
    else:
        token = handle.peek(1)

    if token == '>':
        return 'fasta'
    return 'fastq'


def guess_header_format(handle):
    """
    Guess the header format.

    :arg stream handle: Open readable handle to an NGS data file.

    :return str: Either 'normal', 'x' or 'unknown'.
    """
    if handle.name != '<stdin>':
        line = handle.readline().strip('\n')
        handle.seek(0)
    else:
        line = handle.peek(1024).split('\n')[0]

    if line.count('#') == 1 and line.split('#')[1].count('/') == 1:
        return 'normal'
    if line.count(' ') == 1 and line.split(' ')[1].count(':') == 3:
        return 'x'
    return 'unknown'


def _write_seq(handle, seq, name, file_format='fasta'):
    record = SeqRecord(Seq.Seq(seq), name, '', '')
    SeqIO.write(record, handle, file_format)


def _edits_read(handle):
    """
    Parse a FASTA file that contains edits.

    :arg stream input_handle: Open readable handle to a FASTA file.

    :returns dict: A list of edits (ranges and replacements) per chromosome.
    """
    records = defaultdict(list)

    for record in SeqIO.parse(handle, 'fasta'):
         chrom, start, end = re.split(':|_', record.description.split()[-1])
         records[chrom].append([int(start), int(end), record.seq])
    for reference in records:
        records[reference].sort(reverse=True)
    return records


def _find_motif(record, motif):
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


def batch_iterator(self,iterator, batch_size):
    """
    Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    :param self:
    :param iterator:
    :param batch_size:
    :return:
    """
    entry = True
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)

        if batch:
            yield batch

























