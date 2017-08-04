import re
from collections import defaultdict
from Bio import Seq, SeqIO
from Bio.SeqRecord import SeqRecord


class UMIExtractor(object):

    def __init__(self, input_handle, output_handle, location, bases_start, bases_end, file_format ):

        self.input_handle = input_handle
        self.output_handle = output_handle
        self.location = location
        self.bases_start = bases_start
        self.bases_end = bases_end
        self.file_format = file_format

        self._value_pair_checker()

    def extractor(self):

        record_iter = SeqIO.parse(open(self.input_handle.name), self.file_format)
        # batching makes it faster
        for i, batch in enumerate(batch_iterator(record_iter, 10000)):

            for record in batch:

                if self.location == 'start':
                    seq_identifier, fixed_seq, fixed_qual = self.extract_from_start(record)
                elif self.location == 'end':
                    seq_identifier, fixed_seq, fixed_qual = self.extract_from_end(record)
                elif self.location == 'both_ends':
                    seq_identifier, fixed_seq, fixed_qual = self.extract_from_both_ends(record)
                else:
                    seq_identifier = self.extract_from_index2(record)
                    fixed_seq = record.seq
                    fixed_qual =  record.letter_annotations["phred_quality"]

                record.letter_annotations = {}
                record.seq = fixed_seq
                record.letter_annotations["phred_quality"] = fixed_qual
                fixeid = self.append_umi_identifier_to_read_header(record.description,seq_identifier)
                record.description = ""
                record.id = str(fixeid)

                SeqIO.write(record, self.output_handle, self.file_format)


    def append_umi_identifier_to_read_header(self,description,identifier):

        header_format = guess_header_format(self.input_handle)

        if header_format == "x":
            return description.split(" ")[0] + "_" + identifier + " " + description.split(" ")[1]
        elif header_format == "normal":
            return description.split("#")[0] + "_" + identifier + "#" + description.split("#")[1]
        else:
            raise RuntimeError("Not Valid Header Format")


    def extract_from_index2(self,record):

        header_format = guess_header_format(self.input_handle)

        if header_format == "x":
           try:
               return  record.description.split("+")[1]
           except:
               raise RuntimeError("No second index found in fasta header: "+ record.description)

        elif header_format == "normal":
            try:
                # not sure though if second index supported here. COuld not find something for this format
                return record.description.split('#')[1].split("+")[1].split("/")[0]

            except:
                raise RuntimeError("No second index found in fasta header: "+ record.description)
        else:
            return RuntimeError("Uknown Header Format")


    def extract_from_start(self, record):
        return record.seq[0:self.bases_start], record.seq[self.bases_start:len(record.seq)], record.letter_annotations["phred_quality"][self.bases_start:len(record.seq)]

    def extract_from_end(self, record):
        return record.seq[len(record.seq) - self.bases_end:len(record.seq)], record.seq[0:len(record.seq) - self.bases_end], record.letter_annotations["phred_quality"][0:len(record.seq) - self.bases_end]

    def extract_from_both_ends(self, record):
            return record.seq[0:self.bases_start] + record.seq[len(record.seq) - self.bases_end:len(record.seq)], record.seq[self.bases_start:len(record.seq) - self.bases_end], record.letter_annotations["phred_quality"][self.bases_start:len(record.seq) - self.bases_end]


    def _value_pair_checker(self):

        if self.location == 'start':

            if self.bases_start is None:
                raise Exception("Ivalid Value-Pair:"+self.location+"-"+str(self.bases_start))

        if self.location == 'end':

            if self.bases_end is None:
                raise Exception("Ivalid Value-Pair:"+self.location+"-"+ str(self.bases_end))

        if self.location == 'both_ends':
            if self.bases_start is None or self.bases_end is None:
                raise Exception("Ivalid Value-Pair:"+self.location+"-"+ str(self.bases_start) +":" + str(self.bases_end))


class UMIAppender(object):

    def __init__(self, input_handle, output_handle,file_format ):

        self.input_handle = input_handle
        self.output_handle = output_handle
        self.file_format = file_format

    def get_umi_identifier_from_read_header(self,description):

        header_format = guess_header_format(self.input_handle)

        if header_format == "x":
            return  description.split("_")[1].split(" ")[0]
        elif header_format == "normal":
            return  description.split("_")[1].split("#")[0]
        else:
            raise RuntimeError("Not Valid Header Format")

    def umi_appender(self):


        record_iter = SeqIO.parse(open(self.input_handle.name), self.file_format)

        for i, batch in enumerate(batch_iterator(record_iter, 10000)):
            for record in batch:

                umi = self.get_umi_identifier_from_read_header(record.description)
                seq = record.seq
                fixed_seq = umi + seq
                fixed_qual = [40]*len(umi)+record.letter_annotations["phred_quality"]

                record.letter_annotations = {}
                record.seq = fixed_seq
                record.letter_annotations["phred_quality"] = fixed_qual


                SeqIO.write(record,self.output_handle, self.file_format)

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


def batch_iterator(iterator, batch_size):
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
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)

        if batch:
            yield batch

























