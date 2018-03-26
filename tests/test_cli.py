"""Tests for the fastools CLI."""
import StringIO

from Bio import SeqIO

from fastools import fastools

from shared import md5_check

class TestCLI(object):
    def setup(self):
        self._sanitised_fa = open('data/sanitised.fa')
        self._sanitised_fq = open('data/sanitised.fq')
        self._length_fa = open('data/length.fa')
        self._output = StringIO.StringIO()
        self._null = open('/dev/null', 'w')

    def _md5_check(self, md5sum):
        return md5_check(self._output.getvalue(), md5sum)

    def test_sanitise(self):
        fastools.sanitise(open('data/unsanitised.fa'), self._output)
        assert self._md5_check('95c6eb5c9ae6949bbc4b4153cf976de4')

    def test_fa2fq(self):
        fastools.fa2fq(self._sanitised_fa, self._output, 30)
        assert self._md5_check('75d8abea91a1ec4d24c0e10235a28a7f')

    def test_fq2fa(self):
        fastools.fq2fa(self._sanitised_fq, self._output)
        assert self._md5_check('33d5334b5e210f681f5a23b9c03806be')

    def test_add(self):
        fastools.add(self._sanitised_fq, self._output, 'ACGT', 31)
        assert self._md5_check('a9931faafa17c50b7f9737dac07ca95d')

    def test_aln(self):
        assert fastools.aln(
            [self._sanitised_fa, open('data/sanitised_ed.fa')])[0][2] == 2

    def test_maln(self):
        pass # FIXME.

    def test_length(self):
        assert fastools.length(self._length_fa) == [1, 2, 3]

    def test_list_enzymes(self):
        assert 'EcoRI' in fastools.list_enzymes()

    def test_restrict(self):
        assert fastools.restrict(
            self._sanitised_fa, ['BssMI', 'AgsI']) == [11, 16, 74]

    def test_collapse_fasta(self):
        fastools.collapse_fasta(self._sanitised_fa, self._output, 2)
        assert self._md5_check('d2b7778038d8b61beeca74f0a3306e35')

    def test_s2i(self):
        fastools.s2i(self._sanitised_fq, self._output)
        assert self._md5_check('aa58abde3d3042a1eee9ac5f4064f29b')

    def test_count_tags(self):
        assert fastools.count_tags(self._length_fa, 'ACG', 1) == 2

    def test_select(self):
        fastools.select(self._sanitised_fa, self._output, 5, 10)
        assert self._md5_check('8a3a5bfd1de1c054e6274aa6cfcf93b0')

    def test_rselect(self):
        fastools.rselect(self._sanitised_fa, self._output, 'sanitised', 5, 10)
        assert self._md5_check('8a3a5bfd1de1c054e6274aa6cfcf93b0')

    def test_fa2gb(self):
        fastools.fa2gb(self._sanitised_fa, self._output, 'NM_000000.0')
        assert self._md5_check('a58ca021a538f76737fe641209451f09')

    def test_gb2fa(self):
        fastools.gb2fa(open('data/sanitised.gb'), self._output)
        assert self._md5_check('3c44cba3e9269ca729a6cd6292e4c05e')

    def test_mangle(self):
        fastools.mangle(self._sanitised_fa, self._output)
        assert self._md5_check('b427c19cf161cf23256cd76a044977d0')

    def test_generate_dna(self):
        fastools.generate_dna(10, self._output, 'name', 'description')
        self._output.seek(0)
        record = SeqIO.parse(self._output, 'fasta').next()
        assert (
            len(record.seq) == 10 and record.description == 'name description')

    def test_get_reference(self):
        pass # FIXME.

    def test_cat(self):
        assert fastools.cat(self._sanitised_fa).next().startswith('TGAGCGGAAC')

    def test_raw2fa(self):
        fastools.raw2fa(
            StringIO.StringIO('ACGT'), self._output, 'name', 'description')
        assert self._md5_check('6361fecba38154e9f3563d13c521154d')

    def test_descr(self):
        assert fastools.descr(self._sanitised_fa).next() == 'sanitised'

    def test_seq_split_1(self):
        fastools.seq_split(
            self._length_fa, [self._output, self._null], 'C')
        assert self._md5_check('df04208746b395c26bee2589798ed084')

    def test_seq_split_2(self):
        fastools.seq_split(
            self._length_fa, [self._null, self._output], 'C')
        assert self._md5_check('329a7a950ee7f48e65f784900158d0f8')

    def test_length_split_1(self):
        fastools.length_split(
            self._length_fa, [self._output, self._null], 2)
        assert self._md5_check('df04208746b395c26bee2589798ed084')

    def test_length_split_2(self):
        fastools.length_split(
            self._length_fa, [self._null, self._output], 2)
        assert self._md5_check('329a7a950ee7f48e65f784900158d0f8')

    def test_reverse(self):
        fastools.reverse(self._sanitised_fa, self._output)
        assert self._md5_check('89a3e1c61aecfff9f15a6702661d7170')

    def test_merge(self):
        fastools.merge(
            [self._sanitised_fa, open('data/sanitised.fa')], self._output, 3)
        assert self._md5_check('075b6720e3bd94f777fbb2b7ffa25ada')

    def test_fa_motif2bed(self):
        fastools.fa_motif2bed(self._sanitised_fa, self._output, 'AC')
        assert self._md5_check('d2b0dec731be890350bca49357d753f4')

    def test_csv2fa2_1(self):
        fastools.csv2fa2(
            open('data/primers.csv'), [self._output, self._null], True)
        assert self._md5_check('58b4f8a3832fdbbe7c91a8547c39c472')

    def test_csv2fa2_2(self):
        fastools.csv2fa2(
            open('data/primers.csv'), [self._null, self._output], True)
        assert self._md5_check('33259a9357b75e7baca2909bfc254e61')

    def test_edit(self):
        fastools.edit(self._sanitised_fa, open('data/edits.fa'), self._output)
        assert self._md5_check('b865c2069b8900df35d7733abd8c39e0')

    def test_dna2rna(self):
        fastools.dna2rna(self._sanitised_fa, self._output)
        assert self._md5_check('df76f76e77d5d5f785837c91a73f1a72')

    def test_rna2dna(self):
        fastools.rna2dna(open('data/sanitised_rna.fa'), self._output)
        assert self._md5_check('33d5334b5e210f681f5a23b9c03806be')
