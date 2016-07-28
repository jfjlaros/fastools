"""
Tests for demultiplex.
"""
import StringIO

from fastools import split_fasta

from shared import FakeOpen, md5_check


class TestSplitFasta(object):
    def setup(self):
        fake_open = FakeOpen()
        self._handles = fake_open.handles
        split_fasta.open = fake_open.open

        self._input = open('data/split_fasta.fa')
        self._library = open('data/library.csv')
        self._output = fake_open.open('stdout')

    def _md5_check(self, name, md5sum):
        return md5_check(self._handles[name].getvalue(), md5sum)

    def test_split(self):
        split_fasta.split_fasta(self._input, self._library, self._output)
        assert len(self._handles) == 7
        assert self._md5_check('stdout', 'e1c685ef32bc0e5eff44b4471d428c62')
        assert self._md5_check(
            'one_counted.txt', '32989fa6c7577611c81d18479712589d')
        assert self._md5_check(
            'Unrecognised.txt', '7950997464325678f3f8c1f87d6511ec')
