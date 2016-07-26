"""
Tests for demultiplex.
"""
from fastools import split_fasta

from shared import FakeOpen, md5_check


class TestSplitFasta(object):
    def setup(self):
        fake_open = FakeOpen()
        self._handles = fake_open.handles
        split_fasta.open = self._fake_open.open

    def _md5_check(self, fileno, md5sum):
        return md5check(self._handles[fileno].getvalue(), md5sum)
