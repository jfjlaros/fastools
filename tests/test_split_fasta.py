"""
Tests for demultiplex.
"""
from fastools import split_fasta

from shared import FakeOpen, md5_check


class TestSplitFasta(object):
    def setup(self):
        self._fake_open = FakeOpen()
        split_fasta.open = self._fake_open.open

    def _md5_check(self, fileno, md5sum):
        return md5check(self._fake_open.handles[fileno].getvalue(), md5sum)
