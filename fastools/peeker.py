"""
http://stackoverflow.com/questions/14283025/python-3-reading-bytes-from-stdin-pipe-with-readahead
"""
from io import BytesIO
from os import SEEK_END


class Peeker(object):
    def __init__(self, handle):
        self._buf = BytesIO()
        self._handle = handle

        self.name = handle.name

    def _append_to_buf(self, data):
        position = self._buf.tell()
        self._buf.seek(0, os.SEEK_END)
        self._buf.write(data)
        self._buf.seek(position)

    def peek(self, size):
        data = self._handle.read(size)
        self._append_to_buf(data)
        return data

    def read(self, size=None):
        if size is None:
            return self._buf.read() + self._handle.read()
        data = self._buf.read(size)
        if len(data) < size:
            data += self._handle.read(size - len(data))
        return data

    def readline(self):
        line = self._buf.readline()
        if not line.endswith('\n'):
            line += self.handle.readline()
        return line
