from importlib.metadata import PackageNotFoundError, metadata
from re import split
from typing import Callable

from .fastools import *
from .peeker import Peeker
from .utils import guess_file_format, guess_header_format


def _extract(key: str, delim: str = r'[^\s\S]', index: int = 0) -> str:
    try:
        value = metadata(__package__).get(key, '')
    except PackageNotFoundError:
        return '<NO DATA>'
    return split(delim, value)[index]


def doc_split(func: Callable) -> str:
    return func.__doc__.split('\n\n')[0]


_project = _extract('Name')
_version = _extract('Version')
_year = 2013-2026
_author = _extract('Author-email', r'"', 1)
_email = _extract('Author-email', r'<|>', 1)
_description = _extract('Summary')
_copyright = f'Copyright (c) {_year} by {_author} <{_email}>'
_url = _extract('Project-URL', r', ', 1)
_info = f'{_project} version {_version}\n\n{_copyright}\nHomepage: {_url}'


#_copyright_notice = 'Copyright (c) {} <{}>'.format(
#    _get_metadata('Author'), _get_metadata('Author-email'))
#
#usage = [_get_metadata('Summary'), _copyright_notice]
#
#
#def version(name):
#    return '{} version {}\n\n{}\nHomepage: {}'.format(
#        _get_metadata('Name'), _get_metadata('Version'), _copyright_notice,
#        _get_metadata('Home-page'))
