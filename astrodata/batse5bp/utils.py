"""
Utilities for the batse5bp package.

Created 2012-05-08 by Tom Loredo
2019:  Converted to Python 3
"""

from os.path import exists, join, split
import urllib.request
import gzip
from subprocess import check_call


def retrieve_gzip(url, cache):
    """
    Return an open file accessing data from the file at the http URL `url`
    that may be persistantly stored as a gzipped file named from the URL tail
    and stored in the directory `cache`.

    If the cached file exists, simply return an opened file object accessing
    it.

    If it does not exist, retrieve the file using the URL, cache a gzipped
    version, and return an opened file object accessing it.
    """
    head, tail = split(url[6:])
    path = join(cache, tail)
    gzpath = path + '.gz'
    if not exists(gzpath):
        print('Accessing', tail, 'at CGRO SSC...')
        name, hdrs = urllib.request.urlretrieve(url, path)
        if name != path:
            raise ValueError('URL target/name mismatch!')
        check_call(['gzip', path])
        return gzip.open(gzpath, 'rt')
    return gzip.open(gzpath, 'rt')


def write_seq(fob, seq, format, per_line, label=None):
    """
    Write the contents of the sequence object seq to the open file
    object fob, writing each element according to the passed format,
    with per_line elements per line.  label gives an optional label to
    start each line.
    """

    n = len(seq)
    full_lines, nlast = divmod(n, per_line)
    if label:
        full_format = label + per_line*format + '\n'
        last_format = label + nlast*format + '\n'
    else:
        full_format = per_line*format + '\n'
        last_format = nlast*format + '\n'
    start = 0
    for i in range(full_lines):
        s = full_format % tuple(seq[start:start+per_line])
        fob.write(s)
        start = start + per_line
    s = last_format % tuple(seq[start:])
    fob.write(s)
