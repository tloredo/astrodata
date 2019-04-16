"""
Read data from the "current" BATSE catalog (dubbed 5Bp here, with "p" for
"preliminary," since an official 5B successor to the 4B catalog has not yet
been released).  Provide access to catalog data and other GRB data via
a GRBCollection instance providing access to its individual GRB elements
in three ways:

* as an OrderedDict indexed by BATSE trigger number
* via attributes of the form .t#, with # = BATSE trigger number
* via attributes of the form .b#, with # = YYMMDD burst designation;
  this returns a list of triggers matching the designation (there will
  be >1 if BATSE detected multiple bursts on the specified date)

Only one function in this module is intended for users:  load().

This module was created to access the data as released in Jul-Sep 2000.

Created 2012-05-06 by Tom Loredo
2019:  Converted to Python 3
"""

from os.path import abspath, exists, join
from os import mkdir
import pickle

from .grb import GRB, GRBCollection
from .locations import *
from .utils import retrieve_gzip

__all__ = ['load_catalog']


# TODO:  get_grb_classes is presently unused; is there a use case?  May
# only be useful if the pickled files are unpickled outside this module.

def get_grb_classes(modname, classname):
    """
    Return class objects from the "grb" module that have the specified
    `classname`.

    This function is for identifying classes encountered when unpickling
    BATSE 5Bp data; it satisfies the cPickle "find_global" interface.
    """
    # print 'Module:', modname, ' -- Class:', classname
    if classname == 'GRB':
        return GRB
    elif classname == 'GRBCollection':
        return GRBCollection
    else:
        raise RuntimeError('Unrecognized class in pickled data: %s, %s' % (modname, classname))


def read_summaries():
    """
    Read GRB summary information from a pre-existing pickled data file.
    """
    try:
        sfile = open(join(root, summaries), 'rb')
    except FileNotFoundError:
        raise IOError('Summary data file does not exist!')

    # Define an unpickler that will recognize grb classes even when unpickling
    # is done elsewhere elsewhere (in which case grb classes may not be on the
    # top level, thwarting normal unpickling).
    loader = pickle.Unpickler(sfile)
    # loader.find_global = get_grb_classes
    GRBs = loader.load()
    sfile.close()
    print('Loaded summary data for', len(GRBs), 'GRBs comprising the 5Bp catalog.')
    return GRBs


def get_grb_bright(bfile):
    """
    Read the brightness data for a single GRB from the brightness data file
    `bfile`; return the trigger number and a list of data entries (strings).
    """
    line = bfile.readline()
    if line == '':
        return None, None
    data = []
    # 1:  trigger, ch1 fluence & err, ch2 fluence & err
    words = line.strip().split()
    trig = int(words[0])
    data.extend(words[1:])
    # 2:  ch3 fluence & err, ch4 fluence & err
    words = bfile.readline().strip().split()
    data.extend(words)
    # 3:  64ms peak flux, err, time
    words = bfile.readline().strip().split()
    data.extend(words)
    # 4:  256ms peak flux, err, time
    words = bfile.readline().strip().split()
    data.extend(words)
    # 5:  1024ms peak flux, err, time
    words = bfile.readline().strip().split()
    data.extend(words)
    return trig, data


def fetch_summaries():
    """
    Fetch GRB summary information from the CGRO SSC; return it in a
    GRBCollection instance.
    """

    # Get access to the raw data files, either cached or fetched from the SSC.
    cache = join(root, raw_cache)
    basic = retrieve_gzip(basic_url, cache)
    bright4 = retrieve_gzip(bright_url4, cache)
    bright5 = retrieve_gzip(bright_url5, cache)
    durn = retrieve_gzip(durn_url, cache)
    comments = retrieve_gzip(comments_url4, cache)

    # Read basic data, defining the GRB objects.  Add the trigger data path.
    GRBs = GRBCollection()
    ncomp = 0  # count complete GRBs (not overwritten by subsequent GRB)
    for line in basic:
        if not line:  # in case of empty lines at end
            break
        grb = GRB(line)
        if grb.trigger in GRBs:
            raise ValueError('Duplicate entries for trigger %i !' % grb.trigger)
        GRBs.add(grb)
        if not grb.incomplete:
            ncomp += 1
    basic.close()
    print('Read data for', len(GRBs), 'triggers from basic table,', ncomp,
          'complete...')
    print()

    # Add brightness (flux, fluence) data.
    nf = 0
    extra = []  # collect triggers in flux table but not basic table
    while True:
        trigger, data = get_grb_bright(bright4)
        if trigger is None:
            break
        if trigger in GRBs:
            GRBs[trigger].set_bright(trigger, data)
            nf += 1
        else:
            extra.append(trigger)
    bright4.close()
    while True:
        trigger, data = get_grb_bright(bright5)
        if trigger is None:
            break
        if trigger in GRBs:
            GRBs[trigger].set_bright(trigger, data)
            nf += 1
        else:
            extra.append(trigger)
    bright5.close()
    print('Read flux data for', nf, 'basic table triggers.')
    print('Extraneous flux data for:', extra)
    if extra:
        print('***** Data for these GRBs were ignored!!! *****')
    print()

    # Add duration data.
    ndur = 0
    extra = []
    for line in durn:
        if not line:
            break
        data = line.strip().split()
        trigger = int(data[0])
        if trigger in GRBs:
            GRBs[trigger].set_durn(trigger, data[1:])
            ndur += 1
        else:  #
            extra.append(trigger)
    durn.close()
    print('Read duration data for', ndur, 'basic table triggers.')
    print('Extraneous data for:', extra)
    if extra:
        print('***** Data for these GRBs were ignored!!! *****')
    print()

    # Add comments.
    ncom = 0
    extra = []
    for line in comments:
        if not line:
            break
        if line[0] == '#':  # header
            continue
        trigger = int(line[:6].strip())
        flag = line[11]
        com = line[14:].strip()
        if trigger in GRBs:
            GRBs[trigger].comments.append((flag, com))
            ncom += 1
        else:  #
            extra.append(trigger)
    durn.close()
    print('Read comment data for', ncom, 'basic table triggers.')
    print('Extraneous data for:', extra)
    if extra:
        print('***** Data for these GRBs were ignored!!! *****')

    return GRBs


def load_catalog(root_dir=root):
    """
    Establish access to GRB data from the BATSE '5B' catalog, stored in the
    `root_dir` directory.  Return a GRBCollection providing burst-by-burst
    access keyed by trigger number and via trigger and YYMMDD (date)
    attributes.

    If no catalog has yet been established, the directory is created and
    summary data for all GRBs are fetched from the CGRO SSC and stored
    locally for future use.

    Detailed data for specific bursts is fetched, parsed, and cached
    lazily as requested.
    """
    # TODO:  Probably a better way to handle this than with a global;
    # maybe via a module as a singleton....
    global root

    # Make sure root directory exists.
    root = abspath(root_dir)  # assigns full path throughout package
    rc_dir = join(root, raw_cache)
    if not exists(root):
        mkdir(root)
    if not exists(rc_dir):
        mkdir(rc_dir)

    try:
        GRBs = read_summaries()
    except IOError:
        GRBs = fetch_summaries()
        sfile = open(join(root,summaries), 'wb')
        pickle.dump(GRBs, sfile)
        sfile.close()

    return GRBs
