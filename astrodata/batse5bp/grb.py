"""
Define a GRB object holding BATSE data associated with a single GRB trigger,
and a GRBCollection to hold GRB instances for many GRBs.

Created 2012-05-06 by Tom Loredo
(adapted from an earlier version for the 4B catalog)
2019:  Converted to Python 3
"""

from collections import OrderedDict
from os.path import join, exists
from os import mkdir
import urllib.request

from PIL import Image

from .locations import root, trigger_paths, trigger_url, raw_cache
from .locations import ascii64_paths, ascii64_url
from .drm import DRMs_DISCSC
from .ascii64 import ASCII64


class GRB:
    """
    Provide access to BATSE data for a single GRB, both through persistant
    local storage and via access to the CGRO SSC data archive.
    """

    # In order of columns (after trigger column) in the flux table:
    flux_attributes = [
        'F1', 'F1_err', 'F2', 'F2_err', 'F3', 'F3_err', 'F4', 'F4_err',  # 0-7
        'F64ms', 'F64ms_err', 't64ms',        # 8-10
        'F256ms', 'F256ms_err', 't256ms',     # 11-13
        'F1024ms', 'F1024ms_err', 't1024ms',  # 14-16
    ]

    # In order of columns (after trigger column) in the duration table:
    duration_attributes = [
        'T50', 'T50_err', 'T50_start',  # 0-2
        'T90', 'T90_err', 'T90_start',  # 3-5
    ]

    def __init__(self, basic_record):
        """Define a burst from its record in the basic table."""
        # Break record into the trigger/name, Ulysses flag, and rest,
        # since there may be no space between the name and Uflag.
        # (The Ulysses flag is only in the published 4B catalog.)
        start = basic_record[:18].strip().split()
        # Uflag = basic_record[16]
        cols = basic_record[18:].strip().split()

        self.trigger = int(start[0])
        self.name = start[1] + '_' + start[2]
        self.desig = start[2]
        self.in_4B = start[1] == '4B'

        self.TJDate = int(cols[0])  # integer date
        self.secs = float(cols[1])  # seconds on that day
        self.TJD = self.TJDate + self.secs/86400.
        # All angles are in degrees.
        self.RA = float(cols[2])
        self.dec = float(cols[3])
        self.long = float(cols[4])
        self.lat = float(cols[5])
        self.drxn_err = float(cols[6])  # error circle radius
        self.geo = float(cols[7])  # angle from nadir
        if cols[8] == 'Y':  # overwrote previous trigger
            self.overwrote = True
        else:
            self.overwrote = False
        if cols[9] == 'Y':  # overwritten by later trigger
            self.incomplete = True
        else:
            self.incomplete = False
        self.has_flux = False  # bursts with gaps will have no flux data
        self.has_durn = False
        self.comments = []  # store (flag, comment) tuples

        # Store path elements for local and remote data access.
        self.group, self.dir, tail = trigger_paths(self.trigger)
        self.local_dir = join(self.group, self.dir)
        self.remote_dir = trigger_url + tail

        self.a64_group, self.a64_rfname, tail = ascii64_paths(self.trigger)
        self.a64_remote = ascii64_url + tail

    def set_bright(self, trigger, cols):
        """Add peak flux & fluence data from a brightness table record."""
        if trigger != self.trigger:
            raise ValueError('Flux data for wrong trigger!')
        for i, name in enumerate(self.flux_attributes):
            setattr(self, name, float(cols[i]))
        self.has_flux = True

    def set_durn(self, trigger, cols):
        """Add duration data from a duration table record."""
        if trigger != self.trigger:
            raise ValueError('Duration data for wrong trigger!')
        for i, name in enumerate(self.duration_attributes):
            setattr(self, name, float(cols[i]))
        self.has_durn = True

    # Replaced below.
    # def grb_dir(self):
    #     """
    #     Return the path to the directory in the local database holding data
    #     associated with this GRB, creating the directory if necessary.
    #     """
    #     group = join(root, self.group)
    #     if not exists(group):
    #         mkdir(group)
    #     dir_path = join(root, self.local_dir)
    #     if not exists(dir_path):
    #         mkdir(dir_path)
    #     return dir_path

    def set_grb_dir(self):
        """
        Set the path to this GRB's local data directory, creating the directory
        and its containing group if needed.

        This is not done in __init__ so the group and burst directories are
        created only when detailed data are loaded for a burst.
        """
        group = join(root, self.group)
        if not exists(group):
            mkdir(group)
        self.grb_dir = join(root, self.local_dir)
        if not exists(self.grb_dir):
            mkdir(self.grb_dir)

    def file_check(self, fname):
        """
        Check to see if the named file exists in this bursts data directory.
        """
        path = join(self.grb_dir, fname)
        return exists(path)

    def cached_path(self, fname):
        """
        Return the path to a persistant local copy of a remote file.

        Point to a local copy of the resource if present; otherwise, fetch it
        remotely and save a copy locally, then return the local path.
        """
        resource = join(self.grb_dir, fname)
        if not exists(resource):
            remote = self.remote_dir + fname
            urllib.request.urlretrieve(remote, resource)
        return resource

    def raw_cached_path(self, fname):
        """
        Return the path to a local copy of a remote file; it is stored in the
        raw_cache directory for temporary use.

        Point to an existing local copy of the resource if present; otherwise,
        fetch it remotely and save a copy locally, then return the local path.
        """
        resource = join(root, raw_cache, fname)
        if not exists(resource):
            remote = self.remote_dir + fname
            print('*** Copying:', remote)
            urllib.request.urlretrieve(remote, resource)
        return resource

    def open_cached(self, fname):
        raise NotImplementedError()

    def clear_cached(self, fname):
        raise NotImplementedError()

    def save_pickle(self, obj, fname):
        """
        Save an object associated with this GRB as a pickle in the database.
        """
        raise NotImplementedError()

    def load_pickle(self, fname):
        """
        Retrive an object associated with this GRB from a pickle in the
        database.
        """

    def show_lc(self):
        """
        Show light curve GIF files in the default image browser, merging 4ch
        and summed light curves.
        """
        ch_path = '%d_4ch.gif' % self.trigger
        sum_path = '%d_sum.gif' % self.trigger
        ch_path = self.cached_path(ch_path)
        sum_path = self.cached_path(sum_path)
        ch_im = Image.open(ch_path)
        # Some ch plots are higher than others, so note the height.
        w, h = ch_im.size
        sum_im = Image.open(sum_path)
        both = Image.new("RGB", (1610, max(600, h)))
        both.paste(ch_im, (0,0))
        both.paste(sum_im, (810,0))
        both.show()

        # This version takes less space, but looks horrible due to interpolation!
        # both = Image.new("RGB", (1210, 600))
        # both.paste(ch_im, (0,0))
        # sum_im.thumbnail((400,300), Image.BILINEAR)
        # both.paste(sum_im, (810,150))
        # both.show()

    def load_discsc_drms(self):
        """
        Load DISCSC DRM data for the triggered detectors.
        """
        self.discsc_drms = DRMs_DISCSC(self)

    def load_ascii64(self):
        """
        Load 64 ms count data from COSSC ASCII files.

        These files compile DISCLA, PREB, and DISCSC data over an interval
        determined by examination of trigger and pre-trigger data by the
        BATSE team.

        NOTE:  The DISCLA data is in 1024 ms bins; the counts are divided by
        16 to distribute across the tabulated 64 ms bins.  The DISCLA bins
        can overlap with the 1st PREB bin, in which case the PREB overlap
        is subtracted and the remaining counts divided among earlier bins.

        For full documentation, see:

        ftp://legacy.gsfc.nasa.gov/compton/data/batse/ascii_data/64ms/README
        """
        self.ascii64 = ASCII64(self)

    def __getattr__(self, name):
        """
        Implement auto-loading of some data by catching attributes not yet
        assigned in the instance dict.
        """
        if name == 'grb_dir':
            self.set_grb_dir()
            return self.grb_dir
        elif name == 'discsc_drms':
            self.load_discsc_drms()
            return self.discsc_drms
        elif name == 'ascii64':
            self.load_ascii64()
            return self.ascii64
        else:
            raise AttributeError('No attribute "%s"!' % name)

    def __str__(self):
        s = 'Trigger:  %i\n' % self.trigger
        s += 'Name:  %s\n' % self.name
        s += 'TJD:   %f\n' % self.TJD
        s += '(RA, Dec) (long, lat) (deg):    (%6.2f  %6.2f)  (%6.2f  %6.2f)\n' % \
             (self.RA, self.dec, self.lat, self.long)
        s += 'Drxn err (deg):   %5.2f\n' % self.drxn_err
        s += 'Complete?  %s' % (not self.incomplete)
        if self.has_flux:
            s += '\nCh 1-4 fluence (ergs):  %7.1e  %7.1e  %7.1e  %7.1e\n' % \
                 (self.F1, self.F2, self.F3, self.F4)
            s += 'Ch 1-4 S/N:             %7.1e  %7.1e  %7.1e  %7.1e\n' % \
                 (self.F1/self.F1_err, self.F2/self.F2_err, self.F3/self.F3_err, self.F4/self.F4_err)
            s += '64ms F_pk (cts/cm^2/s), S/N, t:  %6.2f  %6.2f  %7.2f' % \
                 (self.F64ms, self.F64ms/self.F64ms_err, self.t64ms)
        if self.has_durn:
            s += '\nT50, T90:  %6.2f  %6.2f' % (self.T50, self.T90)
        if self.comments:
            for entry in self.comments:
                s += '\n%s:  %s' % entry
        return s


class GRBCollection(OrderedDict):
    """
    Store GRB objects in an ordered dict accessed by BATSE trigger number
    (with the key being the integer trigger number, not a string).

    GRB data may also be accessed as attributes in two ways:

      .t#    where # is the trigger number; returns a GRB instance

      .b#    where # is YYMMDD from the traditional GRBYY... designation;
             returns a list of GRB instances

    Note that multiple bursts may occur on the same day, so the .b# value is a
    list of GRBs that occurred on that day.  The 4B catalog distinguished
    such bursts with designators that add a letter suffix (B, C, D after
    the first burst), nominally ordered by intensity.  The 5Bp bursts (even
    those from 4B) are not so labeled.
    """

    def __init__(self, *args, **kwds):
        # Note:  "*args, **kwds" are needed to support re-creating the instance
        # after pickling, via OrderedDict's pickling interface.
        OrderedDict.__init__(self, *args, **kwds)
        # These dicts map from 't' and 'b' attributes to trigger numbers.
        self.t_attrs = {}
        self.b_attrs = {}

    def add(self, grb):
        """
        Add a GRB instance to the collection.
        """
        self[grb.trigger] = grb
        t_attr = 't%i' % grb.trigger
        b_attr = 'b%s' % grb.desig[0:6]
        self.t_attrs[t_attr] = grb.trigger
        if not (b_attr in self.b_attrs):
            self.b_attrs[b_attr] = [grb]
        else:
            self.b_attrs[b_attr].append(grb)

    def __getattr__(self, name):
        """
        Return values associated with .t# or .b# attributes.
        """
        if name[0] == 't':
            tnum = self.t_attrs[name]
            return self[tnum]
        elif name[0] == 'b':
            return [grb for grb in self.b_attrs[name]]
        else:
            raise AttributeError('Illegal attribute!')

# These are experiments to support unpickling; the solution was to provide
# *args and **kwds arguments to the superclass initializer.

#     def __reduce__(self):
#         """
#         Method used to support pickling of an OrderedDict subclass.
#         """
#         t = OrderedDict.__reduce__(self)
#         return (t[0], ()) + t[2:]

#     def __reduce__(self):
#         """Return state information for pickling."""
#         items = [[k, self[k]] for k in self]
#         inst_dict = vars(self).copy()
#         for k in vars(OrderedDict()):
#             inst_dict.pop(k, None)
#         if inst_dict:
#             return (self.__class__, (items,), inst_dict)
#         return self.__class__, (items,)

#     def __reduce__(self):
#         return self.__class__, (OrderedDict(self),)
