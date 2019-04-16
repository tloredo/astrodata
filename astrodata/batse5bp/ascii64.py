"""
Objects providing access to 64 ms count data as compiled in the ASCII data
files at the COSSC.

These files compile DISCLA, PREB, and DISCSC data over an interval determined by
examination of trigger and pre-trigger data by the BATSE team.

NOTE:  The DISCLA data is in 1024 ms bins; the counts are divided by 16 to
distribute across the tabulated 64 ms bins.  The DISCLA bins can overlap with
the 1st PREB bin, in which case the PREB overlap is subtracted and the remaining
counts divided among earlier bins.

For full documentation, see:

ftp://legacy.gsfc.nasa.gov/compton/data/batse/ascii_data/64ms/README

That text is provided here as the `docn` string.

Created 2012-10-22 by Tom Loredo
2019:  Converted to Python 3
"""

from os.path import join, exists
import pickle
import urllib.request

import numpy as np
from numpy import loadtxt, empty, linspace, array

from .locations import root, raw_cache

# This docstring contains the text of the ascii_data README file at
# the CGRO SSC:
# ftp://legacy.gsfc.nasa.gov/compton/data/batse/ascii_data//64ms//README

docn = """Concatenated 64-ms Burst Data in ASCII Format

This data type is a concatenation of three standard BATSE data types, DISCLA,
PREB, and DISCSC.  All three data types are derived from the on-board data
stream of BATSE's eight Large Area Detectors (LADs), and all three data types
have four energy channels, with approximate channel boundaries:  25-55 keV,
55-110 keV, 110-320 keV, and >320 keV.

DISCLA data is a continuous stream with 1.024-s resolution, and is independent
of burst occurrence.  PREB data covers the interval 2.048 s just prior to a
burst trigger.  DISCSC nominally covers an interval of ~4 minutes starting
just after  trigger (early in the GRO mission, up to trigger #2099), or ~10-11
minutes (trigger #2101 and later).  The adjustment to longer time was made to
accommodate solar flares.  Nearly all cosmic gamma-ray bursts are contained
within the earlier  DISCSC allotment of ~4 minutes.  PREB and DISCSC have
64-ms resolution.

To realize a uniformly binned data stream that covers a sufficient interval
prior to a burst (for purpose of fitting  background), we constructed the
concatenated data type for all cosmic bursts where DISCLA, PREB, and DISCSC
exist, and where DISCLA data overlaps PREB (no gap allowed at that juncture),
as follows:

    The counts in each 1.024-s DISCLA sample were divided equally into 16
    64-ms   samples, with the last DISCLA sample, which can overlap the PREB
    data, being   truncated.  The truncation was accomplished by subtracting
    the counts in that   portion of the PREB data which overlaps the last
    DISCLA sample, then equally  dividing the remaining DISCLA counts into 16
    samples or fewer (once in ~16 occurrences there is no overlap with
    PREB).  Whereas DISCLA and PREB exist for all eight LADs, DISCSC (as
    transmitted to the ground) is summed only over  those detectors which
    triggered on a burst.  Therefore, only counts from the triggered
    detectors for the DISCLA and PREB types are summed and concatenated  with
    DISCSC.  Several count-summation and timing checks were made per burst to
    determine that the correct temporal correspondence was realized.

In the earlier part of the mission when only ~4 minutes of DISCSC data was
available, additional DISCLA data was concatenated after the end of DISCSC
data to attempt to cover the whole burst and provide some background interval,
pre and post burst.  In any case, the maximum number of 64-ms samples
concatenated is 8187 (524 s), save for nine exceptions for very long bursts
(triggers 148, 2287, 3448, 3567, 3639, 3930, 5446, 5693, 6125) which can use
up to 16379 samples.

Infrequently, gaps occur in some portion of the data stream; gaps are filled
with zeroes in all four channels.

The concatenated data type files are written in ASCII format for ease of
interpretation, with a 2-line header:

trig#  npts   nlasc  1preb    followed by 4-chan count rates (64-ms bins)
  105  2496    1856      0

The first line (always the same) describes four important parameters, example
values of which are listed in the second line (ASCII format = 4I6):

    1st parm (trig#): unique BATSE trigger number
    2nd parm (npts): total number of samples to follow, per energy channel
    3rd parm (nlasc): total number of DISCLA 64-ms samples concatenated prior
        to PREB
    4th parm (1preb): first PREB sample number after last 1.024-s DISCLA
        sample

Thus the file for burst #105 has 2496 samples in each of four energy channels,
1856 of these samples were concatenated  prior to PREB, and the last 1.024-s
DISCLA sample used did not overlap the first 64-ms PREB sample (1preb = 0).

Following the 2-line header for burst #105 are 2496 lines with 4 values per
line corresponding to counts per 64-ms sample in channels 1, 2, 3 and 4
(ASCII format = 4E16.8).
"""


# TODO:  Should probably fetch and cache the compressed group directories
# to save bandwidth and speed access of large parts of the catalog
# (perhaps optionally).


class ASCII64:
    """
    Provide access to the 64 ms count data compiled in the ASCII data files
    hosted by the COSSC.
    """

    # Attributes to pickle.
    pkl_attrs = ['n_asc64', 'n_discla', 'olap', 'n_bins',
                 'n_early', 'n_mid', 'n_late', 't_early', 't_mid', 't_late',
                 'counts', 'times']

    def __init__(self, grb):
        """
        Load ASCII 64 ms count data for the GRB instance `grb`.

        If the data is not already locally cached, it is fetched from the SSC,
        loaded, and cached for future use.  Otherwise, the cached copy is
        loaded.
        """
        self.grb = grb
        # Note local pkl file name is remote file name + '.pkl'
        a64_path = join(grb.grb_dir, grb.a64_rfname+'.pkl')
        if exists(a64_path):
            # print 'Reading existing a64 pickle:', a64_path
            ifile = open(a64_path, 'rb')
            data = pickle.load(ifile)
            ifile.close()
            for attr in self.pkl_attrs:
                setattr(self, attr, data[attr])
        else:
            a64_raw = join(root, raw_cache, grb.a64_rfname)
            if not exists(a64_raw):
                # print 'Fetching a64 data from SSC...'
                urllib.request.urlretrieve(grb.a64_remote, a64_raw)
            dfile = open(a64_raw, 'r')
            dfile.readline()  # header
            meta = dfile.readline().strip().split()
            try:
                assert int(meta[0]) == grb.trigger
            except AssertionError:
                raise ValueError('Trigger \# mismatch for ASCII64 data!')
            self.n_asc64 = int(meta[1])
            self.n_discla = int(meta[2])
            self.olap = int(meta[3])
            # Data are indexed as [ch, time].
            raw = loadtxt(dfile, unpack=True)
            self.raw = raw  # useful for debugging

            # Break up the data into the (early) DISCLA bins, the (late) 64 ms
            # bins (PREB + DISCSC), and a possible (mid) truncated DISCLA bin.
            self.n_early = self.n_discla // 16
            early = empty((4, self.n_early), int)
            for ch in range(0,4):
                for i in range(self.n_early):
                    j = 16 * i
                    early[ch,i] = int(round(16*raw[ch,j]))
            self.early = early

            n_trail = self.n_discla - 16*self.n_early
            if n_trail == 0:
                self.n_mid = 0
                mid = None
            else:
                self.n_mid = 1
                mid = empty((4,1), int)  # must be 4x1 for concat with others
                j = self.n_early*16
                for ch in range(0,4):
                    mid[ch] = int(round(n_trail*raw[ch,j]))
            self.mid = mid

            self.n_late = self.n_asc64 - self.n_discla
            late = empty((4,self.n_late), int)
            for ch in range(4):
                for j in range(self.n_late):
                    late[ch,j] = int(raw[ch,self.n_discla+j])
            self.late = late

            self.n_bins = self.n_early + self.n_mid + self.n_late

            # Create vectors of bin time interval boundaries.
            self.t_early = linspace(0, 1.024*self.n_early, self.n_early+1)
            last = self.t_early[-1]
            if self.n_mid == 0:
                self.t_mid = array([last, last])
            else:
                self.t_mid = array([last, last+n_trail*0.064])
                last = self.t_mid[-1]
            self.t_late = linspace(last, last + self.n_late*0.064, self.n_late+1)

            # Shift the times to be relative to the trigger time, which is
            # 2.048 s after the start of the PREB data (the first 64 ms data).
            t_trig = self.t_late[0] + 2.048
            self.t_early = self.t_early - t_trig
            self.t_mid = self.t_mid - t_trig
            self.t_late = self.t_late - t_trig

            # Create a concatenated set of counts and times.
            if self.n_mid == 0:
                self.counts = np.concatenate((self.early, self.late), axis=1)
                self.times = np.append(self.t_early, self.t_late[1:])
            else:
                self.counts = np.concatenate((self.early, self.mid, self.late),
                                             axis=1)
                self.times = np.append(self.t_early, self.t_late)

            # Archive this instance as a pickle.
            data = {}
            for attr in self.pkl_attrs:
                data[attr] = getattr(self, attr)
            ofile = open(a64_path, 'wb')
            pickle.dump(data, ofile, -1)  # use highest protocol
            ofile.close()

    def rates(self):
        """
        Return time bin centers & widths (vectors) & count rates (4-row matrix).
        These values are useful for plotting.
        """
        # TODO:  Add options to report root-cts or bin-wise HPD region error
        # estimates.
        centers = 0.5*(self.times[:-1] + self.times[1:])
        widths = self.times[1:] - self.times[:-1]
        rates = self.counts / widths
        return centers, widths, rates
