"""
A class (TTE) and support functions for reading and writing BATSE TTE data
from tte_ibdb FITS files, loosely based on IDL code by Jay Norris.

NOTE:  Jay Norris has found that a few tte_ibdb files have corrupt
data.  There is no simple way to identify these automatically, but
they are obvious when one plots a histogram of the TTE data.  The
problem appears to be due to misidentified triggered detectors.
In addition, Adam Kruger has found more complicated corruption in
a few files; this is reported in ApJ v. 576, 532 (2002).

FTP users take note!  If you are downloading a FITS file by ftp from an
archive, make sure you download it in BINARY format.  Some ftp programs
(e.g. 'Fetch' on the Mac) work by default in an 'automatic' mode that
attempts to distinguish text files from binary files automatically.  FITS
files mix ASCII and binary data and can confuse these programs.  An
improperly downloaded file will appear to read correctly, but the data
will be corrupted.

Created:  25 May 2000 by Tom Loredo
Modified: 01 Jun 2000 TL
          14 Dec 2006 TL - Convert to PyFITS
          2012-05-09 TL - Adapted for batse5bp
          2019:  Converted to Python 3
"""

# import pyfits
from astropy.io import fits
from numpy import not_equal, array

from .utils import write_seq


# Start with utilities for unpacking packet data into standard types.

# A collection of 1-bit masks for the 16 least significant bits in a Python int.
masks = [0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100,
         0x0080, 0x0040, 0x0020, 0x0010, 0x0008, 0x0004, 0x0002, 0x0001]
masks.reverse()  # Reverse them so they are ordered lsb to msb.


# Pre-list comprehension versions ("_" suffix); not currently used:

def bytes2bits_(d):
    """Convert a 2-byte datum into a list of 16 bits ordered from
    lsb to msb.  Uses logical operations.  This approach is transparent
    to whether the datum has been read from a signed or unsigned
    representation."""

    bits = []
    for mask in masks:
        bits.append((d & mask) != 0)
    return bits


def byte2bits_(d):
    """Convert a 1-byte datum into a list of 8 bits ordered from
    lsb to msb."""

    bits = []
    for mask in masks[0:8]:
        bits.append((d & mask) != 0)
    return bits


def bytes2bits(d):
    """Convert a 2-byte datum into a list of 16 bits ordered from
    lsb to msb.  Uses logical operations.  This approach is transparent
    to whether the datum has been read from a signed or unsigned
    representation."""

    # Note that the comprehension returns a list of bools; use array
    # to cast it to ints.  The bools work without the cast since Pyton
    # treats True == 1, but just to play it safe....
    return array([(d & mask) != 0 for mask in masks], dtype=int)


def byte2bits(d):
    """Convert a 1-byte datum into a list of 8 bits ordered from
    lsb to msb."""

    return array([(d & mask) != 0 for mask in masks[0:8]], dtype=int)


def unpack_TTESD(d):
    """Unpack a TTE SCIENCE_DATA datum from its 2 byte packed format
    into (detector #, channel #'s, # of 2us clock cycles), returned
    as a tuple.  Since mutiple channels can be triggered in a 2us
    interval, the channel value is itself a tuple."""

    # Get the bits from the 2-byte datum.
    bits = bytes2bits(d)

    # The time is in the 9 lsb's.
    fac = 1
    t = 0
    for i in range(9):
        t = t + fac*bits[i]
        fac = 2*fac

    # The next 4 bits give the energy channel, with each bit assigned
    # to a channel.  Compile a list of channels triggered.
    # *** Some packets don't have a channel bit set; what does this mean?
    ch = []
    for i in range(1,5):
        if bits[i+8]:
            ch.append(i)

    # The last 3 bits give the detector number; this is 0-based to match
    # DSELB/DSELH.
    fac = 1
    det = 0
    for i in range(3):
        det = det + fac*bits[i+13]
        fac = 2*fac
    return (det, tuple(ch), t)


# Some test data from trigger 830.
if 0:
    print('Bits are LSB first...')
    print('127 in bits:', bytes2bits(127))
    print('-1 in bits: ', bytes2bits(-1))
    print('2**16-1 in bits: ', bytes2bits(65535))
    sd830 = [601, 18532, -7578, -15254, -31616, 2242, 17125, 757, 1271, -15600, 8976,
             9490, -31964, 17708, -31933, 20804, -3731, -24177,
             -22127, 9115, 18863, -11855, -3660, 17336, -11844,
             25537, 9689, -7192, -22037, -7186, -23055, 26107,
             -24575, 18434, 9226, -7142, -24019]
    print('\nsign, detector, (channels), tics, usec')
    for d in sd830:
        det,c,t = unpack_TTESD(d)
        s = '-'
        if d >= 0:
            s = '+'
        print(s, det, c, 2*t, 2*t+1024)


class TTE(object):
    """
    Instances of this class read and store TTE IDBD data from a BATSE
    tte_ibdb FITS file.  TTE.primary stores the primary header data.  Each
    IBDB field value that is constant (from packet to packet) is stored as
    a TTE attribute (e.g. TTE.BSTNO).  Data from 'good' events that were
    detected in a single energy channel in one of the 4 maximum-count-rate
    detectors are stored as follows:

        TTE.nevents # of events (int)
        TTE.time    integer # of us for each event (list of int)
        TTE.det     detector # for each event (list of int)
        TTE.chan    channel # for each event (list of int)

    Other events are either 'bad' events (multiple channels hit in 2us) or
    'background' events (detected with detectors other than the 'good'
    ones).  These are collected in lists of tuples:

        TTE.background  (t, det, ch) for each bg event
        TTE.bad     (t, det, ch) for each bad event

    In these tuples, ch is a tuple of the channels for that event.

    The initializer does a lot of internal consistency checking of the FITS
    file contents as it reads and parses the data; it takes ~10 s per file on
    a 2 GHz processor.  Store the instance for fast repeated use of the data.

    FWIW: I don't think multi-channel 'bad' events should be ignored, unless
    this is a discriminator error!
    """

    # TODO:  Reconsider 'bad' event handling?
    # TODO:  Can the rollover check be wrong (e.g., for low rates)?

    # These 14 fields should be the same for every packet (row) for a given
    # burst.  They will be made attributes of the TTE instance, with access
    # returning the (constant) field value associated with the burst.
    nconst = 14
    const_fields = ['SCHED_PTR', 'LOOP_CTR', 'ITER_NUM', 'BSTST', 'BSTLT',
                    'BSTNO', 'DSELB', 'NSCALE', 'BERFGS', 'TTRIG', 'SOLFLG',
                    'TIME', 'DSELH', 'SPARE']
    const_sizes = [1, 1, 1, 3, 3,
                   1, 1, 1, 1, 1, 1,
                   3, 4, 1]

    def __init__(self, fname):
        """Open the indicated TTE FITS file, read in the primary header
        and the IBDB TTE binary extension, and close the file.  Convert
        the binary data into a more usable form."""

        # Open the file and read and save the primary header.
        hdus = fits.open(fname)
        self.primary = hdus[0]

        # Read the IBDB binary extension, close the file, and gather
        # some basic info (field names and sizes).
        exten = hdus[1]
        table = exten.data
        hdus.close()
        if exten.name != 'BURST_AUX_IBDB':
            raise ValueError("Wrong ext'n name:  File is not a TTE IBDB FITS file!")
        tot_fields = exten.header['TFIELDS']
        if tot_fields != 19:
            raise ValueError('Not 19 fields; not valid IBDB file!')
        fields = []
        sizes = []
        self.nfields = 0
        # Read in the 19 field names and sizes.
        for n in range(tot_fields):
            field_name = exten.header['TTYPE%d' % (n+1)]
            key = 'TDIM%d' % (n+1)
            if key in exten.header:
                size = int(exten.header[key][1:-1])
            else:
                size = 1
            fields.append(field_name)
            sizes.append(size)
            if field_name not in TTE.const_fields:
                self.nfields = self.nfields + 1
            elif TTE.const_sizes[TTE.const_fields.index(field_name)] != size:
                raise ValueError('Field size does not match expectations!')

        # Make a list of rows, each row being a dictionary of the non-constant
        # IBDB data from the corresponding row in the BURST_AUX_IBDB binary
        # extension. Along the way, gather the values for the special fields
        # holding values that should be constant for a given burst.  Store
        # those constant values as attributes of this TTE instance.  Verify
        # they are indeed equal for every packet (row/record).
        self.nrows = len(table)
        self.rows = []
        self.time = []
        self.chan = []
        self.det = []
        self.nevents = 0
        self.background = []
        self.bad = []
        n512 = 0  # keeps track of event time rollovers beyond the 10th bit
        prev = None
        for i in range(self.nrows):
            if (i+1) % 32 == 0:
                print('Processing row ', i+1, '...')
            row = {}
            for n in range(tot_fields):
                field_name = fields[n]
                datum = table[i][field_name]

                # Handle constant packet info specially: save it the 1st time;
                # verify it the rest.
                if field_name in TTE.const_fields:
                    if i == 0:
                        setattr(self, field_name, datum)
                    else:
                        try:  # for non-array values
                            if datum != getattr(self, field_name):
                                err = '%s: %s %s' % (field_name, getattr(self,field_name), datum)
                                raise ValueError('Packet inconsistency for '+err)
                        except ValueError:  # compare arrays element-wise
                            if not_equal(datum, getattr(self, field_name)).all():
                                err = '%s: %s %s' % (field_name, getattr(self,field_name), datum)
                                raise ValueError('Packet inconsistency for '+err)

                # Save SCIENCE_DATA for special processing.
                elif field_name == 'SCIENCE_DATA':
                    event_list = datum

                # Otherwise, this is a non-constant field; store its value.
                else:
                    row[field_name] = datum
            self.rows.append(row)

            # First time thru, parse the DSELB byte into a list of
            # triggered detectors.
            if i == 0:
                bits = byte2bits(self.DSELB)
                self.tdlist = []
                for j in range(8):
                    if bits[j]:
                        self.tdlist.append(j)

            # Update event data from the SCIENCE_DATA.
            for event in event_list:
                det, ch, t = unpack_TTESD(event)
                # Watch for 9-bit rollover.
                if prev is not None:
                    if t < prev:
                        n512 = n512 + 1
                prev = t
                t = 2*(t + n512*512)  # presumes a 2us clock
                if len(ch) == 1:
                    if det not in self.tdlist:
                        self.background.append((t, det, ch))
                    else:
                        self.time.append(t)
                        self.chan.append(ch[0])
                        self.det.append(det)
                        self.nevents = self.nevents + 1
                else:
                    self.bad.append((t, det, ch))

    def __str__(self):
        """Return the constant fields, and a description of the data."""
        s = ''
        for field_name in TTE.const_fields:
            ss = '%-8s= %s\n' % (field_name, getattr(self,field_name))
            s = s + ss
        s = s + 'Triggered detectors: ' + str(self.tdlist) + '\n'
        s = s + 'Row data includes: ' + str(list(self.rows[0].keys())) + '\n'
        ss = 'Events (good, bg, bad): %d %d %d\n' % \
            (self.nevents, len(self.background), len(self.bad))
        s = s + ss
        return s

    def write_ASCII(self, fname=None):
        """Write the good event data in the ASCII TTE format adopted by
        the COSSC.  This method opens the indicated file for writing,
        and closes it on completion."""

        if fname is None:
            fname = 'TTE%05i.ASCII' % self.BSTNO
        ofile = open(fname, 'w')
        ofile.write('  BATSE burst trigger = %05i\n' % self.BSTNO)
        ofile.write('  npts =        %d\n' % self.nevents)
        ofile.write('  1st array:  "npts" times (unit = 1 microsec)\n')
        ofile.write('  2nd array:  associated energy chan #s (1,2,3 or 4)\n')
        ofile.write('  3rd array:  associated detector #s (0=>7)\n')
        write_seq(ofile, self.time, '%8i', 10)
        write_seq(ofile, self.chan, '%2i', 40)
        write_seq(ofile, self.det, '%2i', 40)
        ofile.close()

#-------------------------------------------------------------------------------


# A test case using trigger 830.
if 0:
    tte830 = TTE('tte_ibdb_830.fits')
    # Total # of events should be 2**15 (TTE buffer size):
    print(len(tte830.time)+len(tte830.background)+len(tte830.bad), 2**15)
    print(tte830.bad[0:4])
    print('DSELH: ', tte830.DSELH, tte830.det[-10:])
    print('times: ', tte830.time[0:10])
    print('chans: ', tte830.chan[0:10])
    print('dets:  ', tte830.det[0:10])
    print(tte830)
    tte830.write_ASCII('830.dat')
    print('Completed test!')

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    sys.exit()

    # Read in a trigger # from the command line, read in the associated
    # tte_ibdb FITS file, and write a TTE.ASCII file with data from the good
    # events.

    print(__name__)
    print('Processing file from command line...')
    try:
        tnum = int(sys.argv[1])
    except:
        raise ValueError('Give an integer trigger number to identify the FITS file!')
    fname = 'tte_ibdb_%i.fits' % tnum
    print('Getting TTE data from %s...' % fname)
    tte_tnum = TTE(fname)
    print(tte_tnum)
    tte_tnum.write_ASCII()
