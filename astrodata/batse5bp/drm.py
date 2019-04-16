"""
Provide access to DRM data stored in FITS files at COSSC.

Documentation of the compressed DRM format from the DRM FITS header
(for 1-based array indexing):

The detector response matrix is stored in a compressed format:        

Elements are listed in order of row by column (n_e_bin by n_e_chan), with
all elements equal to 0 at the top of the column ommited.

The row index of the first non-zero element for each column is given by
n_zeros.

Thus, for the ith column: insert n_zeros(i) - 1 0s, followed by the next
n_e_bin - n_zeros(i) + 1 elements of the drm (either drm_dir, drm_sct or
drm_sum).

The last column should exhaust the elements of the compressed drm.

Created 2012-10-16 by Tom Loredo
2019:  Converted to Python 3
"""

from os.path import join, exists
import pickle

from numpy import zeros, minimum, array, array_equal
# import pyfits
from astropy.io import fits

from .prodquad import ProdQuadRule, CompositeQuad


# This list of the binary table header fields was obtained from a DISCSC DRM
# file by examining the binary table header __dict__.  They are from the
# _cards attribute of the header.
hdr_fields = [
    ('TTYPE1', 'DET_NUM', 'the active detector for one row'),
    ('TFORM1', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE2', 'MAT_TYPE', '0=Direct,1=Scattered, 2= Both, 3=Summed'),
    ('TFORM2', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE3', 'TIME', 'time drm is calculated for'),
    ('TFORM3', '1D', 'data format of the field: DOUBLE PRECISION'),
    ('TUNIT3', 'TJD', 'physical unit field'),
    ('TTYPE4', 'NUMEBINS', 'number of rows in PHT_EDGE (N_E_BINS+1)'),
    ('TFORM4', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE5', 'NUMCHAN', 'number of columns in E_EDGES (N_E_CHAN+1)'),
    ('TFORM5', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE6', 'NUMZERO', 'number of of zeroes top of each column'),
    ('TFORM6', '1I', 'data format of the field: 2-byte INTEGER'),
    ('TTYPE7', 'DIRDRM', 'number of rows in the DIR drm'),
    ('TFORM7', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE8', 'SCTDRM', 'number of rows in the  SCT drm'),
    ('TFORM8', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE9', 'SUMDRM', 'number of rows in the  SUM drm'),
    ('TFORM9', '1J', 'data format of the field: 4-byte INTEGER'),
    ('TTYPE10', 'PHT_EDGE', 'input bin energy edges'),
    ('TFORM10', 'PE(276)', 'data format of the field'),
    ('TUNIT10', 'keV', 'physical unit field'),
    ('TTYPE11', 'E_EDGES', 'channel energy edges'),
    ('TFORM11', 'PE(253)', 'data format of the field'),
    ('TUNIT11', 'keV', 'physical unit field'),
    ('TTYPE12', 'N_ZEROS', 'list of how many 0 elements in each column'),
    ('TFORM12', 'PI(252)', 'data format of the field'),
    ('TTYPE13', 'DRM_DIR', 'detector response matrix'),
    ('TFORM13', 'PE(69552)', 'data format of the field'),
    ('TUNIT13', 'cm**2', 'physical unit field'),
    ('TTYPE14', 'DRM_SCT', 'detector response matrix'),
    ('TFORM14', 'PE(69552)', 'data format of the field'),
    ('TUNIT14', 'cm**2', 'physical unit field'),
    ('TTYPE15', 'DRM_SUM', 'detector response matrix'),
    ('TFORM15', 'PE(69552)', 'data format of the field'),
    ('TUNIT15', 'cm**2', 'physical unit field'),
    ('TTYPE16', 'CAL_NAME', 'Name of energy calibration method'),
    ('TFORM16', '16A', 'data format of the field')]


# Collect information for the TTYPE fields in a dict.
ttype_fields = {}
for attr, val, comment in hdr_fields:
    if attr[:5] == 'TTYPE':
        col = int(attr[5:]) - 1
        ttype_fields[col] = (val, comment)


class DRM_DISCSC:
    """
    Provide access to the detector response matrix (DRM) for a single BATSE
    detector's discriminator science (DISCSC) data.

    The DRM units according to the FITS header is area in cm^2.
    """

    fields = ttype_fields

    def __init__(self, row_data, n_ch, n_E):
        """
        Retrieve DRM info from a single row of parsed binary FITS table data.
        """
        # Copy FITS TTYPE data as attributes.
        for col, info in self.fields.items():
            # info[0] is the column name; info[1] is the FITS comment.
            setattr(self, info[0], row_data[col])

        self.det_num = self.DET_NUM

        # Discriminator bin boundaries (nominal energy loss units):
        self.ch_bins = self.E_EDGES  # n_ch+1 values

        # Incident photon energy "edges;" there should be n_E of these but
        # they are treated as if binned somehow....
        self.E_bins = self.PHT_EDGE

        # Unpack the zero-shortened format.
        self.start = self.N_ZEROS - 1  # N_ZEROS is 1-based index of 1st non-0 element
        self.drm = zeros((n_ch, n_E))
        run = 0
        for i in range(n_ch):
            num_nz = n_E - self.start[i]
            # Use the summed (direct + atmospheric scattering) matrix; typically
            # the direct and scattering entries are empty.
            self.drm[i,self.start[i]:] = self.DRM_SUM[run:run+num_nz]
            run += num_nz

        # Make views providing detector channel response vectors.
        self.crv = []
        for i in range(n_ch):
            self.crv.append(self.drm[i,:])


class DRMs_DISCSC:
    """
    Provide access to detector response matrix (DRM) data from a BATSE DISCSC DRM
    FITS file for all triggered detectors associated with a BATSE trigger; also
    calculate the summed response.
    """

    # Metadata attributes stored in the "-meta" file.
    meta_attrs = ['burst_id', 'file_id', 'n_ch', 'n_E', 'alpha', 'n_det']

    def __init__(self, grb):
        """
        Load DRM data for the DISCSC data associated with GRB instance `grb`.

        If the DRM data is not locally cached, it is fetched from the SSC,
        loaded, and cached for future use.  Otherwise, the cached copy is
        loaded.
        """
        self.grb = grb

        # Get DRMs for the triggered detectors.
        drm_path = 'discsc_drm_%i.pkl' % grb.trigger  # stores detector DRMs
        drm_path = join(self.grb.grb_dir, drm_path)
        meta_path = 'discsc_drm_%i-meta.pkl' % grb.trigger  # stores metadata
        meta_path = join(self.grb.grb_dir, meta_path)
        if exists(drm_path) and exists(meta_path):  # already parsed FITS
            ifile = open(meta_path, 'rb')
            metadata = pickle.load(ifile)
            ifile.close()
            for attr in self.meta_attrs:
                setattr(self, attr, metadata[attr])
            ifile = open(drm_path, 'rb')
            self.detectors = pickle.load(ifile)
            ifile.close()
            try:
                assert self.n_det == len(self.detectors)
            except AssertionError:
                raise RuntimeError('Data mismatch in archive DRM files!')
        else:  # must parse FITS
            fits_name = 'discsc_drm_%i.fits.gz' % grb.trigger
            fits_path = grb.raw_cached_path(fits_name)
            self.load_from_fits(fits_path, drm_path, meta_path)

        # The full-instrument DRM sums the triggered detector DRMs.
        self.drm = zeros((self.n_ch, self.n_E))
        for i, d in enumerate(self.detectors):
            self.drm += d.drm
            # Keep track of lowest non-zero entry for each summed channel.
            if i == 0:
                self.start = d.start[:]
            else:
                self.start = minimum(self.start, d.start)

        # Make views providing detector channel response vectors.
        self.crv = []
        for i in range(self.n_ch):
            self.crv.append(self.drm[i,:])

        # Get energy node info for sum DRM from one of the detectors.
        # Discriminator bin boundaries (nominal energy loss units), n_ch+1 vals:
        self.ch_bins = self.detectors[0].ch_bins

        # Incident photon energy "edges;" there should be n_E of these but
        # they are treated as if binned somehow....
        self.E_bins = self.detectors[0].E_bins

        # Do a quick check of consistency between detectors.
        # TODO:  Verify more thoroughly with array_equal?
        # TODO:  Do this only when reading the FITS file.
        for n, detector in enumerate(self.detectors[1:]):
            for i in [0,1]:  # just check 1st 2 bins
                if detector.E_bins[i] != self.E_bins[i]:
                    raise ValueError('E_bins mismatch between detectors for (det, bin) = (%i, %i)!' % (n,i))

    def load_from_fits(self, fits_path, drm_path, meta_path):
        """
        Read DRM data from a FITS file.  The path `fname` may be to a gzipped
        version.  Archive the data at the DRM file paths.
        """
        hdus = fits.open(fits_path)
        primary = hdus[0].header

        # Verify some key primary header elements.
        try:
            assert primary['filetype'] == 'BATSE_DRM'
            assert primary['det_mode'] == 'LAD'
        except:
            raise ValueError('Primary header is not as expected for DRM FITS file!')

        self.burst_id = primary['object']  # GRByymmdd designator
        self.file_id = primary['file-id']  # contains trigger #

        self.n_ch = primary['n_e_chan']  # number of discriminator channels
        # number of of incident photon energies (should *not* be bins):
        self.n_E = primary['n_e_bins']
        self.alpha = primary['alpha']  # weighting across input bins for direct matrix

        # Read the DRM binary table extension, close the file, and gather
        # some basic info (field names and sizes).
        exten = hdus[1].header
        try:
            assert exten['exttype'] == 'BATSEDRM'
        except:
            raise ValueError('Table header is not as expected for DRM FITS file!')

        # These are handy for interactive exploring/checking of FITS data:
        # self.primary = primary
        # self.exten = exten
        # self.data = hdus[1].data

        # Each row contains the DRM for one of the triggered detectors.
        self. n_det = len(hdus[1].data)  # number of triggered detectors
        self.detectors = []
        for row in hdus[1].data:
            self.detectors.append(DRM_DISCSC(row, self.n_ch, self.n_E))
        hdus.close()

        # Store the data in a format for easier reloading.
        metadata = {}
        for attr in self.meta_attrs:
            metadata[attr] = getattr(self, attr)
        ofile = open(meta_path, 'wb')
        pickle.dump(metadata, ofile)
        ofile.close()
        ofile = open(drm_path, 'wb')
        pickle.dump(self.detectors, ofile, -1)  # use highest protocol
        ofile.close()

    def set_sum_prodquad(self, m=1, n=2):
        """
        Set up a composite interpolatory product quadtrature rule for the
        detector-summed DRM.

        Currently only a composite (1,2) rule is available, with a separate
        rule for each DRM interval.
        """
        # ***Caution***  It's not clear what the energies should be; should
        # check this against RMFIT.
        # Use the centers of the intervals defined by the "edges" as the
        # energies associated with the DRM values.

        # TODO:  Do bookkeeping to keep track of nonzero DRM start entries.
        # Do this in a way so that the spectrum can be evaluated on one
        # grid for all channels; coordinate with chan_quad.

        # Build & store a composite rule for each channel.

        if m != 1 and n != 2:
            raise ValueError('Requested rule unavailable!')

        self.quad_rules = []
        for ch in range(self.n_ch):
            rules = []
            # ns = self.start[ch]  # see TODO above
            ns = 0
            E_vals = (self.E_bins[ns:-1] + self.E_bins[1:]) / 2.
            for i in range(len(E_vals)-1):
                l, u = E_vals[i], E_vals[i+1]
                fvals = self.drm[ch,i:i+2]
                # pq = ProdQuad12(l, u, l, u, fvals=fvals)  # use Gauss-Legendre g nodes
                pq = ProdQuadRule([l, u], 3, l, u, f=fvals)  # use Gauss-Legendre g nodes
                rules.append(pq.quad_object())
            self.quad_rules.append(CompositeQuad(*rules))

        # Pull out the nodes for easy access.
        self.quad_nodes = self.quad_rules[0].nodes.copy()
        # The nodes (but not wts) should be the same for all channels.
        for ch in range(1, self.n_ch):
            if not array_equal(self.quad_nodes, self.quad_rules[ch].nodes):
                raise ValueError('Quadrature nodes do not match across channels!')

        # Record the integration range; it will be somewhat larger than the
        # span of quad_nodes because an open rule is used.
        self.E_l = E_vals[0]
        self.E_u = E_vals[-1]

    def chan_quad(self, ch, spec, *args):
        """
        Return the quadrature of the product of `spec` with the response function
        for channel `ch` as defined by the DRM.

        `spec` may be either a function evaluating the incident spectrum,
        or an array of values of the spectrum pre-computed on the quadrature
        nodes.  If it is a function, it should have the signature
        g(E, arg1, ...), where the first argument is the incident photon
        energy, and any optional remaining arguments get set via *args.

        Since the channels share common nodes, it will be most efficient
        to pass an array of spectrum values, unless only one of the channels
        is of interest.
        """
        if callable(spec):
            svals = array([spec(E, *args) for E in self.quad_nodes])
        elif len(spec) == self.quad_nodes.shape[0]:
            svals = spec
        else:
            raise ValueError('Argument must be callable or array of values!')
        return self.quad_rules[ch].quad(svals)

    def quad_nodes_wts(self, l=None, u=None):
        """
        Return the nodes and weights for quadrature rules integrating DRMs
        times an incident spectrum.

        `l` and `u` specify the incident energy range for the quadrature.  If
        either value is None, the corresponding limit is set equal to the
        DRM limit.  If either limit is beyond the range of the DRM, it is
        truncated.

        A tuple is returned:  (`nodes`, `wts`) where `wts` is a list of
        weight vectors, `wts[c]` being the vector of weights for channel `c`.

        If `svals` holds incident spectrum values evaluated on the nodes,
        the expected count rate for channel `c` is:

            sum(svals*wts[c])
        """
        if l is not None and l <= self.E_bins[0]:
            l = None
        if u is not None and u >= self.E_bins[-1]:
            u = None

        # Full-range case:
        if l is None and u is None:
            wts = []
            for ch in range(self.n_ch):
                wts.append(self.quad_rules[ch].wts)
            return self.quad_nodes, wts

        # Otherwise we have to modify the rules to fit the range.
        wts = []
        for ch in range(self.n_ch):
            print('*** Ch:', ch)
            nodes, ch_wts = self.quad_rules[ch].range_nodes_wts(l, u)
            wts.append(ch_wts)
        return nodes, wts
