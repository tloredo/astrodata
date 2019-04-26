"""
Container classes for TESS data products.
"""

import string
import re
import datetime

import requests
import webbrowser
from lxml import html
import jdcal

import numpy as np
# import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import *

from astropy.io import fits


__all__ = ['TIDData']


sec2day = 1./(24*3600)  # for timedelta conversion


def datetime2mjd2k(dt):
    """
    Calculate the modified Julian date (MJD) for a datetime object, expressed
    in days since midnight starting 2000-01-01.

    Note that the conventional MJD is an integer a *noon*, not midnight.
    The calculated MJD2k is shifted by 0.5 to make the origin midnight.
    Subtract 0.5 to get the noon-based MJD.
    """
    # This is slightly complicated because jdcal only works for y/m/d.
    # Start with date, and delta = time b/t input and the date.
    date = datetime.datetime(dt.year, dt.month, dt.day)
    delta = dt - date  # a timedelta object
    # Convert delta to days.
    delta = delta.days + sec2day*(delta.seconds + 1.e-6*delta.microseconds)
    # MJD for date, in units of days since noon JD2000:
    mjd2k = jdcal.gcal2jd(date.year, date.month, date.day)[1] - jdcal.MJD_JD2000
    # Add fractional day, & shift to midnight boundary.
    return mjd2k + delta + 0.5


def read_target_metadata(ctnr, phdr):
    """
    Read target parameters from a primary FITS headers and store 
    it in the provided container instance.
    """
    # TID as an int, and TIC version:
    ctnr.ticid = phdr['TICID']
    ctnr.ticver = phdr['TICVER']

    # Target star info:
    ctnr.equinox = phdr['EQUINOX']
    ctnr.ra = phdr['RA_OBJ']  # deg
    ctnr.dec = phdr['DEC_OBJ']  # deg
    ctnr.pm_ra = phdr['PMRA']  # proper motion, mas/yr
    ctnr.pm_dec = phdr['PMDEC']
    ctnr.pm_tot = phdr['PMTOTAL']
    ctnr.mag = phdr['TESSMAG']
    ctnr.T_eff = phdr['TEFF']  # K
    ctnr.logg = phdr['LOGG']  # cm/s^2]
    ctnr.metal = phdr['MH']  # log10([M/H])
    ctnr.radius = phdr['radius']  # solar radii


def read_sector_metadata(ctnr, phdr, dhdr):
    """
    Read sector-level parameters from FITS headers (primary & data) and store 
    it in the provided container instance.
    """
    read_target_metadata(ctnr, phdr)

    # Sector and camera info:
    ctnr.sector = phdr['SECTOR']
    ctnr.camera = phdr['CAMERA']
    ctnr.ccd = phdr['CCD']

    # Observation span:
    ctnr.t_start = dhdr["TSTART"]  # BJD
    ctnr.t_stop = dhdr["TSTOP"]  # BJD
    ctnr.t_elapse = dhdr["TELAPSE"]  # days
    ctnr.deadc = dhdr["DEADC"]  # deadtime correction factor
    ctnr.t_live = dhdr["LIVETIME"]  # elapsed * deadc

    ctnr.num_frm = dhdr["NUM_FRM"]  # number of frames

    # Per-frame photon integration time in days:
    ctnr.int_time = dhdr["INT_TIME"]

    # TODO:  Understand deadtime corrections; are these in int_time?

    # Total frame time (integration + readout) in days:
    ctnr.framtim = dhdr["FRAMETIM"]

    # Number of frames and total exposure:
    ctnr.num_frm = dhdr["NUM_FRM"]
    # Note `texp` is not the same as the FITS EXPOSURE value, which is the
    # time on-source.
    ctnr.texp = sec2day * ctnr.int_time * ctnr.num_frm


class SectorImageData:
    """
    Container for image data associated with a TID from a particular sector.
    """

    def __init__(self, tid, sector, tpf_url, cache_fits=False):
        """
        Load image data from the TPF file for a particular TID and sector.
        """
        self.tid, self.sector = tid, sector
        self.tpf_url = tpf_url
        with fits.open(tpf_url, cache=cache_fits) as hdus:
            # Primary HDU has target info, incl. stellar params.
            self.phdr = hdus[0].header
            # Next HDU is a binary table with the image data.
            self.dhdr = hdus[1].header
            self.data = hdus[1].data

        read_sector_metadata(self, self.phdr, self.dhdr)

        # Epochs for the images:
        # Time is adjusted barycentric Julian date (days): BJD - 2457000 (BTJD)
        self.time = self.data["TIME"]

        # Flux images and errors, masked:
        image = self.data["FLUX"]
        self.flux = image  # all fluxes, including anomalous values
        # QUALITY != 0 indicates an "anomaly" for the datum at that epoch;
        # see SDPDD Sxn 9, Tbl 28.
        # Omit anomalies and images with any inf of NaN.
        mask = np.any(np.isfinite(image), axis=(1, 2)) & (self.data["QUALITY"] == 0)
        self.time_mask = mask
        self.mtime = self.time_mask[mask]
        self.mflux = np.ascontiguousarray(image[mask], dtype=np.float64)

        image = self.data["FLUX_ERR"]
        self.flux_err = image  # unmasked
        self.mflux_err = np.ascontiguousarray(image[mask], dtype=np.float64)

        # Time relative to midpoint:
        ref_time = 0.5 * (np.min(self.time)+np.max(self.time))
        self.rtime = self.time - ref_time
        ref_time = 0.5 * (np.min(self.mtime)+np.max(self.mtime))
        self.rmtime = self.mtime - ref_time

        # Mean, median images for masked data:
        self.mflux_mean = np.mean(self.flux, axis=0)
        self.mflux_median = np.median(self.flux, axis=0)

    def plot_mean(self):
        plt.figure()
        plt.imshow(self.mflux_mean.T, cmap="gray_r")
        plt.title('Mean TPF image, {}:s{}'.format(self.tid, self.sector))
        plt.xticks([])
        plt.yticks([])

    def plot_median(self):
        plt.figure()
        plt.imshow(self.mflux_median.T, cmap="gray_r")
        plt.title('Median TPF image, {}:s{}'.format(self.tid, self.sector))
        plt.xticks([])
        plt.yticks([])

    def plot_mmm(self):
        plt.figure()
        plt.imshow(self.mflux_median.T-self.mflux_mean.T, cmap="gray_r")
        plt.title('Median-mean TPF image, {}:s{}'.format(self.tid, self.sector))
        plt.xticks([])
        plt.yticks([])


class SectorLCData:
    """
    Container for light curve data associated with a TID from a particular sector.
    """

    def __init__(self, tid, sector, lc_url, cache_fits=False):
        """
        Load light curve data from the LC file for a particular TID and sector.
        """
        self.tid, self.sector = tid, sector
        self.lc_url = lc_url
        with fits.open(lc_url, cache=cache_fits) as hdus:
            # Primary HDU has target info, incl. stellar params.
            self.phdr = hdus[0].header
            # Next HDU is a binary table with the image data.
            self.dhdr = hdus[1].header
            self.data = hdus[1].data

        read_sector_metadata(self, self.phdr, self.dhdr)

        # Epochs for the photometry:
        # Time is adjusted barycentric Julian date (days): BJD - 2457000 (BTJD)
        self.time = self.data["TIME"]

        # Simple aperture photometry (SAP) fluxes and errors, masked:
        self.sap = self.data["SAP_FLUX"]
        # QUALITY != 0 indicates an "anomaly" for the datum at that epoch;
        # see SDPDD Sxn 9, Tbl 28.  For LC data, this is "aperture photometry
        # quality.""
        # Omit anomalies and data with inf or NaN values.
        mask = np.isfinite(self.sap) & (self.data["QUALITY"] == 0)
        self.mask = mask
        self.mtime = self.time[mask]
        self.msap = np.ascontiguousarray(self.sap[mask], dtype=np.float64)

        self.sap_err = self.data["SAP_FLUX_ERR"]
        self.msap_err = np.ascontiguousarray(self.sap_err[mask], dtype=np.float64)

        # PDC-corrected SAP photometry:
        self.sap_pdc = self.data["PDCSAP_FLUX"]
        self.msap_pdc = np.ascontiguousarray(self.sap_pdc[mask], dtype=np.float64)
        self.sap_pdc_err = self.data["PDCSAP_FLUX_ERR"]
        self.msap_pdc_err = np.ascontiguousarray(self.sap_pdc_err[mask],
                                                 dtype=np.float64)

        # Centroid location & err from PSF fitting, pxl units:
        self.x = self.data["PSF_CENTR1"]
        self.mx = np.ascontiguousarray(self.x[mask], dtype=np.float64)
        self.x_err = self.data["PSF_CENTR1_ERR"]
        self.mx_err = np.ascontiguousarray(self.x_err[mask], dtype=np.float64)
        self.y = self.data["PSF_CENTR2"]
        self.my = np.ascontiguousarray(self.y[mask], dtype=np.float64)
        self.y_err = self.data["PSF_CENTR2_ERR"]
        self.my_err = np.ascontiguousarray(self.y_err[mask], dtype=np.float64)

        # Time relative to midpoint:
        ref_time = 0.5 * (np.min(self.time)+np.max(self.time))
        self.rtime = self.time - ref_time
        ref_time = 0.5 * (np.min(self.mtime)+np.max(self.mtime))
        self.rmtime = self.mtime - ref_time


class SectorData:
    """
    Container class for storing data from a single TESS sector.
    """

    # regexp for name of target pixel file (TPF), a FITS file:
    tpf_re = re.compile(r'''
        tess                        # start of file name
        (?P<yr>\d{4})(?P<od>\d{3})             # start of timestamp: yyyy ddd
        (?P<hr>\d{2})(?P<min>\d{2})(?P<sec>\d{2})   # timestamp hh mm ss
        -s(?P<sctr>\d{4})           # sector
        -\d{16}                     # TIC
        -(?P<scid>\d{4})            # spacecraft config id
        -(?P<crmp>[a-z])            # cosmic ray mitigation procedure
        _tp.fits''', re.VERBOSE)

    # regexp for name of light curve file, a FITS file:
    lc_re = re.compile(r'''
        tess                        # start of file name
        (?P<yr>\d{4})(?P<od>\d{3})             # start of timestamp: yyyy ddd
        (?P<hr>\d{2})(?P<min>\d{2})(?P<sec>\d{2})   # timestamp hh mm ss
        -s(?P<sctr>\d{4})           # sector
        -\d{16}                     # TIC
        -(?P<scid>\d{4})            # spacecraft config id
        -(?P<crmp>[a-z])            # cosmic ray mitigation procedure
        _lc.fits''', re.VERBOSE)

    # regexp for valid data file names in archive folder:
    fname_re = re.compile('^[\w-]*.\w*$')

    # regexps for PDF files:
    dvr_re = re.compile('^[\w-]*_dvr.pdf$')
    dvs_re = re.compile('^[\w-]*_dvs.pdf$')

    def __init__(self, tid, sector, folder, folder_url,
                 get_lc=True, get_im=False, cache_fits=False):
        self.tid = tid
        self.sector = sector
        self.folder = folder
        self.folder_url = folder_url

        # Find the URLs for the various files by scraping the archive's HTML
        # via xpath, looking for href values in anchors to the files.
        r = requests.get(folder_url)
        content = html.fromstring(r.text)
        hrefs = content.xpath('//a/@href')
        # hrefs includes links to column, parent dir; keep just data files.
        self.file_list = [href for href in hrefs if self.fname_re.match(href)]
        # We assume there is only one TPF file...
        self.tpf_name = None
        for name in self.file_list:
            if self.dvr_re.match(name):
                self.dvr_name = name
            if self.dvs_re.match(name):
                self.dvs_name = name
            tpfm = self.tpf_re.match(name)
            # print(name, tpfm, sep='\n', end='\n\n')
            if tpfm:
                m = tpfm
                self.tpf_name = name
            lcm = self.lc_re.match(name)
            if lcm:
                self.lc_name = name

        # Return if no TPF; user should check self.tpf != None.
        if self.tpf_name is None:
            print('WARNING: No TPF file found for {}:s{}!'.format(self.tid, self.sector))
            return

        # Verify the sector matches.
        if sector != int(m.group('sctr')):
            raise ValueError('Sector mismatch in file name!')

        # Convert file name content into attributes.
        # Note that we get this info from the TPF file name; we are assuming
        # it matches for the LC file.
        # TODO:  Verify the TPF-LC match?

        # Spacecraft configuration map ID:
        self.scid = m.group('scid')  # 4-digit ID as a string

        # Cosmic ray mitigation procedure:
        self.crmp = m.group('crmp')  # just a string flag

        # Timestamp: year + ordinal day + h/m/s -> datetime object:
        s = '{} {} {} {} {}'.format(
            m.group('yr'), m.group('od'),
            m.group('hr'), m.group('min'), m.group('sec'))
        self.timestamp = datetime.datetime.strptime(s, '%Y %j %H %M %S')

        # MJD version (a simple float instead of a datetime object):
        self.mjd2k = datetime2mjd2k(self.timestamp)

        # Collect the sector light curve data.
        if get_lc:
            lc_url = self.folder_url + self.lc_name
            self.lc = SectorLCData(self.tid, self.sector, lc_url, cache_fits)

        # Collect the sector image data.
        if get_im:
            tpf_url = self.folder_url + self.tpf_name
            self.im = SectorImageData(self.tid, self.sector, tpf_url, cache_fits)

    def show_dvr(self):
        """
        Show the data validation report as a PDF in the default browser.

        This is an extensive report (30+ pages).
        """
        webbrowser.open(self.folder_url + self.dvr_name)

    def show_dvrs(self):
        """
        Show the data validation report summary as a PDF in the default browser.

        This is a single-page document, comprising mostly plots.
        """
        webbrowser.open(self.folder_url + self.dvs_name)

    def show_folder(self):
        """
        Show the MAST archive folder with data for this TID.
        """
        webbrowser.open(self.folder_url)


# TODO:  Max sectors is hard-wired to 10; fix this.

class TIDData:
    """
    Provide access to TESS TPF and LC data for a specified TIC target ID.
    """

    url_base = 'https://archive.stsci.edu/missions/tess/tid/'
    folder_tmp = string.Template('s${sector}/${tid1}/${tid2}/${tid3}/${tid4}/')

    def __init__(self, tid, get_lc=True, get_im=False, cache_fits=False):
        """
        Provide access to TESS data for the specified target ID.
        """
        self.tid = tid
        self.tid16 = '{:0>16d}'.format(int(tid))
        print(self.tid16)

        # Find valid sectors by checking for a valid HTTP header.
        # TODO:  What's the proper upper limit on # of sectors?
        all_sectors = range(26)
        self.sectors = {}
        for sector in all_sectors:
            folder = self.sector2folder(sector)
            url = self.url_base + folder
            r = requests.head(url)
            if r.status_code < 400:
                self.sectors[sector] = SectorData(tid, sector, folder, url,
                                                  get_lc, get_im, cache_fits)
            # print(sector, url, r.status_code, sep='\n', end='\n\n')
        self.sector_list = list(self.sectors.keys())

        # Use the FITS primary header from the 1st available sector to get
        # target metadata.
        s1 = self.sectors[self.sector_list[0]]
        try:
            read_target_metadata(self, s1.lc.phdr)
        except NameError:
            read_target_metadata(self, s1.im.phdr)

        # Collect data for valid sectors.
        # for sector, sdata in self.sectors.items():
        #     url = sdata.folder_url

    def sector2folder(self, sector):
        """
        Return the folder path for a given sector.
        """
        ss = '{:0>4d}'.format(sector)
        return self.folder_tmp.substitute(
            sector=ss,
            tid1=self.tid16[0:4], tid2=self.tid16[4:8],
            tid3=self.tid16[8:12], tid4=self.tid16[12:16])
