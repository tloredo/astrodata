"""
tesser - TESS archive intERaction

Provide access to TESS archive data via TIC TID.

TIC = Target Input Catalog

TID = Target ID

TPF = Target Pixel File
https://heasarc.gsfc.nasa.gov/docs/tess/data-access.html

The data exposed include:
* Target data, including stellar parameters
* Observational parameters (sectors, camera, CCD, times)
* TPF image time series (fluxes and errors)
* Simple Aperture Photometry (SAP) light curve (fluxes and errors)
* SAP light curve with Presearch Data Conditioning (PDC)

Flux calibration and photometry are from the TESS Science Processing Operations
Center (SPOC) pipeline; for info about the pipeline see:

The TESS science processing operations center (2016)
Jenkins et al., SPIE Proceedings
https://spie.org/Publications/Proceedings/Paper/10.1117/12.2233418
https://heasarc.gsfc.nasa.gov/docs/tess/docs/jenkinsSPIE2016-copyright.pdf

Technical information about the FITS file contents are available in the
TESS Science Data Products Description Document (SDPDD), available at the
URL accessed by `show_products()`.

Example URL for TPF FITS file:
https://archive.stsci.edu/missions/tess/tid/s0001/0000/0002/6113/6679/tess2018206045859-s0001-0000000261136679-0120-s_tp.fits

Directory path: tid/s{sctr}/{tid1}/{tid2}/{tid3}/{tid4}/
File name: tess{date-time}-s{sctr}-{tid}-{scid}-{cr}_tp.fits

Notes on directory path:
{sctr} = A zero-padded, four-digit integer indicating the sector in which
  the data were collected, starting with Sector 1
{tid1} = A zero-padded, four-digit integer consisting of digits 1-4 of the
  full, zero-padded TIC ID.
{tid2} = A zero-padded, four-digit integer consisting of digits 5-8 of the
  full, zero-padded TIC ID.
{tid3} = A zero-padded, four-digit integer consisting of digits 9-12 of the
  full, zero-padded TIC ID.
{tid4} = A zero-padded, four-digit integer consisting of digits 13-16 of the
   full, zero-padded TIC ID.

Notes on file name:
{date-time} = The timestamp associated with this file, in the yyyydddhhmmss
  format.
{sctr} = A zero-padded, four-digit integer indicating the sector in which
  the data were collected, starting with Sector 1
{tid} = A zero-padded, 16-digit target identifier that refers to an object in
  the TESS Input Catalog.
{scid} = A zero-padded, four-digit identifier of the spacecraft configuration
  map used to process this data.
{cr} = A string character that denotes the cosmic ray mitigation procedure.
  Possible values are:
      'x': No mitigation performed at the SPOC.
      's': Mitigation performed on the spacecraft.
      'a': A SPOC mitigation algorithm was used.
      'b': Both a SPOC and onboard spacecraft algorithm was used.


This module was built for the community of students participating in
Cornell University's monthly TESS Hack Days, organized by Nikole Lewis
in Spring 2019.  Archive access is loosely based on Thea Kozakis's find_tp
function.  Data masking is based on Dan Foreman-Mackey's `exoplanet` package
tutorial on TESS LC fitting.

Created 2019-01-19 by Tom Loredo
"""

from .docn import *
from .containers import TIDData
