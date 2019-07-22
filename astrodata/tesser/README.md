# `tesser`: Pythonic access to TESS data by TIC TID

The `tesser` module provides a single class, `TIDData`, for accessing
TESS data for objects from the TESS Input Catalog (TIC) identified by
Target ID (TID).  It also provides convenience functions for accessing
documentation at the MAST TESS archive, both for the TESS mission, and
for specific targets.

The data exposed include:
* Target metadata, including stellar parameters
* Observational parameters (sectors, camera ID, CCD ID, times)
* TPF image time series (pixel fluxes and errors)
* Simple Aperture Photometry (SAP) light curve (fluxes and errors)
* SAP light curve with Presearch Data Conditioning (PDC)

Acronyms used here:

* TIC: TESS Input Catalog (see: [TIC Documentation](https://archive.stsci.edu/tess/tic_doc.html), which also covers the CTL)
* CTL: Candidate Target List (a ranked subset of TIC targets, used to select targets for TESS
  2-min cadence observations)
* TID: Target ID
* TPF: Target Pixel File (containing a time series of images)
* LC: Light curve
* SAP: Simple Aperture Photometry
* PDC: Presearch Data Conditioning (removal of instrumental trends and other artifacts)


## TESS documentation access

The following functions will open TESS documentation web pages in your
default browser:

`show_tess_manual()`: 
Display the [*TESS Archive Manual*](https://outerspace.stsci.edu/display/TESS/TESS+Archive+Manual).

`show_products()`: 
Display the TESS [Data Product Summary page](https://archive.stsci.edu/tess/all_products.html), which provides a table describing the file name conventions and locations for various data products, and links to relevant sections of the [TESS Science Data Products Description Document](https://archive.stsci.edu/missions/tess/doc/EXP-TESS-ARC-ICD-TM-0014.pdf) for each data product.

`show_alerts_table()`: 
Display the [Data Products From TESS Data Alerts](https://archive.stsci.edu/prepds/tess-data-alerts/) page providing TPFs, LCs, and data validation reports for *planet candidates*.  The table is organized by TID and includes brief comments for each TID that triggered an alert.  Note that this table links to LC and TPF FITS files that are named differently from those accessed by this module—for alerts, the data are additionally archived as high-level science products (HLSPs).

Access to target-specific data validation reports is provided via the `TIDData` class; see below.


## Data access

To access data for TID 261136679, for example, create an instance of the
`TIDData` class:
```python
from astrodata import tesser

td = tesser.TIDData(261136679)
```
By default, this will download only pre-processed LC data (not the more voluminous TPF data), and the FITS files
containing the data will not be locally cached.  You can alter the cache 
behavior via the `cache_fits` argument.  You can control whether or not
LC or TPF data are downloaded via `get_lc` and `get_im` arguments.  For example,
this call gets only the TPF (image) data, and caches it:

```python
td = tesser.TIDData(261136679, get_lc=False, get_im=True, cache_fits=True)
```
Note that the FITS files for TPF data are large (~100 MB).  If you cache the FITS files, monitor your cache!  Its location is platform-dependent but is likely to be in a `.astropy` folder in your home directory. (To locate it, call `astropy.utils.data._get_download_cache_locs()`).

The `TIDData` class queries the MAST archive for data from all possible sectors, and provides sector-specific access to the data if finds via a `sectors` attribute (a dictionary, indexed by sector number):

```python
In [1]: td.sectors
Out[1]:
{1: <astrodata.tesser.containers.SectorData at 0x1042e83c8>,
 4: <astrodata.tesser.containers.SectorData at 0x1c1e31beb8>,
 8: <astrodata.tesser.containers.SectorData at 0x1c2ac4fb70>}
```

For convenience, the `sector_list` attribute is a list of sector numbers identifying the sectors with data for the TID:

```python
In [2]: td.sector_list
Out[2]: [1, 4, 8]
```

Each item in `sectors` is an instance of the  `SectorData` class, which collects LC and TPF data for a single sector.

Each LC and TPF FITS file includes metadata for the target in its primary FITS header. `TIDData` uses the primary header from one of the sectors to load metadata as attributs of the `TIDData` instance.

### Light curve data

The `lc` atribute collects LC data for a sector, as produced by the TESS pipeline.

