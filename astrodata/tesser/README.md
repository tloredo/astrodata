# `astrodata.tesser`: Pythonic access to TESS data by TIC TID

The `tesser` module provides a single class, `TIDData`, for accessing
TESS data for objects from the TESS Input Catalog (TIC) by
Target ID (TID).  It also provides convenience functions for accessing
documentation at the MAST TESS archive, both for the TESS mission, and
for specific sources.

The data exposed include:
* Target metadata, including stellar parameters
* Observational parameters (sectors, camera ID, CCD ID, times)
* TPF image time series (fluxes and errors)
* Simple Aperture Photometry (SAP) light curve (fluxes and errors)
* SAP light curve with Presearch Data Conditioning (PDC)

Acronyms used here:

* TIC: TESS Input Catalog
* TID: Target ID
* TPF: Target Pixel File (containing a time series of images)
* LC: Light curve
* SAP: Simple Aperture Photometry
* PDC: Presearch Data Conditioning (removal of instrumental trends and other artifacts)


## Documentation access

The following functions will open TESS documentation web pages in your
default browser:

`show_tess_manual()`: 
Display the TESS archive manual.

`show_products()`: 
Display the TESS data product summary page.

`show_alerts_table()`: 
Display the TESS alerts data access page.  The table is organized by TID and includes brief comments for each TID that triggered an alert.  Note that this table links to LC and TPF FITS files that are named differently from those accessed by this module---they are archived as high-level science products (HLSPs).


## Data access

To access data for TID 261136679, for example, create an instance of the
`TIDData` class:
```python
from astrodata import tesser

td = tesser.TIDData(261136679)
```
By default, this will download only pre-processed LC data, and the FITS file
containing the data will not be locally cached.  You can alter the cache 
behavior via the `cache_fits` argument.  You can control whether or not
LC or TPF data are downloaded via `get_lc` and `get_im` arguments.  For example,
this call gets only the TPF (image) data, and caches it:
```python
td = tesser.TIDData(261136679, get_lc=False, get_im=True, cache_fits=True)
```
Note that the TPF FITS files are large (~100 MB).  Monitor your cache!  Its location is platform-dependent but is likely to be in a `.astropy` folder in your home directory. (To locat it, call `astropy.utils.data._get_download_cache_locs()`).

The `TIDData` class queries the MAST archive for data from all possible sectors, and provides access to the data by sector, via a `sectors` attribute (an ordered dictionary, indexed by sector number):

```python
In [1]: td.sectors
Out[1]:
{1: <astrodata.tesser.containers.SectorData at 0x1042e83c8>,
 4: <astrodata.tesser.containers.SectorData at 0x1c1e31beb8>,
 8: <astrodata.tesser.containers.SectorData at 0x1c2ac4fb70>}
```

For convenience, the `sector_list` attribute is a list of sector numbers with data for the TID:

```python
In [2]: td.sector_list
Out[2]: [1, 4, 8]
```

Each item in `sectors` is an instance of the  `SectorData` class.