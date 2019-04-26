astrodata.tesser: Pythonic access to TESS data by TIC TID
=========================================================

The `tesser` module provides a single class, `TIDData`, for accessing
TESS data for objects from the TESS Input Catalog (TIC) by
Target ID (TID).  It also provides convenience functions for accessing
documentation at the MAST TESS archive, both for the TESS mission, and
for specific sources.

Acronyms used here:

* TIC: TESS Input Catalog
* TID: Target ID
* TPF: Target Pixel File (containing a time series of images)
* LC: Light curve


Documentation access
--------------------

The following functions will open TESS documentation web pages in your
default browser:

`show_tess_manual()`:
    Display the TESS archive manual

`show_products()`:
    Display the TESS data product summary page

`show_alerts_table()`:
    Display the TESS alerts data access page.

    The table is organized by TID and includes brief comments for each TID
    that triggered an alert.

    Note that this table links to LC and TPF FITS files that are named
    differently from those accessed by this module---they are archived
    as high-level science products (HLSPs).


Data access
-----------

To access data for TID 261136679, for example, create an instance of the
`TIDData` class:

    from astrodata import tesser

    td = tesser.TIDData(261136679)

By default, this will download only pre-processed LC data, and the FITS file
containing the data will not be locally cached.  You can alter the cache 
behavior via the `cache_fits` argument.  You can control whether or not
LC or TPF data are downloaded via `get_lc` and `get_im` arguments.  For example,
this call gets only the TPF (image) data, and caches it:

    td = tesser.TIDData(261136679, get_lc=False, get_im=True, cache_fits=True)

Note that the TPF FITS files are large (~100 MB).  Monitor your cache!  Its location is platform-dependent but is likely to be in a `.astropy` folder in your home directory. (To locat it, call `astropy.utils.data._get_download_cache_locs()`).
