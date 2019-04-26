astrodata.tesser: Pythonic access to TESS data by TIC TID
=========================================================

The `tesser` module provides a single class, `TIDData`, for accessing
TESS data for objects from the TESS Input Catalog (TIC) by
Target ID (TID).  It also provides convenience functions for accessing
documentation at the MAST TESS archive, both for the TESS mission, and
for specific sources.


Documentation access
--------------------

The following functions will open TESS documentation web pages in your
default browser:

* `show_tess_manual()`:  Display the TESS archive manual
* `show_products()`:  Display the TESS data product summary page
* `show_alerts_table()`:  Display the TESS alerts data access page.
    The table is organized by TID and includes brief comments for each TID
    that triggered an alert.
    Note that this table links to LC and TPF FITS files that are named
    differently from those accessed by this module---they are archived
    as high-level science products (HLSPs).

