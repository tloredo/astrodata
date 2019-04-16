"""
Access to online documentation for TESS data products.
"""

import webbrowser

__all__ = ['show_tess_manual', 'show_products', 'show_alerts_table']


def show_tess_manual():
    """
    Show the TESS archive manual web site @ MAST in the default browser.
    """
    webbrowser.open('https://outerspace.stsci.edu/'
                    'display/TESS/TESS+Archive+Manual')


def show_products():
    """
    Show the TESS data product summary page @ MAST in the default browser.

    This page links to the current TESS Science Data Products Description
    Document (SDPDD).
    """
    webbrowser.open('https://archive.stsci.edu/tess/all_products.html')


def show_alerts_table():
    """
    Show the TESS alerts data access page @ MAST in the default browser.

    The table is organized by TID and includes brief comments for each TID
    that triggered an alert.

    Note that this table links to LC and TPF FITS files that are named
    differently from those accessed by this module---they are archived
    as high-level science products (HLSPs).
    """
    # TODO:  Verify that the HLSP LC and TPF files are the same as those
    # accessed here, for TIDs that triggered an alert.
    webbrowser.open('https://archive.stsci.edu/prepds/tess-data-alerts/')
