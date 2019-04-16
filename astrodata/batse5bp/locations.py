"""
Specify locations, locally and the COSSC, for accessing GRB data.

No modules are imported here; 'from locations import *' should be safe.

Created 2012-05-07 by Tom Loredo
2019:  Converted to Python 3
"""

# Local locations for persistant data:
root = 'BATSE_5Bp_Data'  # default name in CWD; can change via load_catalog()
raw_cache = 'raw_cache'  # subdir for raw data files
summaries = 'summaries.pkl'


# CGRO SSC URLs for 5Bp catalog summary data:

basic_url = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/current/tables/basic_table.txt'
durn_url = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/current/tables/duration_table.txt'

# Brightnesses are in two tables; 4B and new in 5Bp.
bright_url4 = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/current/tables/flux_table_4b.txt'
bright_url5 = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/current/tables/flux_table.txt'

# Comments are only available for 4B bursts.
comments_url4 = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/4b/tables/4br_grossc.comments'


# COSSC FTP site for trigger data:
trigger_url = 'ftp://legacy.gsfc.nasa.gov/compton/data/batse/trigger/'

# TODO:  ASCII data should be used only for getting the BATSE team's
# subjective start/end times for GRBs; the actual count data should
# come from DISCLA, PREB, and DISCSC files; it will not all be in 64 ms bins.

# COSSC FTP site for ASCII 64ms count data:
ascii64_url = 'ftp://legacy.gsfc.nasa.gov/compton/data/batse/ascii_data/64ms/'


def trigger_paths(tnum):
    """
    Return components of the path to the data for trigger `tnum`, both at
    the COSSC FTP site and locally:  (group, trigger, remote).

    group       name of the group directory containing the trigger directory
    trigger     name of the trigger directory, "<tnum>_burst"
    remote      tail of the ftp path @ COSSC, "<group>/<trigger>"
    """
    # Trigger data is in groups spanning a trigger range of 200.
    fac = tnum//200
    l, u = fac*200+1, (fac+1)*200
    group = '%05d_%05d' % (l, u)
    trigger = '%05d_burst' % tnum
    remote = group + '/' + trigger + '/'  # don't use join; not a local path
    return group, trigger, remote


def ascii64_paths(tnum):
    """
    Return components of the path to the data for trigger `tnum`, both at
    the COSSC FTP site and locally:  (group, fname, remote).

    group       name of the group directory containing the trigger directory
    rfname      name of the remote ASCII data file, "cat64ms.<tnum>"
    remote      tail of the ftp path @ COSSC, "<group>/<fname>"
    """
    # ASCII data is in groups spanning a trigger range of 1000.
    fac = tnum//1000
    u = fac*1000
    group = 'trig%05d' % u
    rfname = 'cat64ms.%05d' % tnum
    remote = group + '/' + rfname  # don't use join; not a local path
    return group, rfname, remote
