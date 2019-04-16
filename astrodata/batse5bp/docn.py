"""
Documentation supporting the 5Bp data.

Created 2012-05-07 by Tom Loredo
2019:  Converted to Python 3
"""

import webbrowser

__all__ = ['show_ssc', 'show_archive', 'show_tech', 'show_sw',
           'show_4B', 'show_5Bp', 'show_problems']


# Access to CGRO SSC documents, incl. 4B, 5Bp ("Current") catalog web pages:
url_ssc = 'http://heasarc.gsfc.nasa.gov/docs/cgro/'
url_archive = 'http://heasarc.gsfc.nasa.gov/docs/journal/cgro7.html'
url_tech = 'http://heasarc.gsfc.nasa.gov/docs/cgro/cossc/nra/appendix_g.html#V.%20BATSE%20GUEST%20INVESTIGATOR%20PROGRAM'
url_sw = 'ftp://legacy.gsfc.nasa.gov/compton/software/'
url_4B = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/4b/'
url_5Bp = 'http://gammaray.msfc.nasa.gov/batse/grb/catalog/current/'
url_problems = 'ftp://legacy.gsfc.nasa.gov/compton/data/batse/trigger/data_problems/'


def show_ssc():
    """
    Show the COSSC web site in the user's default web browser.
    """
    webbrowser.open(url_ssc)


def show_archive():
    """
    Show the *Legacy* article with a CGRO archive web site overview in the
    user's default web browser.
    """
    webbrowser.open(url_archive)


def show_4B():
    """
    Show the 4B catalog web site in the user's default web browser.
    """
    webbrowser.open(url_4B)


def show_5Bp():
    """
    Show the 5Bp ("current") catalog web site in the user's default web browser.
    """
    webbrowser.open(url_5Bp)


def show_problems():
    """
    Show the BATSE data problems archive in the user's default web browser.
    """
    webbrowser.open(url_problems)


def show_tech():
    """
    Display the BATSE technical description (from the CGRO GI program RA) in the
    user's default web browser.
    """
    webbrowser.open(url_tech)


def show_sw():
    """
    Display the SSC software directory in the user's default web browser.
    """
    webbrowser.open(url_sw)


# Basic table description:
basic = """
There are twelve columns in the Basic Table file:

 1. The BATSE trigger number
 2. The BATSE Catalog burst name
 3. The truncated Julian Date (TJD) of the trigger; TJD = JD - 2440000.5
 4. The time in decimal seconds of day (UT) of the trigger
 5. Right ascension (J2000) in decimal degrees
 6. Declination (J2000) in decimal degrees
 7. Galactic longitude in decimal degrees
 8. Galactic latitude in decimal degrees
 9. Radius in decimal degrees of positional error box
10. Angle in decimal degrees of geocenter (the angle between the burst and the nadir, as measured from the satellite)
11. Overwrite flag: Y(true) if this burst overwrote an earlier, weaker trigger. N(false) otherwise
12. Overwritten flag: Y(true) if this burst was overwritten by a later, more intense trigger. N(false) otherwise
"""

# Flux table description, from 4B catalog:
flux = """\
The FLUX Table contains the fluence and peak flux values for the BATSE
gamma-ray bursts between 19 April, 1991 and 29 August, 1996. There are 1292
bursts, each specified by the BATSE trigger number. This table contains 18
columns. All fluences and their errors have units of of ergs/cm^2. All peak
fluxes and their errors have units of photons/cm^2/sec. The errors are one
sigma statistical errors. The peak flux times are expressed in decimal seconds
relative to the burst trigger time for the end of the interval in which the
flux was calculated. The channel 1,2,3 and 4 fluences cover the energy ranges
20-50 keV, 50-100 keV, 100-300 keV, and E > 300 keV respectively. The peak
flux energy range is 50-300 keV, coinciding with the energy range of the
nominal BATSE on-board burst trigger.

Since channel 4 is an integral channel, fluences given for this channel are
quite sensitive to the assumed spectral form. Spectral analyses in this energy
range should be performed with higher resolution data types.

Many of the bursts between March 1992 and March 1993 have significant gaps in
the data and are not included in the table.
"""

# Duration table description, from 4B catalog:
durn = """\
The DURATION Table contains values for T90 and T50, quantities related to
burst duration, for 1234 gamma-ray bursts that triggered the BATSE LAD
detectors between April 1991 and 29 August 1996. T90 measures the duration of
the time interval during which 90% of the total observed counts have been
detected. The start of the T90 interval is defined by the time at which 5% of
the total counts have been detected, and the end of the T90 interval is
defined by the time at which 95% of the total counts have been detected. T50
is similarly defined using the times at which 25% and 75% of the counts have
been detected. T90 and T50 are calculated using data summed over the 4 LAD
discriminator channels and using data summed over only those detectors that
satisfied the BATSE trigger criteria.

Users must note that T90 and T50 are not available for those bursts which
suffer from data gaps during the event; the integration procedure inherently
fails in these cases. However, visual estimates of the burst duration are
provided in the BATSE Comments table for those bursts with sufficient data
coverage. Users may also find other pertinent comments concerning the
calculated value of T90 and T50 in the BATSE COMMENTS table, and it is highly
recommended that the COMMENTS table be consulted before any distribution
selected on T90 or T50 is used.
"""

# Comments description; only for 4B bursts:
comments = """\
The COMMENTS Table contains comments relevant to BATSE gamma-ray burst data
found in the GROSSC BATSE burst catalog. Each category of comment is sorted by
a flag identification and ordered by BATSE trigger number. A gamma-ray burst
may have more than one entry or may have no entry.

Flag  Definition
Q     comments on data quality
A     additional observations by other instruments
O     general comments
L     comments on the gamma-ray burst coordinates
T     comments on the gamma-ray burst duration
"""

# Reclassification information as of 2005-08-01; tab-delimited:
reclass = """\
288  GRB. Reclassified from solar flare in May 1995. Changed from 00288_flare to 00288_burst

1530     Probably GRB; very little data. Not available at COSSC

2154     GRB. No data except max rates. Not available at COSSC

2162     Probably GRB. No data except max rates. Not available at COSSC

2411     Thought to be a particle event. Changed from 02411_burst to 02411_particles

2548     Commanded trigger to get data on a transient source. Not available at COSSC

2777     Solar flare

2778     Solar flare

3452     Unknown - possibly a GRB. 03452_burst

3720     Triggered on a Cyg X-1 fluctuation in dets. 4,5; a GRB occurred at approx. T+25 s, primarily in dets. 7,5,3,1. 03720_burst

3904     The data in the TTE file seem to be partially corrupted. The reason is not understood but may be a hardware problem. The TTE data seem to change from the trigger detectors to two other detectors. The trigger detectors (DSELB) were 2 and 6. The first quarter of TTE (first half of packet ID 33Hex), which is pre-load time, has all detectors as it should. The first two packets and part of the third packet of the second quarter of the TTE data (post-load time) have detectors 2 and 6 as they should. However, part of the way through the third packet (packet sequence number 67 of packet ID 33Hex) and for the remainder of the TTE file, only detectors 0 and 1 appear.

3931     Terrestrial Gamma-ray Flash. Changed from 03931_burst to 03931_tgf TTE & STTE data

4327     During flight software revisions, the TTE data in the second board became bad for a limited time. Therefore, triggers no. 4317 - 4336 bad TTE data in the second board (second quarter of TTE data). (The TTE data for trigger no. 4317 should be all bad for related but different reasons.)

5339     TTE data are OK but STTE data are for trigger 5338 because of overwriting.
"""
