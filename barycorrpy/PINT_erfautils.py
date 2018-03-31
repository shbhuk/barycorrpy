# This is taken from PINT - https://github.com/nanograv/PINT/blob/master/pint/erfautils.py. For licensing see license file on PINT github repository. 
# https://github.com/nanograv/PINT/blob/ea1c3b9f71c06a73d6613d647770b99924f7c9bf/LICENSE.md
from __future__ import division
import numpy as np
import astropy.units as u
import datetime
try:
    import astropy.erfa as erfa
except ImportError:
    import astropy._erfa as erfa
import astropy.table as table
from astropy.time import Time
from astropy.utils.iers import conf,IERS_A, IERS_A_URL, IERS_B, IERS_B_URL, IERS
from astropy.utils.data import download_file

SECS_PER_DAY = erfa.DAYSEC
conf.auto_max_age=15

iers_b_file = download_file(IERS_B_URL, cache=True)
iers_b = IERS_B.open(iers_b_file)
IERS.iers_table = iers_b
iers_tab = IERS.iers_table

# Earth rotation rate in radians per UT1 second
#
# This is from Capitaine, Guinot, McCarthy, 2000 and is
# in IAU Resolution B1.8 on the Earth Rotation Angle (ERA)
# and the relation of it to UT1.  The number 1.00273781191135448
# below is a defining constant.  See here:
# http://iau-c31.nict.go.jp/pdf/2009_IAUGA_JD6/JD06_capitaine_wallace.pdf
OM = 1.00273781191135448 * 2.0 * np.pi / SECS_PER_DAY

# arcsec to radians
asec2rad = 4.84813681109536e-06


def gcrs_posvel_from_itrf(loc, toas,tts):
    """Return a list of PosVel instances for the observatory at the TOA times.
    Observatory location should be given in the loc argument as an astropy
    EarthLocation object. This location will be in the ITRF frame (i.e.
    co-rotating with the Earth).

    This routine returns a list of PosVel instances, containing the
    positions (m) and velocities (m / s) at the times of the toas and
    referenced to the Earth-centered Inertial (ECI, aka GCRS) coordinates.
    This routine is basically SOFA's pvtob() [Position and velocity of
    a terrestrial observing station] with an extra rotation from c2ixys()
    [Form the celestial to intermediate-frame-of-date matrix given the CIP
    X,Y and the CIO locator s].
    
    
    EDIT: 
        Shubham Kanodia 1st December 2017
        Use TT times from Barycorrpy leap second management routine (utc_tdb) and not Astropy.
        
        tts : Time Object in TT format. Generated by utc_tdb.JDUTC_to_JDTDB()
    
    
    """
    # If the input is a single TOA (i.e. a row from the table),
    # then put it into a list
    if type(toas) == table.row.Row:
        ttoas = Time([toas['mjd']])
        N = 1
    elif type(toas) == table.table.Table:
        N = len(toas)
        ttoas = toas['mjd']
    else:
        if toas.isscalar:
            ttoas = Time([toas])
        else:
            ttoas = toas
        N = len(ttoas)
        


    # Get various times from the TOAs as arrays
    #tts = np.asarray([(t.jd1, t.jd2) for t in ttoas.tt]).T
    tts = np.asarray([tts.jd1,tts.jd2])  # Use Leap second management routine from Barycorrpy

    ut1s = np.asarray([(t.jd1, t.jd2) for t in ttoas.ut1]).T
    mjds = np.asarray(ttoas.mjd)

    # Get x, y coords of Celestial Intermediate Pole and CIO locator s
    X, Y, S = erfa.xys00a(*tts)

    # Get dX and dY from IERS B in arcsec and convert to radians
    dX = np.interp(mjds, iers_tab['MJD'], iers_tab['dX_2000A']) * asec2rad
    dY = np.interp(mjds, iers_tab['MJD'], iers_tab['dY_2000A']) * asec2rad

    # Get GCRS to CIRS matrices
    rc2i = erfa.c2ixys(X+dX, Y+dY, S)

    # Gets the TIO locator s'
    sp = erfa.sp00(*tts)

    # Get X and Y from IERS A in arcsec and convert to radians
    #xp = np.interp(mjds, iers_tab['MJD'], iers_tab['PM_X_B']) * asec2rad
    #yp = np.interp(mjds, iers_tab['MJD'], iers_tab['PM_Y_B']) * asec2rad
    # Get X and Y from IERS B in arcsec and convert to radians
    xp = np.interp(mjds, iers_tab['MJD'], iers_tab['PM_x']) * asec2rad
    yp = np.interp(mjds, iers_tab['MJD'], iers_tab['PM_y']) * asec2rad

    # Get the polar motion matrices
    rpm = erfa.pom00(xp, yp, sp)

    # Observatory geocentric coords in m
    xyzm = np.array([a.to(u.m).value for a in loc.geocentric])
    x, y, z = np.dot(xyzm, rpm).T

    # Functions of Earth Rotation Angle
    theta = erfa.era00(*ut1s)
    s, c = np.sin(theta), np.cos(theta)
    sx, cx = s * x, c * x
    sy, cy = s * y, c * y

    # Initial positions and velocities
    iposs = np.asarray([cx - sy, sx + cy, z]).T
    ivels = np.asarray([OM * (-sx - cy), OM * (cx - sy), \
                        np.zeros_like(x)]).T
    # There is probably a way to do this with np.einsum or something...
    # and here it is .
    poss = np.empty((N, 3), dtype=np.float64)
    vels = np.empty((N, 3), dtype=np.float64)
    poss = np.einsum('ij,ijk->ik', iposs, rc2i)
    vels = np.einsum('ij,ijk->ik', ivels, rc2i)
    return poss,vels


