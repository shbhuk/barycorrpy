from __future__ import division
from __future__ import print_function
#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric
from astropy.time import Time
import math
import astropy.constants as ac
import numpy as np
import os


from . import find_hip
from . import PINT_erfautils as PINT
from . import utc_tdb
from .utils import flux_weighting,get_stellar_data


AU = ac.au.value # [m]
c = ac.c.value # Speed of light [m/s]
pctoau = 3600*180/np.pi # No. of AU in one parsec
year = 365.25*3600*24 # [s]
kmstoauyr = 1000 * year/(AU)


def SolarBarycentricCorrection(JDUTC, loc, zmeas, ephemeris, leap_dir, leap_update):
    """
    Perform barycentric correction for Solar observations.
    INPUTS:
        JDUTC: Astropy Time object in UTC scale with JD format
        loc: Astropy EarthLocation Object) 
        zmeas : Measured redshift (e.g., the result of cross correlation with template spectrum). Default is 0.
                The redshift measured by the spectrograph before any barycentric correction. Therefore zmeas includes the barycentric 
                velocity of the observatory.
        ephemeris: Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                    For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                    ['de432s','de430',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir: Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is 
        script directory.
        leap_update: If True, when the leap second file is more than 6 months old will attempt to download a new one.
                      If False, then will just give a warning message. Default is True.
    OUTPUTS:
        v_true: The true radial velocity of the Sun.
                If zmeas is included to show the measured absolute redshift for the Sun as measured by an instrument, 
                then in this formulation, ztrue will show the synodic period of Jupiter
                as seen from Earth.
        v_predicted: Ideal redshift measured for the Sun from Earth for given location and time.
                This output returns the theoretical prediction for the redshift which includes the barycentric component.
        
        The formula used is ztrue = ((1.+zb)*(1.+zmeas)-1.)
        Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) times c.
        
    
    """

    # Convert times to obtain TDB and TT

    JDTDB, JDTT, warning, error = utc_tdb.JDUTC_to_JDTDB(JDUTC)

    ##### NUTATION, PRECESSION, ETC. #####

    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    ##### EPHEMERIDES #####

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris) # [km]
    r_earth = r_eci + earth_geo[0].xyz.value*1000. # [m]
    v_geo = earth_geo[1].xyz.value*1000./86400.  # [m/s]

    # Relativistic Addition of Velocities
    v_earth = (v_eci+v_geo) / (1.+v_eci*v_geo/c**2) # [m/s]
    beta_earth = v_earth / c

    gamma_earth = 1. / np.sqrt(1.-sum(beta_earth**2))



    ##ignoring retarded time for now##

    target_ephem = get_body_barycentric_posvel('sun', JDTDB, ephemeris=ephemeris)

    r_target = target_ephem[0].xyz.value*1000. #[m]
    v_target = target_ephem[1].xyz.value*1000./86400.  # [m/s]
    beta_target = v_target / c

    gamma_target = 1. / np.sqrt(1.-sum(beta_target**2))

    Pos_vector = r_target - r_earth
    Pos_mag = np.sqrt(sum(Pos_vector**2)) #[m]
    Pos_hat = Pos_vector/Pos_mag



    zGREarth =  ac.G.value * ac.M_earth.value / ((ac.c.value**2)*(np.sqrt(np.sum(r_eci**2))))\
                + ac.G.value * ac.M_sun.value / ((ac.c.value**2)*(np.sqrt(np.sum(Pos_vector**2))))

    zGRSun =  - (ac.G.value * ac.M_sun.value) / ((ac.c.value**2) * (np.sqrt(np.sum(ac.R_sun.value**2))))

    zpredicted= ( (gamma_target * (1 + np.dot(beta_target,Pos_hat))*(1+zGREarth)) / (gamma_earth * (1 + np.dot(beta_earth,Pos_hat)) * (1+zGRSun))  ) - 1

    zb = ((gamma_earth*(1 + np.dot(beta_earth,Pos_hat))) / (1+zGREarth)) - 1

    v_true = c * ((1.+zb)*(1.+zpredicted)-1.)  # [m/s]
    v_predicted = c * zpredicted # [m/s]

    
    return v_true, v_predicted