from __future__ import division
from __future__ import print_function
#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric
from astropy.time import Time
from astroquery.jplhorizons import Horizons

import astropy.constants as ac
import numpy as np
import os


from . import find_hip
from . import PINT_erfautils as PINT
from . import utc_tdb
from .utils import flux_weighting, get_stellar_data, CalculatePositionVector
from .PhysicalConstants import *


def SolarBarycentricCorrection(JDUTC, loc, zmeas=0, ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True, predictive=False):
    """
    Perform barycentric correction for Solar observations.
    In the solar case, ignores retarded time, i.e. the change in position angle of the Sun wrt the observatory in the 8 minutes it takes for the light to reach Earth.
    This is because the Solar velocity is low enough (~ 10m/s), that in 8 mins the change in position over 1 AU is negligible.
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

        predictive : If True, then instead of returning v_true, returns v_predicted.
        Default: False, and return is v_true from Wright and Eastman (2014)

        See OUTPUTs for description

    OUTPUTS:
            If predictive = False
                v_true: The true radial velocity of the Sun for an observer at the observatory but in an inertial frame not moving.
                    If zmeas is included to show the measured absolute redshift for the Sun as measured by an instrument,
                    then in this formulation, v_true will show the motion of the Sun,
                    which is mostly dominated by the synodic period of Jupiter as seen from Earth.

            Else if predictive = True
                v_predicted: Ideal redshift measured for the Sun from Earth for given location and time.
                    This output returns the theoretical prediction for the measured redshift which includes the barycentric component.
                    This will be the measurement of a noiseless RV instrument observing the Sun.


        The formula used is ztrue = ((1.+ zb)*(1.+ zmeas)-1.)
        Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) x speed of light (c).


    """

    # Convert times to obtain TDB and TT

    JDTDB, JDTT, warning, error = utc_tdb.JDUTC_to_JDTDB(JDUTC)

    ##### NUTATION, PRECESSION, ETC. #####

    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    ##### EPHEMERIDES #####

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris) # [km]
    PosVector_EarthSSB = r_eci + earth_geo[0].xyz.value*1000. # [m]
    v_geo = earth_geo[1].xyz.value*1000./86400.  # [m/s]

    # Relativistic Addition of Velocities
    VelVector_EarthSSB = (v_eci+v_geo) / (1.+ np.sum(v_eci*v_geo)/c**2) # [m/s]
    BetaEarth = VelVector_EarthSSB / c

    GammaEarth = 1. / np.sqrt(1.- np.sum(BetaEarth**2))

    ## Ignoring retarded time for Solar vectors

    solar_ephem = get_body_barycentric_posvel('sun', JDTDB, ephemeris=ephemeris)

    PosVector_SolSSB = solar_ephem[0].xyz.value*1000. #[m]
    VelVector_SolSSB = solar_ephem[1].xyz.value*1000./86400.  # [m/s]
    BetaSolar = VelVector_SolSSB / c

    GammaSolar = 1. / np.sqrt(1.- np.sum(BetaSolar**2))

    PosVector_SolEarth, PosMag_SolEarth, PosHat_SolEarth = CalculatePositionVector(r1=PosVector_SolSSB, r2=PosVector_EarthSSB)

    # To correct for redshift experienced by photon due to Earth's gravity well,
    # also as it travels further in the Sol System from the Sun to the position of the Earth (Solar gravity well)
    # Correcting for this should give a blue shift.
    zGREarth =   - ac.G.value * ac.M_earth.value / ((ac.c.value**2)*(np.sqrt(np.sum(r_eci**2))))\
                - ac.G.value * ac.M_sun.value / ((ac.c.value**2)*(np.sqrt(np.sum(PosVector_SolEarth**2))))

    # To correct for redshift experienced by photon and place it in the SSB
    zGRSun =  - (ac.G.value * ac.M_sun.value) / ((ac.c.value**2) * (np.sqrt(np.sum(ac.R_sun.value**2))))

    zpredicted= ( (GammaSolar * (1 + np.dot(BetaSolar,PosHat_SolEarth))*(1+zGREarth)) /
                (GammaEarth * (1 + np.dot(BetaEarth,PosHat_SolEarth)) * (1+zGRSun))  ) - 1

    zb = ((GammaEarth*(1 + np.dot(BetaEarth,PosHat_SolEarth))) / (1+zGREarth)) - 1

    v_true = c * ((1.+zb)*(1.+ zmeas)-1.)  # [m/s]
    v_predicted = c * zpredicted # [m/s]


    if predictive:
        vel = v_predicted
    else:
        vel = v_true

    return vel, warning, error


def ReflectedLightBarycentricCorrection(SolSystemTarget, JDUTC, loc, zmeas=0, HorizonsID_type='smallbody',
    ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True, predictive=False):
    """
    Computing the barycentric corrections for reflected light observations of a target in the Solar system.

    INPUTS:
        TargetName
        JDUTC: Astropy Time object in UTC scale with JD format
        loc: Astropy EarthLocation Object)
        zmeas : Measured redshift (e.g., the result of cross correlation with template spectrum). Default is 0.
                The redshift measured by the spectrograph before any barycentric correction. Therefore zmeas includes the barycentric
                velocity of the observatory.
        ephemeris: Ephemeris used only for calculating the position and velocity vectors of the Earth wrt the Solar System Barycenter (SSB).
                For reflected light target and the Sun, the HORIZONS query is used from Astroquery.
                Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                ['de432s','de430',
                'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir: Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is
                script directory.
        leap_update: If True, when the leap second file is more than 6 months old will attempt to download a new one.
                If False, then will just give a warning message. Default is True.

        predictive : If True, then instead of returning v_true, returns v_predicted.
        Default: False, and return is v_true from Wright and Eastman (2014)

        See OUTPUTs for description

    OUTPUTS:
            If predictive = False
                v_true: The true radial velocity of the SolSystemTarget for an observer at the observatory but in an inertial frame not moving.
                    If zmeas is included to show the measured absolute redshift for the SolSystemTarget as measured by an instrument,
                    then in this formulation, v_true will show the motion of the SolSystemTarget.
                    Else if predictive = True
                v_predicted: Ideal redshift measured for the SolSystemTarget from Earth for given location and time.
                    This output returns the theoretical prediction for the redshift which includes the barycentric component.
                    This will be the measurement of a noiseless RV instrument observing the SolSystemTarget.


    """


    # Convert times to obtain TDB and TT
    JDTDB, JDTT, warning, error = utc_tdb.JDUTC_to_JDTDB(JDUTC)


    # Need dictionary object for HORIZONS call
    longi = loc.lon.degree
    lat = loc.lat.degree
    alt = loc.height.value
    loc_dict = {'lon': longi,
                'lat': lat,
                'elevation': alt}

    # Reflecting object
    # First find the light travel time for Object with respect to observatory.

    try:
        TargetObj1 = Horizons(id=SolSystemTarget, location=loc_dict, epochs=JDTDB, id_type=HorizonsID_type).vectors()
    except ValueError:
        warning+= ['Unable to use Vector query for Horizons search using exact observatory coordinates, reverting to using Geocenter. ']
        TargetObj1 = Horizons(id=SolSystemTarget, location='399', epochs=JDTDB, id_type=HorizonsID_type).vectors()
    EarthTargetLightTravel = TargetObj1['lighttime'][0] #days

    # Subtract light time and find pos and vel for target wrt SSB
    TargetObjTime = JDTDB - (EarthTargetLightTravel)
    TargetObj2 = Horizons(id=SolSystemTarget, location='@0', epochs=TargetObjTime, id_type=HorizonsID_type)
    TargetVectors = TargetObj2.vectors(refplane='earth')
    TargetSSBLightTravel = TargetVectors['lighttime'][0] #days

    [TargetVectors[i].convert_unit_to('m') for i in ['x','y','z']] # Convert from AU to km
    [TargetVectors[i].convert_unit_to('m/s') for i in ['vx','vy','vz']] # Convert from AU/d to km/d
    PosVector_TargetSSB = TargetVectors['x','y','z'].as_array()[0].view((float,3))
    VelVector_TargetSSB = TargetVectors['vx','vy','vz'].as_array()[0].view((float,3))
    BetaTarget = VelVector_TargetSSB / c

    ## SUN
    # First find the light travel time for Sun with respect to Observatory
    #SolObj1 = Horizons(id='Sun', location='@0', epochs=JDTDB, id_type='majorbody').ephemerides()
    #SolSSBLightTravel = SolObj1['lighttime'][0] # Units of days


    # Ignoring difference in light travel time between the Sun and SSB for finding the Solar vectors
    # Subtract light time and find pos and vel for Sun wrt SSB
    TargetSolLightTravelDelay = EarthTargetLightTravel + TargetSSBLightTravel
    SolObj2 = Horizons(id='Sun', location='@0', epochs=JDTDB-(TargetSolLightTravelDelay), id_type='majorbody')
    SolVectors = SolObj2.vectors(refplane='earth')

    [SolVectors[i].convert_unit_to('m') for i in ['x','y','z']]
    [SolVectors[i].convert_unit_to('m/s') for i in ['vx','vy','vz']]
    PosVector_SolSSB = SolVectors['x','y','z'].as_array()[0].view((float,3))
    VelVector_SolSSB = SolVectors['vx','vy','vz'].as_array()[0].view((float,3))
    BetaSolar = VelVector_SolSSB / c
    GammaSolar = 1. / np.sqrt(1.- np.sum(BetaSolar**2))


    ##### NUTATION, PRECESSION, ETC. #####
    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris) # [km]
    PosVector_EarthSSB = r_eci + earth_geo[0].xyz.value*1000. # [m]
    v_geo = earth_geo[1].xyz.value*1000./86400.  # [m/s]
    VelVector_EarthSSB = (v_eci+v_geo) / (1.+ np.sum(v_eci*v_geo)/c**2) # [m/s]

    BetaEarth = VelVector_EarthSSB / c
    GammaEarth = 1. / np.sqrt(1. - np.sum(BetaEarth**2))


    PosVector_SolEarth, PosMag_SolEarth, PosHat_SolEarth = CalculatePositionVector(r1=PosVector_SolSSB, r2=PosVector_EarthSSB)
    PosVector_SolTarget, PosMag_SolTarget, PosHat_SolTarget = CalculatePositionVector(PosVector_SolSSB, PosVector_TargetSSB)
    PosVector_TargetEarth, PosMag_TargetEarth, PosHat_TargetEarth = CalculatePositionVector(PosVector_TargetSSB, PosVector_EarthSSB)

    # Calculates the the redshift experienced by photon due to Earth's + Sun's gravity well from infinity.
    # Correcting for this should give a blue shift.
    zGREarth =  - ac.G.value * ac.M_earth.value / ((ac.c.value**2)*(np.sqrt(np.sum(r_eci**2))))\
                - ac.G.value * ac.M_sun.value / ((ac.c.value**2)*(np.sqrt(np.sum(PosVector_SolEarth**2))))

    # To correct for redshift experienced by photon and place it in the SSB
    zGRSun =  - (ac.G.value * ac.M_sun.value) / ((ac.c.value**2) * (np.sqrt(np.sum(ac.R_sun.value**2))))

    zGR = (1+zGREarth)/(1+zGRSun)
    zSR = GammaSolar/GammaEarth
    zClassical1 = (1 + np.dot(BetaSolar, PosHat_SolTarget))/(1 + np.dot(BetaTarget, PosHat_SolTarget))
    zClassical2 =  (1 + np.dot(BetaTarget, PosHat_TargetEarth))/(1 + np.dot(BetaEarth, PosHat_TargetEarth))


    zpredicted = zGR*zSR*zClassical1*zClassical2 - 1

    zb = (GammaEarth * (1+np.dot(BetaEarth, PosHat_TargetEarth)) * \
        (1 + np.dot(BetaTarget, PosHat_SolTarget))/(1 + np.dot(BetaTarget, PosHat_TargetEarth)) \
        / (1+zGREarth)) - 1

    v_true = c * ((1.+zb)*(1.+ zmeas)-1.)  # [m/s]
    v_predicted = c * zpredicted # [m/s]

    if predictive:
        vel = v_predicted
    else:
        vel = v_true

    return vel, warning, error
