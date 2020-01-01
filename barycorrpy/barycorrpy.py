from __future__ import division
from __future__ import print_function
#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric
from astropy.time import Time
import math
import numpy as np
import os


from . import find_hip
from . import PINT_erfautils as PINT
from . import utc_tdb
from .SolarSystemBC import SolarBarycentricCorrection, ReflectedLightBarycentricCorrection
from .utils import flux_weighting, get_stellar_data, CalculatePositionVector
from .PhysicalConstants import *

### Need to install jplephem ###
#de430 is 100 MB in size


def get_BC_vel(JDUTC,
       starname= '', hip_id=None, ra=None, dec=None, epoch=None, pmra=None, pmdec=None, px=None, rv=None,
       obsname='', lat=0., longi=0., alt=0., zmeas=0.,
       ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True, predictive=False,
       SolSystemTarget=None, HorizonsID_type='smallbody'):
    '''
    Barycentric Velocity Correction at the 1 cm/s level, as explained in Wright & Eastman, 2014.
    Calling procedure for barycorrpy. Accepts vector time object (i.e., multiple observation JD values).

    Note: Can be used for stellar observations, or for observations of the Sun or other Solar System objects.

    INPUT:
        JDUTC : Can enter multiple times in Astropy Time object or as float. Will loop through and find barycentric velocity correction corresponding to those times.
                In UTC Scale. If using float, be careful about format and scale used.
        starname : Name of target. Will query SIMBAD database.
                                OR / AND
        hip_id : Hipparcos Catalog ID. (Integer). Epoch will be taken to be Catalogue Epoch or J1991.25
                If specified then ra,dec,pmra,pmdec,px, and epoch need not be specified.
                                OR / AND
        ra, dec : RA and Dec of star [degrees].
        epoch : Epoch of coordinates in Julian Date. Default is J2000 or 2451545.
        pmra : Proper motion in RA [mas/year]. Eg. PMRA = d(RA)/dt * cos(dec). Default is 0.
        pmdec : Proper motion in Dec [mas/year]. Default is 0.
        px : Parallax of target [mas]. Default is 0.
        rv : Radial Velocity of Target [m/s]. Default is 0. This is the bulk RV (systemic) of the target at the ~100 km/s precision.
             Can be ignored for most targets but for high velocity stars. This does not include the barycentric velocity and is
             only required to correct for proper motion and secular acceleration.

        obsname : Name of Observatory as defined in Astropy EarthLocation routine. Can check list by EarthLocation.get_site_names().
                  If obsname is not used, then can enter lat,long,alt.
                                OR
        lat : Latitude of observatory in [degrees]. North (+ve) and South (-ve).
        longi : Longitude of observatory [degrees]. East (+ve) and West (-ve).
        alt : Altitude of observatory [m].

        zmeas : Measured redshift (e.g., the result of cross correlation with template spectrum). Default is 0.
                The redshift measured by the spectrograph before any barycentric correction. Therefore zmeas includes the barycentric
                velocity of the observatory.
        ephemeris : Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                    For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                    ['de432s','de430',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir : Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is
                script directory.
        leap_update : If True, when the leap second file is more than 6 months old will attempt to download a new one.
                      If False, then will just give a warning message. Default is True.
        SolSystemTarget : When running barycentric correction for a stellar target, SolSystemTarget = None. Default value = None
                To correct for Solar RV observations set SolSystemTarget = 'Sun', for reflected light observations ...?
                For Solar observations:
                Only inputs required are -
                    1. JDUTC
                    2. loc
                    3. zmeas
                    4. ephemeris
                    5. leap_dir
                    6. leap_update
                    7. predictive

        predictive : If True, then instead of returning v_true, returns v_predicted.
                     Default: False, and return is v_true from Wright and Eastman (2014)

                     See OUTPUTs for description


    OUTPUT:
        FOR STELLAR OBSERVATIONS (not the Sun)
        If SolSystemTarget == None
            If predictive = False
                v_true : The barycenter-corrected RV [m/s] as defined in Wright & Eastman, 2014.
                NOTE: This is not just the barycentric velocity that can be subtracted directly from the measured RV.
                    The measured RV must be entered in the code as zmeas. This is because the relativistic cross product between zbary and zmeas is
                    required. This matters at ~ m/s level and hence must be included.

                The formula used is ztrue = ((1.+ zb)*(1.+ zmeas) - 1.)
                Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) x speed of light (c).
            Else if predictive = True
                v_predicted
                In this predictive mode, the expected RV shift measured by the instrument is output.
                NOTE: zpredicted = zmeas for an ideal noiseless observation, where the RV from the star (ztrue) is zero (or a flat line vs time)


            warning : Warning and Error message from the routine.
            status : Status regarding warning and error message. Returns the following -
                    0 - No warning or error.
                    1 - Warning message.
                    2 - Error message.

        FOR SOLAR OBSERVATIONS
        If SolSystemTarget == 'Sun'
            If predictive = False
                v_true: The true radial velocity of the Sun for an observer at the observatory but in an inertial frame not moving.
                    If zmeas is included to show the measured absolute redshift for the Sun as measured by an instrument,
                    then in this formulation, v_true will show the motion of the Sun,
                    which is mostly dominated by the synodic period of Jupiter as seen from Earth.
                    Our barycentric correction includes the gravitational redshift of the Sun, and therefore the true
                    velocity (or redshift) is centered at 0 m/s.
            Else if predictive = True
                v_predicted: Ideal redshift measured for the Sun from Earth for given location and time.
                    This output returns the theoretical prediction for the redshift which includes the barycentric component.
                    This will be the measurement of a noiseless RV instrument observing the Sun.
        Else:
            Computing the barycentric corrections for reflected light observations of a target in the Solar system
            If predictive = False
                v_true: The true radial velocity of the SolSystemTarget for an observer at the observatory but in an inertial frame not moving.
                    If zmeas is included to show the measured absolute redshift for the SolSystemTarget as measured by an instrument,
                    then in this formulation, v_true will show the motion of the SolSystemTarget.
            Else if predictive = True
                v_predicted: Ideal redshift measured for the SolSystemTarget from Earth for given location and time.
                    This output returns the theoretical prediction for the redshift which includes the barycentric component.
                    This will be the measurement of a noiseless RV instrument observing the SolSystemTarget.

            The formula used is ztrue = ((1.+zb)*(1.+zmeas)-1.)
            Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) x speed of light (c).

            warning : Warning and Error message from the routine.
            status : Status regarding warning and error message. Returns the following -
                    0 - No warning or error.
                    1 - Warning message.
                    2 - Error message.


    Example:
    >>> from astropy.time import Time
    >>> JDUTC = Time(2458000, format='jd', scale='utc')
    >>> get_BC_vel(JDUTC, hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9)
    (array([15403.95089287]), [[], []], 0)

    '''

    warning = []
    error = []
    status = 0

    # Check for JDUTC type
    if type(JDUTC)!=Time:
         warning += [['Warning: Float JDUTC entered. Verify time scale (UTC) and format (JD)']]
         JDUTC = Time(JDUTC, format='jd', scale='utc')

    if JDUTC.isscalar:
        JDUTC = Time([JDUTC])
    else:
       if np.size(zmeas)>1 and np.size(zmeas)!=np.size(JDUTC):
          error+= [['Error: Size mismatch. JDUTC is a vector, zmeas must also be a vector of same length corresponding to those dates']]
          raise IndexError('Error: Size mismatch. JDUTC is a vector, zmeas must be a vector of same length corresponding to those dates')

    if obsname:
        loc = EarthLocation.of_site(obsname)
        lat = loc.lat.value
        longi = loc.lon.value
        alt = loc.height.value
        warning += ['Warning: Taking observatory coordinates from Astropy Observatory database. Verify precision. Latitude = %f  Longitude = %f  Altitude = %f'%(lat,longi,alt)]
    else:
        loc = EarthLocation.from_geodetic(longi, lat, height=alt)

    ## STELLAR TARGET ##
    if SolSystemTarget is None:
        # Running correction for stellar observation
        star_par = {'ra':ra,'dec':dec,'pmra':pmra,'pmdec':pmdec,'px':px,'rv':rv,'epoch':epoch}
        star_simbad = {}
        star_hip = {}
        star_zero = {'ra':0.,'dec':0.,'pmra':0.,'pmdec':0.,'px':0.,'rv':0.,'epoch':2451545.0}
        star_output = {}

        if starname:
            star_simbad,warning1 = get_stellar_data(starname)
            warning += warning1
        if hip_id:
            if starname:
                warning += ['Warning: Querying SIMBAD and Hipparcos Catalogue']
            star_hip = find_hip(hip_id)

        star_output = star_simbad.copy()
        star_output.update({k:star_hip[k] for k in star_hip if star_hip[k] is not None})
        star_output.update({k:star_par[k] for k in star_par if star_par[k] is not None})
        star_output.update({k:star_zero[k] for k in star_output if star_output[k] is None})
        warning+=['Following are the stellar positional parameters being used - ',star_output]
        vel = []

        for jdutc,zm in zip(JDUTC,np.repeat(zmeas,np.size(JDUTC)/np.size(zmeas))):
            a = BCPy(JDUTC=jdutc,
                    zmeas=zm,
                    loc=loc,
                    ephemeris=ephemeris, leap_dir=leap_dir, leap_update=leap_update, predictive=predictive, **star_output)
            vel.append(a[0])
            warning.append(a[1])
            error.append(a[2])

    ## SOLAR OBSERVATIONS ##
    elif SolSystemTarget == 'Sun' or SolSystemTarget == 'Sol':
        vel = []
        for jdutc,zm in zip(JDUTC,np.repeat(zmeas,np.size(JDUTC)/np.size(zmeas))):
            a = SolarBarycentricCorrection(JDUTC=jdutc, loc=loc, zmeas=zm,
                    ephemeris=ephemeris, leap_dir=leap_dir, leap_update=leap_update, predictive=predictive)
            vel.append(a[0])
            warning.append(a[1])
            error.append(a[2])

    else:
        vel = []
        for jdutc,zm in zip(JDUTC,np.repeat(zmeas,np.size(JDUTC)/np.size(zmeas))):
            a = ReflectedLightBarycentricCorrection(SolSystemTarget=SolSystemTarget, JDUTC=jdutc, loc=loc, zmeas=zm, HorizonsID_type=HorizonsID_type,
                    ephemeris=ephemeris, leap_dir=leap_dir, leap_update=leap_update, predictive=predictive)
            vel.append(a[0])
            warning.append(a[1])
            error.append(a[2])
    # Status messages to check for warning or error
    if not all(vel): error += ['Check inputs. Error in code']
    if any(error):   status |= 2
    if any(warning): status |= 1
    # Convert velocity from list to numpy array
    vel = np.array(vel)

    return vel, warning+error, status






def BCPy(JDUTC,
    ra=0.0, dec=0.0, epoch=2451545.0, pmra=0.0, pmdec=0.0, px=0.0, rv=0.0, zmeas=0.0,
    loc=None,
    ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True, predictive=False):
    '''
    Barycentric Velocity Correction at the 1 cm/s level, as explained in Wright & Eastman, 2014.

    See get_BC_vel() for parameter description.

    '''

    # Convert times to obtain TDB and TT
    JDTDB, JDTT, warning, error = utc_tdb.JDUTC_to_JDTDB(JDUTC, fpath=leap_dir, leap_update=leap_update)

    ##### NUTATION, PRECESSION, ETC. #####

    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)

    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]

    ##### EPHEMERIDES #####

    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris) # [km]
    PosVector_EarthSSB = r_eci + earth_geo[0].xyz.value*1000. # [m]
    v_geo = earth_geo[1].xyz.value*1000./86400.  # [m/s]

    # Relativistic Addition of Velocities
    VelVector_EarthSSB = (v_eci+v_geo) / (1.+np.sum(v_eci*v_geo)/c**2) # [m/s]
    BetaEarth = VelVector_EarthSSB / c

    ##### Convert Star RA DEC to R0hat vector #####

    r0hat = np.array([math.cos(ra*np.pi/180.)*math.cos(dec*np.pi/180.),
                      math.sin(ra*np.pi/180.)*math.cos(dec*np.pi/180.),
                                              math.sin(dec*np.pi/180.)])
    # Eq 14 to 17
    up = [0., 0., 1.]
    east = np.cross(up, r0hat)
    east = east / math.sqrt(np.sum(east*east))
    north = np.cross(r0hat, east)
    mu = (pmra*east+pmdec*north)/pctoau/1000 # Divided by 1000 since the Proper motion is in milli-arcseconds.


    ##### Stellar position corrected for motion #####

    epoch0 = 2000. + (epoch-2451545.)/365.25
    yearnow = 2000. + (JDTDB.jd-2451545.)/365.25

    T = yearnow - epoch0                           # [yr]
    vpi = rv/1.e3 * kmstoauyr * (px/1.e3/pctoau)   # [rad/yr]
    vel = mu + vpi*r0hat                           # [rad/yr] (m in AA)
    r = r0hat + vel*T                              # [rad]    (p1 in AA)
    rhat = r / math.sqrt(np.sum(r*r))

    # Parallax correction #
    if px>0:
        rho = 1000.*rhat/px*pctoau - PosVector_EarthSSB/AU # [AU]
        rhohat = rho / math.sqrt(np.sum(rho*rho)) # Unitless
        r0 = 1000./px*pctoau*AU # [m]
        BetaStar = r0*mu/c/year + rv*r0hat/c

        zlighttravel = rv*r0*np.sum(mu*mu)*T/(year*c*c)

    else:
        rhohat = rhat
        BetaStar = [0., 0., 0.]
        zlighttravel = 0.


    ##### Calculate Gravitaional Redshift due to Solar system bodies #####

    Sum_GR = 0.
    zshapiro = 0.

    for ss_body in ss_bodies:
        if ss_body == 'Earth':
            PosVector_SSObject = earth_geo[0].xyz.value*1000. # [m]
        else:
            jplephem = get_body_barycentric(ss_body, JDTDB, ephemeris=ephemeris)
            PosVector_SSObject = jplephem.xyz.value*1000. # [m]

        # Vector from object barycenter to Observatory
        PosVector_EarthSSObject, PosMag_EarthSSObject, PosHat_EarthSSObject = CalculatePositionVector(r1=PosVector_EarthSSB, r2=PosVector_SSObject)

        # Add Shapiro Delay
        a = np.dot((rhohat-np.dot(PosHat_EarthSSObject,rhohat)*PosHat_EarthSSObject), BetaEarth)
        zshapiro += -2.*GM[ss_body]*a / ((c*c)*PosMag_EarthSSObject*(1+np.dot(PosHat_EarthSSObject, rhohat)))   # Eq 27

        if PosMag_EarthSSObject:
            Sum_GR += GM[ss_body] / PosMag_EarthSSObject  # [(m/s)^2]

    zgravity = 1./(1+Sum_GR/(c*c)) - 1


    ##### Determine the Barycentric RV correction (Eq 28) #####

    GammaEarth = 1. / math.sqrt(1.- np.sum(BetaEarth**2))

    zb = -1. - zshapiro - zlighttravel + \
        GammaEarth*(1+np.dot(BetaEarth, rhohat))*(1+np.dot(r0hat, BetaStar))/((1.+np.dot(BetaStar, rhohat))*(1.+zgravity)) # Eq 28

    if not predictive:
        v_final = c * ((1.+zb)*(1.+zmeas)-1.)  # [m/s]
    else:
        v_final = ((1/(1.+zb)) - 1) * c # [m/s]


    return v_final, warning, error


def exposure_meter_BC_vel(JDUTC, expmeterflux,
       starname = '', hip_id=None, ra=None, dec=None, epoch=None, pmra=None, pmdec=None, px=None, rv=None,
       obsname='', lat=0., longi=0., alt=0., zmeas=0.,
       ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True,
       SolSystemTarget=None, HorizonsID_type='smallbody', predictive=False):

    '''
    Calculate Barycentric velocity weighted by flux from exposure meter to account for long exposure time.
    Enter JDUTC and expmeterflux from exposure meter readings to calculate barycentric velocity correction for exposure.

    INPUT:
        JDUTC : Can enter multiple times in Astropy Time object. In UTC Scale.
                Can also accept list of float JDUTC times. Be cautious about scale used.
        expmeterflux : Array or List of exposure meter fluxes corresponding to each JDUTC.
                       The resultant barycentric correction will be calculated at each JDUTC and then weighted by the exposure meter fluxes. (See Landoni 2014)
        starname : Name of target. Will query SIMBAD database.
                                OR / AND
        hip_id : Hipparcos Catalog ID. (Integer) . Epoch will be taken to be Catalogue Epoch or J1991.25
                If specified then ra,dec,pmra,pmdec,px, and epoch need not be specified.
                                OR / AND
        ra, dec : RA and Dec of star [degrees].
        epoch : Epoch of coordinates in Julian Date. Default is J2000 or 2451545.
        pmra : Proper motion in RA [mas/year]. Eg. PMRA = d(RA)/dt * cos(dec). Default is 0.
        pmdec : Proper motion in Dec [mas/year]. Default is 0.
        px : Parallax of target [mas]. Default is 0.
        rv : Radial Velocity of Target [m/s]. Default is 0. This is the bulk RV (systemic) of the target at the ~100 km/s precision.
             Can be ignored for most targets but for high velocity stars. This does not include the barycentric velocity and is
             only required to correct for proper motion and secular acceleration.

        obsname : Name of Observatory as defined in Astropy EarthLocation routine. Can check list by EarthLocation.get_site_names().
                  If obsname is not used, then can enter lat,long,alt.
                                OR
        lat : Latitude of observatory in [degrees]. North (+ve) and South (-ve).
        longi : Longitude of observatory [degrees]. East (+ve) and West (-ve).
        alt : Altitude of observatory [m].


        zmeas : Measured redshift (e.g., the result of cross correlation with template spectrum). Default is 0.
                The redshift measured by the spectrograph before any barycentric correction. Therefore zmeas includes the barycentric
                velocity of the observatory.
        ephemeris : Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430.
                    For first use Astropy will download the Ephemeris ( for DE430 ~100MB). Options for ephemeris inputs are
                    ['de432s','de430',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                    'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        leap_dir : Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'. Default is
            script directory.
        leap_update : If True, when the leap second file is more than 6 months old will attempt to download a new one.
                      If False, then will just give a warning message. Default is True.
        SolSystemTarget : When running barycentric correction for a stellar target, Target = None. Default value = None
                To correct for Solar RV observations set Target = 'Sun', for reflected light observations ...?
                For Solar observations:
                Only inputs required are -
                    1. JDUTC
                    2. loc
                    3. zmeas
                    4. ephemeris
                    5. leap_dir
                    6. leap_update
                    7. predictive

        predictive : If True, then instead of returning v_true, returns v_predicted.
                Default: False, and return is v_true from Wright and Eastman (2014)

    See OUTPUTs for description
        FOR STELLAR OBSERVATIONS (not the Sun)
        If SolSystemTarget == None
            If predictive = False
            weighted_vtrue : The barycenter-corrected RV (m/s) for the exposure as weighted by the flux from the exposure meter as defined in Wright & Eastman, 2014.
                    NOTE: This is not just the barycentric velocity that can be subtracted directly from the measured RV.
                    The measured RV must be entered in the code as zmeas. This is because the relativistic cross product between zbary and zmeas is
                    required. This matters at ~ m/s level and hence must be included.

                The formula used is ztrue = ((1.+zb)*(1.+zmeas)-1.)
                Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) x speed of light (c).

            Else if predictive = True
                weighted_vpredicted : In this predictive mode, the expected RV shift measured by the instrument is output.
                NOTE: zpredicted = zmeas for an ideal noiseless observation, where the RV from the star (ztrue) is zero (or a flat line vs time)

            JDUTCMID : The flux weighted midpoint time is returned.
            warning : Warning and Error message from the routine.
            status : Status regarding warning and error message. Returns the following -
                    0 - No warning or error.
                    1 - Warning message.
                    2 - Error message.

        FOR SOLAR OBSERVATIONS
        If SolSystemTarget == 'Sun'
            If predictive = False
            weighted_v_true: Flux weighted true radial velocity of the Sun for an observer at the observatory but in an inertial frame not moving.
                    If zmeas is included to show the measured absolute redshift for the Sun as measured by an instrument,
                    then in this formulation, v_true will show the motion of the Sun,
                    which is mostly dominated by the synodic period of Jupiter as seen from Earth.
            The formula used is ztrue = ((1.+zb)*(1.+zmeas)-1.)
            Therefore if zmeas is set to 0, then ztrue = zb. The velocities are just the redshift (z) x speed of light (c).

            Else if predictive = True
            weighted_v_predicted: Flux weighted ideal redshift measured for the Sun from Earth for given location and time.
                    This output returns the theoretical prediction for the redshift which includes the barycentric component.
                    This will be the measurement of a noiseless RV instrument observing the Sun.

            JDUTCMID : The flux weighted midpoint time is returned.
            warning : Warning and Error message from the routine.
            status : Status regarding warning and error message. Returns the following -
                    0 - No warning or error.
                    1 - Warning message.
                    2 - Error message.


    '''
    expmeterflux = np.array(expmeterflux)
    error = []

    if isinstance(JDUTC,(int,float)):
        JDUTC = np.array([JDUTC])

    ## Check for size of Flux Array ##
    if len(JDUTC)!=len(expmeterflux):
        print('Error: Size of JDUTC array is not equal to expmeterflux (Flux) array')
        error += ['Error: Size of JDUTC array is not equal to expmeterflux (Flux) array']

    ####################
    ## STELLAR TARGET ##
    if SolSystemTarget is None:
        # Running correction for stellar observation
        ## Calculate barycentric velocity at each instance of exposure meter reading ##
        vel,warning,status = get_BC_vel(JDUTC=JDUTC,
                                    starname=starname,hip_id=hip_id,ra=ra,dec=dec,epoch=epoch,pmra=pmra,pmdec=pmdec,px=px,
                                    obsname=obsname,lat=lat,longi=longi,alt=alt,
                                    rv=rv,zmeas=zmeas,ephemeris=ephemeris,leap_dir=leap_dir,leap_update=leap_update,
                                    SolSystemTarget=SolSystemTarget, HorizonsID_type=HorizonsID_type, predictive=predictive)



    ## Weight it by flux ##

    weighted_vel = flux_weighting(flux = expmeterflux, qty = vel)
    try:
        JDUTC = JDUTC.jd
    except:
        JDUTC = np.array(JDUTC)

    JDUTCMID = flux_weighting(flux=expmeterflux, qty=JDUTC)

    # Error handling
    if error:
        warning.append(error)
        status = 2

    return weighted_vel, JDUTCMID, warning, status
