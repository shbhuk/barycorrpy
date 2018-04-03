from __future__ import division
from __future__ import print_function
#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import get_body_barycentric_posvel, get_body_barycentric
from astropy.time import Time
import math
import scipy.constants as sc
import astropy.constants as ac
import numpy as np
import os


from . import find_hip
from . import PINT_erfautils as PINT
from . import utc_tdb
from .utils import flux_weighting,get_stellar_data

### Need to install jplephem ###
#de430 is 100 MB in size

# Parsing constants #
AU = sc.astronomical_unit # [m]
c = sc.c # Speed of light [m/s]
pctoau = 3600*180/np.pi # No. of AU in one parsec
year = 365.25*3600*24 # [s]
kmstoauyr = year/(1000*AU)

M_moon = 73476730924573500000000 # Mass of the Moon [kg]

# Mass and normalised mass of solar system bodies
ss_bodies = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Moon']
M = dict(zip(ss_bodies, [ac.M_sun.value, 0.3301e24, 4.867e24, ac.M_earth.value, 0.6417e24, ac.M_jup.value, 568.5e24, 86.82e24, 102.4e24, M_moon])) # [kg]

#Reduced number of bodies, to increase speed. Does not affect precision at 1 cm/s level
ss_bodies = ['Sun','Earth','Jupiter','Saturn']
M = dict(zip(ss_bodies, [ac.M_sun.value, ac.M_earth.value,ac.M_jup.value,568.5e24])) # [kg]

GM = {k:sc.G*M[k] for k in ss_bodies}

def get_BC_vel(JDUTC,
       starname = '', hip_id=None, ra=None, dec=None, epoch=None, pmra=None, pmdec=None, px=None, rv=None,
       obsname='', lat=0., longi=0., alt=0., zmeas=0.,
       ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True):
    '''
    Barycentric Velocity Correction at the 1 cm/s level, as explained in Wright & Eastman, 2014.
    Calling procedure for barycorrpy. Accepts vector time object (i.e., multiple observation JD values).
    
    INPUT:
        JDUTC : Can enter multiple times in Astropy Time object or as float. Will loop through and find barycentric velocity correction corresponding to those times. 
                In UTC Scale. If using float, be careful about format and scale used.
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
    
    OUTPUT:
        vel : The barycenter-corrected RV [m/s] as defined in Wright & Eastman, 2014. 
        NOTE: This is not just the barycentric velocity that can be subtracted directly from the measured RV.
            The measured RV must be entered in the code as zmeas. This is because the relativistic cross product between zbary and zmeas is
            required. This matters at ~ m/s level and hence must be included.
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
    
    vel = []
    warning = []
    error = []
    status = 0
    
    # Check for JDUTC type   
    if type(JDUTC)!=Time:
         warning += [['Warning: Float JDUTC entered. Verify time scale (UTC) and format (JD)']]
         JDUTC=Time(JDUTC, format='jd', scale='utc')
       
    if JDUTC.isscalar:
        JDUTC = Time([JDUTC])
    else:
       if np.size(zmeas)>1 and np.size(zmeas)!=np.size(JDUTC):
          error+= [['Error: Size mismatch. JDUTC is a vector, zmeas must also be a vector of same length corresponding to those dates']]
          raise IndexError('Error: Size mismatch. JDUTC is a vector, zmeas must be a vector of same length corresponding to those dates')

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
           
    if obsname:
        loc = EarthLocation.of_site(obsname)
        lat = loc.lat.value
        longi = loc.lon.value
        alt = loc.height.value
        warning += ['Warning: Taking observatory coordinates from Astropy Observatory database. Verify precision. Latitude = %f  Longitude = %f  Altitude = %f'%(lat,longi,alt)]
    else: 
        loc = EarthLocation.from_geodetic(longi, lat, height=alt)
        

    for jdutc,zm in zip(JDUTC,np.repeat(zmeas,np.size(JDUTC)/np.size(zmeas))):
        a = BCPy(JDUTC=jdutc,
                 zmeas=zm, 
                 loc=loc,
                 ephemeris=ephemeris, leap_dir=leap_dir, leap_update=leap_update,**star_output)
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
    ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True):
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
    r_obs = r_eci + earth_geo[0].xyz.value*1000. # [m]
    v_geo = earth_geo[1].xyz.value*1000./86400.  # [m/s]
    
    # Relativistic Addition of Velocities
    v_obs = (v_eci+v_geo) / (1.+v_eci*v_geo/c**2) # [m/s]
    beta_earth = v_obs / c

    ##### Convert Star RA DEC to R0hat vector #####
    
    r0hat = np.array([math.cos(ra*np.pi/180.)*math.cos(dec*np.pi/180.),
                      math.sin(ra*np.pi/180.)*math.cos(dec*np.pi/180.),
                                              math.sin(dec*np.pi/180.)])
    # Eq 14 to 17
    up = [0., 0., 1.]
    east = np.cross(up, r0hat)
    east = east / math.sqrt(sum(east*east))
    north = np.cross(r0hat, east)
    mu = (pmra*east+pmdec*north)/pctoau/1000 # Divided by 1000 since the Proper motion is in milli-arcseconds.
    
    
    ##### Stellar position corrected for motion #####
    
    epoch0 = 2000. + (epoch-2451545.)/365.25
    yearnow = 2000. + (JDTDB.jd-2451545.)/365.25
    
    T = yearnow - epoch0                           # [years]
    vpi = rv/1.e3 * kmstoauyr * (px/1.e3/pctoau)   # [rad/yr]
    vel = mu + vpi*r0hat                           # [rad/yr] (m in AA)
    r = r0hat + vel*T                              # [rad]    (p1 in AA)
    rhat = r / math.sqrt(sum(r*r))
    
    # Parallax correction #
    if px>0:
        rho = 1000.*rhat/px*pctoau - r_obs/AU # [AU]
        rhohat = rho / math.sqrt(sum(rho*rho)) # Unitless
        r0 = 1000./px*pctoau*AU # [m]
        beta_star = r0*mu/c/year + rv*r0hat/c
        
        zlighttravel = rv*r0*sum(mu*mu)*T/(year*c*c)
    
    else:
        rhohat = rhat
        beta_star = [0., 0., 0.]
        zlighttravel = 0.

    
    ##### Calculate Gravitaional Redshift due to Solar system bodies #####
        
    Sum_GR = 0.
    zshapiro = 0.
    
    for ss_body in ss_bodies:
        if ss_body == 'Earth':
            pos_obj = earth_geo[0].xyz.value*1000. # [m]
        else:
            jplephem = get_body_barycentric(ss_body, JDTDB, ephemeris=ephemeris)
            pos_obj = jplephem.xyz.value*1000. # [m]
        
        # Vector from object barycenter to Observatory
        X = np.array(r_obs-pos_obj)
        Xmag = math.sqrt(sum(X*X)) # [m]
        Xhat = X / Xmag # Unitless
        
        # Add Shapiro Delay
        a = np.dot((rhohat-np.dot(Xhat,rhohat)*Xhat), beta_earth)
        zshapiro += -2.*GM[ss_body]*a / ((c*c)*Xmag*(1+np.dot(Xhat, rhohat)))   # Eq 27
        
        if Xmag!=0.:
            Sum_GR += GM[ss_body] / Xmag  # [(m/s)^2]
        
    zgravity = 1./(1+Sum_GR/(c*c)) - 1
    
    
    ##### Determine the Barycentric RV correction (Eq 28) #####
    
    gamma_earth = 1. / math.sqrt(1.-sum(beta_earth**2))
    
    zb = -1. - zshapiro - zlighttravel + gamma_earth*(1+np.dot(beta_earth, rhohat))*(1+np.dot(r0hat, beta_star))/((1.+np.dot(beta_star, rhohat))*(1.+zgravity)) # Eq 28
    v_final = c * ((1.+zb)*(1.+zmeas)-1.)  # [m/s]
    
    
    return v_final, warning, error

 
def exposure_meter_BC_vel(JDUTC,expmeterflux,
       starname = '', hip_id=None, ra=None, dec=None, epoch=None, pmra=None, pmdec=None, px=None, rv=None,
       obsname='', lat=0., longi=0., alt=0., zmeas=0.,
       ephemeris='de430', leap_dir=os.path.join(os.path.dirname(__file__),'data'), leap_update=True):
       
    '''
    Calculate Barycentric velocity weighted by flux from exposure meter to account for long exposure time. 
    Enter JDUTC and expmeterflux from exposure meter readings to calculate barycentric velocity correction for exposure. 
    
    INPUT: 
        JDUTC : Can enter multiple times in Astropy Time object.In UTC Scale. 
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

       
    OUTPUT:
        weighted_vel : The barycenter-corrected RV (m/s) for the exposure as weighted by the flux from the exposure meter as defined in Wright & Eastman, 2014. 
                NOTE: This is not just the barycentric velocity that can be subtracted directly from the measured RV.
                The measured RV must be entered in the code as zmeas. This is because the relativistic cross product between zbary and zmeas is
                required. This matters at ~ m/s level and hence must be included.
        JDUTCMID : The flux weighted midpoint time is returned. 
        warning : Warning and Error message from the routine
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
        error+=['Error: Size of JDUTC array is not equal to expmeterflux (Flux) array']

    
    ## Calculate barycentric velocity at each instance of exposure meter reading ##        
    vel,warning,status = get_BC_vel(JDUTC=JDUTC,
                                starname=starname,hip_id=hip_id,ra=ra,dec=dec,epoch=epoch,pmra=pmra,pmdec=pmdec,px=px,
                                obsname=obsname,lat=lat,longi=longi,alt=alt,
                                rv=rv,zmeas=zmeas,ephemeris=ephemeris,leap_dir=leap_dir,leap_update=leap_update)   

    ## Weight it by flux ##

    weighted_vel = flux_weighting(flux = expmeterflux, qty = vel)
    try:
        JDUTC = JDUTC.jd
    except:
        JDUTC = np.array(JDUTC)
        
    JDUTCMID = flux_weighting(flux = expmeterflux, qty = JDUTC)
    
    # Error handling
    if error:
        warning.append(error)
        status = 2

    return weighted_vel, JDUTCMID, warning, status




