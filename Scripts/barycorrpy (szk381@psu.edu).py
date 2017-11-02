from __future__ import division
#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time
import datetime
import math
import scipy.constants as const
import astropy.constants as u
import numpy as np
import os
import sys
import inspect

from read_HIP import find_hip
import PINT_erfautils as PINT
import utc_tdb

### Need to install jplephem ### 
#de430 is 100 MB in size


#location=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))   
#sys.path.append(location)
from Eastman_applet import bvc


'''
ra, dec - In degrees
JDUTC - Julian Date in UTC, eg 2450000.0

'''
#For Tau Ceti

JDUTC = Time(datetime.datetime.utcnow(),format='datetime',scale='utc')
#leap_dir=os.getcwd()+'/Box Sync/BaryCorr/barycorrpy/'
leap_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'\\'
#JDUTC=Time(2458000,format='jd',scale='utc')

ra=1.734757638888889*15
ra=26.0213645867
dec=-15.9395557246

obsname='CTIO'
lat=-30.169283
longi=-70.806789
alt=2241.9

epoch = 2451545.0 # Default is J2000 - Julian Day January 1st 2000, 12 noon

pmra = -1721.05 # Proper motion in RA in mas/year (scalar). Eg. PMRA = d(RA)/dt * cos(dec)
pmdec = 854.16 # Proper motion in dec, in mas(year)
px = 273.96 # Parallax in mas 
rv = 0.0 # RV of Target in m/s
zmeas=0.0

ephem=['de432s','de430','https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp','https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'
#ephemeris='de430'
hip_id=8102


def call_BCPy(JDUTC,hip_id=0,ra=0.0,dec=0.0,obsname='',lat=0.0,longi=0.0,alt=0.0,epoch=2451545.0,pmra=0.0,pmdec=0.0,px=0.0,rv=0.0,zmeas=0.0,ephemeris='de430',leap_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/', leap_update = True):
    '''
    INPUTS:
        See BCPy()
        JDUTC : Can enter multiple times in Astropy Time object. Will loop through and find barycentric velocity correction corresponding to those times.
        hip_id : Hipparcos Catalog ID. (Integer) 
                 If specified then ra,dec,pmra,pmdec,px, and epoch need not be specified. Epoch will be taken to be Catalogue Epoch or J1991.25
    
    
    '''
    
    if (type(hip_id) == int) and (hip_id > 0):
        _,ra,dec,px,pmra,pmdec,epoch = find_hip(hip_id)
        print 'Reading from Hipparcos Catalogue'
    
    vel=[]
    for i in range(0,np.size(JDUTC)):               
        vel.append(BCPy(JDUTC=JDUTC[i],ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_dir=leap_dir,leap_update=leap_update))

    return vel





def BCPy(JDUTC,ra=0.0,dec=0.0,obsname='',lat=0.0,longi=0.0,alt=0.0,epoch=2451545.0,pmra=0.0,pmdec=0.0,px=0.0,rv=0.0,zmeas=0.0,ephemeris='de430',leap_dir=os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+'/',leap_update = True ) :
    '''
    Barycentric Velocity Correction at the 1 cm/s level, as defined in Wright & Eastman, 2014.
    
    INPUTS:
        JDUTC : Astropy Time Object 
        
        All subsequent inputs are SCALARS
        
        ra , dec : RA and Dec of object in DEGREES
        
        obsname : Name of Observatory as defined in Astropy EarthLocation routine. Can check list by EarthLocation.get_site_names(). 
                  If observatory is not included in Astropy, then can enter lat,long,alt.
                                OR 
        lat : Latitude of observatory in degrees. North (+ve) and South (-ve)
        longi : Longitude of observatory in DEGREES. East (+ve) and West (-ve)
        alt : Altitude of observatory in METERS
        
        epoch : Epoch of coordinates in Julian Date. Default is J2000 or 2451545.0
        pmra : Proper motion in RA, in MAS/YEAR. Eg. PMRA = d(RA)/dt * cos(dec). Default is 0.0
        pmdec : Proper motion in Dec, in MAS/YEAR. Default is 0.0
        px : Parallax of target in MAS. Default is 0.0
        rv : Radial Velocity of Target in M/S. Default is 0.0
        zmeas : Measured redshift (e.g., the result of cross correlation with template spectrum). Default is 0.0
        ephemeris : Name of Ephemeris to be used. List of Ephemeris as queried by jplephem. Default is DE430. 
                    For first use Astropy will involve the Ephemeris (DE430 ~100MB)
        
        ephem=['de432s','de430',
                'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp',
                'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp']
        
        leap_dir : Directory where leap seconds file will be saved and maintained (STRING). Eg. '/Users/abc/home/savehere/'
        leap_update : If True, when the leap second file is more than 6 months old will attempt to download a new one.
                    If False, then will just give a warning message. Default is True.
         
    OUTPUTS:
        The barycenter-corrected RV (M/S) as defined in Wright & Eastman, 2014.
    
    '''
    
    AU=const.astronomical_unit # 1 AU in meters
    c=const.c # Speed of light in m/s
    pctoau=3600*180/np.pi #No. of AU in one parsec
    year = 365.25*3600*24 # 1 year in seconds
    kmstoauyr = year/(1000*AU)
    
    M_moon=73476730924573500000000 # Mass of the Moon in kg
    M_earth=u.M_earth.value
    M_sun=u.M_sun.value
    
    # Sun,  Mercury , Venus, Earth, Mars, Jupiter, Saturn , Uranus , Neptune, Moon
    GM=[const.G*x for x in [M_sun,0.3301e24,4.867e24,M_earth,0.6417e24,u.M_jup.value,568.5e24,86.82e24,102.4e24,M_moon]]
    
    # Convert times to obtain TDB and TT 
    JDTDB,JDTT=utc_tdb.JDUTC_to_JDTDB(JDUTC,fpath=leap_dir,leap_update=leap_update)

    
    ##### OBSERVATORY EUCLIDEAN COORDINATES #####       
    
    if len(obsname)!=0:
        loc=EarthLocation.of_site(obsname)
        lat=loc.lat.value  # Only need for applet. Can remove ########3
        longi=loc.lon.value
        alt=loc.height.value
    else:       
        loc=EarthLocation.from_geodetic(longi,lat,height=alt)
    R_ITRF=loc.value
    
    ##### NUTATION , PRECESSION , ETC. #####
    
    r_pint,v_pint=PINT.gcrs_posvel_from_itrf(loc,JDUTC)
    
    r_eci=r_pint[0]  # meters
    v_eci=v_pint[0]  # meters/second
    
    ##### EPHEMERIDES #####
            
    earth_geo=get_body_barycentric_posvel('earth',JDTDB,ephemeris=ephemeris)
    r_obs=r_eci+earth_geo[0].xyz.value*1000. # meters
    v_geo=earth_geo[1].xyz.value*1000./(86400.)  # meters/second
    
    # Relativistic Addition of Velocities
    v_obs=(v_eci+v_geo)/(1.+v_eci*v_geo/c**2) # m/s
    beta_earth=v_obs/c
   
    
    ##### Convert Star RA DEC to R0hat vector #####
    
    r0hat=np.array([math.cos(ra*np.pi/180.)*math.cos(dec*np.pi/180.),math.sin(ra*np.pi/180.)*math.cos(dec*np.pi/180.),math.sin(dec*np.pi/180.)])
    
    up = [0.0,0.0,1.0] 
    east = np.cross(up,r0hat)
    east = east/math.sqrt(sum(east*east))
    north = np.cross(r0hat,east)
    mu =  ((pmra*east+pmdec*north)/pctoau/1000) # Divided by 1000 since the Proper motion is in milli-arcseconds.
    
    
    ##### Stellar position at each time #####
    
    epoch0 = 2000.0+(epoch-2451545.0)/365.25
    yearnow = 2000.0+(JDTDB.jd-2451545.0)/365.25
    
    T = (yearnow-epoch0)                           #years
    vpi = rv/1.e3 * kmstoauyr * (px/1.e3/pctoau)     # rad/yr
    vel = mu + vpi*r0hat                           # rad/yr (m in AA)
    r = r0hat + vel*T                              # rad    (p1 in AA)
    rhat= r/math.sqrt(sum(r*r))
    
    # Parallax correction 
    if px>0:
        rho=1000.0*rhat/px*pctoau-r_obs/AU # In AU
        rhohat = rho/math.sqrt(sum(rho*rho)) # Unitless
        r0=1000.0/px*pctoau*AU # In meters
        beta_star = r0*mu/c/year + rv*r0hat/c
        
        zlighttravel = rv*r0*sum(mu*mu)*T/(year*c*c)
    
    else:
        rhohat=rhat
        beta_star=[0.0,0.0,0.0]
        zlighttravel = 0.0
    
    
    ##### Calculate Gravitaional Redshift due to Solar system bodies #####
    
    bodies=solar_system_ephemeris.bodies
    
    ss_bodies=[bodies[1],bodies[3],bodies[4],bodies[0],bodies[6],bodies[7],bodies[8],bodies[9],bodies[10],bodies[2]]   # Sun, Mercury , Venus, Earth, Mars, Jupiter, Saturn , Uranus , Neptune, Moon
    
    Sum_GR=0.0
    zshapiro=0.0
    
    for i in range(0,len(ss_bodies)-3):
        jplephem=get_body_barycentric_posvel(ss_bodies[i],JDTDB,ephemeris=ephemeris)
        pos_obj=jplephem[0].xyz.value*1000. # meters
        vel_obj=jplephem[1].xyz.value*1000./(86400.) # meters/second
        
        # Vecgtor from object barycenter to Observatory
        X=np.array(r_obs-pos_obj)
        Xmag=math.sqrt(sum(X*X)) # In meters
        Xhat=X/(Xmag) # Unitless
        
        # Add Shapiro Delay
        a=np.dot((rhohat-np.dot(Xhat,rhohat)*Xhat),beta_earth)
        zshapiro+= -2.0*GM[i]*a/((c*c)*(Xmag*(1+np.dot(Xhat,rhohat))))
        
        if Xmag!=0.0:
            Sum_GR+=GM[i]/Xmag  # (m/s)^2    
        
    zgravity=1.0/(1+Sum_GR/(c*c))-1
    
    
    ##### Determine the Barycentric RV correction (eq 28) #####
    
    gamma_earth=1./math.sqrt(1.-sum(beta_earth**2))
    

    zb = -1.0 - zshapiro - zlighttravel + gamma_earth*(1+np.dot(beta_earth,rhohat))*(1+np.dot(r0hat,beta_star))/((1+np.dot(beta_star,rhohat))*(1+zgravity)) # Eq 28
    v_final=c*((1.0+zb)*(1.0+zmeas)-1.0)

    ##### Call Eastman applet to compare #####
    res = bvc(jd_utc=JDUTC.jd, ra=ra, dec=dec, lat=lat, lon=longi, elevation=alt,pmra=pmra,pmdec=pmdec,parallax=px,rv=rv,zmeas=zmeas, epoch=epoch)
    
    return v_final,res

b=BCPy(JDUTC=JDUTC,ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_dir=leap_dir,leap_update=False)
print b





 
    
    
