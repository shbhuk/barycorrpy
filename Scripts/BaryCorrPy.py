#import de423    #https://pypi.python.org/pypi/jplephem
from astropy.coordinates import EarthLocation
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
from astropy.time import Time
import datetime
from astropy.utils import iers
import math
import scipy.constants as const
import astropy.constants as u
import numpy as np
import os
import sys
#import matplotlib.pyplot as plt




location=os.getcwd()+'\\Box Sync\\BaryCorr\\'
sys.path.append(location)
import utc_tdb

'''
ra, dec - In degrees
JDUTC - Julian Date in UTC, eg 2450000.0

'''
#For Tau Ceti


ra=1.734757638888889*15
dec=-15.93955572

lat=31.9599
longi=-111.5997
alt=2096

R_ITRF =[1814985.3, -5213916.8, -3187738.1] # Needed only if Lat,Long and Alt are not specified. Values for CTIO # In meters

epoch = 2451545.0 # Default is J2000 - Julian Day January 1st 2000, 12 noon
tbase = 0 # Baseline subtracted from JD_UTCs for higher precision times

pmra = -1721.05 # Proper motion in RA in mas/year (scalar). Eg. PMRA = d(RA)/dt * cos(dec)
pmdec = 854.16 # Proper motion in dec, in mas(year)

px = 273.96 # Parallax in mas 
rv = 00 # RV of Target in m/s









AU=const.astronomical_unit # 1 AU in meters
c=const.c # Speed of light in m/s
pctoau=3600*180/np.pi #No. of AU in one parsec
year = 365.25*3600*24 # 1 year in seconds
kmstoauyr = year/(1000*AU)

GMsun=u.M_sun*const.G # G*M-sun   
M_moon=73476730924573500000000 # Mass of the Moon in kg
M_earth=u.M_earth.value
M_sun=u.M_sun.value

# ? , Mercury , Venus, Earth, Mars, Jupiter, Saturn , Uranus , Neptune, Moon, Sun
GM=[const.G*x for x in [0.3301e24,4.867e24,M_earth,0.6417e24,u.M_jup.value,568.5e24,86.82e24,102.4e24,M_moon,M_sun]]


'''
X=328900.56 # M_sun/(M_earth+M_moon)
X=u.M_sun.value/(u.M_earth.value+7.3476E22) # M_sun/(M_earth+M_moon)
Y=81.30818226359628 # M_earth/M_moon
Y=(u.M_earth.value/7.3476E22) # M_earth/M_moon
#GM=[GMsun/x for x in [np.nan,6023600.0, 408523.71,(X*(Y+1.))/Y,3098708.0,1047.3486,3497.898,22902.98,19412.24,(X*(Y+1.)),1]] #Excluding pluto. GM/R
'''


t = Time("2014-09-22 23:22")
t=Time(datetime.datetime.now())

JDUTC = Time(datetime.datetime.utcnow()).jd

JDTDB,TDB,JDTT=utc_tdb.JDUTC_to_JDTDB(JDUTC)
ntimes=np.size(JDUTC)

dpi=np.pi/180.

    
##### OBSERVATORY EUCLIDEAN COORDINATES ######        
                
sites=EarthLocation.get_site_names()


Lat=lat*dpi
Long=(360-longi)*dpi
F=1.0/298.257223563
CC=1/math.sqrt(math.cos(Lat)**2+((1.0-F)*math.sin(Lat))**2)
S=((1.0-F)**2)*CC
A=6378137.0


R_ITRF=[((A*CC+alt)*math.cos(Lat)*math.cos(Long)),((A*CC+alt)*math.cos(Lat)*math.sin(Long)),(A*S+alt)*math.sin(lat)] # meters

EarthLocation.of_site('Kitt Peak')


######### NUTATION , PRECESSION , ETC. #############


iers.IERS_A_README


r_eci=[10,10,10] # meters 
v_eci=[0,0,0] # meters/second



######### EPHEMERIDES ######## 





#eph=Ephemeris(de423)
#print eph.names
#print eph.position_and_velocity('earthmoon', t.jd)  # JD date
#solar_system_ephemeris.set('de423') 
#print str(get_body_barycentric_posvel('earth',Time(TDB),ephemeris='de430'))  + '\t de'+str(430)+' ephemeris \n'  # Velocity in km / d
#print str(get_body_barycentric_posvel('earth',Time(TDB),ephemeris='de432s'))  + '\t de'+str(432)+'s ephemeris \n'
print str(get_body_barycentric_posvel('earth',Time(TDB),ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp')) + '\t de'+str(405)+ ' ephemeris \n'
#print str(get_body_barycentric_posvel('earth',Time(TDB),ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de423_for_mercury_and_venus/de423.bsp')) + '\t de'+str(423)+ ' ephemeris \n'
#print get_body_barycentric_posvel('earth',Time(TDB))  # Velocity in AU/d
    
earth_geo=get_body_barycentric_posvel('earth',Time(TDB),ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp')
r_obs=r_eci+earth_geo[0].xyz.value*1000. # meters
v_geo=earth_geo[1].xyz.value*1000./(86400.)  # meters/second

# Relativistic Addition of Velocities

v_obs=(v_eci+v_geo)/(1.+v_eci*v_geo/c**2)
beta_earth=v_obs/c



######## Convert Star RA DEC to R0hat ####



r0hat= np.tile(np.array([math.cos(ra*np.pi/180.)*math.cos(dec*np.pi/180.),math.sin(ra*np.pi/180.)*math.cos(dec*np.pi/180.),math.sin(dec*np.pi/180.)]),(ntimes,1))[0]

up = [0.0,0.0,1.0] 
east = np.cross(up,r0hat)
east = east/math.sqrt(sum(east*east))
north = np.cross(r0hat,east)
mu =  ((pmra*east+pmdec*north)/pctoau/1000) # Divided by 1000 since the Proper motion is in milli-arcseconds.

#### Stellar position at each time

epoch0 = 2000+(epoch-2451545.0)/365.25
yearnow = 2000+(JDTDB+tbase-2451545.0)/365.25

T = (yearnow-epoch0)                           #years
vpi = rv/1.e3 * kmstoauyr * (px/1e3/pctoau)     # rad/yr
vel = mu + vpi*r0hat                           # rad/yr (m in AA)
r = r0hat + vel*T                              # rad    (p1 in AA)
rhat= r/math.sqrt(sum(r*r))


px=1000    # DUMMY

if px>0:
    rho=1000.0*rhat/px*pctoau-r_obs/AU # In AU
    rhohat = rho/math.sqrt(sum(rho*rho)) # Unitless
    r0=1000.0/px*pctoau*AU # In meters
    beta_star = r0*mu/c/year + rv*r0hat/c
    
    zlighttravel = rv*r0*sum(mu*mu)*T/(c*c)

else:
    rhohat=rhat
    beta_star=0.0
    zlightravel = 0.0


## Calculate Gravitaional Redshift due to SOlar system bodies

bodies=solar_system_ephemeris.bodies

calcndx=[bodies[1],bodies[0],bodies[7],bodies[8],bodies[2],bodies[4],bodies[9],bodies[10],bodies[6],bodies[3]]   # Sun, Earth, Jupiter, Saturn, Moon, Venus, Uranus, Neptune. #Excluding pluto
ss_bodies=[bodies[3],bodies[4],bodies[0],bodies[6],bodies[7],bodies[8],bodies[9],bodies[10],bodies[2],bodies[1]]   # Mercury , Venus, Earth, Mars, Jupiter, Saturn , Uranus , Neptune, Moon, Sun

Sum_GR=0.0
zshapiro=0.0
for i in range(0,len(ss_bodies)-1):
    jplephem=get_body_barycentric_posvel(ss_bodies[i],Time(TDB),ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp')
    pos_obj=jplephem[0].xyz.value*1000. # meters
    vel_obj=jplephem[1].xyz.value*1000./(86400.) # meters/second
    
    
    X=np.array(r_obs-pos_obj)
    Xmag=math.sqrt(sum(X*X)) # In meters
    Xhat=X/(Xmag) # Unitless
    
    # Add Shapiro Delay
    a=np.dot((rhohat-np.dot(Xhat,rhohat)*Xhat),beta_earth)
    zshapiro+= -2.0*GM[i]*(beta_earth*a)/((c*c)*(Xmag*(1+np.dot(Xhat,rhohat))))
    
    if Xmag!=0.0:
        Sum_GR+=GM[i]/Xmag  # (m/s)^2
    
zgravity=1.0/(1+Sum_GR/(c*c))-1

#determine the Barycentric RV correction (eq 28)

gamma_earth=1./math.sqrt(1.-sum(beta_earth**2))
zb = -1.0 - zshapiro - zlighttravel + gamma_earth*(1+np.dot(beta_earth,rhohat))*(1+np.dot(r0hat,beta_star))/((1+np.dot(beta_star,rhohat))*(1+zgravity))


    
    
    
