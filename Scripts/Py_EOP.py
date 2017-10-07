from astropy.coordinates import EarthLocation
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time
import datetime
from astropy.utils import iers
import math
import scipy.constants as const
import astropy.constants as u
import numpy as np
import os
import sys
import datetime
import matplotlib.pyplot as plt
sys.path.append(os.getcwd()+'\\Box Sync\\BaryCorr\\')

import utc_tdb

mjdutc=Time(datetime.datetime.utcnow()).mjd
mjdutc=2458000 - 2400000.5

jdtt=utc_tdb.JDUTC_to_JDTDB(mjdutc+2400000.5)[2]

jason_leap=np.loadtxt(os.getcwd()+'\\Box Sync\\BaryCorr\\'+'leap0')
iau_1980_nut=np.loadtxt(os.getcwd()+'\\Box Sync\\BaryCorr\\'+'nut_iau1980.txt')


def EOPdata(mjd=mjdutc,download=0):
    '''
    Analog of IDL code eopdata.pro
    
    Input:
        If download==0 then will use local Astropy copy of IERS bulletin. If download==1 then will download 
    
    
    Output:
        mjdutc: MJD date being calculated for
        UT1_UTC
        pmx,pmy,dpsi,deps in arcseconds
    
    '''
    
    as2r = np.pi/(3600.*180.)
    
    #iers.IERS_A_README
    if download==0:
        dat_a = iers.IERS_Auto.open()  
    else:    
        dat_a = iers.IERS_A.open(iers.IERS_A_URL)  
    #dat_b = iers.IERS_B.open(iers.IERS_B_URL)  
    
    colnames = dat_a.colnames
    
    MJD=dat_a['MJD'].value
    pmx0=np.deg2rad(dat_a['PM_x_A'].value/3600.)
    pmy0=np.deg2rad(dat_a['PM_y_A'].value/3600.)     
    dpsi0=np.deg2rad(dat_a['dX_2000A_A'].value/3600.)
    deps0=np.deg2rad(dat_a['dY_2000A_A'].value/3600.)
    UT1=dat_a['UT1_UTC_A'].value  
       
    wh=np.append(np.where(np.array([UT1[x]-UT1[x-1] for x in range(1,len(UT1))])>0.8)[0],len(UT1)-1) # Check for discontinuities 
    leap0=np.zeros(len(UT1))
    
    for i in range(0,len(wh)-1):
        leap0[wh[i]+1:wh[i+1]+1] = (leap0[wh[i]])-np.round(UT1[wh[i]+1]-UT1[wh[i]])
    
    
    roi=[np.where((MJD<=(mjdutc+5)) & (MJD>=(mjdutc-5)) & (UT1 != 0 ))][0][0] #Region of interest
    
    MJD1 = MJD[roi]
    
    ii = np.where(MJD1==np.round(mjdutc))[0]
    jj = roi[ii]
       
    
    UT11=(UT1+leap0)  
    
    #INTERPOLATION
    dt = (mjdutc-MJD1[ii])/(MJD1[ii+1]-MJD1[ii])

    UT1_UTC = UT1[jj]+dt*(UT11[jj+1]-UT11[jj])    
        
    dpsi=np.rad2deg(dpsi0[jj]+dt*(dpsi0[jj+1]-dpsi0[jj]))*3600.
    deps=np.rad2deg(deps0[jj]+dt*(deps0[jj+1]-deps0[jj]))*3600.   
    
    pmx=np.rad2deg(pmx0[jj]+dt*(pmx0[jj+1]-pmx0[jj]))*3600.
    pmy=np.rad2deg(pmy0[jj]+dt*(pmy0[jj+1]-pmy0[jj]))*3600.
    
    return mjdutc,UT1_UTC,pmx,pmy,dpsi,deps



'''
argfacts=iau_1980_nut[:,0:4]
psiamps=iau_1980_nut[:,6:8]
epsamps=iau_1980_nut[:,8:10]


# IAU1980
#c1 = mean anomaly of Moon
#c2 = mean anomaly of Sun
#c3 = mean longitude of the Moon minus the mean longitude of Moon's node
#c4 = mean elongation of Moon from Sun
#c5 = mean longitude of ascending node of the Moon
#c6 = mean anomaly of Mercury
#c7 = mean anomaly of Venus
#c8 = mean anomaly of Earth
#c9 = mean anomaly of Mars
#ca = mean anomaly of Jupiter
#cb = mean anomaly of Saturn

# All values in arcseconds
c1 = [134.96298139*3600, 1717915922.6330,  31.310,   0.064]
c2 = [357.52772333*3600,  129596581.2240,  -0.577,  -0.012]
c3 = [ 93.27191028*3600, 1739527263.1370, -13.257,   0.011]
c4 = [297.85036306*3600, 1602961601.3280,  -6.891,   0.019]
c5 = [125.04452222*3600,   -6962890.5390,   7.455,   0.008]
c6 = [252.3       *3600,     149472.7,      0,       0    ]
c7 = [179.9       *3600,      58517.8,      0,       0    ]
c8 = [ 98.4       *3600,      35999.4,      0,       0    ]
c9 = [353.3       *3600,      19140.3,      0,       0    ]
ca = [ 32.3    *3600,       3034.9,      0,       0    ]
cb = [ 48.0    *3600,       1222.1,      0,       0    ]

args=[c1,c2,c3,c4,c5,c6,c7,c8,c9,ca,cb]
#Conver to radians.
arg80=(np.deg2rad(args))/3600.


fepoch=2451545.0

# From epoch of date in centuries from J2000.0

t=(jdtt-2451545.0)/36525.0


# NOT INTEPROLATING JPL EPHEMERIDES OF NUTATION. NOT USING GET_XTECAL


'''



