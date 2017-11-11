from __future__ import print_function
from .barycorrpy import call_BCPy 
from astropy.time import Time
import datetime


JDUTC = Time(2458000,format='jd',scale='utc')

# Observation of Tau Ceti taken from CTIO on JD 2458000. 
# Observatory location manually entered. Stellar positional parameters taken from Hipparcos Catalogue
result  = call_BCPy(JDUTC=JDUTC,hip_id=8102,lat=-30.169283,longi=-70.806789,alt=2241.9,ephemeris='de430')


# Observation of Tau Ceti taken from CTIO on JD 2458000. Observatory location taken from Astropy list. 
# Stellar positional parameters taken from Hipparcos Catalogue

result2  = call_BCPy(JDUTC=JDUTC,hip_id=8102,obsname='CTIO',ephemeris='de430')

# Observation of Tau Ceti taken from CTIO on JD 2458000,2458010,2458020. 
# Observatory and stellar parameters entered by user.
# Use DE405 ephemeris

ra=26.0213645867
dec=-15.9395557246  
obsname=''
lat=-30.169283
longi=-70.806789
alt=2241.9
epoch = 2451545.0  
pmra = -1721.05
pmdec = 854.16 
px = 273.96 
rv = 0.0 
zmeas=0.0
JDUTC=Time([2458000,2458010,2458020],format='jd',scale='utc')
ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'

result3=call_BCPy(JDUTC=JDUTC,ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,
    pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_update=True)
