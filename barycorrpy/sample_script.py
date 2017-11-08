from .barycorrpy import call_BCPy 
from astropy.time import Time
import datetime
import numpy as np

JDUTC = Time(datetime.datetime.utcnow(),format='datetime',scale='utc')

# Sample case for Tau Ceti taken from CTIO.
results  = call_BCPy(JDUTC=JDUTC,hip_id=8102,lat=-30.169283,longi=-70.806789,alt=2241.9,ephemeris='de430')
print results