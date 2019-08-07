from __future__ import print_function
from __future__ import division
from astropy.time import Time
import numpy as np
from .barycorrpy import get_BC_vel , exposure_meter_BC_vel
from . import utc_tdb

def run_sample():
            a=[]
            b=0

            JDUTC = 2458000 # Also accepts float input for JDUTC. Verify scale and format

            # Observation of Tau Ceti taken from CTIO on JD 2458000.
            # Observatory location manually entered. Stellar positional parameters taken from Hipparcos Catalogue
            result = get_BC_vel(JDUTC=JDUTC, hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9, ephemeris='de430', zmeas=0.0)

            if np.isclose(a = result[0], b = 15403.9508, atol = 1e-2, rtol = 0):
                a.append('result')
                b+=1


            # Observation of Tau Ceti taken from CTIO on JD 2458000. Observatory location taken from Astropy list.
            # Stellar positional parameters taken from Hipparcos Catalogue
            JDUTC = Time(2458000, format='jd', scale='utc')
            result2  = get_BC_vel(JDUTC=JDUTC, hip_id=8102, obsname='CTIO', ephemeris='de430')

            if np.isclose(a = result2[0], b = 15403.9608, atol = 1e-2, rtol = 0):
                a.append('result2')
                b+=1

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
            JDUTC=[2458000,2458000.00001,2458000.00002] # Can also enter JDUTC as float instead of Astropy Time Object
            ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'

            result3=get_BC_vel(JDUTC=JDUTC, ra=ra, dec=dec, obsname=obsname, lat=lat, longi=longi, alt=alt, pmra=pmra,
                pmdec=pmdec, px=px, rv=rv, zmeas=zmeas,epoch=epoch, ephemeris=ephemeris, leap_update=True)

            if np.allclose([result3[0][0],result3[0][1],result3[0][2]],[15407.4860,15407.4723,15407.4586],atol = 1e-2, rtol = 0):
                a.append('result3')
                b+=1



            # Exposure meter calculation
            flux = [4.5,1.5,2] # Same number of elements as JDUTC

            result4,JDUTCMID,warning4,status4=exposure_meter_BC_vel(JDUTC=JDUTC, expmeterflux = flux,
                ra=ra, dec=dec, obsname=obsname, lat=lat, longi=longi, alt=alt, pmra=pmra,
                pmdec=pmdec, px=px, rv=rv, zmeas=zmeas, epoch=epoch, ephemeris=ephemeris, leap_update=True)

            if np.isclose(a = result4, b = 15407.4765, atol = 1e-2, rtol = 0):
                a.append('result4')
                b+=1

            # JDUTC to BJDTDB time converter
            corr_time = utc_tdb.JDUTC_to_BJDTDB(JDUTC, hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9)

            if np.allclose([corr_time[0]], [2458000.00505211,  2458000.00506211,  2458000.00507211], atol = 1e-7, rtol = 0):
                a.append('corr_time')
                b+=1


            if b==5:
                print('***********SUCCESS**************\nAll barycentric correction velocities match expected values to 1 cm/s\n')
            else:
                print('{} out of 5 results match. Compare outputs vs those on the github wiki. Check others - \n'.format(b,a))

            return result, result2, result3, result4, JDUTCMID, warning4, status4, corr_time
