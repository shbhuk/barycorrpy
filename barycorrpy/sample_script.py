from __future__ import print_function
from __future__ import division
from .barycorrpy import get_BC_vel , exposure_meter_BC_vel
from . import utc_tdb
from astropy.time import Time


def run_sample():
            JDUTC = 2458000 # Also accepts float input for JDUTC. Verify scale and format
            
            # Observation of Tau Ceti taken from CTIO on JD 2458000. 
            # Observatory location manually entered. Stellar positional parameters taken from Hipparcos Catalogue
            result = get_BC_vel(JDUTC=JDUTC,hip_id=8102,lat=-30.169283,longi=-70.806789,alt=2241.9,ephemeris='de430',zmeas=0.0)
            
            
            # Observation of Tau Ceti taken from CTIO on JD 2458000. Observatory location taken from Astropy list. 
            # Stellar positional parameters taken from Hipparcos Catalogue
            JDUTC = Time(2458000,format='jd',scale='utc')
            result2  = get_BC_vel(JDUTC=JDUTC,hip_id=8102,obsname='CTIO',ephemeris='de430')
            
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
            
            result3=get_BC_vel(JDUTC=JDUTC,ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,
                pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_update=True)
                                
            if result3[2]:
                print('Check ')   
                        
            # Exposure meter calculation
            flux = [4.5,1.5,2] # Same number of elements as JDUTC
                        
            result4,JDUTCMID,warning4,status4=exposure_meter_BC_vel(JDUTC=JDUTC,expmeterflux = flux,
                ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,
                pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_update=True)

            # JDUTC to BJDTDB time converter
            corr_time = utc_tdb.JDUTC_to_BJDTDB(JDUTC,hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9)

            return result,result2,result3,result4,JDUTCMID,warning4,status4,corr_time
