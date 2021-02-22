import unittest
from barycorrpy import get_BC_vel , exposure_meter_BC_vel
import numpy as np
from astropy.time import Time
from barycorrpy import utc_tdb



class Barycorrpy_tests(unittest.TestCase):

    def test_hip_id(self):
        JDUTC = 2458000 # Also accepts float input for JDUTC. Verify scale and format
        result = get_BC_vel(JDUTC=JDUTC,hip_id=8102,lat=-30.169283,longi=-70.806789,alt=2241.9,ephemeris='de430',zmeas=0.0)
        self.assertTrue(np.isclose(a = result[0], b = 15403.9508,atol = 1e-2, rtol = 0), 4)

    def test_hip_id_astropy_obs(self):
        JDUTC = Time(2458000,format='jd',scale='utc')
        result  = get_BC_vel(JDUTC=JDUTC,hip_id=8102,obsname='CTIO',ephemeris='de430')
        self.assertTrue(np.isclose(a = result[0],b = 15403.9608,atol = 1e-2, rtol = 0))

    def test_inputs(self):
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

        result=get_BC_vel(JDUTC=JDUTC,ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,
                pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,ephemeris=ephemeris,leap_update=True)

        self.assertTrue(np.allclose([result[0][0],result[0][1],result[0][2]],[15407.4860,15407.4723,15407.4586],atol = 1e-2, rtol = 0))

    def test_flux_weighting(self):
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
        JDUTC=[2458000,2458000.00001,2458000.00002]

        flux = [4.5,1.5,2] # Same number of elements as JDUTC

        result,JDUTCMID,warning4,status4=exposure_meter_BC_vel(JDUTC=JDUTC,expmeterflux = flux,
                ra=ra,dec=dec,obsname=obsname,lat=lat,longi=longi,alt=alt,pmra=pmra,
                pmdec=pmdec,px=px,rv=rv,zmeas=zmeas,epoch=epoch,leap_update=True)

        self.assertTrue(np.isclose(a = result,b = 15407.4765,atol = 1e-2, rtol = 0))


    def test_JDUTC_BJDTDB(self):
        JDUTC=[2458000,2458000.00001,2458000.00002]
        corr_time = utc_tdb.JDUTC_to_BJDTDB(JDUTC,hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9)

        self.assertTrue(np.allclose([corr_time[0]],[2458000.00505211,  2458000.00506211,  2458000.00507211],atol = 1e-7, rtol = 0))

    def test_stellar_predictive(self):
        JDUTC = 2458000
        result5 = get_BC_vel(JDUTC=JDUTC, hip_id=8102, lat=-30.169283, longi=-70.806789, alt=2241.9, ephemeris='de430', zmeas=0.0, predictive=True)

        self.assertTrue(np.isclose(a = result5[0], b = -15403.15938, atol = 1e-2, rtol = 0))


    def test_SolarBC(self):
        JDUTC = 2458000

        result6 = get_BC_vel(JDUTC=2458000, lat=-30.169138888, longi=-70.805888, alt=2379.5, zmeas=0.0, SolSystemTarget='Sun')
        self.assertTrue(np.isclose(a = result6[0], b = 819.4474, atol = 1e-2, rtol = 0))

    def test_SolarEmissionTDB(self):

        result7 = utc_tdb.JDUTC_to_SolarEmissionTDB(JDUTC=2458000, obsname='KPNO')
        self.assertTrue(np.isclose(a=result7[0], b=2457999.99497543, atol=1e-7, rtol=0))



def main():
    unittest.main()

if __name__ == "__main__":
    main()
