import numpy as np
from astroquery.simbad import Simbad
import astropy.constants as u
from astropy.coordinates import SkyCoord


def flux_weighting(flux,qty):
    '''
    INPUT:
        flux - Numpy Array. Will weight qty array by flux
        qty - Quantity to be normalized / weighted.
    
    OUTPUT:
        Normalized value
    '''
    
    return sum(qty * flux)/sum(flux)
    
def get_stellar_data(name='',use_hip_if_exists = True):
    '''
    Function to query Simbad for following stellar information RA, Dec, PMRA, PMDec, Parallax Epoch
    INPUTS:
        name = Name of source. Example 
        use_hip_if_exists = If source exists in Hipparcos catalog and has a HIP ID, use Hipparcos
    
    
    '''
    customSimbad = Simbad()
    customSimbad.add_votable_fields('pm', 'plx','parallax','rv_value')
    #Simbad.list_votable_fields()
    #customSimbad.remove_votable_fields( 'coordinates')
    #Simbad.get_field_description('orv')
    obj = customSimbad.query_object(name)
    
    pos = SkyCoord(ra=obj['RA'],dec=obj['DEC'],unit=(u.hourangle, u.deg))
    ra = pos.ra.value
    dec = pos.dec.value
    pmra = obj['PMRA'][0]
    pmdec = obj['PMDEC'][0]
    plx = obj['PLX_VALUE'][0]
    rv = obj['RV_VALUE'][0]
  