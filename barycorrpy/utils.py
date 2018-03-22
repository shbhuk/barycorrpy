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
    
    obj = Simbad.query_object(name)
    
    pos = SkyCoord(ra=obj['RA'],dec=obj['DEC'],unit=(u.hourangle, u.deg))
    
    