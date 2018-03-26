import numpy as np
from astroquery.simbad import Simbad
import astropy.units as u
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
    
def get_stellar_data(name=''):
    '''
    Function to query Simbad for following stellar information RA, Dec, PMRA, PMDec, Parallax Epoch
    INPUTS:
        name = Name of source. Example 
    
    
    '''
    warning = []
    
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(2;A;ICRS;J2000)', 'dec(2;D;ICRS;J2000)','pm', 'plx','parallax','rv_value')
    #Simbad.list_votable_fields()
    customSimbad.remove_votable_fields( 'coordinates')
    #Simbad.get_field_description('orv')
    obj = customSimbad.query_object(name)
    if obj is None:
        raise ValueError('ERROR: {} target not found. Check target name or enter RA,Dec,PMRA,PMDec,Plx,RV,Epoch manually\n\n'.format(name))
    else:        
        warning += ['{} queried from SIMBAD.'.format(name)]

    # Check for masked values
    if all([not x for x in [obj.mask[0][i] for i in obj.colnames]])==False:
        warning += ['Masked values present in queried dataset']


    obj = obj.filled(None)
    
    pos = SkyCoord(ra=obj['RA_2_A_ICRS_J2000'],dec=obj['DEC_2_D_ICRS_J2000'],unit=(u.hourangle, u.deg))
    ra = pos.ra.value[0]
    dec = pos.dec.value[0]
    pmra = obj['PMRA'][0]
    pmdec = obj['PMDEC'][0]
    plx = obj['PLX_VALUE'][0]
    rv = obj['RV_VALUE'][0]
    epoch = 2451545.0
    
    star = {'ra':ra,'dec':dec,'pmra':pmra,'pmdec':pmdec,'px':plx,'rv':rv,'epoch':epoch}
    
    # Fill Masked values with None. Again. 
    for i in [k for k in star if star[k]==1e20]:
        star[i]=None
        
    warning += ['Values queried from SIMBAD are {}'.format(star)]
  
    
    return star,warning