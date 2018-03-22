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

    # Check for masked values in the most convoluted way possible    
    if all([not x for x in [obj.mask[0][i] for i in obj.colnames]])==False:
        warning = 'Masked values present in queried dataset'
        print(warning)

    
    obj = obj.filled(0)
    
    
    pos = SkyCoord(ra=obj['RA_2_A_ICRS_J2000'],dec=obj['DEC_2_D_ICRS_J2000'],unit=(u.hourangle, u.deg))
    ra = pos.ra.value
    dec = pos.dec.value
    pmra = obj['PMRA'][0]
    pmdec = obj['PMDEC'][0]
    plx = obj['PLX_VALUE'][0]
    rv = obj['RV_VALUE'][0]
    epoch = 2451545.0
    
    star = {}
    star['ra'] = ra
    star['dec'] = dec
    star['pmra'] = pmra
    star['pmdec'] = pmdec
    star['plx'] = plx
    star['rv'] = rv
    star['epoch'] = epoch
    
    warning = 'Querying object {} from Simbad. Values used are RA:{} Dec:{} PMRA:{} PMDEC:{} Plx:{} RV:{} Epoch:{} '.format(name,ra,dec,pmra,pmdec,plx,rv,epoch)
    print(warning)
    
    
    return ra,dec,pmra,pmdec,plx,rv,epoch
  