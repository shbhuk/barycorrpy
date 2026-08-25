from __future__ import print_function
import numpy as np
try:
    from astroquery.simbad import Simbad
except:
    print('Cannot import astroquery.simbad')
import astropy.units as u
from astropy.coordinates import SkyCoord
import os
# import winsound
# frequency = 2000  # Set Frequency To 2500 Hertz
# duration = 500 

def flux_weighting(flux,qty):
    '''
    INPUT:
        flux - Numpy Array. Will weight qty array by flux
        qty - Quantity to be normalized / weighted.
    
    OUTPUT:
        Normalized value
    '''
    
    return sum(qty * flux)/sum(flux)
    
def _astroquery_simbad_is_new():
    '''
    True if the installed astroquery uses the post-0.4.8 SIMBAD field schema.
    '''
    from astroquery import __version__ as aq_version
    from packaging.version import parse
    return parse(aq_version) >= parse('0.4.8')


def _query_simbad_object(name):
    '''
    Query SIMBAD for `name` across astroquery versions.

    Many pipelines still run on astroquery < 0.4.8, which used a different set
    of votable field names and coordinate units than current releases. All of
    that version-specific knowledge is isolated here so that the caller works
    with a single, uniform set of keys regardless of the installed version.

    Note: the `keys` map the old names to the new SIMBAD names, so this can
    be deprecated when older versions of astroquery are no longer supported.

    OUTPUT:
        obj : the astropy Table returned by SIMBAD (or None if not found)
        keys : dict mapping logical name -> column name in `obj`
        ra_dec_unit : units to feed SkyCoord for the ra/dec columns
    '''
    customSimbad = Simbad()

    if _astroquery_simbad_is_new():
        # astroquery >= 0.4.8 (new SIMBAD field names)
        customSimbad.add_votable_fields('ra', 'dec', 'pmra', 'pmdec', 'plx_value', 'rvz_radvel')
        keys = {'ra': 'ra', 'dec': 'dec', 'pmra': 'pmra', 'pmdec': 'pmdec',
                'plx': 'plx_value', 'rv': 'rvz_radvel'}
        ra_dec_unit = (u.deg, u.deg)
    else:
        # astroquery < 0.4.8 (legacy SIMBAD field names)
        customSimbad.add_votable_fields('ra(2;A;ICRS;J2000)', 'dec(2;D;ICRS;J2000)',
                                        'pm', 'plx', 'parallax', 'rv_value')
        customSimbad.remove_votable_fields('coordinates')
        keys = {'ra': 'RA_2_A_ICRS_J2000', 'dec': 'DEC_2_D_ICRS_J2000',
                'pmra': 'PMRA', 'pmdec': 'PMDEC',
                'plx': 'PLX_VALUE', 'rv': 'RV_VALUE'}
        ra_dec_unit = (u.hourangle, u.deg)

    obj = customSimbad.query_object(name)
    return obj, keys, ra_dec_unit


def get_stellar_data(name=''):
    '''
    Function to query Simbad for following stellar information RA, Dec, PMRA, PMDec, Parallax Epoch
    INPUTS:
        name = Name of source. Example


    '''
    warning = []

    obj, keys, ra_dec_unit = _query_simbad_object(name)
    if obj is None:
        raise ValueError('ERROR: {} target not found. Check target name or enter RA,Dec,PMRA,PMDec,Plx,RV,Epoch manually\n\n'.format(name))
    else:
        warning += ['{} queried from SIMBAD.'.format(name)]

    # Check for masked values
    if all([not x for x in [obj.mask[0][i] for i in obj.colnames]])==False:
        warning += ['Masked values present in queried dataset']


    obj = obj.filled(None)

    pos = SkyCoord(ra=obj[keys['ra']], dec=obj[keys['dec']], unit=ra_dec_unit)
    ra = pos.ra.value[0]
    dec = pos.dec.value[0]
    pmra = obj[keys['pmra']][0]
    pmdec = obj[keys['pmdec']][0]
    plx = obj[keys['plx']][0]
    rv = obj[keys['rv']][0] * 1000 #SIMBAD output is in km/s. Converting to m/s
    epoch = 2451545.0

    star = {'ra':ra,'dec':dec,'pmra':pmra,'pmdec':pmdec,'px':plx,'rv':rv,'epoch':epoch}
    
    # Fill Masked values with None. Again. 
    for i in star:
        if star[i] > 1e10:
            star[i] = None           
       
    warning += ['Values queried from SIMBAD are {}'.format(star)]
  
    
    return star,warning
    
def find_hip(hip_index,cat_dir=os.path.join(*[os.path.dirname(__file__),'data','hip2.dat'])):
    '''
    NOTE: Hipparcos Catalogue Epoch is J1991.25 or JD 2448349.0625
    
    INPUT:
        hip_index : The index of the star that needs to be searched.
        cat_dir : Directory where catalogue is saved. 
    OUTPUT:
        hip_id : Index of Star
        ra : RA of star in degrees
        dec : Declination of star in degrees
        px_mas : Parallax angle in milli-arcseconds
        pmra : Proper motion in RA in milli-arcseconds/year
        pmdec : Proper motion in Declination in milli-arcseconds/year 
        epoch : Epoch of Catalogue - J1991.25 , JD 2448349.0625       
    
    '''
    epoch = 2448349.0625
    
    hip_id=[]
    ra=[]
    dec=[]
    px_mas=[]
    pmra=[]
    pmdec=[]

    
    with open(cat_dir) as f:
        for line in f:
            a = line.split()
            hip_id.append(float(a[0]))
            ra.append((float(a[4])*180.)/np.pi) # Convert to degrees
            dec.append((float(a[5])*180.)/np.pi) # Convert from radians to degrees
            px_mas.append(float(a[6]))  # in mas
            pmra.append(float(a[7])) # in mas/year
            pmdec.append(float(a[8]))  # in mas/year
        
    index=np.where(np.array(hip_id)==hip_index)[0][0]
    
    star = {'ra':ra[index],'dec':dec[index],'pmra':pmra[index],'pmdec':pmdec[index],'px':px_mas[index],'epoch':epoch}
    
    return star

      
def CalculatePositionVector(r1, r2):
    '''
    INPUTS:
        r1, r2 = (x,y,z) positions for Obj 1 and Obj 2.
    OUTPUTS:    
        vec(X) = vec(r1) - vec(r2)
        Xmag, Xhat = Magnitude of the position vector, and its unit vector resp. are also returned.
    
    '''
    # Vector from object barycenter to Observatory
    X = np.array(r1-r2)
    Xmag = np.sqrt(np.sum(X*X)) # [m]
    Xhat = X / Xmag # Unitless
    
    return X, Xmag, Xhat