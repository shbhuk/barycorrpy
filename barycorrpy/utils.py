import numpy as np
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
import os


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