import numpy as np
import inspect
import os

def find_hip(hip_index,cat_dir=os.path.join(os.path.dirname(__file__),'hip2.dat')):
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

    
    #HIP=np.loadtxt(f)
    
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
    
    return hip_id[index],ra[index],dec[index],px_mas[index],pmra[index],pmdec[index],epoch
