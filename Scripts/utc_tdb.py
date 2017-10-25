from __future__ import division
import urllib
import astropy
import os
import datetime
from astropy.time import Time
import numpy as np



def JDUTC_to_JDTDB(utctime,fpath):
    '''
    Enter UTC time as Astropy Time Object
    
    fpath - Path to where the file would be saved.
    
    
    Output:
        JDTDB : Julian Date Barycentric Dynamic Time (Astropy Time object)
        JDTT: Julian Date Terrestrial Dynamic time (Astropy Time object)
    '''
    

    url='http://maia.usno.navy.mil/ser7/tai-utc.dat'

    # Location of saved file
    fpath+='utc_tai.txt'

    

    
    if os.path.isfile(fpath)==False:
        openURL = urllib.URLopener()
        openURL.retrieve(url, fpath)
    else:
        filetime=datetime.datetime.fromtimestamp(os.path.getmtime(fpath))  
        now=datetime.datetime.now()  # Current time 
    
        file_history=now-filetime  # How old is the file?
        if file_history>datetime.timedelta(days=180):
            print 'File more than 180 days old'
            openURL = urllib.URLopener()
            openURL.retrieve(url, fpath)
            print "Downloaded Leap second file"

    f=open(fpath,'r')
    jd=[]
    offset=[]
    
    for line in f:
        a=(line.split())
        jd.append(float(a[4]))
        offset.append(float(a[6]))
    f.close()
    
    if type(utctime)!=astropy.time.core.Time: # Check if it is in Astropy Time Object
        print "Input Time should be as Astropy Time Object"
        raise TypeError("Input Time should be as Astropy Time Object")
    
    else:
        JDUTC=utctime.jd
        check_time=utctime.datetime

    tai_utc=offset[np.max(np.where(JDUTC>=np.array(jd)))] # Leap second offset value to convert UTC to TAI
    new_tai=check_time+datetime.timedelta(seconds=tai_utc)  # Add offset and convert to TAI
    new_tt=new_tai+datetime.timedelta(seconds=32.184)  # Add 32.184 to convert TAI to TT
    g=(357.53+0.9856003*( JDUTC - 2451545.0 ))*np.pi/180.  # Earth's mean anomaly
    TDB = new_tt+datetime.timedelta(seconds=0.001658*np.sin(g)+0.000014*np.sin(2*g)) # TT to TDB
    JDTT=Time(new_tt,scale='tt',format='datetime')
    JDTT.format='jd'
    JDTDB = Time(TDB,scale='tdb',format='datetime')
    JDTDB.format='jd'

   
    return JDTDB,JDTT




