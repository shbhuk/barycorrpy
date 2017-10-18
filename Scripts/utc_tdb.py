import urllib
import numpy
import astropy
import os
import datetime
from astropy.time import Time
import numpy as np
import math


def JDUTC_to_JDTDB(utctime):
    '''
    Enter UTC time as Astropy Time Object
    
    
    Output:
        JDTDB : Julian Date Barycentric Dynamic Time (Astropy Time object)
        TDB: Barycentric Dynamic Time (datetime.datetime object)
        JDTT: Julian Date Terrestrial Dynamic time (Astropy Time object)
    '''
    
    #url='https://hpiers.obspm.fr/eoppc/bul/bulc/UTC-TAI.history'
    url='http://maia.usno.navy.mil/ser7/tai-utc.dat'
    location=os.getcwd()+'/Box Sync/BaryCorr/utc_tai.txt'
    
    utc_tai=datetime.datetime.fromtimestamp(os.path.getmtime(location))  # Location of saved file
    now=datetime.datetime.now()  # Current time 
    
    file_history=now-utc_tai  # How old is the file?
    
    if file_history>datetime.timedelta(days=20):
        print 'File more than 10 day old'
        openURL = urllib.URLopener()
        openURL.retrieve(url, location)
        print "Downloaded Leap second file"


    f=open(location,'r')
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
      
    #print '\nConverting UTC time = '+str(check_time) +' \nJDUTC = '+ (str(JDUTC)) 
    

    tai_utc=offset[np.max(np.where(JDUTC>=np.array(jd)))]  # Leap second offset value to convert UTC to TAI
    
    new_tai=check_time+datetime.timedelta(seconds=tai_utc)  # Add offset and convert to TAI
    new_tt=new_tai+datetime.timedelta(seconds=32.184)  # Add 32.184 to convert TAI to TT

    g=(357.53+0.9856003*( JDUTC - 2451545.0 ))*np.pi/180.  # Earth's mean anomaly
    TDB = new_tt+datetime.timedelta(seconds=0.001658*np.sin(g)+0.000014*np.sin(2*g)) # TT to TDB
    JDTT=Time(new_tt,scale='tt',format='datetime')
    JDTT.format='jd'
    JDTDB = Time(TDB,scale='tdb',format='datetime')
    JDTDB.format='jd'

    #print '\nJDTDB = '+str(JDTDB)
    #print 'TDB = '+str(TDB)    
    #print 'JDTT = '+str(JDTT)
   
    return JDTDB,JDTT




