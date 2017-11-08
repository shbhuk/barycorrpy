from __future__ import division
from __future__ import print_function
import urllib
import astropy
import os
import datetime
from astropy.time import Time
import numpy as np


def staleness_check(file_time,now):
    '''
    Check whether the leap second file needs to be updated.
    
    Leap second updates generally announced on December 31st or June 30th. Thus assuming a month lag for the list to be updated, 
    check if if a February 1st or August 1st has passed between when the file was downloaded and when the check is being run. 
    
    INPUT:
        Both inputs are in datetime format
        file_time : Time the leap second file was last downloaded. Stored in leap_log.txt. 
        now : UTC time right now. 
    
    
    Return 1 if need to update, else return 0. 
       
    '''
    
    if file_time.month > 8 :
        year=file_time.year+1
    else: year=file_time.year
    
    if (file_time < datetime.datetime(year,2,1) <= now) or (file_time < datetime.datetime(year,9,1) <= now):
        return 1
    else: return 0
    
    
def leap_download(ls_fpath,log_fpath):
    '''
    Download the leap second file and update the log file. 
    INPUT: 
        ls_fpath : Path to where the file will be saved. 
        log_fpath : Path to where the log file is created. 
    
    '''
    
    url='http://maia.usno.navy.mil/ser7/tai-utc.dat'
    
    openURL = urllib.URLopener()
    openURL.retrieve(url, ls_fpath)
        
    with open(log_fpath,'w') as f:
        f.write(str(datetime.datetime.utcnow())) # Write date of download in log file
         
   
def leap_manage(utctime,fpath,leap_update):
    '''
    
    Calculates the offset between UTC and TAI from the leap second file.
    Also, checks for 'staleness' (how old is the file) of the leap second file. 
    When the leap second file gets about 6 months old, it will update the leap second file, to check if there has been 
    a new leap second announced. 
    
    
    INPUT: 
        Enter UTC time as Astropy Time Object. In UTC Scale.
        
        fpath : Path to where the file would be saved.
        leap_update : If True, when the leap second file is more than 6 months old it will attempt to download a new one.
                    If False, then will just give a warning message.  
    
    OUTPUT: 
        tai_utc : Leap second offset between TAI and UTC (Scalar)
        
        Offset in seconds, between TAI and UTC.   
    '''
    
    warning=[]    

    # Location of saved leap second file
    ls_fpath=os.path.join(fpath,'leap_sec.txt')
    
    # Location of log file to record how old is the leap second file
    log_fpath=os.path.join(fpath,'leapsec_log.txt')
       
    # If log file does not exist, then download the leap second file and create a log file
    if (os.path.isfile(log_fpath)==False or os.path.isfile(ls_fpath)==False):
        if leap_update==True:
            leap_download(ls_fpath=ls_fpath,log_fpath=log_fpath)
            warning+=['Downloaded leap second file from http://maia.usno.navy.mil/ser7/tai-utc.dat ']            
        else:
            warning+=['ERROR : LEAP SECOND FILE / LOG FILE DOES NOT EXIST. Please set leap_update = True to update file. Corrections may not be accurate ']
            return 0,warning


    # If file exists then check if need to update
    
    else:
        
        with open(log_fpath,'r') as f:
            file_time=datetime.datetime.strptime(f.readline() , '%Y-%m-%d %H:%M:%S.%f') # Read date of download in log file        
         
        now=datetime.datetime.utcnow()  # Current time  

        
        # Check if need to update file
        if staleness_check(file_time,now)==1:
            warning+=['WARNING : Need to update leap second file. Corrections may not be accurate to 1 cm/s with old file']
            if leap_update==True:
                leap_download(ls_fpath=ls_fpath,log_fpath=log_fpath)  
                warning+=['Downloaded leap second file from http://maia.usno.navy.mil/ser7/tai-utc.dat ']            
                    
            else: warning+=["ERROR : Leap second file should be updated. Set leap_update = True to download file"]

    f=open(ls_fpath,'r')
    jd=[]
    offset=[]
    
    for line in f:
        a=(line.split())
        jd.append(float(a[4]))
        offset.append(float(a[6]))
    f.close()
    
    if type(utctime)!=astropy.time.core.Time: # Check if it is in Astropy Time Object
        print("Input Time should be as Astropy Time Object")
        raise TypeError("Input Time should be as Astropy Time Object")
    else:
        JDUTC=utctime.jd

    tai_utc=offset[np.max(np.where(JDUTC>=np.array(jd)))] # Leap second offset value to convert UTC to TAI
        
    return tai_utc  , warning
 

def JDUTC_to_JDTDB(utctime,leap_update,fpath=os.path.dirname(__file__)):
    '''
    Convert JDUTC to JDTDB
    
    utctime : Enter UTC time as Astropy Time Object. In UTC Scale.
    
    fpath : Path to where the file would be saved. Default is script directory.
    leap_update : If True, when the leap second file is more than 6 months old it will attempt to download a new one.
                  If False, then will just give a warning message.     
    
    Output:
        JDTDB : Julian Date Barycentric Dynamic Time (Astropy Time object)
        JDTT: Julian Date Terrestrial Dynamic time (Astropy Time object)
        warning : Error message if any from 
        
        
        

    '''
    JDUTC=utctime.jd
    if JDUTC<2441317.5 : 
        return JDUTC.tdb,JDUTC.tt,['WARNING : Precise leap second history is not maintained here for before 1972. Defaulting to Astropy. Corrections maybe inaccurate']

    # Call function to check leap second file and find offset between UTC and TAI. 
    tai_utc,warning=leap_manage(utctime=utctime,fpath=fpath,leap_update=leap_update)
    
    if tai_utc==0:
        return JDUTC.tdb,JDUTC.tt,warning+['ERROR : Unable to maintain leap second file. Defaulting to AstroPy version']        

    check_time=utctime.datetime
    
    
    # Add leap seconds to convert UTC to TAI   
    new_tai=check_time+datetime.timedelta(seconds=tai_utc)  # Add offset and convert to TAI
    new_tt=new_tai+datetime.timedelta(seconds=32.184)  # Add 32.184 to convert TAI to TT

    g=(357.53+0.9856003*( JDUTC - 2451545.0 ))*np.pi/180.  # Earth's mean anomaly
    
    TDB = new_tt+datetime.timedelta(seconds=0.001658*np.sin(g)+0.000014*np.sin(2*g)) # TT to TDB
    
    JDTT=Time(new_tt,scale='tt',format='datetime')
    JDTT.format='jd'
    JDTDB = Time(TDB,scale='tdb',format='datetime')
    JDTDB.format='jd'

    return JDTDB,JDTT,warning
