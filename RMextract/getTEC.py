import RMextract.PosTools
import RMextract.getIONEX as ionex
import os
import numpy as np
from math import *

ION_HEIGHT=450.e3

def getTEC(MS=None,
           server="ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/",
           prefix='codg',
           ionexPath="IONEXdata/",
           earth_rot=0,ha_limit=-1000,
           **kwargs):
    '''optional arguments are  radec or pointing : [ra,dec] in radians,
    azel : [az,el] in radians, If azel is specified, radec will be ignored
    timestep in s, timerange = [start, end]in MJD seconds, 
    otherwise use start_time/end_time (either MJD or casa timestring (NEEDS PYRAP), eg. 2012/11/21/12:00:00  or 56252.5d  ),    
    stat_names = list of strings per station, 
    stat_positions = list of length 3 numpy arrays, containing station_position in ITRF meters.
    Returns the (timegrid,timestep,TEC) where TEC is a dictionary containing 1 enumpyarray per station in stat_names. 
    If stat_names is not given, the station names will either be extracted from the MS or st1...stN '''

    stat_names=[]
    stat_pos=[PosTools.posCS002]
    #stat_names=['LOFAR_CS002']
    # the following parameter is used to extend the observation range by 120 sec
    # before and after the actual specified time. If we want to correct an
    # actual data set, this is required for the scipy 1d interpolation routine to
    # ensure that the calculated range of data exceeds the time range actually
    # observed - I have used the same value in ALBUS - Tony
    TIME_OFFSET = 0
    if not (MS is None):

        (timerange,timestep,pointing,stat_names,stat_pos)=PosTools.getMSinfo(MS);
    useAzel=False
    for key in kwargs.keys():
        if not useAzel and (key=='radec' or key=='pointing'):
            pointing = kwargs[key]
        if key=='start_time':
            start_time=kwargs[key]
            time_in_sec = False
        if key=='end_time':
            end_time=kwargs[key]
            time_in_sec = False
        if key=='azel':
            pointing = kwargs[key]
            useAzel=True
        if key=='timestep':
            timestep=kwargs[key]
        if key=='timerange':
            timerange=kwargs[key]
        if key=='stat_names':
            stat_names=kwargs[key]
        if key=='TIME_OFFSET':
           TIME_OFFSET=kwargs[key]
        if key=='stat_positions':
            stat_pos=kwargs[key]
            if not stat_names:
                stat_names =['st%d'%(i+1) for i in range(len(stat_pos))]

    if timerange != 0:
      start_time = timerange[0]-TIME_OFFSET
      end_time = timerange[1]+TIME_OFFSET
      time_in_sec = True
      reference_time = start_time
      str_start_time=PosTools.obtain_observation_year_month_day_hms(reference_time)
    else:
        timerange,str_start_time,reference_time=PosTools.get_time_range(start_time,end_time,timestep,time_in_sec,TIME_OFFSET)
    if str_start_time==-1:
       return

    #IONEX files go per day, check if more than one file is  needed.
    times,timerange=PosTools.getIONEXtimerange(timerange,timestep)
    #times.append(np.arange(timerange[0],timerange[1]+timestep,timestep)) #add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
    timegrid=np.array([])
    TECs={};
    ams={};
    flags={}
    for st in stat_names:
        TECs[st]=[]
        ams[st]=[]
        flags[st]=[]
    for time_array in times:
        #get RM per timeslot
        starttime=time_array[0];
        date_parms =  PosTools.obtain_observation_year_month_day_fraction(starttime)
        #get relevant ionex file
        ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath)
        if ionexf==-1:
            print ("error opening ionex data")
            return
       
        tecinfo=ionex.readTEC(ionexf)
        for station,position in  zip(stat_names,stat_pos):
          #print ("getting TEC for",station,"at",position,"time",time_array )
          for time in time_array:
            result =  PosTools.obtain_observation_year_month_day_fraction(time)
            part_of_day= result[3] * 24
            if useAzel:
                  az=pointing[0]
                  #el=0.5*np.pi-pointing[1]
                  el=pointing[1]
                  flags[station].append(1)
            else:

                az,el=PosTools.getAzEl(pointing,time,position,ha_limit)
                if az==-1 and el==-1:
                    return
                if az==999 and el==999:
                    #outsite hadec range
                    TECs[station].append(-999)  #sTEC value
                    ams[station].append(-999)
                    flags[station].append(0)
                    continue
                flags[station].append(1) 
            latpp,lonpp,height,lon,lat,am1=PosTools.getlonlatheight(az,el,position)
            if latpp==-1 and lonpp==-1 and height==-1:
               return
            vTEC=ionex.getTECinterpol(time=part_of_day,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)

            TECs[station].append(vTEC)  #sTEC value
            ams[station].append(am1)
        timegrid=np.concatenate((timegrid,time_array));
    
    for st in stat_names:
            TECs[st]=np.array(TECs[st])
            ams[st]=np.array(ams[st])
            flags[st]=np.array(flags[st])
    return (timegrid,timestep,TECs,ams,flags)
        


