#!/usr/bin/env python
from RMextract import PosTools
from RMextract import getIONEX as ionex;
import os
import numpy as np
from datetime import date

from RMextract.EMM import EMM as EMM

ION_HEIGHT=PosTools.ION_HEIGHT
#####################  main processing function #####################
def getRM(MS=None,
           server="ftp://ftp.aiub.unibe.ch/CODE/",
           prefix='codg',
           ionexPath="IONEXdata/",
           earth_rot=0,
           timerange=0,
           use_azel = False,ha_limit=-1000,use_filter=None,
           **kwargs):
    '''optional arguments are:
    radec or pointing : [ra,dec] in radians, or if use_azel =True, az + el in radians,
    timestep in s, timerange = [start, end]in MJD seconds, 
    otherwise use start_time/end_time (casa timestring,eg. 2012/11/21/12:00:00  or 56252.5d  NEEDS PYRAP) ,
    stat_names = list of strings per station, 
    stat_positions = list of length 3 numpy arrays, containing station_position in ITRF meters.
    useEMM = boolean, use EMM for Earth magnetic field, otherwise WMM cooefficients will be used.
    out_file = string, if given the data points will be written to a text file.
    use_mean = True if you only want report for mean of station positions
    use_filter =standard deviation,or list of standard deviations (time,long, lat) to gaussian filter TEC data 
    TIME_OFFSET = float, offset time at start and end to ensure all needed values are calculated,
    Returns the (timegrid,timestep,TEC) where TEC is a dictionary containing 1 enumpyarray per station in stat_names. 
    If stat_names is not given, the station names will either be extracted from the MS or st1...stN '''

    #print ('earth_rot', earth_rot)
    #print ('timerange ', timerange)
    #print ('use_azel ', use_azel)
    stat_names=[]
    useEMM=False
    use_mean = False
    myobject = ''
    timestep=60
    pointing=[0.,0.5*np.pi]
    out_file=''
    stat_pos=[PosTools.posCS002]
    use_proxy = False
    #stat_names=['LOFAR_CS002']
    if not (MS is None):

      (timerange,timestep,pointing,stat_names,stat_pos)=PosTools.getMSinfo(MS);
    
    for key in kwargs.keys():
        if key=='use_mean':
            use_mean = kwargs[key]
            stat_pos_mean = False
        if key=='radec' or key=='pointing':
            pointing = kwargs[key]
        if key=='object':
            myobject = kwargs[key]
            #out_file = 'RMextract_report_' + myobject
        if key=='start_time':
            start_time=kwargs[key]
            time_in_sec = False
        if key=='end_time':
            end_time=kwargs[key]
            time_in_sec = False
        if key=='timestep':
            timestep=kwargs[key]
        if key=='stat_names':
            stat_names=kwargs[key]
        if key=='out_file':
            out_file=kwargs[key]
        if key=='TIME_OFFSET':
           TIME_OFFSET=kwargs[key]
        if key=='stat_positions':
            stat_pos=kwargs[key]
            if not stat_names:
                stat_names =['st%d'%(i+1) for i in range(len(stat_pos))]
        if key=='useEMM':
            useEMM=kwargs[key]
        if key=="proxy_server":		#Check to see if user wants to use a proxy for downloading IONEX files.
            use_proxy = True
            proxy_server=kwargs["proxy_server"]
            proxy_type=kwargs["proxy_type"]
            proxy_port=kwargs["proxy_port"]
            proxy_user=kwargs["proxy_user"]
            proxy_pass=kwargs["proxy_pass"]
      
    if timerange != 0:
      start_time = timerange[0]
      end_time = timerange[1]
      time_in_sec = True
      reference_time = start_time
      timerange[0] = start_time - timestep
      timerange[1] = end_time + timestep
      str_start_time=PosTools.obtain_observation_year_month_day_hms(reference_time)
    else:
        timerange,str_start_time,reference_time=PosTools.get_time_range(start_time,end_time,timestep,time_in_sec,0)
    if str_start_time==-1:
       return
        
    if useEMM:
        #print ("USING EMM for EarthMagnetic Field")
        emm= EMM.EMM()
    else:
        emm= EMM.WMM()
    
    times,timerange=PosTools.getIONEXtimerange(timerange,timestep)
    if len(times[-1])==0 or times[-1][-1]<timerange[1]:
        timestmp=list(times[-1])
        timestmp.append(timerange[1]) #add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
        times[-1]=np.array(timestmp)
    timegrid=np.array([])
    TECs={};
    Bs={}
    Bpars={}
    air_mass={};
    RMs={};
    azimuth={};
    elevation={};
    flags={}
    if len(out_file) and len(myobject):
        out_file = out_file + '_' + myobject
    else:
      if len(myobject):
        out_file = 'RMextract_report_' + myobject
# print header section of report
    if len(out_file):
       if os.path.exists(out_file):
         os.remove(out_file)
       log = open(out_file, 'a')
       log.write ('Observing %s\n' % myobject)
       if not (MS is None):
           log.write ('Using measurement set %s\n' % MS)
       if use_azel:
         log.write ('observing at a fixed azimuth and elevation \n')
       log.write ('station_positions %s \n' % stat_pos)
       log.write ('\n')

    for st in stat_names:
        azimuth[st] = []
        elevation[st] = []
        air_mass[st] = []
        Bs[st]=[]
        Bpars[st]=[]
        RMs[st]=[]
        TECs[st]=[]
        air_mass[st]=[]
        flags[st]=[]
    for time_array in times:
        #get RM per timeslot
        starttime=time_array[0]
        #print ("getting ionexfile for",starttime)
        date_parms =  PosTools.obtain_observation_year_month_day_fraction(starttime)
        dayofyear = date(date_parms[0],date_parms[1],date_parms[2]).timetuple().tm_yday  
        emm.date=date_parms[0]+float(dayofyear)/365.
        #get relevant ionex file
        if not use_proxy:
            if not "http" in server: #ftp server use ftplib
                ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath)
            else:
                ionexf=ionex.get_urllib_IONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath)
        else:
                ionexf=ionex.get_urllib_IONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath,proxy_server=proxy_server,proxy_type=proxy_type,proxy_port=proxy_port,proxy_user=proxy_user,proxy_pass=proxy_pass)

        assert (ionexf!=-1),"error getting ionex data"
           
        tecinfo=ionex.readTEC(ionexf,use_filter=use_filter)
        if use_mean:
          if not stat_pos_mean:
            stn_mean = stat_pos.mean(0)
            stat_pos = []
            stat_pos.append(stn_mean)
            stat_pos_mean = True
            #print ('stat_pos.mean', stn_mean)
        for station,position in  zip(stat_names,stat_pos):
          #print ('generating data for station ', station)
 
          
          for time in time_array:
            result =  PosTools.obtain_observation_year_month_day_fraction(time)
            part_of_day= result[3] * 24
            if use_azel:
              az = pointing[0]
              el = pointing[1]
              flags[station].append(1)   
            else:
               az,el = PosTools.getAzEl(pointing,time,position,ha_limit)
               if az==-1 and el==-1:
                  return
               if az==999 and el==999:
                    #outsite hadec range
                  air_mass[station].append(-999)
                  azimuth[station].append(-999)
                  elevation[station].append(-999)
                  Bs[station].append([-999,-999,-999])
                  Bpars[station].append(-999)
                  RMs[station].append(-999)
                  TECs[station].append(-999)  #sTEC value
                  flags[station].append(0)
                  continue
            flags[station].append(1)   
            latpp,lonpp,height,lon,lat,am1=PosTools.getlonlatheight(az,el,position)

            if latpp==-1 and lonpp==-1 and height==-1:
               return
            #get VTEC from IONEX interpolation
            vTEC=ionex.getTECinterpol(time=part_of_day,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)

            TECs[station].append(vTEC*am1)  #STEC value
            emm.lon = lonpp
            emm.lat = latpp
            emm.h= ION_HEIGHT / 1.0e3
            Bpar=-1*emm.getProjectedField(lon,lat)# minus sign since the radiation is towards the Earth
            BField=emm.getXYZ()
            Bs[station].append(BField)

            air_mass[station].append(am1)
            azimuth[station].append(az)
            elevation[station].append(el)
            Bpars[station].append(Bpar)
            #calculate RM constant comes from VtEC in TECU,B in nT, RM = 2.62
            RMs[station].append((Bpar*vTEC*am1)*2.62e-6)
          if use_mean:
            break

            
        timegrid=np.concatenate((timegrid,time_array));
    
    for st in stat_names:
      air_mass[station] = np.array(air_mass[st])
      azimuth[station] = np.array(azimuth[st])
      elevation[station] = np.array(elevation[st])
      TECs[st]=np.array(TECs[st])
      Bs[st]=np.array(Bs[st])
      Bpars[st]=np.array(Bpars[st])
      RMs[st]=np.array(RMs[st]) 
      flags[st]=np.array(flags[st])
      if use_mean:
        break
       
    big_dict={}
    big_dict['STEC']=TECs
    big_dict['Bpar']=Bpars
    big_dict['BField']=Bs
    big_dict['AirMass']=air_mass
    big_dict['elev']=elevation
    big_dict['azimuth']=azimuth
    big_dict['RM']=RMs 
    big_dict['times']=timegrid
    big_dict['timestep']=timestep
    big_dict['station_names'] = stat_names
    big_dict['stat_pos'] = stat_pos
    big_dict['flags'] = flags
    big_dict['reference_time'] = reference_time 
    # finish writing computed data to report
    if len(out_file):
       time_range=[big_dict['times'][0],big_dict['times'][-1]]
       log.write ('start and end times %s %s \n' % (time_range[0], time_range[1]))
       log.write ('reference time for rel_time=0: year month day hr min sec %s %s %s %s %s %s \n' % str_start_time)
       if use_azel:
         log.write ('observing at azimuth and elevation %s %s \n' % (pointing[0], pointing[1]))
       else:
         log.write ('observation direction %s %s \n' % (pointing[0], pointing[1]))
       log.write ('\n')
       k = 0
       for key in big_dict['station_names']:
         seq_no = 0
         if use_mean:
           log.write ('data for station mean position at %s\n' % (stat_pos[k]))
         else:
           log.write ('data for station %s  at position %s\n' % (key, stat_pos[k]))
         log.write ('seq  rel_time time_width El         Az         STEC           RM (rad/m2)   VTEC factor  \n')
         for i in range (timegrid.shape[0]):
           el = big_dict['elev'][key][i]
           if el < 0 :
             ok = 1
             stec = 0.0
             rm = 0.0
             vtec_factor = 1.0
           else:
             ok = 0
             stec =big_dict['STEC'][key][i]
             rm = big_dict['RM'][key][i]
             vtec_factor = 1.0 / big_dict['AirMass'][key][i]
           az = big_dict['azimuth'][key][i]
           rel_time = timegrid[i] - reference_time
           if i  == 0:
             time_width = reference_time - timegrid[i] 
           else:
             time_width = timegrid[i] - timegrid[i-1]
           log.write("%s : %s %s %s %s %s %s %s %s\n" % (seq_no, ok, rel_time, time_width, el, az, stec, rm, vtec_factor))
           seq_no = seq_no + 1
         k = k + 1
         if use_mean:
           break
         log.write (' \n')
       log.close()
       print ('****** finished ionosphere predictions report: ', out_file)
    else:
      print ('*********** finished ionosphere predictions ***************')


    return big_dict


