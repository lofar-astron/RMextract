#!/usr/bin/env python
import PosTools
import getIONEX as ionex;
import os
import numpy as np
from math import *
try:
   import cPickle as pickle
except:
   import pickle

HAS_PYRAP = True
try:
  from pyrap import tables as tab
  import pyrap.quanta as qu
  from pyrap.measures import measures
  print ('pyrap will be used to compute positions')
except:
  print ('We will need PyEphem to perform calculations!')
  print ('the accuracy of results might decease a bit')
  HAS_PYRAP = False

#HAS_PYRAP = False

print ('using pyrap?', HAS_PYRAP)

HAS_EPHEM = True
try:
  import ephem
  if not HAS_PYRAP:
    print ('PyEphem will be used to perform calculations!')
except:
  HAS_EPHEM = False

import math

# the following parameter is used to extend the observation range by 120 sec
# before and after the actual specified time. If we want to correct an
# actual data set, this is required for the scipy 1d interpolation routine to
# ensure that the calculated range of data exceeds the time range actually
# observed - I have used the same value in ALBUS - Tony
TIME_OFFSET = 120.0

# height of thin layer ionosphere
ION_HEIGHT=450.e3;

########################################
sm_a = 6378137.0
invf = 298.257223563
f = 1.0 / invf

# The following two functions were provided by Emil Lenc
# Convert latitude (radians), longitude (radians) and elevation (metres) 
# to ITRF XYZ and vice versa
def WGS84ToITRF(lat, lon, h): # WGS-84 to ITRF (input in radians)
	SINK = math.sin(lat)
	COSK = math.cos(lat)
	e2 = 2.0 * f - f * f
	v = sm_a / math.sqrt(1.0 - e2 * SINK * SINK)
	x = (v + h) * COSK * math.cos(lon)
	y = (v + h) * COSK * math.sin(lon)
	z = ((1 - e2) * v + h) * SINK
	return x, y, z

def ITRFToWGS84(x, y, z):
	e2 = 2.0 * f - f * f
	E = e2 / (1.0 - e2)
	b = sm_a * (1.0 - f)
	p = math.sqrt(x * x + y * y)
	q = math.atan2(z * sm_a, (p * b))
	lat = math.atan2((z + E * b * math.sin(q) * math.sin(q) * math.sin(q)), (p - e2 * sm_a * math.cos(q) * math.cos(q) * math.cos(q)))
	v = sm_a / math.sqrt(1.0 - e2 * math.sin(lat) * math.sin(lat))
	lon = math.atan2(y, x)
	h = (p / math.cos(lat)) - v
        lat = math.degrees(lat)
        lon = math.degrees(lon)
	return lat, lon, h       # output in degrees here

########################################
# convert geodetic latitude to geocentric latitude
# input and output in radians
def GeodeticToGeocentricLat(geodetic_lat, height):
  l_sin = math.sin(geodetic_lat)
  e2 = 2.0 * f - f * f
  div = math.sqrt(1 - e2 * l_sin**2)
  rn =  sm_a / div 
  rn_div = rn + height
  ratio = 1 - e2 * rn / rn_div
  tan_geocentric_lat = ratio * math.tan(geodetic_lat) 
  geocentric_lat = math.atan(tan_geocentric_lat)
  return geocentric_lat

########################################
# see http://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt

# input expected in degrees here
def aer2ecef(azimuthDeg, elevationDeg, slantRange, obs_lat, obs_lon, obs_alt):
    obs_lat_r = math.radians(obs_lat)
    obs_lon_r = math.radians(obs_lon)
    sitex, sitey, sitez = WGS84ToITRF(obs_lat_r,obs_lon_r,obs_alt)

    #some needed calculations
    slat = math.sin(math.radians(obs_lat))
    slon = math.sin(math.radians(obs_lon))
    clat = math.cos(math.radians(obs_lat))
    clon = math.cos(math.radians(obs_lon))

    azRad = math.radians(azimuthDeg)
    elRad = math.radians(elevationDeg)

    # az,el,range to sez convertion
    south  = -slantRange * math.cos(elRad) * math.cos(azRad)
    east   =  slantRange * math.cos(elRad) * math.sin(azRad)
    zenith =  slantRange * math.sin(elRad)


    x = ( slat * clon * south) + (-slon * east) + (clat * clon * zenith) + sitex
    y = ( slat * slon * south) + ( clon * east) + (clat * slon * zenith) + sitey
    z = (-clat *        south) + ( slat * zenith) + sitez
    lat,lon,h = ITRFToWGS84(x, y, z)
    return lat, lon, h  # data are in units of degrees

################ from JMA_tools in ALBUS ########################
def get_ymdf_from_JD(JD):
    """get the year, month, day, day_fraction from an MJD

Taken from _Astronomical Algorithms_, Meeus, 1991
    """
    JD2 = JD + 0.5
    Z = int(JD2)
    F = JD2 - Z
    A = Z
    if(Z>= 2299161):
        alpha = int((Z-1867216.25)/36524.25)
        A = Z + 1 + alpha - int(alpha/4)
    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)
    day_total = B - D - int(30.6001*E) + F
    month = E - 1
    if(E >= 14): month = E - 13
    year = C - 4716
    if(month <= 2): year += 1
    day = int(day_total)
    day_fraction = day_total - day
    return year, month, day, day_fraction

################################################################################
def get_hms_from_frac(day_fraction):
    """get hours, minues, seconds from a fractional day.

Does not worry about leap seconds.
"""
    h = day_fraction * 24.0
    hour = int(h+2E-13)
    m = (h - hour) * 60.0
    minute = int(m+1E-11)
    second = (m - minute) * 60.0
    return hour, minute, second


################################################################################
def get_ymdh_from_JD(JD):
    """get hours, minues, seconds from a fractional day.
"""
    year, month, day, day_fraction = get_ymdf_from_JD(JD)
    hour, minute, second = get_hms_from_frac(day_fraction)
    return year, month, day, hour, minute, second

################################################################################
def obtain_observation_year_month_day_fraction(start_time):
    julian_day = (start_time / 86400.0) + 2400000.5
    result = get_ymdf_from_JD(julian_day)
    return result

################################################################################
# Get the day of year from the Year, month, day for the start of observations
def obtain_observation_year_month_day_hms(start_time):
    if HAS_PYRAP:
      date_list = qu.quantity(str(start_time)+'s').formatted("YMD").split("/")
      year = int(date_list[0])
      month = int(date_list[1])
      day = int(date_list[2])
      time_list=date_list[3].split(":")
      return (year, month, day,int(time_list[0]),int(time_list[1]),float(time_list[2]))
    else:
      julian_day = (start_time / 86400.0) + 2400000.5
      year, month, day, hour, minute, second  = get_ymdh_from_JD(julian_day)
      return (year, month, day,hour, minute, second)

#####################  main processing function #####################
def getRM(MS=None,
           server="ftp://ftp.unibe.ch/aiub/CODE/",
           prefix='CODG',
           ionexPath="IONEXdata/",
           earth_rot=0,
           timerange=0,
           use_azel = False,
           **kwargs):
    '''optional arguments are:
    radec or pointing : [ra,dec] in radians,
    timestep in s, timerange = [start, end]in MJD seconds (sorry I shall write the converter soon...), 
    stat_names = list of stings per station, 
    stat_positions = list of length 3 numpy arrays, containing station_position in ITRF meters.
    useEMM = boolean, use EMM for Earth magnetic field, you should specify the path to the coefficients via
    EMMpath = path. If the path is given it will be automatically assumed that you want to use EMM for Earthm agnetic field 
    useWMM = boolean, use WMM for Earth magnetic Field (see above). You can specify the path to coefficients with
    WMMpath= path.
    Returns the (timegrid,timestep,TEC) where TEC is a dictionary containing 1 enumpyarray per station in stat_names. 
    If stat_names is not given, the station names will either be extracted from the MS or st1...stN '''

    print ('earth_rot', earth_rot)
    print ('timerange ', timerange)
    print ('use_azel ', use_azel)
    stat_names=[]
    useWMM=False
    useEMM=False
#   EMMpath='/Users/elenc/local/src/RMextract/EMM/'
    EMMpath = './'
    object = ''
    timestep=60
    out_file = 'RMextract_report'
    pointing=[0.,np.pi]
    if not (MS is None):

      (timerange,timestep,pointing,stat_names,stat_pos)=PosTools.getMSinfo(MS);
    
    for key in kwargs.keys():
        if key=='radec' or key=='pointing':
            pointing = kwargs[key]
        if key=='object':
            object = kwargs[key]
            out_file = 'RMextract_report_' + object
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
        if key=='stat_positions':
            stat_pos=kwargs[key]
            if not stat_names:
                stat_names =['st%d'%(i+1) for i in range(len(stat_pos))]
        if key=='EMMpath':
            EMMpath=kwargs[key]
            if not useWMM:
                useEMM=True
        if key=='WMMpath':
            EMMpath=kwargs[key]
            if not useEMM:
                useWMM=True
        if key=='useWMM':
            useWMM=kwargs[key]
        if key=='useEMM':
            useEMM=kwargs[key]
    
    if timerange != 0:
      start_time = timerange[0]
      end_time = timerange[1]
      time_in_sec = True

    if timerange ==0 and HAS_PYRAP:
      try:
        start_time = qu.quantity(start_time).get_value()
        end_time = qu.quantity(end_time).get_value()
        print ('**** specified start and end time ', start_time, end_time)
        reference_time = start_time * 86400.0 - TIME_OFFSET
        st = reference_time - timestep 
        et = end_time * 86400.0 +  timestep + TIME_OFFSET
        str_start_time =  obtain_observation_year_month_day_hms(reference_time)
        timerange= [st, et]
      except:
        print ('no time range given')
        print ('exiting')
        return
    elif HAS_EPHEM:
      if time_in_sec:
        dublin_start = start_time / 86400.0 -15019.5
        dublin_end = end_time / 86400.0 -15019.5
        start_time = ephem.julian_date(ephem.Date(dublin_start)) - 2400000.5
        end_time = ephem.julian_date(ephem.Date(dublin_end)) - 2400000.5
      else:
        start_time = ephem.julian_date(ephem.Date(start_time)) - 2400000.5
        end_time = ephem.julian_date(ephem.Date(end_time)) - 2400000.5
      print ('ephem start and end time ', start_time, end_time)
      reference_time = start_time * 86400.0 - TIME_OFFSET
      st = reference_time - timestep 
      et = end_time * 86400.0 + timestep + TIME_OFFSET
      str_start_time =  obtain_observation_year_month_day_hms(reference_time)
      timerange= [st, et]
    else:
      print ('unable to get time range so exiting!')
      return

# print header section of report
    if os.path.exists(out_file):
      os.remove(out_file)
    log = open(out_file, 'a')
    log.write ('Observing %s\n' % object)
    log.write ('station_positions %s \n' % stat_pos)
    if use_azel:
      log.write ('observing at a fixed azimuth and elevation \n')
    log.write ('start and end times %s %s \n' % (st, et))
    log.write ('reference time for rel_time=0: year month day hr min sec %s %s %s %s %s %s \n' % str_start_time)
    log.write ('\n')

    if useWMM  or useEMM :
        import EMM.EMM as EMM
        
        if useEMM:
            print ("USING EMM for EarthMagnetic Field, remember to have set the path to the coeffiecients (EMMpath) correctly.")
            emm=EMM.EMM(cof=EMMpath)
        else:
            print ("USING WMM for EarthMagnetic Field, remember to have set the path to the coeffiecients (WMMpath) correctly.")
            emm=EMM.WMM(cof=EMMpath)
    
    else:
        print ("USING obsolete casacore libraries for getting Earthmagnetic field!")
        import getEarthMagnetic as EM


     
    #IONEX files go per day, check if more than one file is  needed.
    times=[];
    
    while timerange[0] + 86400.0 <= timerange[1]:
        result =  obtain_observation_year_month_day_fraction(timerange[0])
        part_of_day = result[3] * 24.0 * 60 * 60
        nr_remaining_seconds=(24.*60.*60.-part_of_day);
        times.append(np.arange(timerange[0],timerange[0]+nr_remaining_seconds,timestep));
        timerange[0]+=len(times[-1])*timestep;
    print('final timerange ', timerange)
    times.append(np.arange(timerange[0],timerange[1]+timestep,timestep)) #add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
    timegrid=np.array([])
    TECs={};
    Bs={};
    air_mass={};
    RMs={};
    azimuth={};
    elevation={};
    for st in stat_names:
        azimuth[st] = []
        elevation[st] = []
        air_mass[st] = []
        Bs[st]=[]
        RMs[st]=[]
        TECs[st]=[]
        air_mass[st]=[]
    for time_array in times:
        #get RM per timeslot
        starttime=time_array[0];
        date_parms =  obtain_observation_year_month_day_fraction(starttime)
        #get relevant ionex file
        ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath);
        tecinfo=ionex.readTEC(ionexf)
        for station,position in  zip(stat_names,stat_pos):
          print ('generating data for station ', station)
          if not HAS_PYRAP and HAS_EPHEM:
            location_lat, location_lon, location_height = ITRFToWGS84(position[0], position[1], position[2])
            location = ephem.Observer()
            # convert geodetic latitude to geocentric
            # flattening, f, defined above for WGS84 stuff
            geodet_lat = math.radians(location_lat)
            tan_geocentric_latitude =  math.tan(geodet_lat) * (1 - f) **2
            geocentric_latitude = GeodeticToGeocentricLat(geodet_lat, location_height)
            location.lat = geocentric_latitude
            location.lon = math.radians(location_lon)
            location.elevation = location_height
            location.pressure = 0.0
          for time in time_array:
            result =  obtain_observation_year_month_day_fraction(time)
            part_of_day= result[3] * 24
            if use_azel:
              az = radians(180.0)
              el = radians(90.0)
            elif HAS_PYRAP:
                azel=PosTools.radec2azel(pointing[0],pointing[1],time=str(time)+'s',pos=position);
                az=azel['m0']['value'];
                el=azel['m1']['value'];
#               print ('pyrap derived azel = ', math.degrees(az), math.degrees(el))
            elif HAS_EPHEM:
              # convert to Dublin Julian Date for PyEphem
              location.date =  time/86400.0 - 15019.5
              lst = location.sidereal_time()
              equatorial = ephem.Equatorial(str(12.0 * math.degrees(pointing[0])/180),str(math.degrees(pointing[1])))
              body = ephem.FixedBody()
              body._ra = equatorial.ra
              body._dec = equatorial.dec
              body._epoch = equatorial.epoch
              body.compute(location)
              az = math.degrees(body.az)   * math.pi / 180.0
              el = math.degrees(body.alt) * math.pi / 180.0
            else: 
              print ('failure to get azimuth and elevation! Exiting!')
              return
            if HAS_PYRAP:
              lonlat=PosTools.getLonLatStation(az,el,pos=position);

              lon=lonlat['m0']['value'];
              lat=lonlat['m1']['value'];
              # convert to itrf coordinates on sphere with radius 1
              diritrf=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
              # calculate piercepoint in xyz(code from Bas vd Tol)
            
              pos_pp = (ppx1,ppy1,ppz1,am1)=PosTools.getPP(h=ION_HEIGHT,mPosition=position,direction=diritrf)
    
              #get pp in lon,lat h
              me=measures();
              pp1position=me.position("ITRF",str(ppx1)+'m',str(ppy1)+'m',str(ppz1)+'m')
            
              lonpp = degrees(pp1position['m0']['value']);
              latpp = degrees(pp1position['m1']['value']);
              height = pp1position['m2']['value'] / 1.0e3
            elif HAS_EPHEM:
              slantRange = 3.0e20    # seem to need a large value to
                                       # get near the WNB measures equivalent
#             lat,lon,ht = aer2ecef(math.degrees(body.az), math.degrees(body.alt), slantRange, math.degrees(geocentric_latitude), location_lon, ION_HEIGHT)
              lat,lon,ht = aer2ecef(math.degrees(az), math.degrees(el), slantRange, location_lat, location_lon, ION_HEIGHT)
              # convert to itrf coordinates on sphere with radius 1
              lat = math.radians(lat)
              lon = math.radians(lon)
              diritrf=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
              # calculate piercepoint in xyz(code from Bas vd Tol)
            
              pos_pp = (ppx1,ppy1,ppz1,am1)=PosTools.getPP(h=ION_HEIGHT,mPosition=position,direction=diritrf)
              #get pp in lon,lat h
              latpp, lonpp, height =  ITRFToWGS84(ppx1,ppy1,ppz1)
            else:
              print ('unable to compute position parameters - exiting!')
              return

            #get VTEC from IONEX interpolation
            vTEC=ionex.getTECinterpol(time=part_of_day,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)

            TECs[station].append(vTEC*am1)  #STEC value
            if useEMM or useWMM:
                emm.lon = lonpp
                emm.lat = latpp
                emm.h= ION_HEIGHT / 1.0e3
                Bpar=-1*emm.getProjectedField(lon,lat)# minus sign since the radiation is towards the Earth
            else:
                # get IGRF BField and project along LOS
                BField=EM.getField([ppx1,ppy1,ppz1],time=str(time)+'s');
                Bpar=-1*EM.ProjectField(BField,lon,lat); # minus sign since the radiation is towards the Earth
            air_mass[station].append(am1)
            azimuth[station].append(az)
            elevation[station].append(el)
            Bs[station].append(Bpar)
            #calculate RM constant comes from VtEC in TECU,B in nT, RM = 2.62
            RMs[station].append((Bpar*vTEC*am1)*2.62e-6)

            
        timegrid=np.concatenate((timegrid,time_array));
    
    for st in stat_names:
      air_mass[station] = np.array(air_mass[st])
      azimuth[station] = np.array(azimuth[st])
      elevation[station] = np.array(elevation[st])
      TECs[st]=np.array(TECs[st])
      Bs[st]=np.array(Bs[st])
      RMs[st]=np.array(RMs[st]) 

        
    big_dict={}
    big_dict['STEC']=TECs
    big_dict['Bpar']=Bs
    big_dict['AirMass']=air_mass
    big_dict['elev']=elevation
    big_dict['azimuth']=azimuth
    big_dict['RM']=RMs 
    big_dict['times']=timegrid
    big_dict['timestep']=timestep
    big_dict['station_names'] = stat_names

    # finish writing computed data to report
    k = 0
    for key in big_dict['station_names']:
      seq_no = 0
      log.write ('data for station %s  at position %s\n' % (key, stat_pos[k]))
      log.write ('seq  rel_time time_width El         Az         STEC           RM (rad/m2)   VTEC factor  \n')
      for i in range (timegrid.shape[0]):
        el = big_dict['elev'][key][i]
        if el < 0 :
          ok = 1
        else:
          ok = 0
        az = big_dict['azimuth'][key][i]
        stec =big_dict['STEC'][key][i]
        rm = big_dict['RM'][key][i]
        vtec_factor = 1.0 / big_dict['AirMass'][key][i]
        rel_time = timegrid[i] - reference_time
        if i  == 0:
          time_width = reference_time - timegrid[i] 
        else:
          time_width = timegrid[i] - timegrid[i-1]
        log.write("%s : %s %s %s %s %s %s %s %s\n" % (seq_no, ok, rel_time, time_width, el, az, stec, rm, vtec_factor))
        seq_no = seq_no + 1
      k = k + 1
      log.write (' \n')
      
    print ('*********** finished ionosphere predictions ***************')

#   # save in pickle file for other uses
#   pickle_name = 'RMextract_data_for_' + object
#   outputx = open(pickle_name, 'wb')
#   pickle.dump(result, outputx)
#   outputx.close()

    return big_dict


#############################

if __name__ == "__main__":


#Observing  B0329+54_ARO
#ALBUS data processing option: RI_G03
#station position [[  918014.88108001 -4346162.2500943   4561997.40582786]]
#Potential number of receivers for GPS fit  177
#Final number of receivers for GPS fit  50
#set ground station position to 918014.88108 -4346162.25009 4561997.40583
#reference year,month,day,hr,min,sec  2014 7 31 23 58 0.0
#start and end times  4913567880.0 4913740919.0
#observation direction  0.934198390438 0.953421252894
  START_TIME="2014/08/01 00:00:00"
  END_TIME="2014/08/01 23:59:59"


# result = getRM(timestep=300.0,timerange=[4913567880.0, 4913740919.0], radec=[0.934198390438, 0.953421252894], stat_positions= [[918014.88108001, -4346162.2500943,   4561997.40582786]],stat_names=['1'], useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )

# result = getRM(timestep=300.0,timerange=[4913567880.0, 4913740919.0], radec=[0.934198390438, 0.953421252894], stat_positions= [[918014.88108001, -4346162.2500943,   4561997.40582786]],stat_names=['1'] )

# result = getRM(timestep=300.0,start_time = START_TIME, end_time=END_TIME, radec=[0.934198390438, 0.953421252894], stat_positions= [[918014.88108001, -4346162.2500943,   4561997.40582786]],stat_names=['1'], useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )


#Observing  B0329+54_DRAO
  START_TIME="2015/02/08 00:00:00"
  END_TIME="2015/02/11 23:59:59"
# result = getRM(timestep=300.0,start_time=START_TIME, end_time=END_TIME, radec=[0.934198390438, 0.953421252894], stat_positions= [[-2059164.91026054, -3621108.42669587,  4814432.26620199]],stat_names=['1'], useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )

#Observing  3C196

# result = getRM(use_azel=True,object='3C196',earth_rot=0,timestep=3600.0,timerange=[4861550345.0, 4861579375.0], radec=[2.1379861549, 0.844194330425],stat_positions= [[3826578.16006, 461022.741786, 5064892.01151]],stat_names=['1'], useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )

# Observing  Bob Sault object at ASKAP
  START_TIME="2014/12/09 00:00:00"
  END_TIME="2014/12/09 23:59:59"

# print (' ******* STARTING UP SAULT_IONO_TESTT******** ')
# result = getRM(object='B_Sault',start_time=START_TIME,end_time=END_TIME,timestep=300.0,radec=[ 5.15199451295,-1.11138596143],stat_positions= [[-2558271.67083283,5095662.32694705,-2849025.79668766]],stat_names=['1'],useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )

  OBJECT="EoR0"
  START_TIME="2013/11/28 00:00:00"
  END_TIME="2013/11/28 01:59:59"
  TIME_STEP = 300.0

  MWA_antennas = np.array([[-2559314.23084924,5095535.90961438,-2848889.57667157],
   [-2559293.10717106,5095498.79164383,-2848974.05801863],
   [-2559156.42269442,5095418.83230373,-2849233.34162414],
   [-2559094.63804600,5095526.84526066,-2849097.63284488],
   [-2559042.54106487,5095566.88538445,-2849072.42535023],
   [-2559068.53757350,5095654.59288871,-2848892.60844473],
   [-2559161.70851932,5095607.73286033,-2848894.91011893],
   [-2559173.78330034,5095643.10464650,-2848820.20086397],
   [-2559782.26851932,5095304.54438001,-2848875.78669410],
   [-2559644.22829760,5095369.93521424,-2848885.21756417],
   [-2559507.77695003,5095490.23883646,-2848797.39833977],
   [-2559467.43177484,5095508.01973328,-2848802.30233654],
   [-2559460.59086333,5095515.74910944,-2848794.76371318],
   [-2559491.68457220,5095527.67486954,-2848745.44773601],
   [-2559603.60609646,5095563.73884050,-2848579.72258876],
   [-2559631.28428317,5095541.41922988,-2848594.98830325],
   [-2559113.92023486,5095854.59042124,-2848492.05455485],
   [-2559133.51844911,5095831.00304170,-2848517.14718873],
   [-2559018.96896708,5095793.67783611,-2848686.69023686],
   [-2558906.48396095,5095592.28259425,-2849148.93562390],
   [-2558894.77687225,5095720.00453191,-2848930.82517056],
   [-2558880.58102582,5095762.06255238,-2848868.27661380],
   [-2558503.88881043,5095891.11710898,-2848981.31195756],
   [-2558648.85477276,5096060.47633611,-2848544.49069260],
   [-2558998.73468649,5095390.06352995,-2849423.09595365],
   [-2559238.04568324,5095263.75775157,-2849432.88470164],
   [-2558856.49159020,5095257.96516587,-2849788.57821277],
   [-2558761.92575271,5095281.91134845,-2849829.99130606],
   [-2558719.21221208,5095416.28342253,-2849628.99110746],
   [-2558836.79342206,5095555.42415917,-2849277.33903756],
   [-2558850.45931999,5095586.71918979,-2849209.71070222],
   [-2558890.31919482,5095521.92810583,-2849288.42518348]])
  result = getRM(start_time=START_TIME,end_time=END_TIME,object=OBJECT, timestep=300.0,stat_positions= MWA_antennas,useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )
# result = getRM(use_azel= True, timerange=[4892313480.0, 4892572919.0],object=OBJECT, timestep=300.0,stat_positions= MWA_antennas,useEMM=True,EMMpath='/home/twillis/ASKAP_related_ionosphere/RMextract/EMM/' )

# process_ionosphere(telescope_pos=MWA_antennas,Ra=RA,Dec=DEC,start_time=START_TIME,end_time=END_TIME,processing_option="RI_G03",do_serial=1,num_processors=6, gps_data_directory="MWA_Beta2_agw_data_store",object=OBJECT)
# process_ionosphere(telescope_pos=MWA_antennas,start_time=START_TIME,end_time=END_TIME,El=EL,processing_option="RI_G03",do_serial=1,gps_data_directory="/home/twillis1/MWA_Beta2_emil_gw_data_store",use_global_data=0,object=OBJECT)


# set time relative to zero
  timegrid = result['times']
  timegrid = (timegrid - timegrid[0]) / 3600.0
  RMs = result['RM'] 
  TECs = result['STEC']
  times = result['times']
  selected_key = False
  for key in TECs.keys():
    if not selected_key:
      use_key = key
      selected_key = True
  plot_RMs = False
  STEC = True

  from pylab import *
  if plot_RMs:
    plot(timegrid, RMs[use_key],'ro')
  else:
    plot(timegrid, TECs[use_key],'ro')
  xlabel('relative time (hours)')
  if plot_RMs:
    ylabel('RM (radians/m^2)')
    title_string ='RMextract: RM as a function of time'
    plot_file =  'RMextract_RM_plot'
  else:
    if STEC:
      ylabel('STEC (TEC unit)')
      title_string ='CODE: STEC as a function of time'
      plot_file =  'STEC_plot'
    else:
      ylabel('VTEC (TEC unit)')
      title_string ='CODE: VTEC as a function of time'
      plot_file =  'VTEC_plot'
  title(title_string)
  grid(True)

  savefig(plot_file)
  show()


