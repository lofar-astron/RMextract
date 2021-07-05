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

HAS_EPHEM = True
try:
  import ephem
  if not HAS_PYRAP:
    print ('PyEphem will be used to perform calculations!')
except:
  HAS_EPHEM = False

from math import *
import math
import os
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz,Angle,ITRS,FK5
from scipy import interpolate

# height of thin layer ionosphere
ION_HEIGHT=450.e3;


R_earth=6364.62e3;
R_earthkm=6364.62;
earth_ellipsoid_a = 6378137.0;
earth_ellipsoid_a2 = earth_ellipsoid_a*earth_ellipsoid_a;
earth_ellipsoid_b = 6356752.3142;
earth_ellipsoid_b2 = earth_ellipsoid_b*earth_ellipsoid_b;
earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2;
posCS002=[3826577.1095  ,461022.900196, 5064892.758]
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
def get_day_of_year(year, month, day):
    """Get the day-of-year integer from the year/month/day

    year   I  YYYY
    month  I  MM starting from January == 1
    day    I  DD starting from 1 == 1
    """
    day_of_year_list = [0,31,59,90,120,151,181,212,243,273,304,334]
    doy = day_of_year_list[month-1] + day
    if(month>2):
        if((year&0x3)==0):
            if((year % 100 != 0) or (year % 400 == 0)):
                doy = doy+1
    return doy
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
      #print ("getting time", str(start_time)+'s')
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





def getMSinfo(MS=None):
    if MS is None:
        print ("No measurement set given")
        return
    if not HAS_PYRAP:
        print ("Install pyrap to be able to extract info from MS")
        return

        
    if os.path.isdir(MS):
        myMS=tab.table(MS)
    else:
        print ("Do not understand the format of MS",MS,"bailing out")
        return;
    timerange=[np.amin(myMS.getcol('TIME_CENTROID')),np.amax(myMS.getcol('TIME_CENTROID'))]
    timestep=myMS.getcell('INTERVAL',0)
    
    pointing= tab.table(myMS.getkeyword('FIELD')).getcell('PHASE_DIR',0);    
    stations = tab.table(myMS.getkeyword('ANTENNA')).getcol('NAME')
    station_pos = tab.table(myMS.getkeyword('ANTENNA')).getcol('POSITION')

    return (timerange,timestep,pointing.flatten(),stations,station_pos)


def getPPsimple(height=[ION_HEIGHT,],mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
    '''get piercepoints for antenna position mPosition in m, direction ITRF in m on unit sphere and for array of heights, assuming a spherical Earth'''
    height=np.array(height)
    stX=mPosition[0]
    stY=mPosition[1]
    stZ=mPosition[2]
    x=np.divide(stX,(R_earth+height))
    y=np.divide(stY,(R_earth+height))
    z=np.divide(stZ,(R_earth+height))
    
    c = x*x + y*y + z*z - 1.0;
   
    dx=np.divide(direction[0],(R_earth+height))
    dy=np.divide(direction[1],(R_earth+height))
    dz=np.divide(direction[2],(R_earth+height))

    a = dx*dx + dy*dy + dz*dz;
    b = x*dx + y*dy  + z*dz;

    alpha = (-b + np.sqrt(b*b - a*c))/a;


    pp=np.zeros(height.shape+(3,))
    pp[:,0]=stX+alpha*direction[0]
    pp[:,1]=stY+alpha*direction[1]
    pp[:,2]=stZ+alpha*direction[2]

    am=np.divide(1.,pp[:,0]*dx+pp[:,1]*dy+pp[:,2]*dz)
    return pp,am

def getPPsimpleAngle(height=[ION_HEIGHT,],mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
    '''get (lon,lat,h values) of piercepoints for antenna position mPosition in m, direction ITRF in m on unit sphere and for array of heights, assuming a spherical Earth'''
    height=np.array(height)
    stX=mPosition[0]
    stY=mPosition[1]
    stZ=mPosition[2]
    x=np.divide(stX,(R_earth+height))
    y=np.divide(stY,(R_earth+height))
    z=np.divide(stZ,(R_earth+height))
    
    c = x*x + y*y + z*z - 1.0;
   
    dx=np.divide(direction[0],(R_earth+height))
    dy=np.divide(direction[1],(R_earth+height))
    dz=np.divide(direction[2],(R_earth+height))

    a = dx*dx + dy*dy + dz*dz;
    b = x*dx + y*dy  + z*dz;

    alpha = (-b + np.sqrt(b*b - a*c))/a;


    pp=np.zeros(height.shape+(3,))
    pp[:,0]=stX+alpha*direction[0]
    pp[:,1]=stY+alpha*direction[1]
    pp[:,2]=stZ+alpha*direction[2]

    am=np.divide(1.,pp[:,0]*dx+pp[:,1]*dy+pp[:,2]*dz)

    ppl=np.zeros(height.shape+(3,))
    ppl[:,0]=np.arctan2(pp[:,1],pp[:,0])
    ppl[:,1]=np.arctan2(pp[:,2],np.sqrt(pp[:,0]*pp[:,0]+pp[:,1]*pp[:,1]))
    ppl[:,2]=height

    return ppl,am
    

def getPP(h=ION_HEIGHT,mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
    stationX = mPosition[0];
    stationY = mPosition[1];
    stationZ = mPosition[2];

    ion_ellipsoid_a = earth_ellipsoid_a + h;
    ion_ellipsoid_a2_inv = 1.0 / (ion_ellipsoid_a * ion_ellipsoid_a);
    ion_ellipsoid_b = earth_ellipsoid_b + h;
    ion_ellipsoid_b2_inv = 1.0 / (ion_ellipsoid_b * ion_ellipsoid_b);
    
    x = stationX/ion_ellipsoid_a;
    y = stationY/ion_ellipsoid_a;
    z = stationZ/ion_ellipsoid_b;
    c = x*x + y*y + z*z - 1.0;

    dx = direction [0]/ ion_ellipsoid_a;
    dy = direction [1] / ion_ellipsoid_a;
    dz = direction [2] / ion_ellipsoid_b;

    a = dx*dx + dy*dy + dz*dz;
    b = x*dx + y*dy  + z*dz;
    alpha = (-b + sqrt(b*b - a*c))/a;
    pp_x = stationX + alpha*direction[0];
    pp_y = stationY + alpha*direction[1]
    pp_z = stationZ + alpha*direction[2];

    normal_x = pp_x * ion_ellipsoid_a2_inv;
    normal_y = pp_y * ion_ellipsoid_a2_inv;
    normal_z = pp_z * ion_ellipsoid_b2_inv;
    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z;
    norm_normal = sqrt(norm_normal2);
    sin_lat2 = normal_z*normal_z / norm_normal2;

 
    g = 1.0 - earth_ellipsoid_e2*sin_lat2;
    sqrt_g = sqrt(g);

    M = earth_ellipsoid_b2 / ( earth_ellipsoid_a * g * sqrt_g );
    N = earth_ellipsoid_a / sqrt_g;

    local_ion_ellipsoid_e2 = (M-N) / ((M+h)*sin_lat2 - N - h);
    local_ion_ellipsoid_a = (N+h) * sqrt(1.0 - local_ion_ellipsoid_e2*sin_lat2);
    local_ion_ellipsoid_b = local_ion_ellipsoid_a*sqrt(1.0 - local_ion_ellipsoid_e2);

    z_offset = ((1.0-earth_ellipsoid_e2)*N + h - (1.0-local_ion_ellipsoid_e2)*(N+h)) * sqrt(sin_lat2);

    x1 = stationX/local_ion_ellipsoid_a;
    y1 = stationY/local_ion_ellipsoid_a;
    z1 = (stationZ-z_offset)/local_ion_ellipsoid_b;
    c1 = x1*x1 + y1*y1 + z1*z1 - 1.0;

    dx = direction[0] / local_ion_ellipsoid_a;
    dy = direction[1] / local_ion_ellipsoid_a;
    dz = direction[2] / local_ion_ellipsoid_b;
    a = dx*dx + dy*dy + dz*dz;
    b = x1*dx + y1*dy  + z1*dz;
    alpha = (-b + sqrt(b*b - a*c1))/a;

    pp_x = stationX + alpha*direction[0];
    pp_y = stationY + alpha*direction[1]
    pp_z = stationZ + alpha*direction[2];

    normal_x = pp_x / (local_ion_ellipsoid_a * local_ion_ellipsoid_a);
    normal_y = pp_y / (local_ion_ellipsoid_a * local_ion_ellipsoid_a);
    normal_z = (pp_z-z_offset) / (local_ion_ellipsoid_b * local_ion_ellipsoid_b);

    norm_normal2 = normal_x*normal_x + normal_y*normal_y + normal_z*normal_z;
    norm_normal = sqrt(norm_normal2);
    
    pp_airmass = norm_normal / (direction[0]*normal_x + direction[1]*normal_y + direction[2]*normal_z);

    return (pp_x,pp_y,pp_z,pp_airmass)


def getLonLat(pos):
    #converts ITRF pos in xyz to lon lat 
    me=measures()
    a=me.measure(me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m'),"ITRF");
    return (a['m0']['value'],a['m1']['value'])

def getLonLatStation(az=0,el=0,pos=posCS002):
    #gets converts local station direction to ITRF lon/lat
    if not isinstance(az,str):
        az=str(az)+'rad';
    if not isinstance(el,str):
        el=str(el)+'rad';
    me=measures()
    me.do_frame(me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m'))
    #me.do_frame(me.epoch('utc', 'today'))
    direction=me.direction("AZELGEO",az,el)
    return me.measure(direction,"ITRF");


def radec2azel(ra,dec,time, pos):
    me=measures();
    if type(ra)!=str:
        ra=str(ra)+'rad';
    if type(dec)!=str:
        dec=str(dec)+'rad';
    phasedir=me.direction('J2000',ra,dec)
    t=me.epoch("UTC",qu.quantity(time));
    me.do_frame(t);

    p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(p);
    #print ("input radec2azel",phasedir,ra,dec,p,t)
    azel = me.measure(phasedir,'AZELGEO');
    return azel;

def azel2radec(az,el,time, pos):
    me=measures();
    if type(az)!=str:
        az=str(az)+'rad';
    if type(el)!=str:
        el=str(el)+'rad';
    phasedir=me.direction('AZELGEO',az,el)
    t=me.epoch("UTC",qu.quantity(time));
    me.do_frame(t);

    p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(p);

    radec = me.measure(phasedir,'RADEC');
    return radec;


def getuvw(ra,dec,time, pos1,pos2):
    me=measures();
    if type(ra)!=str:
        ra=str(ra)+'rad';
    if type(dec)!=str:
        dec=str(dec)+'rad';
    phasedir=me.direction('J2000',ra,dec)
    me.do_frame(phasedir);
    t=me.epoch("UTC",qu.quantity(time));
    me.do_frame(t);
    
    p = me.position('ITRF',str(pos1[0])+'m',str(pos1[1])+'m',str(pos1[2])+'m')
    bl = me.baseline('ITRF',str(pos2[0]-pos1[0])+'m',str(pos2[1]-pos1[1])+'m',str(pos2[2]-pos1[2])+'m')
    #print (bl)
    me.do_frame(p);
    #print (me.to_uvw(bl)['xyz'])
    #return uvw;
     
def getIONEXtimerange(timerange,timestep):
    #IONEX files go per day, check if more than one file is  needed.
    times=[];
    oldtimerange=-100
    while timerange[0]< timerange[1] and timerange[0]>oldtimerange:
      oldtimerange=timerange[0]
      #print (timerange)
      result =  obtain_observation_year_month_day_fraction(timerange[0])
      part_of_day = result[3]
      result2 =  obtain_observation_year_month_day_fraction(timerange[1])
      if result2[2]==result[2]:  #sameday
        times.append(np.arange(timerange[0],timerange[1]+timestep,timestep)) #make sure to include the last timestep
      else:
        nr_remaining_seconds=(1.-part_of_day) * 24.0 * 60. * 60.;
        #print ("new day",nr_remaining_seconds)
        times.append(np.arange(timerange[0],timerange[0]+nr_remaining_seconds,timestep));
      #print ("in postools",len(times),times[-1],timerange,nr_remaining_seconds,result2,result,times[-1][-1],part_of_day)
      if len(times[-1]):
        timerange[0]=times[-1][-1]+timestep
    return times,timerange

def getlonlatheight(az,el,position,h=ION_HEIGHT):
    if HAS_PYRAP:
      lonlat=getLonLatStation(az,el,pos=position);

      lon=lonlat['m0']['value'];
      lat=lonlat['m1']['value'];
      # convert to itrf coordinates on sphere with radius 1
      diritrf=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
      # calculate piercepoint in xyz(code from Bas vd Tol)

      pos_pp = (ppx1,ppy1,ppz1,am1)=getPP(h=h,mPosition=position,direction=diritrf)

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
      location_lat, location_lon, location_height = ITRFToWGS84(position[0], position[1], position[2])
      lat,lon,ht = aer2ecef(math.degrees(az), math.degrees(el), slantRange, location_lat, location_lon, ION_HEIGHT)
      # convert to itrf coordinates on sphere with radius 1
      lat = math.radians(lat)
      lon = math.radians(lon)
      diritrf=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
      # calculate piercepoint in xyz(code from Bas vd Tol)

      pos_pp = (ppx1,ppy1,ppz1,am1)=getPP(h=ION_HEIGHT,mPosition=position,direction=diritrf)
      #get pp in lon,lat h
      latpp, lonpp, height =  ITRFToWGS84(ppx1,ppy1,ppz1)
    else:
      print ('unable to compute position parameters - exiting!')
      return -1,-1,-1
    return latpp, lonpp, height,lon,lat,am1

def far_to_near(altaz_far, distance, obstime=None):
    """Add distance to AltAz instance"""
    if obstime is None:
        obstime = altaz_far.obstime
    return AltAz(
        location=altaz_far.location,
        obstime=obstime,
        az=altaz_far.az,
        alt=altaz_far.alt,
        distance=distance,
    )

def distance_from_height(obs_altaz, height_start=50., height_end = 20200., height_unit = u.km ):
    """From an AltAz pointing, find the distance to where it crosses a plane 100 km above Earth"""
    try_distances = np.linspace(height_start*0.7,2*obs_altaz.secz*height_end) * height_unit
    dummytime = Time.now()
    try_altaz = far_to_near(obs_altaz, try_distances, obstime=dummytime)
    try_heights = EarthLocation(
        *(try_altaz.transform_to(ITRS(obstime=dummytime)).cartesian.xyz)
    )
    f_height_to_distance = interpolate.interp1d(
        try_heights.height.to(u.km), try_distances.to(u.km)
    )
    f_dist_to_latpp = interpolate.interp1d(
        try_distances.to(u.km), try_heights.lat.to(u.deg)
    )
    f_dist_to_lonpp = interpolate.interp1d(
        try_distances.to(u.km), try_heights.lon.to(u.deg)
    )
    f_dist_to_Rlocal = interpolate.interp1d(
        try_distances.to(u.km), try_heights.to(u.km).itrs.cartesian.norm()
    )
    return f_height_to_distance,f_dist_to_latpp,f_dist_to_lonpp, f_dist_to_Rlocal

def getProfile(source_pos, stat_pos, time, h = np.arange(60,20200,10)*u.km):
  '''Get profile of piercepoints through the atmosphere, given a station and source direction
  
  Args: 
  source_pos : SkyCoord or list(2) Ra,Dec J200
  stat_pos : EarthLocation or array(3) in m
  time : Time or float (mjd) (can be an array)
  h : array of heights

  Returns:
  Tuple[ np.array, np.array, np.array, float, float, np.array ]:
  

  '''
  # first calculate Az,El
  if not type(stat_pos) is EarthLocation:
    stat_pos = EarthLocation.from_geocentric(*list(stat_pos),unit=u.m)
  if not type(source_pos) is SkyCoord:
    source_pos = SkyCoord(source_pos[0],source_pos[1],unit=(u.hourangle,u.deg),frame = FK5)
  if not type(time) is Time:
    time = Time(time,format = 'mjd')
  aa = AltAz(location = stat_pos, obstime = time)
  azel = source_pos.transform_to(aa)
  direction = azel.transform_to(ITRS)
  londir = np.arctan2(direction.y.value,direction.x.value)
  latdir = np.arcsin(direction.z.value)
  
  fh,flat,flon,fRloc = distance_from_height(azel,h[0].to(u.km).value,h[-1].to(u.km).value, u.km)
  hdist = fh(h.to(u.km))
  latpp = flat(hdist)
  lonpp = flon(hdist)
  R_local = stat_pos.to(u.km).itrs.cartesian.norm().value
  Rpp = fRloc(hdist)
  #am = 1./((R_earthkm**2-hdist**2-(R_earthkm+h.to(u.km).value)**2)/(-2.*hdist*(R_earthkm+h.to(u.km).value)))
  am = 1./((R_local**2-hdist**2-Rpp**2)/(-2.*hdist*Rpp))
  return latpp, lonpp, latdir, londir, am
  
  

def getAzEl(pointing,time,position,ha_limit=-1000):

    if HAS_PYRAP:
        if ha_limit==-1000:
            azel=radec2azel(pointing[0],pointing[1],time=str(time)+'s',pos=position);
            az=azel['m0']['value']
            el=azel['m1']['value']
        else:
            me=measures()
            p=me.position("ITRF",str(position[0])+'m',str(position[1])+'m',str(position[2])+'m')
            t=me.epoch("UTC",qu.quantity(str(time)+'s'))
            phasedir=me.direction('J2000',str(pointing[0])+'rad',str(pointing[1])+'rad')
            me.doframe(p)
            me.doframe(t)
            hadec=me.measure(phasedir,"HADEC")
            if abs(hadec['m0']['value'])>ha_limit:
                print ("below horizon",tab.taql('calc ctod($time s)')[0],degrees(hadec['m0']['value']),degrees(hadec['m1']['value']))
                return 999,999
            else:
                azel=me.measure(phasedir,"AZELGEO")
  
                az=azel['m0']['value'];
                el=azel['m1']['value'];
    elif HAS_EPHEM:
        if ha_limit!=-1000:
            print ("limiting on HA/DEC not implemented for PyEphem yet, ignoring")
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
        #  convert to Dublin Julian Date for PyEphem
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
      return -1,-1
    return az,el


def get_time_range(start_time,end_time,timestep,time_in_sec,TIME_OFFSET=0):
    if HAS_PYRAP:
      try:
        start_time = qu.quantity(start_time).get_value()
        end_time = qu.quantity(end_time).get_value()
        print ('**** specified start and end time ', start_time, end_time)
        reference_time = start_time * 86400.0 - TIME_OFFSET
        st = reference_time - timestep 
        et = end_time * 86400.0 +  timestep + TIME_OFFSET
        #print ("getting string",reference_time)
        str_start_time =  obtain_observation_year_month_day_hms(reference_time)
        timerange= [st, et]
      except:
        print ('no time range given')
        print ('exiting')
        return -1,-1,-1
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
      return -1,-1,-1
    return timerange,str_start_time,reference_time
