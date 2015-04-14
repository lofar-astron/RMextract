from pyrap.measures import measures
from math import *
import pyrap.quanta as qa;
import os
from pyrap import tables as tab
import numpy as np

R_earth=6364.62e3;
earth_ellipsoid_a = 6378137.0;
earth_ellipsoid_a2 = earth_ellipsoid_a*earth_ellipsoid_a;
earth_ellipsoid_b = 6356752.3142;
earth_ellipsoid_b2 = earth_ellipsoid_b*earth_ellipsoid_b;
earth_ellipsoid_e2 = (earth_ellipsoid_a2 - earth_ellipsoid_b2) / earth_ellipsoid_a2;
posCS002=[3826577.1095  ,461022.900196, 5064892.758]

def getMSinfo(MS=None):
    if MS is None:
        print "No measurement set given"
        return
    if os.path.isdir(MS):
        myMS=tab.table(MS)
    else:
        print "Do not understand the format of MS",MS,"bailing out"
        return;
    timerange=[np.amin(myMS.getcol('TIME_CENTROID')),np.amax(myMS.getcol('TIME_CENTROID'))]
    timestep=myMS.getcell('INTERVAL',0)
    
    pointing= tab.table(myMS.getkeyword('FIELD')).getcell('PHASE_DIR',0);    
    stations = tab.table(myMS.getkeyword('ANTENNA')).getcol('NAME')
    station_pos = tab.table(myMS.getkeyword('ANTENNA')).getcol('POSITION')

    return (timerange,timestep,pointing.flatten(),stations,station_pos)


def getPPsimple(height=[450.e3,],mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
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

def getPPsimpleAngle(height=[450.e3,],mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
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
    ppl[:,0]=np.atan2(pp[:,1],pp[:0])
    ppl[:,1]=np.atan2(pp[:,2],np.sqrt(pp[:0]*pp[:,0]+pp[:1]*pp[:,1]))
    ppl[:,2]=heigth

    return ppl,am
    

def getPP(h=450e3,mPosition=[0.,0.,0.],direction=[0.,0.,0.]):
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
    direction=me.direction("AZEL",az,el)
    return me.measure(direction,"ITRF");


def radec2azel(ra,dec,time, pos):
    me=measures();
    if type(ra)!=str:
        ra=str(ra)+'rad';
    if type(dec)!=str:
        dec=str(dec)+'rad';
    phasedir=me.direction('J2000',ra,dec)
    t=me.epoch("UTC",qa.quantity(time));
    me.do_frame(t);

    p = me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(p);

    azel = me.measure(phasedir,'azel');
    return azel;

def azel2radec(az,el,time, pos):
    me=measures();
    if type(az)!=str:
        az=str(az)+'rad';
    if type(el)!=str:
        el=str(el)+'rad';
    phasedir=me.direction('AZEL',az,el)
    t=me.epoch("UTC",qa.quantity(time));
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
    t=me.epoch("UTC",qa.quantity(time));
    me.do_frame(t);
    
    p = me.position('ITRF',str(pos1[0])+'m',str(pos1[1])+'m',str(pos1[2])+'m')
    bl = me.baseline('ITRF',str(pos2[0]-pos1[0])+'m',str(pos2[1]-pos1[1])+'m',str(pos2[2]-pos1[2])+'m')
    print bl
    me.do_frame(p);
    print me.to_uvw(bl)['xyz']
    #return uvw;
