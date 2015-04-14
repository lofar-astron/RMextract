from pyrap.measures import measures
from math import *


def getField(pos=[3826577.1095  ,   461022.900196, 5064892.758],time='today'):
    me=measures()
    position=pos;
    if isinstance(position,list):
        position=me.position('ITRF',str(pos[0])+'m',str(pos[1])+'m',str(pos[2])+'m')
    me.do_frame(position);
    me.do_frame(me.epoch('utc', time))
    magn=me.earthmagnetic()
    field=me.measure(magn,"ITRF");
    field['m0']['value']*=-1
    field['m1']['value']*=-1
    field['m2']['value']*=-1
    # sign of field is flipped in casacore
    return [field['m0']['value'],field['m1']['value'],field['m2']['value']]


def ProjectField(field,lon=0,lat=0):
    Tx=field[0]
    Ty=field[1]
    Tz=field[2]
    return sin(lat)*Tz+cos(lat)*cos(lon)*Tx +cos(lat)*sin(lon)*Ty;
 



