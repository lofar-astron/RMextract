from __future__ import absolute_import
import RMextract.EMM.EMM_Model as emm
from math import *
import numpy as np
from pkg_resources import resource_filename
import os

class WMM:
    def __init__(self,cof=resource_filename(__name__,'WMM.COF'),date=2010.,lon=0.,lat=0.,h=0.):
        try:
            cof=cof.replace("WMM.COF","")
            if not os.path.isfile(cof+"WMM.COF"):
                print ("initializing model failed. Did you specify correct location for coefficient files")
                return
            self.WMM_Model=emm.WMM_Model(cof+"WMM",date,lon,lat,h);
        except (RuntimeError, TypeError, NameError):
            print ("initializing mdoel failed. Did you specify correct location for coefficient files")
        self.date=date;
        self.lon=lon;
        self.lat=lat;
        self.h=h;
        self.changed=False

    def getXYZ(self):
        '''get BField, convert to (ITRF) X,Y,Z coordinates'''
        N,E,D=self.getNED()
        lat=radians(self.lat)
        lon=radians(self.lon)
        X=-1*cos(lat)*cos(lon)*D-sin(lat)*cos(lon)*N-sin(lon)*E;
        Y=-1*cos(lat)*sin(lon)*D-sin(lat)*sin(lon)*N+cos(lon)*E;
        Z=-1*sin(lat)*D+N*cos(lat)
        return X,Y,Z

    def getNED(self):
        #check if something has changed
        if self.WMM_Model.getDate() != self.date:
            self.changed=True;
            self.WMM_Model.setDate(self.date);
        if self.WMM_Model.getLon() != self.lon or self.WMM_Model.getLat() != self.lat :
            self.changed=True;
            self.WMM_Model.setLonLat(self.lon,self.lat);
        if self.WMM_Model.getHeight() != self.h:
            self.changed=True;
            self.WMM_Model.setHeight(self.h);
        if self.changed:
            self.WMM_Model.setEM()
            self.changed=False
            
        return self.WMM_Model.getX(),self.WMM_Model.getY(),self.WMM_Model.getZ()
        #North, East,Down
    
    def getProjectedField(self,lon,lat):
        '''get fieldstrength  in direction lon,lat'''
        x,y,z=self.getXYZ();
        return sin(lat)*z+cos(lat)*cos(lon)*x +cos(lat)*sin(lon)*y;

    def getProjectedFieldVector(self,lon,lat):
        '''get field  with one axis projected in direction lon,lat'''
        x,y,z=self.getXYZ();
        return [sin(lat)*z+cos(lat)*cos(lon)*x +cos(lat)*sin(lon)*y,
                -1*sin(lon)*x +cos(lon)*y,
                cos(lat)*z-sin(lat)*cos(lon)*x -sin(lat)*sin(lon)*y]

    def getProjectedFieldArray(self,lon_array,lat_array,h_array,los_dir):
        '''returns numpy array with projected BField along LOS (given in lon,lat) for ray tracing'''
        result=np.zeros((len(lon_array),3))
        
        for idx,self.lon,self.lat,self.h in zip(range(len(lon_array)),lon_array,lat_array,h_array):
            result[idx]=self.getProjectedFieldVector(los_dir[0],los_dir[1])
        return result


class EMM(WMM):
    def __init__(self,cof=resource_filename(__name__,'EMM2000.COF'),date=2010.,lon=0.,lat=0.,h=0.):
        try:
            cof=cof.replace("EMM2000.COF","")
            if not os.path.isfile(cof+"EMM2000.COF"):
                print ("initializing model",cof+"EMM2000.COF","failed. Did you specify correct location for coefficient files")
                return
            
            self.WMM_Model=emm.EMM_Model(cof+"EMM",date,lon,lat,h);
        except (RuntimeError, TypeError, NameError):
            print ("initializing mdoel failed. Did you specify correct location for coefficient files")
        self.date=date;
        self.lon=lon;
        self.lat=lat;
        self.h=h;
        self.changed=False
