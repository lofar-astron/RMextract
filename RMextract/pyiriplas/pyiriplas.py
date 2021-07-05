from . import _iriplas
import numpy as np
from pkg_resources import resource_filename
from scipy.interpolate import interp1d


class pyiriplas:
    def __init__(self,year=2016,month=5,day=12,hour=1.,lon=0.,lat=0.):
        self.lon=lon
        self.lat=lat
        self.yr=year
        self.month=month
        self.day=day
        self.hour=hour
        self.jmag=0 #geographic (0) or geomagnetic (1) coordinates
        self.datapath=resource_filename(__name__,'ig_rz.dat').replace('ig_rz.dat','')
        

    def get_profile(self,heights):
        '''get interpolated density profile at heights (between 80 and 20200km)'''
        #res=_iri.iri_sub(self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        res=_iriplas.iri_plas_main(self.datapath,alati=self.lat,along=self.lon,jmag=self.jmag,iyyyy=self.yr,mmdd=self.month*100+self.day,hours=self.hour)
        h = res[0]
        hidx = h>0
        NE = res[1][hidx]
        h=h[hidx]
        hinterpol = interp1d(h,NE,fill_value='extrapolate')
        
        return hinterpol(heights)

