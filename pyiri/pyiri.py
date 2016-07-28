import _iri 
import numpy as np
from pkg_resources import resource_filename


class pyiri:
    def __init__(self,year=2016,month=5,day=12,hour=1.,lon=0.,lat=0.):
        self.lon=lon
        self.lat=lat
        self.yr=year
        self.month=month
        self.day=day
        self.hour=hour
        self.flags=np.ones(50,dtype=bool)
        self.jmag=0 #geographic (0) or geomagnetic (1) coordinates
        self.datapath=resource_filename(__name__,'if_rz.dat').replace('if_rz.dat','')
        

    def get_profile(self,hstart,hend,hstep):
        #res=_iri.iri_sub(self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        res=_iri.iri_get(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        myshape=np.arange(hstart,hend+hstep,hstep).shape[0]
        return res[:myshape]

