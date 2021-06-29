from . import _iri
import numpy as np
from pkg_resources import resource_filename
from scipy.interpolate import interp1d



class pyiri:
    def __init__(self,year=2016,month=5,day=12,hour=1.,lon=0.,lat=0.):
        self.lon=lon
        self.lat=lat
        self.yr=year
        self.month=month
        self.day=day
        self.hour=hour
        self.flags=np.ones(50,dtype=bool)  #to do, check which flags to set
        self.jmag=0 #geographic (0) or geomagnetic (1) coordinates
        self.datapath=resource_filename(__name__,'ig_rz.dat').replace('ig_rz.dat','')
        

    def get_profile(self,hstart,hend,hstep):
        #res=_iri.iri_sub(self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        res=_iri.iri_get_sub(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        myshape=np.arange(hstart,hend+hstep,hstep).shape[0]
        return res[:myshape]

    def get_tec(self):
        #res=_iri.iri_sub(self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        res=_iri.iri_get_tec(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=50,heiend=1500)
        return res

    def get_sub(self,hstart=399,hend=400,hstep=1):
        res=_iri.iri_get_sub(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)
        #res=_iri.iri_get_tec(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=50,heiend=1500)
        return res

    def get_hprofile(self,h):
        #res=_iri.iri_sub(self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=hstart,heiend=hend,heistp=hstep)[0][0]
        heibeg=h[0]
        heiend=min(2000,h[-1])
        heistp=(heiend-heibeg)/999.
        res=_iri.iri_get_sub(self.datapath,self.flags,jmag=self.jmag,alati=self.lat,along=self.lon,iyyyy=self.yr,mmdd=self.month*100+self.day,dhour=self.hour,heibeg=heibeg,heiend=heiend,heistp=heistp)[0][0]
        heights = np.arange(heibeg,heiend+0.5*heistp,heistp)
        hinterpol = interp1d(heights,res,fill_value='extrapolate')
        
        return hinterpol(h)
      
