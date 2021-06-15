from RMextract.EMM import EMM as EMM
import RMextract.PosTools as PosTools
import numpy as np
import RMextract.getIONEX as ionex
import pyiri.pyiri as pyiri
from datetime import date

def getPParray(pointing,time,position,height_array):
    az,el = PosTools.getAzEl(pointing,time,position)
    lonlat=PosTools.getLonLatStation(az,el,pos=position)
    los_dir=[lonlat['m0']['value'],lonlat['m1']['value']]
    lat=los_dir[1]
    lon=los_dir[0]
    itrfdir=[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
    pp,am=PosTools.getPPsimpleAngle(height=height_array,mPosition=position,direction=itrfdir)
    return pp,am

def getBarray(pointing,time,position,height_array):
    pp,am=getPParray(pointing,time,position,height_array)
    year, month, day, myfrac =  PosTools.obtain_observation_year_month_day_fraction(time)
    dayofyear = date(year,month,day).timetuple().tm_yday
    emm= EMM.WMM(date=year + float(dayofyear) / 365.)
    BField=emm.getProjectedFieldArray(np.degrees(pp[:,0]),np.degrees(pp[:,1]),pp[:,2]/1e3,los_dir)
    return pp,am,BField


def getTECarray(pp,time,server="ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/",prefix="igsg",earth_rot=0,ionexPath="/home/users/mevius/IONEXdata/"):
    date_parms =  PosTools.obtain_observation_year_month_day_fraction(time)
    part_of_day= date_parms[3] * 24
        #get relevant ionex file
    ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=ionexPath)
    tecinfo=ionex.readTEC(ionexf)
    vtecs=[]
    for latpp,lonpp in zip(np.degrees(pp[:,1]),np.degrees(pp[:,0])):
        vTEC=ionex.getTECinterpol(time=part_of_day,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)
        vtecs.append(vTEC)
    return np.array(vtecs)


def getRM_raytrace(time,pointing,position):
    heights=np.arange(50e3,1000e3,10e3)
    pp,am,BField=getBarray(pointing,time,position,heights)
    year, month, day, myfrac =  PosTools.obtain_observation_year_month_day_fraction(time)

    vtecs=getTECarray(pp,time)
    myiri=pyiri.pyiri(year=year,month=month,day=day,hour=myfrac)
    dprofile=[]
    for lon,lat,h in zip(np.degrees(pp[:,0]),np.degrees(pp[:,1]),pp[:,2]/1.e3):
        myiri.lon=lon
        myiri.lat=lat
        dprofile.append(myiri.get_profile(hstart=h,hend=h,hstep=1))
    #dprofile=myiri.iri_sub(flags,jmag=0,alati=52.,along=6.,iyyyy=2000,mmdd=101,dhour=1.5,heibeg=100.,heiend=1000.,heistp=10.)
    return pp,am,BField,vtecs,dprofile

def getTEC_iri(time,pointing,position):
    heights=np.array([300e3])
    pp,am=getPParray(pointing,time,position,heights)
    year, month, day, myfrac =  PosTools.obtain_observation_year_month_day_fraction(time)
    print ("getting data for",year, month, day, myfrac)
    vtecs=getTECarray(pp,time)
    myiri=pyiri.pyiri(year=year,month=month,day=day,hour=myfrac)
    iri_tec=[]
    for lon,lat,h in zip(np.degrees(pp[:,0]),np.degrees(pp[:,1]),pp[:,2]/1.e3):
        myiri.lon=lon
        myiri.lat=lat
        iri_tec.append(myiri.get_tec()[0])
    #dprofile=myiri.iri_sub(flags,jmag=0,alati=52.,along=6.,iyyyy=2000,mmdd=101,dhour=1.5,heibeg=100.,heiend=1000.,heistp=10.)
    return pp,am,vtecs,np.array(iri_tec)
