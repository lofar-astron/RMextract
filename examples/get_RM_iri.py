import RMextract.EMM.EMM as EMM
import RMextract.getIONEX as ionex
import RMextract.PosTools as pt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from astropy.coordinates import Angle
from astropy.coordinates import Latitude
from astropy.time import Time
import astropy.units as u  
import RMextract.pyiriplas.pyiriplas as pyiriplas
from lofarantpos import db
from lofarantpos import geo

mydb=db.LofarAntennaDatabase()
nenufarcentre_geo = (np.deg2rad(2.192400), np.deg2rad(47.376511), 150)
nenufarcentre_xyz = geo.xyz_from_geographic(*nenufarcentre_geo)
nenufarcentre_xyz = [4324017.054,165545.160 ,4670271.072]
#ionex_server="ftp://ftp.aiub.unibe.ch/CODE/"
#
ionex_server = "ftp://cddis.nasa.gov/gnss/products/ionex/"
ionex_server = "ftp://gssc.esa.int/gnss/products/ionex/"
#ionex_prefix='UQRG'
ionex_prefix='CODG'
ionexPath="/home/mevius/IONEX_DATA/"
earth_rot = 0
light_speed=299792458.



def init_Meerkat(src = "00:34:08.8703 -07:21:53.409", times =[58774.9430718395,
             59034.0376808637,
             59064.9741027126,
             59090.8953866498] ): 
    XYZ_Meerkat = [5109360.133,2006852.586,-3238948.127]
    #statpos = EarthLocation.from_geodetic(Angle((21,19,48.00),u.deg),Latitude((-30,49,48.0),u.deg))
    statpos = EarthLocation.from_geocentric(*XYZ_Meerkat,u.m)
    statposarray = np.array([[statpos.x.value,statpos.y.value,statpos.z.value]]*len(times) )
    srcpos = SkyCoord(src, unit=(u.hourangle,u.deg),frame="fk5")
    atimes =  [Time(i,format='mjd') for i in times]

    altazdir = [ srcpos.transform_to(AltAz(obstime = itime, location = statpos)) for itime in atimes ] 
    itrfdir = [i.transform_to('itrs') for i  in altazdir]
    direction = np.array([ [i.x,i.y,i.z] for i in itrfdir]).T
    return statposarray,srcpos,direction,atimes
    


def init_from_file(fname,src ="05:28:52.264 +22:00:04"):
    data = np.loadtxt(fname,usecols=[0,1,2,3,4,5,6])
    stations = np.loadtxt(fname,usecols=[7],dtype=str)
    statpos = [ "CS002LBA" if i=="LOFAR" else i+"HBA" for i in stations]
    statpos=[EarthLocation.from_geocentric(*mydb.phase_centres[i],unit=u.m) for i in statpos]
    statposarray = np.array([[i.x.value,i.y.value,i.z.value]  for i in statpos])
    srcpos = SkyCoord(src, unit=(u.hourangle,u.deg),frame="fk5")
    atimes = Time(data[:,0],format='mjd') 
    altazdir = [ srcpos.transform_to(AltAz(obstime = itime, location = pos)) for pos,itime in zip(statpos,atimes) ]
    itrfdir = [i.transform_to('itrs') for i in altazdir]
    direction = np.array([ [i.x,i.y,i.z] for i in itrfdir]).T
    return statposarray,srcpos,direction,atimes

def init_from_time_file(fname,src = "10:22:57.9992 +10:01:52.78"):
    tms = np.loadtxt(fname,usecols=[1])
    stations = [fname.split("/")[-1][:5],]
    statpos = [ "CS002LBA" if i=="LOFAR" else i+"HBA" for i in stations]*tms.shape[0]
    statpos=[EarthLocation.from_geocentric(*mydb.phase_centres[i],unit=u.m) for i in statpos]
    statposarray = np.array([[i.x.value,i.y.value,i.z.value]  for i in statpos])
    srcpos = SkyCoord(src, unit=(u.hourangle,u.deg),frame="fk5")
    atimes = Time(tms,format='mjd') 
    altazdir = [ srcpos.transform_to(AltAz(obstime = itime, location = pos)) for pos,itime in zip(statpos,atimes) ]
    itrfdir = [i.transform_to('itrs') for i in altazdir]
    direction = np.array([ [i.x,i.y,i.z] for i in itrfdir]).T
    return statposarray,srcpos,direction,atimes
    
def init_from_nat_file(fname,src = "10:22:57.9992 +10:01:52.78",skiprows = 0,stations = None):
    tms = np.loadtxt(fname,usecols=[0],skiprows = skiprows)
    if stations is None:
        stations = [fname.split("/")[-2][:5],]
    statpos = [ "CS002LBA" if i=="LOFAR" else i+"HBA" for i in stations]*tms.shape[0]
    statpos=[EarthLocation.from_geocentric(*mydb.phase_centres[i],unit=u.m) for i in statpos]
    statposarray = np.array([[i.x.value,i.y.value,i.z.value]  for i in statpos])
    srcpos = SkyCoord(src, unit=(u.hourangle,u.deg),frame="fk5")
    atimes = Time(tms,format='mjd') 
    altazdir = [ srcpos.transform_to(AltAz(obstime = itime, location = pos)) for pos,itime in zip(statpos,atimes) ]
    itrfdir = [i.transform_to('itrs') for i in altazdir]
    direction = np.array([ [i.x,i.y,i.z] for i in itrfdir]).T
    return statposarray,srcpos,direction,atimes

def init_from_louis_file(fname,srqc = "10:22:57.9992 +10:01:52.78"):
    tms = np.loadtxt(fname)[0]
    stations = [nenufarcentre_xyz ]*tms.shape[0]
    statpos=[EarthLocation.from_geocentric(*i,unit=u.m) for i in stations]
    statposarray = np.array([[i.x.value,i.y.value,i.z.value]  for i in statpos])
    srcpos = SkyCoord(src, unit=(u.hourangle,u.deg),frame="fk5")
    atimes = Time(tms,format='mjd') 
    altazdir = [ srcpos.transform_to(AltAz(obstime = itime, location = pos)) for pos,itime in zip(statpos,atimes) ]
    itrfdir = [i.transform_to('itrs') for i in altazdir]
    direction = np.array([ [i.x,i.y,i.z] for i in itrfdir]).T
    return statposarray,srcpos,direction,atimes

hstart=60
hmiddle=800
hend=20200
hstep1=10
hstep2=100
hstep1=hstep2=60
#heights = 450.e3*np.ones(mydir.shape[1:])
h =  np.arange(hstart*1e3,hmiddle*1e3+0.5*hstep1*1e3,hstep1*1e3)
idx_hbot = h.shape[0]
h = np.concatenate((h, np.arange(h[-1]+hstep2*1e3,hend*1e3+0.5*hstep2*1e3,hstep2*1e3)))
hweights = np.concatenate((hstep1*np.ones(idx_hbot), hstep2*np.ones(h.shape[0]-idx_hbot)))


def getRM_iri(times, statposarray,directions,heights = h,ionex_server = ionex_server ,ionex_prefix = "uqrg",ionexPath = "/home/mevius/IONEX_DATA/",return_short=True,alt_ionex_prefix='uprg'): 
    h = heights
    hidx = np.argmin(np.abs(h-250e3))
    londir = np.arctan2(directions[1],directions[0])
    latdir = np.arcsin(directions[2])
    pps = []
    vtecs = []
    Bpars = []
    ams = []
    profiles = []
    plasprofiles = []
    topprofiles = []
    
    for itime,(time,statpos) in enumerate(zip(times,statposarray)):
        '''get (lon,lat,h values) of piercepoints for antenna position mPosition in m, direction ITRF in m on unit sphere and for array of heights, assuming a spherical Earth'''
        print (itime,time)
        pp,am = pt.getPPsimpleAngle(height=h,mPosition=statpos,direction=directions[:,itime])
        try:
            print("trying  here",ionex_prefix)
            ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix=ionex_prefix,outpath=ionexPath)
            tecinfo=ionex.readTEC(ionexf)
        except:
            try:
                print("failed trying",alt_ionex_prefix)
                ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix=alt_ionex_prefix,outpath=ionexPath)
                tecinfo=ionex.readTEC(ionexf)

            except:
                try:
                    print("failed trying igsg")
                    ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix="igsg",outpath=ionexPath)
                    tecinfo=ionex.readTEC(ionexf)
                except:
                    print("last trying codg")
                    ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix="codg",outpath=ionexPath)
                    tecinfo=ionex.readTEC(ionexf)

        print(ionexf)
        lonpp = np.degrees(pp[:,0])
        latpp = np.degrees(pp[:,1])
        vTEC=ionex.getTECinterpol(time=np.resize((time.mjd - int(time.mjd))*24,latpp.shape),lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)
        vtecs.append(vTEC)

        emm=EMM.WMM(date=time.decimalyear,lon=lonpp[hidx],lat=latpp[hidx],h=h[hidx]*1e-3)

        BField = emm.getProjectedFieldArray(lonpp,latpp,h*1e-3,los_dir=[londir[itime],latdir[itime]])
        Bpars.append(BField)
    #    myiri=pyiri.pyiri(year=time.ymdhms.year,month=time.ymdhms.month,day=time.ymdhms.day,hour=time.ymdhms.hour + time.ymdhms.minute/60.)

    #    myiri.lon=lonpp[hidx]
    #    myiri.lat=latpp[hidx]
    #    profiles.append(myiri.get_hprofile(h*1e-3))
        myiri=pyiriplas.pyiriplas(year=time.ymdhms.year,month=time.ymdhms.month,day=time.ymdhms.day,hour=time.ymdhms.hour + time.ymdhms.minute/60.)

        myiri.lon=lonpp[hidx]
        myiri.lat=latpp[hidx]
        #plasprofiles.append(myiri.get_profile(h*1e-3))
        plasprofiles.append(np.concatenate((myiri.get_profile(h[:idx_hbot]*1e-3),myiri.get_profile(h[idx_hbot:]*1e-3))))
        print(idx_hbot,plasprofiles[0].shape)
        topprofiles.append([])
        for ih  in range(0,lonpp.shape[0]):
            myiri.lon=lonpp[ih]
            myiri.lat=latpp[ih]
            topprofiles[-1].append(np.concatenate((myiri.get_profile(h[:idx_hbot]*1e-3),myiri.get_profile(h[idx_hbot:]*1e-3))))

        ams.append(am)
        pps.append(pp)
    vtecs = np.array(vtecs)
    Bpars = np.array(Bpars) 
    ams = np.array(ams)
    pps = np.array(pps)
    plasprofiles = np.array(plasprofiles)
    print (plasprofiles.shape,hweights.shape)
    plasprofiles *= hweights
    nplas = plasprofiles/np.sum(plasprofiles,axis=1)[:,np.newaxis]
    topprofiles = np.array(topprofiles) *np.diag(hweights)
    ntop = topprofiles/np.sum(np.diagonal(topprofiles,0,1,2),axis=1)[:,np.newaxis,np.newaxis]
    finalp = []
    for k in range(ntop.shape[0]):
        finalp.append(np.diag(ntop[k]))
    finalp = np.array(finalp)


    RM = np.sum(-1*vtecs*ams*Bpars[:,:,0]*nplas*2.62e-6,axis=1)
    RM2 = np.sum(-1*vtecs*ams*Bpars[:,:,0]*finalp*2.62e-6,axis=1)
    RM3 = -1*vtecs[:,hidx]*ams[:,hidx]*Bpars[:,hidx,0]*2.62e-6

    if return_short:
        return RM2
    else:
        return vtecs,Bpars,ams,finalp,pps,nplas,RM2


def getRM_iri_new(times, statposarray,srcpos,heights = h,ionex_server = ionex_server ,ionex_prefix = "uqrg",ionexPath = "/home/mevius/IONEX_DATA/",return_short=True,alt_ionex_prefix='uprg'): 
    h = heights
    hidx = np.argmin(np.abs(h-250e3))
    hkm = heights*u.m
    pps = []
    vtecs = []
    Bpars = []
    ams = []
    profiles = []
    plasprofiles = []
    topprofiles = []
    
    for itime,(time,statpos) in enumerate(zip(times,statposarray)):
        '''get (lon,lat,h values) of piercepoints for antenna position mPosition in m, sourcepos radec or SkyCoord for array of heights'''
        print (itime,time)
        #pp,am = pt.getPPsimpleAngle(height=h,mPosition=statpos,direction=directions[:,itime])
        latpp, lonpp, latdir, londir, am = pt.getProfile(srcpos,statpos,time, hkm)
        try:
            print("trying  here",ionex_prefix)
            ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix=ionex_prefix,outpath=ionexPath)
            tecinfo=ionex.readTEC(ionexf)
        except:
            try:
                print("failed trying",alt_ionex_prefix)
                ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix=alt_ionex_prefix,outpath=ionexPath)
                tecinfo=ionex.readTEC(ionexf)

            except:
                try:
                    print("failed trying igsg")
                    ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix="igsg",outpath=ionexPath)
                    tecinfo=ionex.readTEC(ionexf)
                except:
                    print("last trying codg")
                    ionexf=ionex.getIONEXfile(time=time.ymdhms,server=ionex_server,prefix="codg",outpath=ionexPath)
                    tecinfo=ionex.readTEC(ionexf)

        print(ionexf)
        vTEC=ionex.getTECinterpol(time=np.resize((time.mjd - int(time.mjd))*24,latpp.shape),lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)
        vtecs.append(vTEC)

        emm=EMM.WMM(date=time.decimalyear,lon=lonpp[hidx],lat=latpp[hidx],h=h[hidx]*1e-3)

        BField = emm.getProjectedFieldArray(lonpp,latpp,h*1e-3,los_dir=[londir,latdir])
        Bpars.append(BField)
    #    myiri=pyiri.pyiri(year=time.ymdhms.year,month=time.ymdhms.month,day=time.ymdhms.day,hour=time.ymdhms.hour + time.ymdhms.minute/60.)

    #    myiri.lon=lonpp[hidx]
    #    myiri.lat=latpp[hidx]
    #    profiles.append(myiri.get_hprofile(h*1e-3))
        myiri=pyiriplas.pyiriplas(year=time.ymdhms.year,month=time.ymdhms.month,day=time.ymdhms.day,hour=time.ymdhms.hour + time.ymdhms.minute/60.)

        myiri.lon=lonpp[hidx]
        myiri.lat=latpp[hidx]
        #plasprofiles.append(myiri.get_profile(h*1e-3))
        plasprofiles.append(np.concatenate((myiri.get_profile(h[:idx_hbot]*1e-3),myiri.get_profile(h[idx_hbot:]*1e-3))))
        print(idx_hbot,plasprofiles[0].shape)
        topprofiles.append([])
        for ih  in range(0,lonpp.shape[0]):
            myiri.lon=lonpp[ih]
            myiri.lat=latpp[ih]
            topprofiles[-1].append(np.concatenate((myiri.get_profile(h[:idx_hbot]*1e-3),myiri.get_profile(h[idx_hbot:]*1e-3))))

        ams.append(am)
        pps.append([lonpp,latpp,h])
    vtecs = np.array(vtecs)
    Bpars = np.array(Bpars) 
    ams = np.array(ams)
    pps = np.array(pps)
    plasprofiles = np.array(plasprofiles)
    print (plasprofiles.shape,hweights.shape)
    plasprofiles *= hweights
    nplas = plasprofiles/np.sum(plasprofiles,axis=1)[:,np.newaxis]
    topprofiles = np.array(topprofiles) *np.diag(hweights)
    ntop = topprofiles/np.sum(np.diagonal(topprofiles,0,1,2),axis=1)[:,np.newaxis,np.newaxis]
    finalp = []
    for k in range(ntop.shape[0]):
        finalp.append(np.diag(ntop[k]))
    finalp = np.array(finalp)


    RM = np.sum(-1*vtecs*ams*Bpars[:,:,0]*nplas*2.62e-6,axis=1)
    RM2 = np.sum(-1*vtecs*ams*Bpars[:,:,0]*finalp*2.62e-6,axis=1)
    RM3 = -1*vtecs[:,hidx]*ams[:,hidx]*Bpars[:,hidx,0]*2.62e-6

    if return_short:
        return RM2
    else:
        return vtecs,Bpars,ams,finalp,pps,nplas,RM2
