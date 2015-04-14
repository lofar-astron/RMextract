import PosTools
import getIONEX as ionex;
from pyrap import tables as tab
import os
import numpy as np
from math import *
from pyrap.measures import measures
import pyrap.quanta as qa;


ION_HEIGHT=450.e3;

def getRM(MS=None,
           server="ftp://ftp.unibe.ch/aiub/CODE/",
           prefix='CODG',
           ionexPath="IONEXdata/",
           earth_rot=0,ha_limit=0.5*np.pi,
           **kwargs):
    '''optional arguments are:
    radec or pointing : [ra,dec] in radians,
    azel : [az,el] in radians, If azel is specified, radec will be ignored,
    timestep in s, timerange = [start, end]in MJD seconds (sorry I shall write the converter soon...), 
    stat_names = list of stings per station, 
    stat_positions = list of length 3 numpy arrays, containing station_position in ITRF meters.
    useEMM = boolean, use EMM for Earth magnetic field, you should specify the path to the coefficients via
    EMMpath = path. If the path is given it will be automatically assumed that you want to use EMM for Earthm agnetic field 
    useWMM = boolean, use WMM for Earth magnetic Field (see above). You can specify the path to coefficients with
    WMMpath= path.
    Returns the (timegrid,timestep,TEC) where TEC is a dictionary containing 1 enumpyarray per station in stat_names. 
    If stat_names is not given, the station names will either be extracted from the MS or st1...stN '''

    me=measures();
    stat_names=[]
    useWMM=False
    useEMM=False
    EMMpath='./'
    timestep=10
    pointing=[0.,np.pi]
    if not (MS is None):

        (timerange,timestep,pointing,stat_names,stat_pos)=PosTools.getMSinfo(MS);
    useAzel=False
    
    for key in kwargs.keys():
        if not useAzel and (key=='radec' or key=='pointing'):
            pointing = kwargs[key]
        if key=='azel':
            pointing = kwargs[key]
            useAzel=True
        if key.lower()=='timestep':
            timestep=kwargs[key]
        if key.lower()=='timerange':
            timerange=kwargs[key][:]
        if key.lower()=='stat_names':
            stat_names=kwargs[key][:]
        if key.lower()=='stat_positions':
            stat_pos=kwargs[key][:]
            if not stat_names:
                stat_names =['st%d'%(i+1) for i in range(len(stat_pos))]
        if key.lower()=='emmpath':
            EMMpath=kwargs[key]
            if not useWMM:
                useEMM=True
        if key.lower()=='wmmpath':
            EMMpath=kwargs[key]
            if not useEMM:
                useWMM=True
        if key.lower()=='usewmm':
            useWMM=kwargs[key]
        if key.lower()=='useemm':
            useEMM=kwargs[key]
    
    if useWMM  or useEMM :
        import EMM.EMM as EMM
        #import EMM
        if useEMM:
            print "USING EMM for EarthMagnetic Field, remember to have set the path to the coeffiecients (EMMpath) correctly."
            emm=EMM.EMM(cof=EMMpath)
        else:
            print "USING WMM for EarthMagnetic Field, remember to have set the path to the coeffiecients (WMMpath) correctly."
            emm=EMM.WMM(cof=EMMpath)
    
    else:
        print "USING obsolete casacore libraries for getting Earthmagnetic field!"
        import getEarthMagnetic as EM


     
    #IONEX files go per day, check if more than one file is  needed.
    times=[];
    while tab.taql('calc cdate($timerange[0] s)')[0] != tab.taql('calc cdate($timerange[1] s)')[0]:
        
        part_of_day = tab.taql('calc time($timerange[0] s)')[0]*24.*60.*60./(2*np.pi);
        nr_remaining_seconds=(24.*60.*60.-part_of_day);
        times.append(np.arange(timerange[0],timerange[0]+nr_remaining_seconds,timestep));
        timerange[0]+=len(times[-1])*timestep;
    times.append(np.arange(timerange[0],timerange[1]+timestep,timestep)) #add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
    
    timegrid=np.array([])
    TECs={};
    Bs={};
    Bpars={};
    ams={};
    flags={}
    RMs={};
    for st in stat_names:
        Bs[st]=[]
        Bpars[st]=[]
        RMs[st]=[]
        TECs[st]=[]
        ams[st]=[]
        flags[st]=[]
    for time_array in times:
        #get RM per timeslot
        starttime=time_array[0];
        starttime=tab.taql('calc ctod($starttime s)')[0]
        #get relevant ionex file
        ionexf=ionex.getIONEXfile(time=starttime,server=server,prefix=prefix,outpath=ionexPath);
        tecinfo=ionex.readTEC(ionexf)
        for station,position in  zip(stat_names,stat_pos):
          print "getting RM for",station,"at",position,"time",time_array.shape
          for time in time_array:
            part_of_day=tab.taql('calc time($time s)')[0]*24./(2*np.pi)
            if useAzel:
                  az=pointing[0]
                  #el=0.5*np.pi-pointing[0]
                  el=pointing[1]
                  flags[station].append(1)
            else:
                  p=me.position("ITRF",str(position[0])+'m',str(position[1])+'m',str(position[2])+'m')
                  t=me.epoch("UTC",qa.quantity(str(time)+'s'))
                  phasedir=me.direction('J2000',str(pointing[0])+'rad',str(pointing[1])+'rad')
                  me.doframe(p)
                  me.doframe(t)
                  hadec=me.measure(phasedir,"HADEC")
                  if abs(hadec['m0']['value'])>ha_limit:
                      print "below horizon",tab.taql('calc ctod($time s)')[0],degrees(hadec['m0']['value']),degrees(hadec['m1']['value'])

                      flags[station].append(0)
                  else:
                      flags[station].append(1)
                  azel=me.measure(phasedir,"AZEL")
                  #                azel=PosTools.radec2azel(pointing[0],pointing[1],time=str(time)+'s',pos=position);
                  az=azel['m0']['value'];
                  el=azel['m1']['value'];
            lonlat=PosTools.getLonLatStation(az,el,pos=position);

            lon=lonlat['m0']['value'];
            lat=lonlat['m1']['value'];
            # convert to itrf coordinates on sphere with radius 1
            diritrf=[cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat)]
            # calculate piercepoint in xyz(code from Bas vd Tol)
            
            (ppx1,ppy1,ppz1,am1)=PosTools.getPP(h=ION_HEIGHT,mPosition=position,direction=diritrf)

            #get pp in lon,lat h
            pp1position=me.position("ITRF",str(ppx1)+'m',str(ppy1)+'m',str(ppz1)+'m')
            
            lonpp = degrees(pp1position['m0']['value']);
            latpp = degrees(pp1position['m1']['value']);
            #get VTEC from IONEX interpolation
            vTEC=ionex.getTECinterpol(time=part_of_day,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)

            TECs[station].append(vTEC*am1)  #STEC value
            if useEMM or useWMM:
                            # get EMM BField and project along LOS  
                pp1wgsposition=me.measure(pp1position,"WGS84")
           
            
                lonpp = degrees(pp1wgsposition['m0']['value']);
                latpp = degrees(pp1wgsposition['m1']['value']);
                emm.lon=lonpp
                emm.lat=latpp
                emm.h=pp1wgsposition['m2']['value']/1.e3
                BField=emm.getXYZ()
                Bpar=-1*emm.getProjectedField(lon,lat)# minus sign since the radiation is towards the Earth
            else:
                # get IGRF BField and project along LOS
                BField=EM.getField([ppx1,ppy1,ppz1],time=str(time)+'s');
                Bpar=-1*EM.ProjectField(BField,lon,lat); # minus sign since the radiation is towards the Earth
            Bs[station].append(BField)
            Bpars[station].append(Bpar)
            ams[station].append(am1)
            #calculate RM constant comes from VtEC in TECU,B in nT, RM = 2.62
            RMs[station].append((Bpar*vTEC*am1)*2.62e-6)

            
        timegrid=np.concatenate((timegrid,time_array));
    
    for st in stat_names:
            TECs[st]=np.array(TECs[st])
            Bs[st]=np.array(Bs[st])
            Bpars[st]=np.array(Bpars[st])
            ams[st]=np.array(ams[st])
            RMs[st]=np.array(RMs[st])
            flags[st]=np.array(flags[st])

        
    big_dict={}
    big_dict['VTEC']=TECs
    big_dict['Bpar']=Bpars
    big_dict['BField']=Bs
    big_dict['AirMass']=ams
    big_dict['RM']=RMs
    big_dict['times']=timegrid
    big_dict['timestep']=timestep
    big_dict['flags']=flags
    
    return big_dict


