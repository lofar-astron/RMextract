import RMextract.getRM as gt
import lofar.parmdb as parmdb
import pyrap.tables as tab
import numpy as np


def createRMParmdb(MS,parmdbname,create=True,patchnames=['phase_center']):
    myParmdb=parmdb.parmdb(parmdbname,create=create)
    if create:
        myParmdb.addDefValues("Gain:0:0:Ampl",1.e-4)
        myParmdb.addDefValues("DirectionalGain:0:0:Ampl",1.)
        myParmdb.addDefValues("DirectionalGain:1:1:Ampl",1.)
        myParmdb.addDefValues("Gain:0:0:Real",1.e-4)
        myParmdb.addDefValues("Gain:1:1:Real",1.e-4)
        myParmdb.addDefValues("DirectionalGain:0:0:Real",1.)
        myParmdb.addDefValues("DirectionalGain:1:1:Real",1.)
    myMS=tab.table(MS)
    stations=tab.table(myMS.getkeyword('ANTENNA')).getcol('NAME')
    stat_pos=tab.table(myMS.getkeyword('ANTENNA')).getcol('POSITION')
    result=gt.getRM(MS=MS,earth_rot=1,ionexPath='IONEXdata/',useWMM=True,WMMpath='EMM/',timestep=900.,stat_names=stations[2:3],stat_positions=stat_pos[2:3])
    RM=result['RM']    
    for st in stations:
        print "storing station",st,(st in RM.keys())
        if not (st in RM.keys()):
            stname=RM.keys()[0]
        RM[stname]=RM[stname].reshape(RM[stname].shape[:1]+(1,))
        myValue=myParmdb.makeValue(values=RM[stname], sfreq=1e10, efreq=2e10, stime=result['times'], etime=np.ones(result['times'].shape,dtype=float)*result['timestep'], asStartEnd=False)
        for patchname in patchnames:
            valuename = "RotationMeasure:%s:%s"%(st,patchname)
            myParmdb.deleteValues(valuename)            
            myParmdb.addValues(valuename,myValue)
