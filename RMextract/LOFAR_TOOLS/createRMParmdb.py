import RMextract.getRM as gt
import lofar.parmdb as parmdb
import pyrap.tables as tab
import numpy as np
import argparse


def createRMParmdb(MS,parmdbname,create=True,patchname='',
                   server="ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/",
                   prefix='codg',
                   ionexPath="IONEXdata/",
                   earth_rot=0,
                   timestep=900.,
                   stat_names=None):
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
    if not stat_names==None:
        if not(stat_names=='all'):
            stat_pos=[stat_pos[stations.index(name)] for name in stat_names] 
        else:
            stat_names=stations[:]
    else:
        stat_names=stations[2:3]
        stat_pos=stat_pos[2:3]
    result=gt.getRM(MS=MS,server=server,ionexPath=ionexPath,prefix=prefix,earth_rot=earth_rot,timestep=timestep,stat_names=stat_names,stat_positions=stat_pos)
    RM=result['RM']    
    for st in stations:
        #print "storing station",st,(st in RM.keys())
        if not (st in RM.keys()):
            stname=RM.keys()[0]
        else:
            stname=st
        RM[stname]=RM[stname].reshape(RM[stname].shape[:1]+(1,))
        myValue=myParmdb.makeValue(values=RM[stname], sfreq=1e10, efreq=2e10, stime=result['times'], etime=np.ones(result['times'].shape,dtype=float)*result['timestep'], asStartEnd=False)
        if patchname:
            valuename = "RotationMeasure:%s:%s"%(st,patchname)
        else:
            valuename = "RotationMeasure:%s"%(st)
        myParmdb.deleteValues(valuename)            
        myParmdb.addValues(valuename,myValue)

def main():
    descriptiontext = "Create a parmDB with the Ionospheric RM predicted from data in IONEX files.\n" + \
                      "Default is to create a parmDB that is compatible with the \"rotationmeasure\" " + \
                      "correction in NDPPP ApplyCal.\n"
    parser = argparse.ArgumentParser(description=descriptiontext)
   
    parser.add_argument('MS',help='Measurement-Set for which the parmDB is to be created.')
    parser.add_argument('-o','--out', help='name of the parmdb',dest="parmdbName",default='RMParmdb');
    parser.add_argument('-p','--patch',dest="patch",default="",type=str,
                        help='If given: create parameter-names that include a patch-name as needed for BBS. '
                        'The default is to create parameter-names without a patch-name as they are needed for '
                        'the "rotationmeasure" correction in NDPPP ApplyCal.')
    #parser.add_argument('-s','--sky', help='name of the skymodel, if no source/patch name is given first will be selected. If you do not set the skymodel, the phasecenter of the MS will be assumed as your direction. If option --use_phase_center is used the direction in the sky model are ignored. However, it could still be useful to automatically select a patch name, which is needed for direction dependent parmareters.',dest="sky")
    #parser.add_argument('-f','--use_phase_center', help='force using of phase center of MS, even if skymodel is supplied',action='store_true',dest="usePhaseCenter");
    parser.add_argument('--IONprefix', help='prefix of IONEX files, either CODG or ROBR',dest="prefix",default='codg');
    parser.add_argument('--IONserver', help='server of IONEX files',dest="server",default='ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/')
    parser.add_argument('--IONpath', help='location of IONEX files',dest="ionexPath",default='./')
    parser.add_argument('--all','-a', help='calculate RM per station (default calculates only for CS002LBA)',action='store_true',dest="allStations")
    parser.add_argument('-t','--timestep', help='timestep in seconds. for values <=0 (default) the timegrid of the MS is used ',dest="timestep",type=float,default=900.)
    parser.add_argument('-e','--smart_interpol', help='float parameter describing how much of earth rotation is taken in to account in interpolation of the IONEX files. 1.0 means time interpolation assumes ionosphere rotates in opposite direction of the Earth. 0.0 (default) means no rotation applied',dest="earth_rot",type=float,default=0)
    parser.add_argument('-r','--overwrite',help='overwrite existing parmdb',action='store_true',dest="overwrite")

    args = parser.parse_args()
        
    stat_names=None
    if args.allStations:
        stat_names='all'
    if not args.overwrite and not os.path.isdir(args.parmdbName):
        args.overwrite=True
        
    createRMParmdb(args.MS,
                   parmdbname=args.parmdbName,
                   create=args.overwrite,
                   patchname=args.patch,
                   server=args.server,
                   prefix=args.prefix,
                   ionexPath=args.ionexPath,
                   earth_rot=args.earth_rot,
                   timestep=args.timestep,
                   stat_names=stat_names)
