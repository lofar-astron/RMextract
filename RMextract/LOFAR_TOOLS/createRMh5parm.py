#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract RM values for a LOFAR observation. Fill h5parm

Created on Tue Aug 7 11:46:57 2018

@author: mevius
"""
from losoto.h5parm import h5parm
from RMextract import getRM
from RMextract import PosTools
import pyrap.tables as pt
import os
import numpy as np
import sys
import logging


def makesolset(MS, data, solset_name):
    solset = data.makeSolset(solset_name)    

    antennaFile = MS+"/ANTENNA"
    logging.info('Collecting information from the ANTENNA table.')
    antennaTable = pt.table(antennaFile, ack=False)
    antennaNames = antennaTable.getcol('NAME')
    antennaPositions = antennaTable.getcol('POSITION')
    antennaTable.close()
    antennaTable = solset.obj._f_get_child('antenna')
    antennaTable.append(zip(*(antennaNames,antennaPositions)))
    
    fieldFile = MS + "/FIELD"
    logging.info('Collecting information from the FIELD table.')
    fieldTable = pt.table(fieldFile, ack=False)
    phaseDir = fieldTable.getcol('PHASE_DIR')
    pointing = phaseDir[0, 0, :]
    fieldTable.close()

    sourceTable = solset.obj._f_get_child('source')
    # add the field centre, that is also the direction for Gain and Common*
    sourceTable.append([('pointing',pointing)])




def main(MSfiles, h5parmdb, solset_name = "sol000",all_stations=False,timestep=300,
                ionex_server="ftp://ftp.aiub.unibe.ch/CODE/",
                ionex_prefix='CODG',ionexPath="./",earth_rot=0):
    '''Add rotation measure to existing h5parmdb
    
    Args:
        MSfiles (string) :  path + filename of Measurement Set 
        h5parmdb (string) : name of existing h5parm
        solset_name (string) : optional name of solset in h5parmdb, 
                            if not set, first one will be chosen
        all_stations (bool) : optional calculate RM values for all stations in the MS,'
                            default only for position of CS002LBA 
        timestep (float) : timestep in seconds
        ionex_server (string) : ftp server for IONEX files
        ionex_prefix (string) : prefix of IONEX files
        ionexPath (string) : location where IONEX files are stored
        earth_rot (float) : parameter to determine how much earth rotation is taken \
        into account when interpolating IONEX data. (0 is none, 1 is full)
    '''
    
    mslist = MSfiles.lstrip('[').rstrip(']').replace(' ','').replace("'","").split(',')
    
    if len(mslist) == 0:
        raise ValueError("Did not find any existing directory in input MS list!")
        pass
    else:
        MS = mslist[0]
        pass
    
    if not all_stations:
        rmdict = getRM.getRM(MS, 
                             server=ionex_server, 
                             prefix=ionex_prefix, 
                             ionexPath=ionexPath, 
                             timestep=timestep,
                             stat_names = ["st0"],
                             stat_pos=[PosTools.posCS002],
                             earth_rot=earth_rot)

    else:
        rmdict = getRM.getRM(MS, 
                             server=ionex_server, 
                             prefix=ionex_prefix, 
                             ionexPath=ionexPath, 
                             timestep=timestep,
                             earth_rot=earth_rot)
    if not rmdict:
        if not server:
            raise ValueError("One or more IONEX files is not found on disk and download is disabled!\n"
                                 "(You can run \"bin/download_IONEX.py\" outside the pipeline if needed.)")
        else:
            raise ValueError("Couldn't get RM information from RMextract! (But I don't know why.)")
 
    data          = h5parm(h5parmdb, readonly=False)
    if not solset_name in data.getSolsetNames():
        makesolset(MS,data,solset_name)
    solset        = data.getSolset(solset_name)
    station_names = sorted(solset.getAnt().keys())

 
    logging.info('Adding rotation measure values to: ' + solset_name + ' of ' + h5parmdb)
    if all_stations:
        rm_vals=np.array([rmdict["RM"][stat] for stat in station_names])[:,:,0]
    else:
        rm_vals=np.ones((len(station_names),rmdict['RM']['st0'].shape[0]),dtype=float)
        rm_vals+=rmdict['RM']['st0'][:,0][np.newaxis]
    new_soltab = solset.makeSoltab(soltype='rotationmeasure', soltabName='RMextract',
                                   axesNames=['ant', 'time'], axesVals=[station_names, rmdict['times']],
                                   vals=rm_vals,
                                   weights=np.ones_like(rm_vals, dtype=np.float16))


    
    ########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Adds CommonRotationAngle to an H5parm from RMextract.')

    parser.add_argument('MSfiles', type=str,
                        help='MS for which the parmdb should be created.')
    parser.add_argument('h5parm', type=str,
                        help='H5parm to which the results of the CommonRotationAngle is added.')
    parser.add_argument('--server', type=str, default='ftp://ftp.aiub.unibe.ch/CODE/',
                        help='URL of the server to use. (default: ftp://ftp.aiub.unibe.ch/CODE/)')
    parser.add_argument('--prefix', type=str, default='CODG',
                        help='Prefix of the IONEX files. (default: \"CODG\")')
    parser.add_argument('--ionexpath', '--path', type=str, default='IONEXdata/',
                        help='Directory where to store the IONEX files. (default: \"IONEXdata/\")')
    parser.add_argument('--solsetName', '--solset', type=str, default='sol000',
                        help='Name of the h5parm solution set (default: sol000)')
    parser.add_argument('--all','-a', help=
                        'calculate RM per station (default calculates only for CS002LBA)',
                        action='store_true',dest="allStations")
    parser.add_argument('-t','--timestep', help=
                        'timestep in seconds. for values <=0 (default) the timegrid of the MS is used ',
                        dest="timestep",type=float,default=300.)
    parser.add_argument('-e','--smart_interpol', help=
                        'float parameter describing how much of Earth rotation is taken in to account \
                        in interpolation of the IONEX files. 1.0 means time interpolation assumes \
                        ionosphere rotates in opposite direction of the Earth. 0.0 (default) means \
                        no rotation applied',dest="earth_rot",type=float,default=0)

    args = parser.parse_args()

    MS = args.MSfiles
    h5parmdb = args.h5parm
    logging.info("Working on:", MS, h5parmdb)
    main(MS, h5parmdb, ionex_server=args.server, ionex_prefix=args.prefix, 
                 ionexPath=args.ionexpath, solset_name=args.solsetName, 
                 all_stations=args.allStations, timestep=args.timestep,
                 earth_rot=args.earth_rot)
    
