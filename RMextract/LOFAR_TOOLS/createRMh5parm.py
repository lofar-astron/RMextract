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
import os
import numpy as np
import sys
import logging


def createh5Parm(ms, h5parmdb, solset_name = "sol000",all_stations=False,timestep=300,
                server="ftp://ftp.aiub.unibe.ch/CODE/",prefix='CODG',ionexPath="./",earth_rot=0):
    '''Add rotation measure to existing h5parmdb
    
    Args:
        ms (string) :  path + filename of Measurement Set 
        h5parmdb (string) : name of existing h5parm
        solset_name (string) : optional name of solset in h5parmdb, 
                            if not set, first one will be chosen
        all_stations (bool) : optional calculate RM values for all stations in the MS,'
                            default only for position of CS002LBA 
        timestep (float) : timestep in seconds
        server (string) : ftp server for IONEX files
        prefix (string) : prefix of IONEX files
        ionexPath (string) : location where IONEX files are stored
        earth_rot (float) : parameter to determine how much earth rotation is taken \
        into account when interpolating IONEX data. (0 is none, 1 is full)
    '''
    
    data          = h5parm(h5parmdb, readonly=False)
    solset        = data.getSolset(solset_name)
    station_names = solset.getAnt().keys()
    if not all_stations:
        rmdict = getRM.getRM(MS, 
                             server=server, 
                             prefix=prefix, 
                             ionexPath=ionexPath, 
                             timestep=timestep,
                             stat_pos=PosTools.posCS002,
                             earth_rot=earth_rot)
        if not rmdict:
            if not server:
                raise ValueError("One or more IONEX files is not found on disk and download is disabled!\n"
                                 "(You can run \"bin/download_IONEX.py\" outside the pipeline if needed.)")
            else:
                raise ValueError("Couldn't get RM information from RMextract! (But I don't know why.)")

        rm_vals=np.ones((len(station_names),rmdict['RM']['st0'].shape[0]),dtype=float)
        rm_vals+=rmdict['RM']['st0'][np.newaxis]
    else:
        rmdict = getRM.getRM(MS, 
                             server=server, 
                             prefix=prefix, 
                             ionexPath=ionexPath, 
                             timestep=timestep,
                             earth_rot=earth_rot)
        if not rmdict:
            if not server:
                raise ValueError("One or more IONEX files is not found on disk and download is disabled!\n"
                                 "(You can run \"bin/download_IONEX.py\" outside the pipeline if needed.)")
            else:
                raise ValueError("Couldn't get RM information from RMextract! (But I don't know why.)")

        rm_vals=np.array([rmdict["RM"][stat] for stat in station_names])

    logging.info('Adding rotation measure values to: ' + solset_name + ' of ' + h5parmdb)
    try:
        new_soltab = solset.getSoltab('RMextract')
        new_soltab.delete()
    except:
        pass
    new_soltab = solset.makeSoltab(soltype='rotationmeasure', soltabName='RMextract',
                                   axesNames=['ant', 'time'], axesVals=[station_names, rmdict['times']],
                                   vals=rm_vals,
                                   weights=np.ones_like(rm_vals, dtype=np.float16))
    new_soltab.addHistory('CREATE (by createRMh5parm (RMextract) script)')


    
    ########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Adds CommonRotationAngle to an H5parm from RMextract.')

    parser.add_argument('MSfile', type=str,
                        help='MS for which the parmdb should be created.')
    parser.add_argument('h5parm', type=str,
                        help='H5parm to which the results of the CommonRotationAngle is added.')
    parser.add_argument('--server', type=str, default='None',
                        help='URL of the server to use. (default: None)')
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

    MS = args.MSfile
    h5parmdb = args.h5parm
    logging.info("Working on:", MS, h5parmdb)
    createh5Parm(MS, h5parmdb, ionex_server=args.server, ionex_prefix=args.prefix, 
                 ionexPath=args.ionexpath, solset_name=args.solsetName, 
                 all_stations=args.allStations, timestep=args.timestep,
                 earth_rot=args.earth_rot)
    
    
    
    