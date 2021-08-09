import sys
import RMextract.getIONEX as ionex
from RMextract import PosTools
import numpy as np
import logging
import argparse


def main():
    parser = argparse.ArgumentParser(description='Downloads relevant IONEX files for a given MS')
    parser.add_argument('MSfiles', type=str,
                            help='MS for which the parmdb should be created.')
    parser.add_argument('--server', type=str, default='ftp://ftp.aiub.unibe.ch/CODE/',
                            help='URL of the server to use. (default: ftp://ftp.aiub.unibe.ch/CODE/)')
    parser.add_argument('--prefix', type=str, default='CODG',
                            help='Prefix of the IONEX files. (default: \"CODG\")')
    parser.add_argument('--ionexpath', '--path', type=str, default='IONEXdata/',
                            help='Directory where to store the IONEX files. (default: \"IONEXdata/\")')
    args = parser.parse_args()

    MS=args.MSfiles
    server=args.server
    prefix=args.prefix
    (timerange,timestep,pointing,stat_names,stat_pos)=PosTools.getMSinfo(MS)
    
    start_time = timerange[0]
    end_time = timerange[1]
    timerange[0] = start_time - timestep
    timerange[1] = end_time + timestep
    times,timerange=PosTools.getIONEXtimerange(timerange,timestep)
    if len(times[-1])==0 or times[-1][-1]<timerange[1]:
        timestmp=list(times[-1])
        timestmp.append(timerange[1]) #add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
        times[-1]=np.array(timestmp)

    for time_array in times[::-1]: #check in reverse order, since datamight not be available for latest days
        starttime=time_array[0]
        date_parms =  PosTools.obtain_observation_year_month_day_fraction(starttime)
        if not "http" in server: #ftp server use ftplib
            ionexf=ionex.getIONEXfile(time=date_parms,server=server,prefix=prefix,outpath=args.ionexpath)
        else:
            ionexf=ionex.get_urllib_IONEXfile(time=date_parms,server=server,prefix=prefix,outpath=args.ionexpath)
        if ionexf==-1:
            if not "igsiono.uwm.edu.pl" in server:
                logging.info("cannot get IONEX data, try fast product server instead")
                server="https://igsiono.uwm.edu.pl"
                prefix="igrg"
                ionexf=ionex.get_urllib_IONEXfile(time=date_parms,server=server,prefix=prefix,outpath=args.ionexpath)
            if ionexf==-1:
                logging.error("IONEX data not available, even not from fast product server")
                return False
        
    return True

