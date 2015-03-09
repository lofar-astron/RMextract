#!/usr/bin/env python

#------------------------------------------------------
# Extract TEC values from an IONEX file
# given a specific time and geographic coordinate.

import numpy as np
import math
from datetime import date
import os
import time as systime




DEFAULT_TIMEOUT=100;

def readTEC(filename): 

	# Opening and reading the IONEX file into memory
	linestring = open(filename, 'r').read()
	LongList = linestring.split('\n')
	# creating a new array without the header and only
	# with the TEC maps
	exponent=0.1 #Default
	for i,line in enumerate(LongList):
		splitted=line.split()
		#print "evaluating",splitted;
		if splitted[-1] == 'DESCRIPTION' or splitted[-1] == 'COMMENT':
			continue;
		if splitted[-1] == 'FILE':
			if splitted[-2] == 'IN':
				NumberOfMaps = float(splitted[0])
				continue;
		if splitted[-1] == 'DHGT':
			IonH = float(splitted[0])
			continue;
		if splitted[-1] == 'EXPONENT':
			exponent = pow(10,float(splitted[0]))
			continue;
		if splitted[-1] == 'DLAT':
			startLat = float(splitted[0])
			endLat = float(splitted[1])
			stepLat = float(splitted[2])
			continue;
		if splitted[-1] == 'DLON':
			startLon = float(splitted[0])
			endLon = float(splitted[1])
			stepLon = float(splitted[2])
			continue;
		if splitted[-1]=='MAP' and (splitted[-4]+splitted[-2] == 'EPOCHFIRST'):
			startYear = float(splitted[0])
			startMonth = float(splitted[1])
			startDay = float(splitted[2])
			date=startYear*366.+startMonth*31.+startDay;
			continue;
		if splitted[0] == 'END':
			if splitted[2] == 'HEADER':
				break;
	NewLongList=LongList[i+1:];
	# Variables that indicate the number of points in Lat. and Lon.
	# 3D array that will contain TEC/RMS values only
	lonarray=np.arange(startLon,endLon+stepLon,stepLon);
	latarray=np.arange(startLat,endLat+stepLat,stepLat);
	pointsLon = lonarray.shape[0]
	pointsLat = latarray.shape[0]
	times=np.zeros(NumberOfMaps,dtype='float32');
	tecdata = np.zeros((NumberOfMaps,pointsLat,pointsLon))
	rmsdata = np.zeros((NumberOfMaps,pointsLat,pointsLon))
	start_fill=False;
	for line in NewLongList:
		splitted =line.split();
		#print "evaluating",splitted;
		if splitted[0]=='END' and splitted[2]=='FILE':
			break;
		if splitted[-1] == 'MAP' and splitted[-4] == 'START':
			start_fill=True;
			# found map start filling
			if splitted[-2] == 'TEC':
				fillarray=tecdata
			else:
				if splitted[-2] == 'RMS':
					fillarray=rmsdata
				else:
					start_fill=False;
					# something else
					continue;
			mapnr=int(splitted[0])-1;
			continue;
		if start_fill:
			if splitted[-1] == 'MAP' and splitted[1] == 'END':
				start_fill=False;
				continue;
			if splitted[-1] == 'MAP' and splitted[-4] == 'EPOCH':
				
				times[mapnr]=float(splitted[3])+float(splitted[4])/60.+float(splitted[5])/3600.;
				if (float(splitted[0])*366+float(splitted[1])*31+float(splitted[2]))>date: #next day
					times[mapnr]+=24.;
				continue;
			if splitted[-1] == 'LAT/LON1/LON2/DLON/H':
				latidx=np.argmin(np.absolute(latarray-float(line[:8])));
				lonidx=0;
				#print "latidx",latidx,float(line[:8])
				continue
			datalength=len(splitted)
			fillarray[mapnr,latidx,lonidx:lonidx+datalength]=\
			    np.array([float(i)*exponent for i in splitted]);
			lonidx+=datalength;

	return (tecdata,rmsdata,lonarray,latarray,times);




def getTECinterpol(time,lat,lon,tecinfo,apply_earth_rotation=False):
	'''Derive interpolated (in lon,lat,time) vTEC values. Lat should be angle in degrees between -90 and 90,Lon angle between -180,180. Time in hour of the day (eg. 23.5 for half past 11PM). With apply_earth_rotation you can specify (with a number between 0 and 1) how much of the earth rotaion is taken in to account in the interpolation step. This is assuming that the TEC maps move according to the rotation Earth (following method 3 of interpolation described in the IONEX document). Experiments with ROB data show that this is not really the case, resulting in strange wavelike structures when applying this smart interpolation. Negative values of this parameter would result in an ionosphere that moves with the roation of the earth''' 
	
	tecdata=tecinfo[0];  #TEC in TECU
	rmsdata=tecinfo[1];  # not used here
	lonarray=tecinfo[2]; # longitude in degrees from West to East (- to +)
	latarray=tecinfo[3]; # lattitude in degrees
	times=tecinfo[4];    # times in hour of the day (eg. 23.5 for half past 11PM)
	exactTime=False;

	#get indices of nearest time,lon,lat
	#assume time is  sorted from early to late
	timeIdx2=timeIdx1=np.argmin(np.absolute(times-time));
	if times[timeIdx1]>time:
		timeIdx1-=1;
		timeIdx2=timeIdx1;
	elif times[timeIdx1]==time:
		exactTime=True;
	# get rotation  angle for method 3 (rotate maps before interpolating)
	rot1=rot2=0;
	wt1=wt2=0.5;
	lonstep = abs(lonarray[1]-lonarray[0]);
	#assume lon is sorted from west to east BUG:: Should be checked! but ok for ODG and ROB files
	if not exactTime:
		# rotation in degrees of the earlier map
		rot1 = ((time-times[timeIdx1])*360./24.)*apply_earth_rotation;
		# rotation in degrees of the later map
		rot2 = ((times[timeIdx1+1]-time)*360./24.)*apply_earth_rotation;
		lonarray1=lonarray-rot1  #create new rotated longitude array
		lonarray2=lonarray+rot2   #create new rotated longitude array
		timeIdx2+=1
		timestep=times[timeIdx2]-times[timeIdx1];
		wt1=(times[timeIdx2]-time)/timestep; #weight of earlier map
		wt2=(time-times[timeIdx1])/timestep; #weight of later map
	else:
		lonarray1=lonarray
		lonarray2=lonarray
	exactLon=False;
	# index of longitude1 in earlier map
	lonIdx1=np.argmin(np.absolute(np.fmod(lonarray1-lon+540.,360.)-180.)); #make sure it lies in [-180,180]
	# index of longitude1 in later map
	lonIdx2=np.argmin(np.absolute(np.fmod(lonarray2-lon+540.,360.)-180.)); #make sure it lies in [-180,180]
	if lonarray1[lonIdx1]>lon or (lonIdx1==0 and (lonarray1[lonIdx1]+lonstep)<lon) :
		lonIdx1-=1;
		if lonIdx1 <0:
			lonIdx1+=lonarray.shape[0]; #BUG: this only works if you have coordinates of a full rotation (-180,180)
	if lonarray2[lonIdx2]>lon or (lonIdx2==0 and (lonarray2[lonIdx2]+lonstep)<lon) :
		lonIdx2-=1;
		if lonIdx2 <0:
			lonIdx2+=lonarray.shape[0];#BUG: this only works if you have coordinates of a full rotation (-180,180)
	elif lonarray1[lonIdx1]==lon and lonarray2[lonIdx2]==lon:
		exactLon=True;

	#indices used for method 3 interpol, remainder correctly treats negative indices, only works correctly if data has all longitudes between -180,180. Otherwise you get into trouble at the edges.
	lon12=lon11=np.remainder(lonIdx1,lonarray.shape[0]); #coordinates of earlier map
	lon22=lon21=np.remainder(lonIdx2,lonarray.shape[0]); #coordinates of later map
	wtlon11=wtlon21=wtlon12=wtlon22=0.5;
	if not exactLon: 
		lon12=np.remainder(lonIdx1+1,lonarray.shape[0]);
		lon22=np.remainder(lonIdx2+1,lonarray.shape[0]);
		wtlon11=np.remainder(lonarray1[lon12]-lon,360.)/lonstep;
		wtlon12=np.remainder(lon-lonarray1[lon11],360.)/lonstep;
		wtlon21=np.remainder(lonarray2[lon22]-lon,360.)/lonstep;
		wtlon22=np.remainder(lon-lonarray2[lon21],360.)/lonstep;
		
		

	exactLat=False;
	latIdx1=np.argmin(np.absolute(latarray-lat));
	#find orientation
	if(latarray[0]>latarray[1]):
		NSorient=True
	else:
		NSorient=False

	if NSorient:
		if latarray[latIdx1]<lat:
			latIdx1-=1;
		elif latarray[latIdx1]==lat or latIdx1==(latarray.shape[0]-1):
			exactLat=True;
	else:
		if latarray[latIdx1]>lat:
			latIdx1-=1;
		elif latarray[latIdx1]==lat or latIdx1==(latarray.shape[0]-1):
			exactLat=True;
	lat1=lat2=latIdx1;
	wtlat1=wtlat2=0.5;
		
	if not exactLat:
		lat2+=1;
		latstep = abs(latarray[lat1]-latarray[lat2]);
		wtlat1=abs(lat-latarray[lat2])/latstep;
		wtlat2=abs(latarray[lat1]-lat)/latstep;
	#now get all needed maps
	timeinterpols=[wt1*tecdata[timeIdx1,lat1,lon11]*wtlon11+wt2*tecdata[timeIdx2,lat1,lon21]*wtlon21,
		       wt1*tecdata[timeIdx1,lat1,lon12]*wtlon12+wt2*tecdata[timeIdx2,lat1,lon22]*wtlon22,
		       wt1*tecdata[timeIdx1,lat2,lon11]*wtlon11+wt2*tecdata[timeIdx2,lat2,lon21]*wtlon21,
		       wt1*tecdata[timeIdx1,lat2,lon12]*wtlon12+wt2*tecdata[timeIdx2,lat2,lon22]*wtlon22]#weighted points
	

	tecvalue=timeinterpols[0]*wtlat1 + timeinterpols[1]*wtlat1 + \
	    timeinterpols[2]*wtlat2 + timeinterpols[3]*wtlat2
	return tecvalue;


def run_command_timeout(command, args, timeout):
    """run some command with a timeout limit

INPUTS:
command   I  The name of the program to run
args      I  The full argument list tot he command.  Note that this should
             be an array, with element 0 the name of the program
timeout   I  The timeout, in s.  If the program has not finished within
             timeout seconds, kill it

OUTPUTS: retval
retval    O  The return code of the command.  Note that if the process times out,
             this function raises a Command_Timeout_Error


    """
    pid = os.fork()
    print "my pid",pid
    if(pid == 0):
        # We are the child, execute the command
	print "executing",command,args
        os.execvp(command, args)
        # If we are here, something bad happened.  It must be the
        # user's fault.  :)
        os._exit(2)
    elif(pid < 0):
        # fork failed
        raise OSError("Hey, I couldn't fork!")
    for i in xrange(int(timeout)):
        print "waiting",i,timeout
        systime.sleep(1)
        status = os.waitpid(pid,os.WNOHANG)
        if(status == (0,0)):
            # still waiting
            continue
        if(os.WIFEXITED(status[1])):
            return os.WEXITSTATUS(status[1])
        if(os.WIFSIGNALED(status[1])):
            # someone killed the process.  This is probably bad, and means the
            # user got tired of waiting.  Try calling this a timeout
            raise Command_Timeout_Error("Command '%s' signalled with %d"%(command, os.WTERMSIG(status[1])))
    os.kill(pid,signal.SIGKILL)
    raise Command_Timeout_Error("Command '%s' timed out after %d seconds"%(command, timeout))

def gunzip_some_file(compressed_file,
                     uncompressed_file,
                     delete_file = 1):
    # make sure there is a file
    if not os.path.isfile(compressed_file):
        if(compressed_file[-1] == 'Z'):
            # try another form
            warnings.warn("No .Z form of compressed file '%s', trying .gz"%compressed_file)
            new_compressed = compressed_file[:-1] + "gz"
            return gunzip_some_file(new_compressed, uncompressed_file,
                                    delete_file)
        raise RuntimeError, "No such file '%s' to uncompress"%compressed_file
    command = "gunzip -dc %s > %s"%(compressed_file,uncompressed_file)
    retcode = os.system(command)
    if(retcode):
        raise RuntimeError, "Could not run '%s'"%command
    if(delete_file):
        os.remove(compressed_file)
    return

		
def getIONEXfile(time="2012/03/23/02:20:10.01",server="ftp://ftp.unibe.ch/aiub/CODE/",prefix="CODG",outpath='',overwrite=False):
	yy=time[2:4];
	year=int(time[:4])
	month=int(time[5:7])
	day=int(time[8:10])
	dayofyear = date(year,month,day).timetuple().tm_yday;
	if prefix!='ROBR':
		filename=str(year)+"/"+prefix+"%03d0.%sI.Z"%(dayofyear,yy);
		backupfilename = str(year)+"/%03d/"%(dayofyear)+prefix+"%03d0.%si.Z"%(dayofyear,yy)
	else:
		filename=str(year)+"/"+prefix+"%03d.%sI.Z"%(dayofyear,yy);
	#print 'file needed:', filename;
	#print "checking",outpath+filename[5:-2],os.path.isfile(outpath+filename[5:-2])
	if not overwrite and os.path.isfile(outpath+filename[5:-2]):
		print "file exists",outpath+filename[5:-2];
		return outpath+filename[5:-2];
	systime.sleep(2)
	print "retreiving",run_command_timeout("URL_download.py",
					       ["URL_download.py", server+filename, outpath+filename[5:],
						"%d"%DEFAULT_TIMEOUT],
					       DEFAULT_TIMEOUT)
	if not os.path.isfile(outpath+filename[5:]):
		print "trying",server+backupfilename
		#try backup name
		systime.sleep(2)
		run_command_timeout("URL_download.py",
				    ["URL_download.py", server+backupfilename, outpath+filename[5:],
				     "%d"%DEFAULT_TIMEOUT],
				    DEFAULT_TIMEOUT)
	count=0
	while count<10 and (not os.path.isfile(outpath+filename[5:]) and not os.path.isfile(outpath+filename[5:-2])):
		print "waiting..."
		systime.sleep(1)
		count+=1
		
	if os.path.isfile(outpath+filename[5:]):

		gunzip_some_file(outpath+filename[5:],outpath+filename[5:-2]);
		
		return outpath+filename[5:-2];
	elif os.path.isfile(outpath+filename[5:-2]):
		#no need to unzip
		return outpath+filename[5:-2]
	else:
		print "file",filename,"not found on server",server;
		return -1;
