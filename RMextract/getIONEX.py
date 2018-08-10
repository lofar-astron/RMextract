#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract TEC values from an IONEX file given a specific time and geographic coordinate.

Created on Tue Apr 24 11:46:57 2018

@author: mevius
"""
import numpy as np
import datetime
import scipy.ndimage.filters as myfilter
import logging
import os
import ftplib

logging.basicConfig(level=logging.DEBUG)


def _read_ionex_header(filep):
    """reads header from ionex file. returns data shape and position of first
    data in the file.

    Args:
        filep (filepointer) : pointer to opened ionex file.

    Returns:
        Tuple[float, np.array, np.array, np.array]:
            multiplication factor,lonarray,latarray,timearray
    """

    filep.seek(0)
    for line in filep:
        if "END OF HEADER" in line:
            break
        stripped = line.strip()
        if stripped.endswith("EPOCH OF FIRST MAP"):
            starttime = datetime.datetime(
                *(int(i) for i in
                  stripped.strip("EPOCH OF FIRST MAP").split()))
        if stripped.endswith("EPOCH OF LAST MAP"):
            endtime = datetime.datetime(
                *(int(i) for i in
                  stripped.strip("EPOCH OF LAST MAP").split()))
        if stripped.endswith("INTERVAL"):
            timestep = float(stripped.split()[0]) / 3600.
        if stripped.endswith("EXPONENT"):
            exponent = pow(10, float(stripped.split()[0]))
        if stripped.endswith("DLON"):
            start_lon, end_lon, step_lon = \
                (float(i) for i in stripped.split()[:3])
        if stripped.endswith("DLAT"):
            start_lat, end_lat, step_lat = \
                (float(i) for i in stripped.split()[:3])
        if stripped.endswith("OF MAPS IN FILE"):
            ntimes = int(stripped.split()[0])

    lonarray = np.arange(start_lon, end_lon + step_lon, step_lon)
    latarray = np.arange(start_lat, end_lat + step_lat, step_lat)
    dtime = endtime - starttime
    dtimef = dtime.days * 24. + dtime.seconds / 3600.
    logging.debug("timerange %f hours. step = %f ", dtimef, timestep)
    timearray = np.arange(0,
                          dtimef + timestep,
                          timestep)
    if timearray.shape[0] < ntimes:
        # bug in ILTF files,last time in header is incorrect
        extratimes = np.arange(timearray[-1] + timestep,
                               timearray[-1]
                               + (ntimes -
                                  timearray.shape[0] + 0.5) * timestep,
                               timestep)
        timearray = np.concatenate((timearray, extratimes))
    timearray += starttime.hour\
        + starttime.minute/60.\
        + starttime.second/3600.

    return exponent, lonarray, latarray, timearray


def read_tec(filename, _use_filter=None):
    """ returns TEC, RMS longitude, lattitude and time read from an IONEX file.

    Args:
        filename (string) : the full path to the IONEXfile
        _use_filter (float) : optional filter the data in space and time
                             with a gaussian filter with sigma use_filter.
                             calls scipy.ndimage.filter.gaussian_filter(tec,
                             use_filter)

    Returns:
        Tuple[np.array, np.array, np.array, np.array, np.array]:
            3D-arrays (time,lat,lon) of (optionally filtered) TEC and RMS +
            longitude, latitude and time array
    """
    ionex_file = open(filename, "r")
    exponent, lonarray, latarray, timearray = _read_ionex_header(ionex_file)
    logging.info("reading data with shapes %d  x %d x %d",
                 timearray.shape[0],
                 latarray.shape[0],
                 lonarray.shape[0])
    tecarray = np.zeros(timearray.shape
                        + latarray.shape + lonarray.shape, dtype=float)
    rmsarray = np.zeros_like(tecarray)
    timeidx = 0
    lonidx = 0
    latidx = 0
    tecdata = False
    rmsdata = False
    readdata = False
    for line in ionex_file:
        if "START OF TEC MAP" in line:
            tecdata = True
            rmsdata = False
            timeidx = int(line.strip().split()[0]) - 1
            continue
        if "START OF RMS MAP" in line:
            rmsdata = True
            tecdata = False
            timeidx = int(line.strip().split()[0]) - 1
            continue
        if "LAT/LON1/LON2/DLON/H" in line:
            readdata = True
            latstr = line.strip().strip("LAT/LON1/LON2/DLON/H")
            lat = np.fromstring(" -".join(latstr.split("-")), sep=" ")
            latidx = np.argmin(np.abs(latarray - lat[0]))
            lonidx = 0
            continue
        if tecdata and ("END OF TEC MAP" in line):
            readdata = False
            continue
        if rmsdata and ("END OF RMS MAP" in line):
            readdata = False
            continue
        if readdata:
            data = np.fromstring(" -".join(line.strip().split("-")),
                                 sep=" ") * exponent
            if tecdata:
                tecarray[timeidx, latidx, lonidx:lonidx + data.shape[0]] = data
            elif rmsdata:
                rmsarray[timeidx, latidx, lonidx:lonidx + data.shape[0]] = data
            lonidx += data.shape[0]
    if not _use_filter is None:
        tecarray = myfilter.gaussian_filter(
            tecarray, _use_filter, mode='nearest')
    return tecarray, rmsarray, lonarray, latarray, timearray


def readTEC(filename, use_filter=None):
    """oldfunction name for compatibility. Use read_tec."""

    logging.warning("function readTEC obsolete, use read_tec instead")
    return read_tec(filename, _use_filter=use_filter)


def _compute_index_and_weights(maparray, mapvalues):
    '''helper function  to get indices and weights for interpolating tecmaps


    Args:
        maparray (np.array) : array to get indices in
        mapvalues (Union[float,np.array]) :  values to get indices for
    Returns:
        Tuple[np.array, np.array, np.array]: idx1,idx2 and weights for idx2,
                                             idx2 is always >= idx1

    '''
    is_reverse = maparray[1] < maparray[0]
    idx1 = np.argmin(np.absolute(maparray[np.newaxis]
                                 - mapvalues[:, np.newaxis]), axis=1)
    idx2 = idx1.copy()
    if not is_reverse:
        idx1[maparray[idx1] > mapvalues] -= 1
        idx2[maparray[idx2] < mapvalues] += 1
    else:
        idx1[maparray[idx1] < mapvalues] -= 1
        idx2[maparray[idx2] > mapvalues] += 1
    idx1[idx1 < 0] = 0
    idx2[idx2 < 0] = 0
    idx1[idx1 >= maparray.shape[0]] = maparray.shape[0] - 1
    idx2[idx2 >= maparray.shape[0]] = maparray.shape[0] - 1
    _steps = np.absolute(maparray[idx2] - maparray[idx1])
    weights = np.absolute(mapvalues - maparray[idx1])
    _steps[_steps == 0] = weights[_steps == 0]
    weights = weights / _steps
    return idx1, idx2, weights


def compute_tec_interpol(times, lats, lons, tecinfo, apply_earth_rotation=0):
    '''Get interpolated TEC for array of times/lats and lons

    Derive interpolated (4 point in lon,lat,2 point in time) vTEC values,
    optionally correcting for earth rotation.

    Args:
        lats (np.array) : angles in degrees between -90 and 90
        lons (np.array) : angles in degrees between -180,180
        times (np.array) : times in decimal hour of the day (eg. 23.5
        for half past 11PM)
        tecinfo (tuple) : tuple with the return values of read_tec function
        apply_earth_rotation(float) : specify (with a number between 0 and 1)
        how much of the earth rotaion is taken in to account in the
        interpolation step.
        This is assuming that the TEC maps move according to the rotation
        Earth (following method 3 of interpolation described in the IONEX
        document). Experiments with high time resolution ROB data show that
        this is not really the case, resulting in strange wavelike structures
        when applying this smart interpolation. Negative values of this
        parameter would result in an ionosphere that moves with the rotation
        of the earth
    Returns:
        np.array : interpolated tecvalues
    '''
    assert times.shape == lats.shape,\
        "times and lats should be array with same shape"
    assert times.shape == lons.shape, \
        "times and lons should be array with same shape"
    tecdata = tecinfo[0]  # TEC in TECU
    lonarray = tecinfo[2]  # longitude in degrees from West to East (- to +)
    latarray = tecinfo[3]  # lattitude in degrees
    # times in hour of the day (eg. 23.5 for half past 11PM)
    maptimes = tecinfo[4]

    # get indices of nearest 2 time frames + inverse distance weights
    # assume time is  sorted from early to late
    timeidx1, timeidx2, time_weights = _compute_index_and_weights(
        maptimes, times)

    # latarray is sorted small to large
    latidx1, latidx2, lat_weights = _compute_index_and_weights(
        latarray, lats)

    # for getting lon idx take into account earth_rotation
    # if longitudes cover all range between -180 and 180 you can modulate the
    # indices, otherwise we have to take the edge of the map.
    lonstep = lonarray[1] - lonarray[0]
    full_circle = np.remainder(lonarray[0] - lonarray[-1], 360.)\
        <= 1.1 * lonstep

    rot1 = ((times - maptimes[timeidx1]) * 360. / 24.) * apply_earth_rotation
    rot2 = ((times - maptimes[timeidx2]) * 360. / 24.) * apply_earth_rotation

    if not full_circle:
        lonidx11, lonidx12, lon_weights1 = _compute_index_and_weights(
            lonarray, lons + rot1)
        lonidx21, lonidx22, lon_weights2 = _compute_index_and_weights(
            lonarray, lons + rot2)
    else:
        lonidx11 = np.argmin(np.absolute(np.remainder(lonarray[np.newaxis]
                                                      - rot1[:, np.newaxis]
                                                      - lons[:, np.newaxis]
                                                      + 180., 360.)
                                         - 180.), axis=1)
        lonidx12 = lonidx11.copy()
        lonidx11[np.remainder(lonarray[lonidx11] - rot1 - lons + 180., 360.)
                 - 180. > 0] -= 1
        lonidx12[np.remainder(lonarray[lonidx12] - rot1 - lons + 180., 360.)
                 - 180. < 0] += 1
        lonidx11[lonidx11 < 0] += lonarray.shape[0]
        lonidx12[lonidx12 < 0] += lonarray.shape[0]
        lonidx11[lonidx11 >= lonarray.shape[0]] -= lonarray.shape[0]
        lonidx12[lonidx12 >= lonarray.shape[0]] -= lonarray.shape[0]
        lon_weights1 = np.absolute(np.remainder(lonarray[lonidx11]
                                                - rot1
                                                - lons + 180., 360.)
                                   - 180.) / lonstep

        lonidx21 = np.argmin(np.absolute(np.remainder(lonarray[np.newaxis]
                                                      - rot2[:, np.newaxis]
                                                      - lons[:, np.newaxis]
                                                      + 180., 360.)
                                         - 180.), axis=1)
        lonidx22 = lonidx21.copy()
        lonidx21[np.remainder(lonarray[lonidx21] - rot2 - lons + 180., 360.)
                 - 180. > 0] -= 1
        lonidx22[np.remainder(lonarray[lonidx22] - rot2 - lons + 180., 360.)
                 - 180. < 0] += 1
        lonidx21[lonidx21 < 0] += lonarray.shape[0]
        lonidx22[lonidx22 < 0] += lonarray.shape[0]
        lonidx21[lonidx21 >= lonarray.shape[0]] -= lonarray.shape[0]
        lonidx22[lonidx22 >= lonarray.shape[0]] -= lonarray.shape[0]
        lon_weights2 = np.absolute(np.remainder(lonarray[lonidx21]
                                                - rot2
                                                - lons + 180., 360.)
                                   - 180.) / lonstep
    logging.debug("inidces time %d %d indices lat %d %d indices \
                  lon %d %d %d %d", timeidx1[0], timeidx2[0],
                  latidx1[0], latidx2[0],
                  lonidx11[0], lonidx12[0],
                  lonidx21[0], lonidx22[0])
    logging.debug("weights time %f lat %f lon %f %f",
                  time_weights[0],
                  lat_weights[0],
                  lon_weights1[0],
                  lon_weights2[0])
    tecs = (tecdata[timeidx1, latidx1, lonidx11] * (1. - lon_weights1)
            + tecdata[timeidx1, latidx1, lonidx12] * lon_weights1) \
        * (1. - time_weights)
    tecs += (tecdata[timeidx2, latidx1, lonidx21] * (1. - lon_weights2)
             + tecdata[timeidx2, latidx1, lonidx22] * lon_weights2) \
        * (time_weights)
    tecs *= 1. - lat_weights
    tecs += lat_weights \
        * (tecdata[timeidx1, latidx2, lonidx11] * (1. - lon_weights1)
           + tecdata[timeidx1, latidx2, lonidx12] * lon_weights1) \
        * (1. - time_weights)
    tecs += lat_weights \
        * (tecdata[timeidx2, latidx2, lonidx21] * (1. - lon_weights2)
           + tecdata[timeidx2, latidx2, lonidx22] * lon_weights2) \
        * (time_weights)
    return tecs


def getTECinterpol(time, lat, lon, tecinfo, apply_earth_rotation=0):
    """old function name for compatibility. Use compute_tec_interpol
    instead"""

    logging.warning("obsolete, use compute_tec_interpol instead")
    if np.isscalar(time):
        time = [time]
        lat = [lat]
        lon = [lon]
    return compute_tec_interpol(np.array(time), np.array(lat), np.array(lon),
                                tecinfo,
                                apply_earth_rotation)


def _combine_ionex(outpath, filenames, newfilename):
    """Helper function to combine separate IONEXfiles into 1 single file
    (needed for 15min ROBR data)"""

    if os.path.isfile(outpath + newfilename):
        logging.info("FILE exists: " + outpath + newfilename)
        return outpath + newfilename
    newf = open(outpath + newfilename, 'w')
    filenames = sorted(filenames)
    firstfile = open(filenames[0])
    lastfile = open(filenames[-1])
    for line in lastfile:
        if "EPOCH OF LAST MAP" in line:
            lastepoch = line
            lastfile.close()
            break
    # write header + tec map
    for line in firstfile:
        if "END OF TEC MAP" in line:
            break
        if "EPOCH OF LAST MAP" not in line:
            if "OF MAPS IN FILE" in line:
                newf.write(line.replace('1', str(len(filenames))))
            else:
                newf.write(line)
        else:
            newf.write(lastepoch)
    tecmapnr = 2
    for myfname in filenames[1:]:
        myf = open(myfname)
        end_of_header = False
        for line in myf:
            if not end_of_header and "END OF HEADER" in line:
                end_of_header = True
            else:
                if end_of_header:
                    if "END OF TEC MAP" in line:
                        newf.write(line.replace('1', str(tecmapnr)))
                        break
                    if "START OF TEC MAP" in line:
                        newf.write(line.replace('1', str(tecmapnr)))
                    else:
                        newf.write(line)
        tecmapnr += 1
    newf.write("END OF FILE\n")
    return os.path.join(outpath, newfilename)
    # ignore RMS map for now, since it is filled with zeros anyway


def _gunzip_some_file(compressed_file,
                      uncompressed_file,
                      delete_file=1):
    command = "gunzip -dc %s > %s" % (compressed_file, uncompressed_file)
    retcode = os.system(command)
    if retcode:
        raise RuntimeError("Could not run '%s'" % command)
    if delete_file:
        os.remove(compressed_file)
    return uncompressed_file


def _store_files(ftp, filenames, outpath, overwrite=False):
    """helper function to store files from ftp server to outpath"""

    nfilenames = []
    for myf in filenames:
        #make sure filename is always stored uppercase
        mypath = os.path.join(outpath, myf.upper())
        if not overwrite and os.path.isfile(mypath):
            nfilenames.append(mypath)
            logging.info("file %s exists", mypath)
        elif not overwrite and\
        os.path.isfile(mypath.strip(".Z")):
            nfilenames.append(mypath.strip(".Z"))
            logging.info("file %s exists",
                         mypath.strip(".Z"))
        else:
            myp = open(mypath, "wb")
            ftp.retrbinary("RETR " + myf, myp.write)
            myp.close()
            nfilenames.append(mypath)
    nfilenames_copy = nfilenames[:]
    nfilenames = []
    for myf in nfilenames_copy:
        if myf.endswith(".Z"):
            nfilenames.append(_gunzip_some_file(
                myf,
                myf.strip(".Z")))
        else:
            nfilenames.append(myf)
    return nfilenames


def _get_IONEX_file(time="2012/03/23/02:20:10.01",
                    server="ftp://cddis.gsfc.nasa.gov/gnss/productsionex/",
                    prefix="codg",
                    outpath='./',
                    overwrite=False):
    """Get IONEX file with prefix from server for a given day

    Downloads files with given prefix from the ftp server, unzips and stores
    the data. For prefix ROBR the data is stored on the server in a separate
    IONEX file every 15 minutes, these are automatically combined for
    compatibility.

    Args:
        time (string or list) : date of the observation
        server (string) : ftp server + path to the ionex directories
        prefix (string) : prefix of the IONEX files (case insensitive)
        outpath (string) : path where the data is stored
        overwrite (bool) : Do (not) overwrite existing data
    """

    if outpath[-1] != "/":
        outpath += "/"
    if not os.path.isdir(outpath):
        try:
            os.makedirs(outpath)
        except:
            logging.error("cannot create output dir for IONEXdata: %s",
                          outpath)
            raise
    try:
        yy = time[2:4]
        year = int(time[:4])
        month = int(time[5:7])
        day = int(time[8:10])
    except:
        year = time[0]
        yy = year - 2000
        month = time[1]
        day = time[2]
    mydate = datetime.date(year, month, day)
    dayofyear = mydate.timetuple().tm_yday
    ftpserver = server.strip("ftp:").strip("/").split("/")[0]
    ftppath = "/".join(server.strip("ftp:").strip("/").split("/")[1:])
    ftp = ftplib.FTP(ftpserver)
    try:
        ftp.login()
    except ftplib.error_perm:
        ftp.login("data-out", "Qz8803#mhR4z")
    ftp.cwd(ftppath)
    totpath = ftppath
    myl = []
    ftp.retrlines("NLST", myl.append)
    logging.info("Retrieving data for %d or %02d%03d", year, yy, dayofyear)
    if str(year) in myl:
        ftp.cwd(str(year))
        totpath += "/%d" % (year)
    elif "%02d%03d" % (yy, dayofyear) in myl:
        ftp.cwd("%02d%03d" % (yy, dayofyear))
        totpath += "/%02d%03d" % (yy, dayofyear)
    myl = []
    ftp.retrlines("NLST", myl.append)
    if str(dayofyear) in myl:
        ftp.cwd(str(dayofyear))
        totpath += "/%d" % (dayofyear)
    logging.info("Retrieving data from %s", totpath)
    myl = []
    ftp.retrlines("NLST", myl.append)
    filenames = [i for i in myl if (prefix.lower() in i.lower()) and 
                 ("%03d"%dayofyear in i.lower()) and
                 (i.lower().endswith("i.z") or i.lower().endswith("i"))]
    logging.info(" ".join(filenames))
    assert len(filenames) > 0, "No files found on %s for %s" % (server,
                                                                prefix)
    if prefix.lower() == "robr" and len(filenames) > 1:
        filenames = sorted(filenames)
        filenames = _store_files(ftp, filenames, outpath, overwrite)
        # get data for next day
        nextday = mydate + datetime.timedelta(days=1)
        nyear = nextday.year
        ndayofyear = nextday.timetuple().tm_yday
        ftp.cwd("/" + ftppath)
        myl = []
        ftp.retrlines("NLST", myl.append)
        if str(nyear) in myl:
            ftp.cwd(str(nyear))
        myl = ftp.retrlines("NLST")
        if str(ndayofyear) in myl:
            ftp.cwd(str(ndayofyear))
        myl = ftp.retrlines("NLST")
        nfilenames = [i for i in myl if (prefix.lower() in i.lower()) and
                      (i.lower().endswith("i.z")) and "A00" in i.upper()]
        nfilenames = _store_files(ftp, nfilenames, outpath, overwrite)
        filenames += nfilenames
        _combine_ionex(outpath, filenames,
                       prefix + "%03d0.%sI" % (dayofyear, yy))
        return os.path.join(outpath, prefix + "%03d0.%sI" % (dayofyear, yy))
    else:
        nfilenames = _store_files(ftp, filenames, outpath, overwrite)
        return nfilenames[0]


def getIONEXfile(time="2012/03/23/02:20:10.01",
                 server="ftp://cddis.gsfc.nasa.gov/gnss/productsionex/",
                 prefix="codg",
                 outpath='./',
                 overwrite=False):
    getIONEXfile.__doc__ = _get_IONEX_file.__doc__
    return _get_IONEX_file(time, server, prefix, outpath, overwrite)

def get_TEC_data(times, lonlatpp, server, prefix, outpath, use_filter=None,earth_rot=0.):
    '''Returns vtec for given times and lonlats.
    If times has the same length as the first axis  of lonlatpp, 
    it is assumed that there is a one to one correspondence.
    Else vTEC is calculated for every combination of lonlatpp and times.
    
    Args:
        times (np.array) : float of time in MJD seconds
        lonlatpp (np.array) : array of time X 2 ,longitude lattitude of 
                              piercepoints
        server (string) : ftp server to get IONEX data from
        prefix (string) : prefix of IONEX data 
        outpath (string) : local location of the IONEX files
    Returns:
        np.array : array with shape times.shape x lonlatpp.shape[0], unless both
                   are equal
    '''
    
    date_parms = PosTools.obtain_observation_year_month_day_fraction(times[0])
    ionexf=get_IONEX_file(time=date_parms,server=server,prefix=prefix,outpath=outpath)
    tecinfo=ionex.readTEC(ionexf,use_filter=use_filter)
    latpp = lonlatpp[:, 1]
    lonpp = lonlatpp[:, 0]
    if latpp.shape == times.shape:
        vtec = compute_tec_interpol(times,lat=latpp,lon=lonpp,tecinfo=tecinfo,apply_earth_rotation=earth_rot)
    else:
        vtec=[]
        for itime in range(times.shape[0]):
            vtec.append(compute_tec_interpol(times[itime]*np.ones_like(latpp),
                                             lat=latpp,
                                             lon=lonpp,
                                             tecinfo=tecinfo,
                                             apply_earth_rotation=earth_rot))
    return np.array(vtec)
            
