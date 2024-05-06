#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract TEC values from an IONEX file given a specific time and geographic coordinate.

Created on Tue Apr 24 11:46:57 2018

@author: mevius
"""

import datetime
import ftplib
import os
import socket
from pathlib import Path
from typing import List, Optional, TypeVar, Union
from urllib import request
from urllib.parse import urlparse

import numpy as np
import scipy.ndimage.filters as myfilter
import socks

from RMextract import PosTools
from RMextract.formatters import KNOWN_FORMATTERS, Formatter
from RMextract.logging import logger

PathLike = TypeVar("PathLike", str, Path)


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
                *(
                    int(float(i))
                    for i in stripped.replace("EPOCH OF FIRST MAP", "").split()
                )
            )
        if stripped.endswith("EPOCH OF LAST MAP"):
            endtime = datetime.datetime(
                *(
                    int(float(i))
                    for i in stripped.replace("EPOCH OF LAST MAP", "").split()
                )
            )
        if stripped.endswith("INTERVAL"):
            timestep = float(stripped.split()[0]) / 3600.0
        if stripped.endswith("EXPONENT"):
            exponent = pow(10, float(stripped.split()[0]))
        if stripped.endswith("DLON"):
            start_lon, end_lon, step_lon = (float(i) for i in stripped.split()[:3])
        if stripped.endswith("DLAT"):
            start_lat, end_lat, step_lat = (float(i) for i in stripped.split()[:3])
        if stripped.endswith("OF MAPS IN FILE"):
            ntimes = int(stripped.split()[0])

    lonarray = np.arange(start_lon, end_lon + step_lon, step_lon)
    latarray = np.arange(start_lat, end_lat + step_lat, step_lat)
    dtime = endtime - starttime
    dtimef = dtime.days * 24.0 + dtime.seconds / 3600.0
    logger.debug("timerange %f hours. step = %f ", dtimef, timestep)
    timearray = np.arange(0, dtimef + timestep, timestep)
    if timearray.shape[0] < ntimes:
        # bug in ILTF files,last time in header is incorrect
        extratimes = np.arange(
            timearray[-1] + timestep,
            timearray[-1] + (ntimes - timearray.shape[0] + 0.5) * timestep,
            timestep,
        )
        timearray = np.concatenate((timearray, extratimes))
    timearray += starttime.hour + starttime.minute / 60.0 + starttime.second / 3600.0

    return exponent, lonarray, latarray, timearray


def read_tec(filename, _use_filter=None):
    """returns TEC, RMS longitude, lattitude and time read from an IONEX file.

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
    logger.info(
        "reading data with shapes %d  x %d x %d",
        timearray.shape[0],
        latarray.shape[0],
        lonarray.shape[0],
    )
    tecarray = np.zeros(timearray.shape + latarray.shape + lonarray.shape, dtype=float)
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
            latstr = line.strip().replace("LAT/LON1/LON2/DLON/H", "")
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
            data = np.fromstring(" -".join(line.strip().split("-")), sep=" ") * exponent
            if tecdata:
                tecarray[timeidx, latidx, lonidx : lonidx + data.shape[0]] = data
            elif rmsdata:
                rmsarray[timeidx, latidx, lonidx : lonidx + data.shape[0]] = data
            lonidx += data.shape[0]
    if _use_filter is not None:
        tecarray = myfilter.gaussian_filter(tecarray, _use_filter, mode="nearest")
    return tecarray, rmsarray, lonarray, latarray, timearray


def readTEC(filename, use_filter=None):
    """oldfunction name for compatibility. Use read_tec."""

    logger.warning("function readTEC obsolete, use read_tec instead")
    return read_tec(filename, _use_filter=use_filter)


def _compute_index_and_weights(maparray, mapvalues):
    """helper function  to get indices and weights for interpolating tecmaps


    Args:
        maparray (np.array) : array to get indices in
        mapvalues (Union[float,np.array]) :  values to get indices for
    Returns:
        Tuple[np.array, np.array, np.array]: idx1,idx2 and weights for idx2,
                                             idx2 is always >= idx1

    """
    is_reverse = maparray[1] < maparray[0]
    idx1 = np.argmin(
        np.absolute(maparray[np.newaxis] - mapvalues[:, np.newaxis]), axis=1
    )
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
    weights[_steps == 0] = 1.0
    weights[_steps != 0] = weights[_steps != 0] / _steps[_steps != 0]
    return idx1, idx2, weights


def compute_tec_interpol(times, lats, lons, tecinfo, apply_earth_rotation=0):
    """Get interpolated TEC for array of times/lats and lons

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
    """
    assert times.shape == lats.shape, "times and lats should be array with same shape"
    assert times.shape == lons.shape, "times and lons should be array with same shape"
    tecdata = tecinfo[0]  # TEC in TECU
    lonarray = tecinfo[2]  # longitude in degrees from West to East (- to +)
    latarray = tecinfo[3]  # lattitude in degrees
    # times in hour of the day (eg. 23.5 for half past 11PM)
    maptimes = tecinfo[4]

    # get indices of nearest 2 time frames + inverse distance weights
    # assume time is  sorted from early to late
    timeidx1, timeidx2, time_weights = _compute_index_and_weights(maptimes, times)

    # latarray is sorted small to large
    latidx1, latidx2, lat_weights = _compute_index_and_weights(latarray, lats)

    # for getting lon idx take into account earth_rotation
    # if longitudes cover all range between -180 and 180 you can modulate the
    # indices, otherwise we have to take the edge of the map.
    lonstep = lonarray[1] - lonarray[0]
    full_circle = np.remainder(lonarray[0] - lonarray[-1], 360.0) <= 1.1 * lonstep

    rot1 = ((times - maptimes[timeidx1]) * 360.0 / 24.0) * apply_earth_rotation
    rot2 = ((times - maptimes[timeidx2]) * 360.0 / 24.0) * apply_earth_rotation

    if not full_circle:
        lonidx11, lonidx12, lon_weights1 = _compute_index_and_weights(
            lonarray, lons + rot1
        )
        lonidx21, lonidx22, lon_weights2 = _compute_index_and_weights(
            lonarray, lons + rot2
        )
    else:
        lonidx11 = np.argmin(
            np.absolute(
                np.remainder(
                    lonarray[np.newaxis]
                    - rot1[:, np.newaxis]
                    - lons[:, np.newaxis]
                    + 180.0,
                    360.0,
                )
                - 180.0
            ),
            axis=1,
        )
        lonidx12 = lonidx11.copy()
        lonidx11[
            np.remainder(lonarray[lonidx11] - rot1 - lons + 180.0, 360.0) - 180.0 > 0
        ] -= 1
        lonidx12[
            np.remainder(lonarray[lonidx12] - rot1 - lons + 180.0, 360.0) - 180.0 < 0
        ] += 1
        lonidx11[lonidx11 < 0] += lonarray.shape[0]
        lonidx12[lonidx12 < 0] += lonarray.shape[0]
        lonidx11[lonidx11 >= lonarray.shape[0]] -= lonarray.shape[0]
        lonidx12[lonidx12 >= lonarray.shape[0]] -= lonarray.shape[0]
        lon_weights1 = (
            np.absolute(
                np.remainder(lonarray[lonidx11] - rot1 - lons + 180.0, 360.0) - 180.0
            )
            / lonstep
        )

        lonidx21 = np.argmin(
            np.absolute(
                np.remainder(
                    lonarray[np.newaxis]
                    - rot2[:, np.newaxis]
                    - lons[:, np.newaxis]
                    + 180.0,
                    360.0,
                )
                - 180.0
            ),
            axis=1,
        )
        lonidx22 = lonidx21.copy()
        lonidx21[
            np.remainder(lonarray[lonidx21] - rot2 - lons + 180.0, 360.0) - 180.0 > 0
        ] -= 1
        lonidx22[
            np.remainder(lonarray[lonidx22] - rot2 - lons + 180.0, 360.0) - 180.0 < 0
        ] += 1
        lonidx21[lonidx21 < 0] += lonarray.shape[0]
        lonidx22[lonidx22 < 0] += lonarray.shape[0]
        lonidx21[lonidx21 >= lonarray.shape[0]] -= lonarray.shape[0]
        lonidx22[lonidx22 >= lonarray.shape[0]] -= lonarray.shape[0]
        lon_weights2 = (
            np.absolute(
                np.remainder(lonarray[lonidx21] - rot2 - lons + 180.0, 360.0) - 180.0
            )
            / lonstep
        )
    logger.debug(
        "inidces time %d %d indices lat %d %d indices \
                  lon %d %d %d %d",
        timeidx1[0],
        timeidx2[0],
        latidx1[0],
        latidx2[0],
        lonidx11[0],
        lonidx12[0],
        lonidx21[0],
        lonidx22[0],
    )
    logger.debug(
        "weights time %f lat %f lon %f %f",
        time_weights[0],
        lat_weights[0],
        lon_weights1[0],
        lon_weights2[0],
    )
    tecs = (
        tecdata[timeidx1, latidx1, lonidx11] * (1.0 - lon_weights1)
        + tecdata[timeidx1, latidx1, lonidx12] * lon_weights1
    ) * (1.0 - time_weights)
    tecs += (
        tecdata[timeidx2, latidx1, lonidx21] * (1.0 - lon_weights2)
        + tecdata[timeidx2, latidx1, lonidx22] * lon_weights2
    ) * (time_weights)
    tecs *= 1.0 - lat_weights
    tecs += (
        lat_weights
        * (
            tecdata[timeidx1, latidx2, lonidx11] * (1.0 - lon_weights1)
            + tecdata[timeidx1, latidx2, lonidx12] * lon_weights1
        )
        * (1.0 - time_weights)
    )
    tecs += (
        lat_weights
        * (
            tecdata[timeidx2, latidx2, lonidx21] * (1.0 - lon_weights2)
            + tecdata[timeidx2, latidx2, lonidx22] * lon_weights2
        )
        * (time_weights)
    )
    return tecs


def getTECinterpol(time, lat, lon, tecinfo, apply_earth_rotation=0):
    """old function name for compatibility. Use compute_tec_interpol
    instead"""

    # logger.warning("obsolete, use compute_tec_interpol instead")
    if np.isscalar(time):
        time = [time]
        lat = [lat]
        lon = [lon]
    return compute_tec_interpol(
        np.array(time), np.array(lat), np.array(lon), tecinfo, apply_earth_rotation
    )


def _combine_ionex(outpath, filenames, newfilename, overwrite=False):
    """Helper function to combine separate IONEXfiles into 1 single file
    (needed for 15min ROBR data)"""

    if not overwrite and os.path.isfile(outpath + newfilename):
        logger.info("FILE exists: " + outpath + newfilename)
        return outpath + newfilename
    newf = open(outpath + newfilename, "w")
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
            newf.write(line)
            break
        if "EPOCH OF LAST MAP" not in line:
            if "OF MAPS IN FILE" in line:
                newf.write(line.replace("1", str(len(filenames))))
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
                        newf.write(line.replace("1", str(tecmapnr)))
                        break
                    if "START OF TEC MAP" in line:
                        newf.write(line.replace("1", str(tecmapnr)))
                    else:
                        newf.write(line)
        tecmapnr += 1
    newf.write("END OF FILE\n")
    return os.path.join(outpath, newfilename)
    # ignore RMS map for now, since it is filled with zeros anyway


def _gunzip_some_file(
    compressed_file: PathLike, uncompressed_file: PathLike, delete_file: bool = True
):
    command = "gunzip -dc %s > %s" % (compressed_file, uncompressed_file)
    retcode = os.system(command)
    if retcode:
        raise RuntimeError("Could not run '%s'" % command)
    if delete_file:
        os.remove(compressed_file)
    return uncompressed_file


def _store_files(
    ftp: ftplib.FTP, filenames: List[str], outpath: Path, overwrite=False
) -> List[str]:
    """helper function to store files from ftp server to outpath"""

    npaths: List[Path] = []
    for myf in filenames:
        # make sure filename is always stored uppercase
        mypath = outpath / myf.upper()
        if not overwrite and mypath.exists():
            npaths.append(mypath)
            logger.info("file %s exists", mypath)
        elif not overwrite and mypath.with_suffix("").exists():
            npaths.append(mypath.with_suffix(""))
            logger.info("file %s exists", mypath.with_suffix(""))
        else:
            with open(mypath, "wb") as myp:
                ftp.retrbinary(f"RETR {myf}", myp.write)
            npaths.append(mypath)

    nfilenames: List[str] = []
    for myf in npaths:
        if myf.suffix == ".Z":
            nfilenames.append(_gunzip_some_file(myf, myf.with_suffix("")).as_posix())
        else:
            nfilenames.append(myf.as_posix())
    return nfilenames


def _get_IONEX_file(
    time="2012/03/23/02:20:10.01",
    server="ftp://gssc.esa.int/gnss/products/ionex/",
    prefix="UQRG",
    outpath=Path("./"),
    overwrite=False,
    backupserver="ftp://ftp.aiub.unibe.ch/CODE/",
) -> str:
    """Get IONEX file with prefix from server for a given day

    Downloads files with given prefix from the ftp server, unzips and stores
    the data. For prefix ROBR the data is stored on the server in a separate
    IONEX file every 15 minutes, these are automatically combined for
    compatibility.

    Args:
        time (string or list) : date of the observation
        server (string) : ftp server + path to the ionex directories
        prefix (string) : prefix of the IONEX files (case insensitive)
        outpath (Path) : path where the data is stored
        overwrite (bool) : Do (not) overwrite existing data
    """
    prefix = prefix.upper()

    if not isinstance(outpath, Path):
        outpath = Path(outpath)  # for backward compatibility
    if not outpath.exists():
        try:
            outpath.mkdir(parents=True)
        except Exception as e:
            logger.error(f"cannot create output dir for IONEXdata: {outpath}")
            raise e

    try:
        yy = time[2:4]
        year = int(time[:4])
        month = int(time[5:7])
        day = int(time[8:10])
    except (IndexError, TypeError):
        year = time[0]
        yy = year - 2000
        month = time[1]
        day = time[2]

    mydate = datetime.date(year, month, day)
    dayofyear = mydate.timetuple().tm_yday
    # If file exists just return filename
    for _test_path in (
        outpath / f"{prefix}{dayofyear:03d}0.{yy:02d}I",
        outpath / f"IGRG{dayofyear:03d}0.{yy:02d}I",  # IGRG (fast files) (UGLY!!)
    ):
        if not overwrite and _test_path.exists():
            logger.info(f"FILE exists: {_test_path}")
            return _test_path.as_posix()

    tried_backup = False
    serverfound = False
    while not serverfound:
        url = urlparse(server)
        ftpserver = url.netloc
        ftppath = url.path
        nr_tries = 0
        try_again = True
        while try_again and nr_tries < 10:
            try:
                ftp = ftplib.FTP(ftpserver)
                ftp.login()
                try_again = False
                serverfound = True
            except ftplib.error_perm:
                if "213.184.6.172" in server:
                    ftp.login("data-out", "Qz8803#mhR4z")
                    try_again = False
                    serverfound = True
                else:
                    try_again = True
                    nr_tries += 1
                    if nr_tries >= 10:
                        if tried_backup or server == backupserver:
                            raise Exception("Could not connect to %s" % ftpserver)
                        else:
                            server = backupserver
                            tried_backup = True
                            logger.warning(
                                f"Primary IONEX host '{server}' resolution "
                                f"failure. Trying backup at '{backupserver}'"
                            )
            except socket.gaierror:
                try_again = True
                nr_tries += 1
                if nr_tries >= 10:
                    if tried_backup or server == backupserver:
                        raise Exception("Could not connect to %s" % ftpserver)
                    else:
                        server = backupserver
                        tried_backup = True
                        logger.warning(
                            f"Primary IONEX host '{server}' resolution "
                            f"failure. Trying backup at '{backupserver}'"
                        )
    if serverfound:
        logger.info(f"Successfully contacted IONEX host '{server}'")
    ftp.cwd(ftppath)
    totpath = ftppath
    myl = []
    ftp.retrlines("NLST", myl.append)
    logger.info("Retrieving data for %d or %02d%03d", year, yy, dayofyear)
    if str(year) in myl:
        ftp.cwd(str(year))
        totpath += "/%d" % (year)
    elif "%02d%03d" % (yy, dayofyear) in myl:
        ftp.cwd("%02d%03d" % (yy, dayofyear))
        totpath += "/%02d%03d" % (yy, dayofyear)
    myl = []
    ftp.retrlines("NLST", myl.append)
    if "%03d" % dayofyear in myl:
        ftp.cwd("%03d" % dayofyear)
        totpath += "/%03d" % (dayofyear)
    logger.info("Retrieving data from %s", totpath)
    myl = []
    ftp.retrlines("NLST", myl.append)
    filenames = [
        i
        for i in myl
        if (prefix.lower() in i.lower())
        and ("%03d" % dayofyear in i.lower())
        and (i.lower().endswith("i.z") or i.lower().endswith("i"))
    ]
    logger.info(" ".join(filenames))
    # assert len(filenames) > 0, "No files found on %s for %s" % (server,prefix)
    if len(filenames) <= 0:
        raise FileNotFoundError(f"No files found on {server} for {prefix}")

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
        nfilenames = [
            i
            for i in myl
            if (prefix.lower() in i.lower())
            and (i.lower().endswith("i.z"))
            and "A00" in i.upper()
        ]
        nfilenames = _store_files(ftp, nfilenames, outpath, overwrite)
        filenames += nfilenames
        _combine_ionex(
            outpath,
            filenames,
            prefix + "%03d0.%sI" % (dayofyear, yy),
            overwrite=overwrite,
        )
        ftp.quit()
        return os.path.join(outpath, prefix + "%03d0.%sI" % (dayofyear, yy))
    else:
        nfilenames = _store_files(ftp, filenames, outpath, overwrite)
        ftp.quit()
        return nfilenames[0]


def get_urllib_IONEXfile(
    time="2012/03/23/02:20:10.01",
    server="http://ftp.aiub.unibe.ch/CODE/",
    prefix="codg",
    outpath=Path("./"),
    overwrite=False,
    backupserver="http://ftp.aiub.unibe.ch/CODE/",
    formatter: Optional[Union[Formatter, str]] = None,
    proxy_server=None,
    proxy_type=None,
    proxy_port=None,
    proxy_user=None,
    proxy_pass=None,
) -> str:
    """Get IONEX file with prefix from server for a given day

    Downloads files with given prefix from the ftp server, unzips and stores
    the data. Uses urllib2 instead of ftplib to have the option to use a ftp proxy server.

    Proxy args are optional.

    Args:
        time (string or list) : date of the observation
        server (string) : ftp server + path to the ionex directories
        prefix (string) : prefix of the IONEX files (case insensitive)
        outpath (string) : path where the data is stored
        overwrite (bool) : Do (not) overwrite existing data
        formatter (Optional, Formatter | str):
                    If a string is given, it will be used as as an index in KNOWN_FORMATTERS
                    If a Formatter is given, it will be used to construct the filenames.
                    Must have the following signature:
                        formatter(server,prefix,year,dayofyear) -> str
                    If not given, the function will try to guess the formatter based on the server
        proxy_server (string): address of proxyserver, either url or ip address
        proxy_type (string): socks4 or socks5
        proxy_port (int): port of proxy server
        proxy_user (string): username for proxyserver
        proxy_pass (string): password for proxyserver
    """
    prefix = prefix.upper()
    if not isinstance(outpath, Path):
        outpath = Path(outpath)  # for backward compatibility
    if not outpath.exists():
        try:
            outpath.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.error(f"cannot create output dir for IONEXdata: {outpath}")
            raise e

    try:
        yy = int(time[2:4])
        year = int(time[:4])
        month = int(time[5:7])
        day = int(time[8:10])
    except (IndexError, TypeError):
        year = time[0]
        yy = year - 2000
        month = time[1]
        day = time[2]

    mydate = datetime.date(year, month, day)
    dayofyear = mydate.timetuple().tm_yday
    # If file exists just return filename
    for _test_path in (
        outpath / f"{prefix}{dayofyear:03d}0.{yy:02d}I",
        outpath / f"IGRG{dayofyear:03d}0.{yy:02d}I",  # IGRG (fast files) (UGLY!!)
    ):
        if not overwrite and _test_path.exists():
            logger.info(f"FILE exists: {_test_path}")
            return _test_path

    # If proxy url is given, enable proxy using pysocks
    if proxy_server and ("None" not in proxy_server):
        s = socks.socksocket()
        if proxy_type == "socks4":
            ProxyType = socks.SOCKS4
        if proxy_type == "socks5":
            ProxyType = socks.SOCKS5
        s.set_proxy(
            ProxyType,
            proxy_server,
            proxy_port,
            rdns=True,
            username=proxy_user,
            password=proxy_pass,
        )

    # Don't do connection tests for local files
    if "file://" not in server:
        # try primary url
        try:
            _ = request.urlopen(server, timeout=30)
        except Exception as e:
            logger.error(f"{e}")
            try:
                _ = request.urlopen(backupserver, timeout=30)
                server = backupserver
                logger.warning(
                    f"Primary IONEX host '{server}' resolution "
                    f"failure. Trying backup at '{backupserver}'"
                )
            except Exception as e:
                logger.error(
                    f"Primary and Backup Server not responding: {e}"
                )  # enable in lover environment

    if isinstance(formatter, str):
        try:
            formatter = KNOWN_FORMATTERS[formatter]
        except KeyError:
            raise ValueError(
                f"Unknown formatter {formatter} - please provide a callable"
            )
    if formatter is None:
        # Check known servers
        for known_server in KNOWN_FORMATTERS.keys():
            if known_server in server:
                formatter = KNOWN_FORMATTERS.get(known_server)
                break
        if formatter is None:
            raise ValueError(f"Unknown server {server} - please provide a formatter")

    url = formatter(server=server, prefix=prefix, year=year, dayofyear=dayofyear)

    logger.debug(f"Constructed {url=}.")

    # Download IONEX file, make sure it is always uppercase
    fname = outpath / Path(url).name.upper()
    out_fname = fname.with_suffix("") if fname.suffix == ".Z" else fname

    # First, if the final file already exists, simply return
    if not overwrite and out_fname.exists():
        return out_fname.as_posix()

    # Here we actually download the file. Lets make sure though that this
    # file is not downloaded by another task first
    if not fname.exists():
        logger.info(f"Downloading to {fname=}.")
        try:
            site = request.urlopen(url, timeout=30)
        except Exception as e:
            logger.error(f"No files found on {server} for {fname}")
            raise e

        with open(fname, "wb") as output:
            output.write(site.read())
        ###### gunzip files
        # Now if the fname and out_fname are different we need to extract.
        # Make sure that the out_fname does not already exist.
        if fname != out_fname and not out_fname.exists():
            _gunzip_some_file(fname, out_fname)

    return out_fname.as_posix()


def getIONEXfile(
    time="2012/03/23/02:20:10.01",
    server="ftp://ftp.aiub.unibe.ch/CODE/",
    prefix="codg",
    outpath="./",
    overwrite=False,
):
    getIONEXfile.__doc__ = _get_IONEX_file.__doc__
    return _get_IONEX_file(time, server, prefix, outpath, overwrite)


def get_TEC_data(
    times, lonlatpp, server, prefix, outpath, use_filter=None, earth_rot=0.0
):
    """Returns vtec for given times and lonlats.
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
    """

    date_parms = PosTools.obtain_observation_year_month_day_fraction(times[0])
    ionexf = getIONEXfile(
        time=date_parms, server=server, prefix=prefix, outpath=outpath
    )
    tecinfo = readTEC(ionexf, use_filter=use_filter)
    latpp = lonlatpp[:, 1]
    lonpp = lonlatpp[:, 0]
    if latpp.shape == times.shape:
        vtec = compute_tec_interpol(
            times, lat=latpp, lon=lonpp, tecinfo=tecinfo, apply_earth_rotation=earth_rot
        )
    else:
        vtec = []
        for itime in range(times.shape[0]):
            vtec.append(
                compute_tec_interpol(
                    times[itime] * np.ones_like(latpp),
                    lat=latpp,
                    lon=lonpp,
                    tecinfo=tecinfo,
                    apply_earth_rotation=earth_rot,
                )
            )
    return np.array(vtec)
