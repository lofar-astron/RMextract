#!/usr/bin/env python
import os
from datetime import date
from pathlib import Path
from typing import List, Literal, Optional, Union

import numpy as np

from RMextract import PosTools
from RMextract import getIONEX as ionex
from RMextract.EMM import EMM as EMM
from RMextract.formatters import Formatter
from RMextract.logging import logger

ION_HEIGHT = PosTools.ION_HEIGHT
#####################  main processing function #####################


def getRM(
    MS: Optional[str] = None,
    server="ftp://ftp.aiub.unibe.ch/CODE/",
    prefix="codg",
    ionexPath="IONEXdata/",
    earth_rot: float = 0,
    timerange: Union[List[float], Literal[0]] = 0,
    start_time: Optional[str] = None,
    end_time: Optional[str] = None,
    pointing: List[float] = [0.0, 0.5 * np.pi],
    radec: Optional[List[float]] = None,
    use_azel=False,
    ha_limit=-1000,
    use_filter: Optional[Union[float, List[float]]] = None,
    use_urlib=False,
    formatter: Optional[Formatter] = None,
    use_mean: Optional[bool] = None,
    stat_names=[],
    useEMM=False,
    object="",
    timestep=60,
    out_file="",
    stat_positions=[PosTools.posCS002],
    use_proxy=False,
    proxy_server: Optional[str] = None,
    proxy_type: Optional[str] = None,
    proxy_port: Optional[int] = None,
    proxy_user: Optional[str] = None,
    proxy_pass: Optional[str] = None,
    overwrite=False,
) -> dict:
    """Get ionspheric RM values for a given observation

    Args:
        MS (Optional[str], optional): MeasurementSet to inspect for observation. Defaults to None.
        server (str, optional): Server to obtain ionospheric data. Defaults to "ftp://ftp.aiub.unibe.ch/CODE/".
        prefix (str, optional): Prefix of the IONEX files. Defaults to 'codg'.
        ionexPath (str, optional): Where to store IONEX files. Defaults to "IONEXdata/".
        earth_rot (float, optional): specify (with a number between 0 and 1) how much of the earth rotaion is taken in to account in the interpolation step. Defaults to 0.
        timerange (int, optional): [start, end]in MJD seconds, or 0 to use start_time and end_time. Defaults to 0.
        start_time (Optional[str], optional): casa timestring,eg. 2012/11/21/12:00:00 or 56252.5d, NEEDS PYRAP! Defaults to None.
        end_time (Optional[str], optional): casa timestring,eg. 2012/11/21/12:00:00 or 56252.5d, NEEDS PYRAP! Defaults to None.
        pointing (List[float], optional): [ra,dec] in radians, or if use_azel =True, az + el in radians. Defaults to [0.,0.5*np.pi].
        radec (Optional[List[float]], optional): Overides pointing. Defaults to None.
        use_azel (bool, optional): Pointing is in az/el. Defaults to False.
        ha_limit (int, optional): _description_. Defaults to -1000.
        use_filter (Optional[Union[float, List[float]]], optional): standard deviation,or list of standard deviations (time,long, lat) to gaussian filter TEC data. Defaults to None.
        use_urlib (bool, optional): if True use urllib to download IONEX files (will also occur if 'http' is in server), otherwise use ftplib. Defaults to False.
        formatter (Optional[Callable], optional): Function to format the url to download the files
                Must have the following signature:
                    formatter(server,prefix,year,dayofyear) -> str
                 If not given, the function will try to guess the formatter based on the server. Defaults to None.
        use_mean (Optional[bool], optional): True if you only want report for mean of station positions. Defaults to None.
        stat_names (list, optional): list of strings per station. Defaults to [].
        useEMM (bool, optional): use EMM for Earth magnetic field, otherwise WMM cooefficients will be used.. Defaults to False.
        object (str, optional): Object of interest to get position. Defaults to ''.
        timestep (int, optional): Timestep to use in s. Defaults to 60.
        out_file (str, optional): if given the data points will be written to a text file. Defaults to ''.
        stat_positions (list, optional): list of length 3 numpy arrays, containing station_position in ITRF meters. Defaults to [PosTools.posCS002].
        use_proxy (bool, optional): Use a proxy server (see proxy_ args). Can only be used with urllib. Defaults to False.
        proxy_server (Optional[str], optional): Proxy server name. Defaults to None.
        proxy_type (Optional[str], optional): socks4 or socks5. Defaults to None.
        proxy_port (Optional[str], optional): port of proxy server. Defaults to None.
        proxy_user (Optional[str], optional): username for proxyserver. Defaults to None.
        proxy_pass (Optional[str], optional): password for proxyserver. Defaults to None.
        overwrite (bool, optional): if True overwrite existing IONEX files and download them again. Defaults to False.

        TIME_OFFSET (Optional[float], optional): _description_. Defaults to None.

    Raises:
        ValueError: If start_time and end_time are not given together

    Returns: The (timegrid,timestep,TEC) where TEC is a dictionary containing 1 enumpyarray per station in stat_names.
    If stat_names is not given, the station names will either be extracted from the MS or st1...stN

    """

    if MS:
        (timerange, timestep, pointing, stat_names, stat_positions) = (
            PosTools.getMSinfo(MS)
        )

    if use_mean is not None:
        stat_pos_mean = False

    if radec:
        logger.info("Using radec instead of pointing")
        pointing = radec

    if start_time and end_time:
        time_in_sec = False
    elif (start_time and not end_time) or (end_time and not start_time):
        raise ValueError("start_time and end_time must be given together")

    if not stat_names:
        stat_names = ["st%d" % (i + 1) for i in range(len(stat_positions))]

    if timerange != 0:
        start_time = timerange[0]
        end_time = timerange[1]
        time_in_sec = True
        reference_time = start_time
        timerange[0] = start_time - timestep
        timerange[1] = end_time + timestep
        str_start_time = PosTools.obtain_observation_year_month_day_hms(reference_time)
    else:
        timerange, str_start_time, reference_time = PosTools.get_time_range(
            start_time, end_time, timestep, time_in_sec, 0
        )
    if str_start_time == -1:
        return

    emm = EMM.EMM() if useEMM else EMM.WMM()

    times, timerange = PosTools.getIONEXtimerange(timerange, timestep)
    if len(times[-1]) == 0 or times[-1][-1] < timerange[1]:
        timestmp = list(times[-1])
        timestmp.append(
            timerange[1]
        )  # add one extra step to make sure you have a value for all times in the MS in case timestep hase been changed
        times[-1] = np.array(timestmp)
    timegrid = np.array([])
    TECs = {}
    Bs = {}
    Bpars = {}
    air_mass = {}
    RMs = {}
    azimuth = {}
    elevation = {}
    flags = {}
    if len(out_file) and len(object):
        out_file = out_file + "_" + object
    else:
        if len(object):
            out_file = "RMextract_report_" + object

    if len(out_file):
        if os.path.exists(out_file):
            os.remove(out_file)
        log = open(out_file, "a")
        log.write("Observing %s\n" % object)
        if MS is not None:
            log.write("Using measurement set %s\n" % MS)
        if use_azel:
            log.write("observing at a fixed azimuth and elevation \n")
        log.write("station_positions %s \n" % stat_positions)
        log.write("\n")

    for st in stat_names:
        azimuth[st] = []
        elevation[st] = []
        air_mass[st] = []
        Bs[st] = []
        Bpars[st] = []
        RMs[st] = []
        TECs[st] = []
        air_mass[st] = []
        flags[st] = []
    for time_array in times:
        # get RM per timeslot
        starttime = time_array[0]
        # print ("getting ionexfile for",starttime)
        date_parms = PosTools.obtain_observation_year_month_day_fraction(starttime)
        dayofyear = (
            date(date_parms[0], date_parms[1], date_parms[2]).timetuple().tm_yday
        )
        emm.date = date_parms[0] + float(dayofyear) / 365.0
        # get relevant ionex file
        if not use_proxy:
            if use_urlib or "http" in server:
                ionexf = ionex.get_urllib_IONEXfile(
                    time=date_parms,
                    server=server,
                    prefix=prefix,
                    outpath=Path(ionexPath),
                    overwrite=overwrite,
                    formatter=formatter,
                )
            else:  # ftp server use ftplib
                ionexf = ionex.getIONEXfile(
                    time=date_parms,
                    server=server,
                    prefix=prefix,
                    outpath=Path(ionexPath),
                    overwrite=overwrite,
                )

        else:
            ionexf = ionex.get_urllib_IONEXfile(
                time=date_parms,
                server=server,
                prefix=prefix,
                outpath=Path(ionexPath),
                proxy_server=proxy_server,
                proxy_type=proxy_type,
                proxy_port=proxy_port,
                proxy_user=proxy_user,
                proxy_pass=proxy_pass,
                overwrite=overwrite,
                formatter=formatter,
            )

        logger.info(f"Using IONEX file: {ionexf}")
        tecinfo = ionex.read_tec(ionexf, _use_filter=use_filter)
        if use_mean:
            if not stat_pos_mean:
                stn_mean = stat_positions.mean(0)
                stat_positions = []
                stat_positions.append(stn_mean)
                stat_pos_mean = True
        for station, position in zip(stat_names, stat_positions):
            for time in time_array:
                result = PosTools.obtain_observation_year_month_day_fraction(time)
                part_of_day = result[3] * 24
                if use_azel:
                    az = pointing[0]
                    el = pointing[1]
                    flags[station].append(1)
                else:
                    az, el = PosTools.getAzEl(pointing, time, position, ha_limit)
                    if az == -1 and el == -1:
                        return
                    if az == 999 and el == 999:
                        # outsite hadec range
                        air_mass[station].append(-999)
                        azimuth[station].append(-999)
                        elevation[station].append(-999)
                        Bs[station].append([-999, -999, -999])
                        Bpars[station].append(-999)
                        RMs[station].append(-999)
                        TECs[station].append(-999)  # sTEC value
                        flags[station].append(0)
                        continue
                flags[station].append(1)
                latpp, lonpp, height, lon, lat, am1 = PosTools.getlonlatheight(
                    az, el, position
                )

                if latpp == -1 and lonpp == -1 and height == -1:
                    return
                # get VTEC from IONEX interpolation
                vTEC = ionex.getTECinterpol(
                    time=part_of_day,
                    lat=latpp,
                    lon=lonpp,
                    tecinfo=tecinfo,
                    apply_earth_rotation=earth_rot,
                )

                TECs[station].append(vTEC * am1)  # STEC value
                emm.lon = lonpp
                emm.lat = latpp
                emm.h = ION_HEIGHT / 1.0e3
                Bpar = -1 * emm.getProjectedField(
                    lon, lat
                )  # minus sign since the radiation is towards the Earth
                BField = emm.getXYZ()
                Bs[station].append(BField)

                air_mass[station].append(am1)
                azimuth[station].append(az)
                elevation[station].append(el)
                Bpars[station].append(Bpar)
                # calculate RM constant comes from VtEC in TECU,B in nT, RM = 2.62
                RMs[station].append((Bpar * vTEC * am1) * 2.62e-6)

        if use_mean:
            break

        timegrid = np.concatenate((timegrid, time_array))

    for st in stat_names:
        air_mass[station] = np.array(air_mass[st])
        azimuth[station] = np.array(azimuth[st])
        elevation[station] = np.array(elevation[st])
        TECs[st] = np.array(TECs[st])
        Bs[st] = np.array(Bs[st])
        Bpars[st] = np.array(Bpars[st])
        RMs[st] = np.array(RMs[st])
        flags[st] = np.array(flags[st])
        if use_mean:
            break

    big_dict = {}
    big_dict["STEC"] = TECs
    big_dict["Bpar"] = Bpars
    big_dict["BField"] = Bs
    big_dict["AirMass"] = air_mass
    big_dict["elev"] = elevation
    big_dict["azimuth"] = azimuth
    big_dict["RM"] = RMs
    big_dict["times"] = timegrid
    big_dict["timestep"] = timestep
    big_dict["station_names"] = stat_names
    big_dict["stat_positions"] = stat_positions
    big_dict["flags"] = flags
    big_dict["reference_time"] = reference_time
    # finish writing computed data to report
    if len(out_file):
        time_range = [big_dict["times"][0], big_dict["times"][-1]]
        log.write("start and end times %s %s \n" % (time_range[0], time_range[1]))
        log.write(
            "reference time for rel_time=0: year month day hr min sec %s %s %s %s %s %s \n"
            % str_start_time
        )
        if use_azel:
            log.write(
                "observing at azimuth and elevation %s %s \n"
                % (pointing[0], pointing[1])
            )
        else:
            log.write("observation direction %s %s \n" % (pointing[0], pointing[1]))
        log.write("\n")
        k = 0
        for key in big_dict["station_names"]:
            seq_no = 0
            if use_mean:
                log.write(
                    "data for station mean position at %s\n" % (stat_positions[k])
                )
            else:
                log.write(
                    "data for station %s  at position %s\n" % (key, stat_positions[k])
                )
            log.write(
                "seq  rel_time time_width El         Az         STEC           RM (rad/m2)   VTEC factor  \n"
            )
            for i in range(timegrid.shape[0]):
                el = big_dict["elev"][key][i]
                if el < 0:
                    ok = 1
                    stec = 0.0
                    rm = 0.0
                    vtec_factor = 1.0
                else:
                    ok = 0
                    stec = big_dict["STEC"][key][i]
                    rm = big_dict["RM"][key][i]
                    vtec_factor = 1.0 / big_dict["AirMass"][key][i]
                az = big_dict["azimuth"][key][i]
                rel_time = timegrid[i] - reference_time
                if i == 0:
                    time_width = reference_time - timegrid[i]
                else:
                    time_width = timegrid[i] - timegrid[i - 1]
                log.write(
                    "%s : %s %s %s %s %s %s %s %s\n"
                    % (seq_no, ok, rel_time, time_width, el, az, stec, rm, vtec_factor)
                )
                seq_no = seq_no + 1
            k = k + 1
            if use_mean:
                break
            log.write(" \n")
        log.close()
        logger.info(f"****** finished ionosphere predictions report: {out_file}")
    else:
        logger.info("*********** finished ionosphere predictions ***************")

    return big_dict
