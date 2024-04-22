#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Formatters for the different IONEX servers
"""
# Url of the primary server has the syntax 
# "ftp://ftp.aiub.unibe.ch/CODE/YYYY/CODGDOY0.YYI.Z" 
# where DOY is the day of the year, padded with 
# leading zero if <100, and YY is the last two digits of year.

# Url of the backup server has the syntax 
# "ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/YYYY/DOY/codgDOY.YYi.Z"
# where DOY is the day of the year, padded with 
# leading zero if <100, and YY is the last two digits of year.

from typing import Dict, Callable

Formatter = Callable[[str, int, int, str], str]

def ftp_aiub_unibe_ch(server: str, year: int, dayofyear: int, prefix: str) -> str:
    """
    Format the URL for the ftp.aiub.unibe.ch server
    """
    return f"{server}/CODE/{year:4d}/{prefix.upper()}{dayofyear:03d}0.{year%100:02d}I.Z"

def cddis_gsfc_nasa_gov(server: str, year: int, dayofyear: int, prefix: str) -> str:
    return f"{server}/gnss/products/ionex/{year:4d}/{dayofyear:03d}/{prefix.lower()}{dayofyear:03d}0.{year%100:02d}i.Z"

def igsiono_uwm_edu_pl(server: str, year: int, dayofyear: int, prefix: str) -> str:
    return f"{server}/data/ilt/{year:4d}/igrg{dayofyear:03d}0.{year%100:02d}i"



KNOWN_FORMATTERS: Dict[str,Formatter] = { 
    "ftp.aiub.unibe.ch": ftp_aiub_unibe_ch,
    "cddis.gsfc.nasa.gov": cddis_gsfc_nasa_gov,
    "igsiono.uwm.edu.pl": igsiono_uwm_edu_pl,
}