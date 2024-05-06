#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tests for the formatters module
"""

from datetime import datetime, timedelta
import random
from RMextract.formatters import KNOWN_FORMATTERS, ftp_aiub_unibe_ch, cddis_gsfc_nasa_gov, igsiono_uwm_edu_pl, chapman_upc_es

def test_known_formatters():
    """
    Test that the known formatters are correct
    """
    assert KNOWN_FORMATTERS["ftp.aiub.unibe.ch"] == ftp_aiub_unibe_ch
    assert KNOWN_FORMATTERS["cddis.gsfc.nasa.gov"] == cddis_gsfc_nasa_gov
    assert KNOWN_FORMATTERS["igsiono.uwm.edu.pl"] == igsiono_uwm_edu_pl
    assert KNOWN_FORMATTERS["http://chapman.upc.es"] == chapman_upc_es

def test_ftp_aiub_unibe_ch():
    """
    Test the ftp_aiub_unibe_ch formatter
    """
    for server in ["ftp://ftp.aiub.unibe.ch", "file://ftp.aiub.unibe.ch"]:
        year = random.randint(2000, 2024)
        dayofyear = random.randint(1, 366)
        prefix = "CODG"
        url = ftp_aiub_unibe_ch(server, year, dayofyear, prefix)
        expected = f"{server}/CODE/{year:4d}/{prefix.upper()}{dayofyear:03d}0.{year%100:02d}I.Z"
        assert url == expected, f"Expected {expected}, got {url}"

def test_cddis_gsfc_nasa_gov():
    """
    Test the cddis_gsfc_nasa_gov formatter
    """
    for server in ["ftp://cddis.gsfc.nasa.gov", "file://cddis.gsfc.nasa.gov"]:
        year = random.randint(2000, 2024)
        dayofyear = random.randint(1, 366)
        prefix = "CODG"
        url = cddis_gsfc_nasa_gov(server, year, dayofyear, prefix)
        expected = f"{server}/gnss/products/ionex/{year:4d}/{dayofyear:03d}/{prefix.lower()}{dayofyear:03d}0.{year%100:02d}i.Z"
        assert url == expected, f"Expected {expected}, got {url}"

def test_igsiono_uwm_edu_pl():
    """
    Test the igsiono_uwm_edu_pl formatter
    """
    for server in ["ftp://igsiono.uwm.edu.pl", "file://igsiono.uwm.edu.pl"]:
        year = random.randint(2000, 2024)
        dayofyear = random.randint(1, 366)
        prefix = "IGRG"
        url = igsiono_uwm_edu_pl(server, year, dayofyear, prefix)
        expected = f"{server}/data/ilt/{year:4d}/{prefix.lower()}{dayofyear:03d}0.{year%100:02d}i"
        assert url == expected, f"Expected {expected}, got {url}"

def test_chapman_upc_es():
    """
    Test the chapman_upc_es formatter
    """
    for server in ["http://chapman.upc.es", "file://chapman.upc.es"]:
        year = random.randint(2000, 2024)
        dayofyear = random.randint(1, 366)
        prefix = "CODG"
        url = chapman_upc_es(server, year, dayofyear, prefix)
        dt = datetime(year,1,1) + timedelta(dayofyear-1)
        expected = f"{server}/tomion/rapid/{year:4d}/{dayofyear:03d}_{year%100:02d}{dt.month:02d}{dt.day:02d}.15min/{prefix.lower()}{dayofyear:03d}0.{year%100:02d}i.Z"
        assert url == expected, f"Expected {expected}, got {url}"