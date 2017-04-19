# Albus_Coordinates.py
# Python stuff for dealing with coordinate system conversions
# 2006 Jun 22  James M Anderson  --JIVE  start
# 2007 Jan 17  JMA  --add conversion from radians to degrees, min, sec



################################################################################
# some import commands  The User should not need to change this.
################################################################################

################################################################################
# miscellaneous stuff
import copy
# import optparse, os, sys
#import re, string
#import inspect
import warnings
import math










################################################################################
# Global variables
M_DEG2RAD = math.pi/180.0
M_RAD2DEG = 180.0/math.pi

M_AS2RAD = M_DEG2RAD / 3600.0
M_RAD2AS = M_RAD2DEG * 3600.0







################################################################################
def deg_num_to_rad(deg, min, sec):
    """ convert sexagesimal degrees to radians

    This function converts a sexagesimal number sequence to an angle in
    radians.  It works for positive angles (take care of the negative part
    yourself."""
    deg_total = ((sec) / 60.0 + min) / 60.0 + deg
    rad_total = deg_total * M_DEG2RAD
    return rad_total



################################################################################
def hour_num_to_rad(hour, min, sec):
    """ convert sexagesimal hours to radians

    This function converts a sexagesimal number sequence to an angle in
    radians.  It works for positive angles (take care of the negative part
    yourself."""
    deg_total = (((sec) / 60.0 + min) / 60.0 + hour) * 15.0
    rad_total = deg_total * M_DEG2RAD
    return rad_total






################################################################################
def deg_str_to_rad(deg_str):
    """converts a string of degrees minutes seconds to a radian float

deg_str  I  string containing a degree angle in sexagesimal notation
            such as '-0 12 43.5' or '45d18\'19\"' or '+0:3:56'
    """
    deg_str = copy.copy(deg_str.strip())
    sign = +1.0
    if(deg_str[0] == '-'): sign = -1.0
    for i in xrange(len(deg_str)):
        if( (deg_str[i].isdigit()) or (deg_str[i]=='.') ):
            pass
        else:
            deg_str = deg_str.replace(deg_str[i],' ')
    deg_deg = map(float, deg_str.strip().split())
    assert(len(deg_deg) == 3)
    deg_rad = deg_num_to_rad(deg_deg[0], deg_deg[1], deg_deg[2]) * sign    
    return deg_rad





################################################################################
# convert RA and Dec strings to radian values
def radec_str_to_rad2(ra_str, dec_str):
    """convert RA and Dec strings to radian values

ra_str   I  String containing a right ascension value in sexagesimal notation
            such as '10 02 34.5507' or '02h35m06.78s' or '23:12:02'
dec_str  I  String containing a declination value in sexagesimal notation
            such as '-0 12 43.5' or '45d18\'19\"' or '+0:3:56'

OUTPUT as: ra_rad, dec_rad
ra_rad   O  ra in radians
dec_rad  O  dec in radians
            """
    ra_str = copy.copy(ra_str.strip())
    for i in xrange(len(ra_str)):
        if( (ra_str[i].isdigit()) or (ra_str[i]=='.') ):
            pass
        else:
            ra_str = ra_str.replace(ra_str[i],' ')
    ra_hours = map(float, ra_str.strip().split())
    assert(len(ra_hours) == 3)
    ra_rad = hour_num_to_rad(ra_hours[0], ra_hours[1], ra_hours[2])
    dec_str = copy.copy(dec_str.strip())
    sign = +1.0
    if(dec_str[0] == '-'): sign = -1.0
    for i in xrange(len(dec_str)):
        if( (dec_str[i].isdigit()) or (dec_str[i]=='.') ):
            pass
        else:
            dec_str = dec_str.replace(dec_str[i],' ')
    dec_deg = map(float, dec_str.strip().split())
    assert(len(dec_deg) == 3)
    dec_rad = deg_num_to_rad(dec_deg[0], dec_deg[1], dec_deg[2]) * sign    
    print 'observing in direction RA, DEC in radians ', ra_rad, dec_rad
    return ra_rad, dec_rad










################################################################################
def radec_str_to_rad(ra_str, dec_str = None):
    """convert RA and Dec strings to radian values

If dec_str is not None,

ra_str   I  String containing a right ascension value in sexagesimal notation
            such as '10 02 34.5507' or '02h35m06.78s' or '23:12:02'
dec_str  I  String containing a declination value in sexagesimal notation
            such as '-0 12 43.5' or '45d18\'19\"' or '+0:3:56'

Else, ra_str is a string containing both the right ascension and declination
information as one big string, with the sexagesimal notation above.

OUTPUT as: ra_rad, dec_rad
ra_rad   O  ra in radians
dec_rad  O  dec in radians
            """
    if(dec_str is not None): return radec_str_to_rad2(ra_str, dec_str)
    # Else, we do it this way with one big string
    # get the numbers
    str = copy.copy(ra_str.strip())
    for i in xrange(len(str)):
        if( (str[i].isdigit()) or (str[i]=='.') or (str[i]=='-') ):
            pass
        else:
            str = str.replace(str[i],' ')
    str = str.strip().split()
    sign = +1.0
    if(str[3][0] == '-'):
        sign = -1.0
        if(len(str[3]) == 1):
            str[3:] = str[4:]
        else:
            str[3] = str[3].replace('-',' ')
    nums = map(float, str)
    assert(len(nums) == 6)
    ra_rad = hour_num_to_rad(nums[0], nums[1], nums[2])
    dec_rad = deg_num_to_rad(nums[3], nums[4], nums[5]) * sign    
    return ra_rad, dec_rad





################################################################################
def rad_to_dms(rad):
    """ convert radians to sexagesimal degrees

    This function converts a positive angle in radians to a sexagesimal
    angle in degrees, minutes, and seconds.

    It works for positive angles (take care of the negative part
    yourself.

INPUTS:
rad      I  positive angle in radians

OUTPUTS: deg min sec
deg      O  degrees (integer)
min      O  minutes (integer)
sec      O  seconds (float)

    """
    d = math.fabs(rad) * M_RAD2DEG
    deg = int(d+2E-13)
    m = (d-deg) * 60.0
    min = int(m+1E-11)
    sec = (m - min) * 60.0
    return deg, min, sec















################################################################################
def hav(theta):
    """haversine mathematical function, angles in radians"""
    # I am expecting relatively small angles sometimes
    # h = 0.5 - 0.5 * math.cos(theta)
    h = math.sin(0.5 * theta)
    h = h*h
    return h

################################################################################
def ahav(theta):
    """inverse haversine mathematical function, angles in radians"""
    # I am expecting relatively small angles sometimes
    h = 2.0 * math.asin(math.sqrt(theta))
    return h


################################################################################
def angular_separation(ra_1, dec_1, ra_2, dec_2):
    """ calculate the angular separation of two sky directions, in radians

Note, all angles in radians.

This function nominally takes right ascension and declination, but would also
work with azimuth and elevation, galactic coordinates, ecliptical coords,
and so on.



Taken from _Astronomical Algorithms_, Meeus, 1991
    """
    delta_dec = dec_1 - dec_2
    delta_ra = ra_1 - ra_2
    hav_d = hav(delta_dec) + math.cos(dec_1)*math.cos(dec_2)*hav(delta_ra)
    d = ahav(hav_d)
    return d

