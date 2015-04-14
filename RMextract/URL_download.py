#!/usr/bin/env python
# Albus_RINEX_download.py
# Download files for RINEX stuff, using a timeout.  The Python timeout
# stuff for < 2.5 seems to not work or was not implemented, so
# I have made my own simple case
# 2007 Jan 19  James M Anderson  --JIVE  start


import sys
import socket
import urllib







################################################################################
def urlreporthook(block_count, block_size, file_size):
    """show user how far file download has progressed"""
    if(file_size <= 0):
        return
    percent = float(block_count * block_size) / file_size * 100.0
    sys.stdout.write("\b\b\b\b\b\b%5.1f%%"%percent)
    return




################################################################################
def main():
    if(len(sys.argv) != 4):
        print "Error: correct usage is %s inURL, outfilename timeout"%sys.argv[0]
        sys.exit(-2)
    try:
        print "getting",sys.argv[1], sys.argv[2];
        timeout = float(sys.argv[3])
        socket.setdefaulttimeout(timeout)
        urllib.urlretrieve(sys.argv[1], sys.argv[2], urlreporthook)
        urllib.urlcleanup()
    except:
        sys.exit(-3)
    sys.exit(0)

        
    

if __name__ == '__main__':
    main()
