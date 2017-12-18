#!/usr/bin/env python
# AGW - Adapted from Albus_RINEX_download.py

import sys
import os
try:
  import pycurl
  HAS_PYCURL = True
except:
  HAS_PYCURL = False
  try:
    import socket
    import urllib
  except:
    import urllib.request

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
    if(len(sys.argv) < 4):
        print ("Error: correct usage is %s inURL, outfilename timeout <username> <passwd>"%sys.argv[0])
        sys.exit(-2)
    if HAS_PYCURL:
      print ('using PyCurl')
      try:
        print("URL=",sys.argv[1]," File=",sys.argv[2])
        try:
          with open(sys.argv[2], 'wb') as f:
               c = pycurl.Curl()
               c.setopt(c.URL, sys.argv[1])
               c.setopt(c.WRITEDATA, f)
               if len(sys.argv)==6:
                 print ("adding username,pwd",sys.argv[4],sys.argv[5])
                 c.setopt(pycurl.USERPWD, '%s:%s'%(sys.argv[4],sys.argv[5]))

               print ('curl getting data at ',sys.argv[1])
               c.perform()
               print ('curl closing for ', sys.argv[1])
               c.close()
        except:
          print ('curl failure - ', sys.argv[1], ' probably not found')
          os.remove(sys.argv[2])
#         sys.exit(-3)
      except:
        pass

    if not HAS_PYCURL:
      print ('using urllib')
      try:
        timeout = float(sys.argv[3])
        socket.setdefaulttimeout(timeout)
        print("URL=",sys.argv[1]," File=",sys.argv[2])
        if len(sys.argv)==6:
          print ("password manager with urllib stillneeds implementation")
        try:
          urllib.urlretrieve(sys.argv[1], sys.argv[2], urlreporthook)
        except:
          urllib.request.urlretrieve(url=sys.argv[1], filename=sys.argv[2], reporthook=urlreporthook)
        urllib.urlcleanup()
      except:
        sys.exit(-3)
    sys.exit(0)

if __name__ == '__main__':
    main()



