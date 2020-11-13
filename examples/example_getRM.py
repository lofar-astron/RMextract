import RMextract.getRM as gt
from astropy.time import Time

t = Time('2010-01-01T00:00:00',format='isot',scale ='utc')
starttime = t.mjd*24*3600.  # getRM still wants MJD time in seconds (casacore definition)
endtime = starttime + 3600. # one hour of data
statpos = [3826577.1095  ,461022.900196, 5064892.758] #LOFAR CS002LBA (center) , ITRF xyz in meters
pointing=[ 2.15374123,  0.8415521 ] #3C196  Ra, Dec in radians

RMdict = gt.getRM(ionexPath='./IONEXdata/', radec=pointing, timestep=100, timerange = [starttime, endtime], stat_positions=[statpos,])

times=RMdict['times']
RM = RMdict['RM']['st1']

print ("TIME(mjd)     RM (rad/m^2)")
for tm,rm in zip(tm,RM):
    print ("%5.2f        %1.3f"%(tm/(3600*24.),rm))
