import RMextract.getRM as gt
import pyrap.tables as tab
from pylab import *
a=tab.taql('calc MJD("2013/03/02/17:02:54")')[0]*3600*24
b=tab.taql('calc MJD("2013/03/03/01:02:50")')[0]*3600*24

statpos=gt.PosTools.posCS002
pointing=array([ 2.15374123,  0.8415521 ]) #3C196
#pointing=array([3.7146860578925645, 0.9111636804140731  ]) #3C295
tec=gt.getRM(ionexPath='./IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=450, timerange = [a, b],stat_positions=[statpos,],server="ftp://gnss.oma.be/gnss/products/IONEX/",prefix="ROBR")
tec2=gt.getRM(ionexPath='./IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=450, timerange = [a, b],stat_positions=[statpos,])
times=tec['times']
flags=tec['flags']['st1']
timess=[tm/(3600*24.) for tm in times]
dates=tab.taql('calc ctod($timess)')
if 'array' in dates.keys():
    dates=dates['array']
else:
    dates=dates[dates.keys()[0]] #backward compatibility with older casacore vresions

format="%Y/%m/%d/%H:%M:%S.%f"
mydatetimes=[datetime.datetime.strptime(mydate,format) for mydate in dates]
maskeddata=np.ma.array(tec2['RM']['st1'],mask=np.logical_not(flags))
plot_date(mydatetimes,maskeddata,'-')
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
#title("RM variation direction 3C196 2015/07/22")
title("RM variation direction 3C196 2013/03/02")
show()
