import RMextract.getRM as gt
import pyrap.tables as tab
from pylab import *
a=tab.taql('calc MJD("2015/07/22/00:00:00")')[0]*3600*24
b=tab.taql('calc MJD("2015/07/22/23:59:00")')[0]*3600*24

statpos=gt.PosTools.posCS002
pointing=array([ 2.15374123,  0.8415521 ]) #3C196
tec=gt.getRM(ionexPath='./IONEXdata/',earth_rot=0,ha_limit=0.5*np.pi,radec=pointing,timestep=450, timerange = [a, b],stat_positions=[statpos,],server="ftp://gnss.oma.be/gnss/products/IONEX/",prefix="ROBR")
tec2=gt.getRM(ionexPath='./IONEXdata/',earth_rot=0,ha_limit=0.5*np.pi,radec=pointing,timestep=450, timerange = [a, b],stat_positions=[statpos,])
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
maskeddata=np.ma.array(tec['RM']['st1'],mask=np.logical_not(flags))
plot_date(mydatetimes,maskeddata,'-')
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
title("RM variation direction 3C196 2015/07/22")
show()
