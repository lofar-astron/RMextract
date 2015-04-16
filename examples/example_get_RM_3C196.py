import RMextract.getRM as gt
import pyrap.tables as tab
from pylab import *
a=tab.taql('calc MJD("2014/03/01/00:00:00")')[0]*3600*24
b=tab.taql('calc MJD("2014/03/20/00:00:00")')[0]*3600*24

statpos=gt.PosTools.posCS002
pointing=array([ 2.15374123,  0.8415521 ]) #3C196
tec=gt.getRM(ionexPath='./IONEXdata/',earth_rot=1,ha_limit=0.25*np.pi,radec=pointing,timestep=1800, timerange = [a, b],stat_positions=[statpos,])
times=tec['times']
flags=tec['flags']['st1']
timess=[tm/(3600*24.) for tm in times]
dates=tab.taql('calc ctod($timess)')
dates=dates[dates.keys()[0]]

format="%Y/%m/%d/%H:%M:%S.%f"
mydatetimes=[datetime.datetime.strptime(mydate,format) for mydate in dates]
maskeddata=np.ma.array(tec['RM']['st1'],mask=np.logical_not(flags))
plot_date(mydatetimes,maskeddata,'-')
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
title("RM variation winter 2012/2013")
show()
