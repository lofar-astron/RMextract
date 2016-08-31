import RMextract.getTEC as gt
import pyrap.tables as tab
from pylab import *
a=tab.taql('calc MJD("2012/11/01/00:00:00")')[0]*3600*24
b=tab.taql('calc MJD("2012/11/05/00:00:00")')[0]*3600*24

statpos=gt.PosTools.posCS002
ra=0.
dec=0.5*np.pi
tec=gt.getTEC(ionexPath='./IONEXdata',radec=[ra,dec],timestep=1800, timerange = [a, b],stat_positions=[statpos,],ha_limit=0.5*np.pi)

times=tec[0]
timess=[tm/(3600*24.) for tm in times]
dates=tab.taql('calc ctod($timess)')
if 'array' in dates.keys():
    dates=dates['array']
else:
    dates=dates[dates.keys()[0]] #backward compatibility with older casacore vresions
format="%Y/%m/%d/%H:%M:%S.%f"
mydatetimes=[datetime.datetime.strptime(mydate,format) for mydate in dates]
flags=tec[-1]['st1']
flaggedtec=np.ma.array(tec[2]['st1'],mask=np.logical_not (flags))
plot_date(mydatetimes,flaggedtec,'-')
plt.gcf().autofmt_xdate()
ylabel("VTEC (TECU)")
show()
