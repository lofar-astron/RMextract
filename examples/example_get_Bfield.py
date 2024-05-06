import pyrap.tables as tab
from pylab import *

import RMextract.getRM as gt

a = tab.taql('calc MJD("2013/11/01/00:00:00")')[0] * 3600 * 24
b = tab.taql('calc MJD("2013/11/09/00:00:00")')[0] * 3600 * 24

statpos = gt.PosTools.posCS002
pointing = array([2.15374123, 0.8415521])  # 3C196
bigdict = gt.getRM(
    ionexPath="./IONEXdata/",
    earth_rot=0,
    ha_limit=45 * np.pi / 180,
    radec=pointing,
    timestep=1800,
    timerange=[a, b],
    stat_positions=[
        statpos,
    ],
)
flags = np.logical_not(bigdict["flags"]["st1"])
times = bigdict["times"][np.logical_not(flags)]
timess = [tm / (3600 * 24.0) for tm in times]
dates = tab.taql("calc ctod($timess)")
if "array" in dates.keys():
    dates = dates["array"]
else:
    dates = dates[
        dates.keys()[0]
    ]  # backward compatibility with older casacore vresions
format = "%Y/%m/%d/%H:%M:%S.%f"
mydatetimes = [datetime.datetime.strptime(mydate, format) for mydate in dates]

bpar = bigdict["Bpar"]["st1"][np.logical_not(flags)]
bfield = bigdict["BField"]["st1"][np.logical_not(flags)]
abs_field = np.sqrt(np.sum(np.square(bfield), axis=1))
bperp = abs_field - np.abs(bpar)
plot_date(mydatetimes, bperp, "-")
plot_date(mydatetimes, bpar, "-")
plot_date(mydatetimes, bperp)
plot_date(mydatetimes, bpar)
plt.gcf().autofmt_xdate()
ylabel("Bpar (nGaus)")
show()
rm = bigdict["RM"]["st1"][np.logical_not(flags)]
plot_date(mydatetimes, rm, "-")
# plot_date(mydatetimes,rm)
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
show()
