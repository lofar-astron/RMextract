import pyrap.tables as tab
from pylab import *

import RMextract.getRM as gt

rc("lines", lw=3)
rc("font", family="serif", size=16)


# a=tab.taql('calc MJD("2013/03/02/17:02:54")')[0]*3600*24
# b=tab.taql('calc MJD("2013/03/03/01:02:50")')[0]*3600*24
a = tab.taql('calc MJD("2012/12/06/22:46:05")')[0] * 3600 * 24
b = tab.taql('calc MJD("2012/12/07/16:35:55")')[0] * 3600 * 24

statpos = gt.PosTools.posCS002
pointing = array([2.15374123, 0.8415521])  # 3C196
# pointing=array([3.7146860578925645, 0.9111636804140731  ]) #3C295
# tec=gt.getRM(ionexPath='/home/users/mevius/IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=300, timerange = [a, b],stat_positions=[statpos,],server="ftp://gnss.oma.be/gnss/products/IONEX/",prefix="ROBR")
# tec2=gt.getRM(ionexPath='/home/mevius/IONEXdata/',earth_rot=0,ha_limit=1*np.pi,radec=pointing,timestep=300, timerange = [a, b],stat_positions=[statpos,],prefix='iltg')
tec3 = gt.getRM(
    ionexPath="/home/users/mevius/IONEXdata/",
    server="ftp://igs-final.man.olsztyn.pl/pub/gps_data/GPS_IONO/cmpcmb/",
    prefix="igsg",
    earth_rot=1,
    radec=pointing,
    timestep=150,
    timerange=[a, b + 150],
    stat_positions=[
        statpos,
    ],
)
tec2 = gt.getRM(
    ionexPath="/home/users/mevius/IONEXdata/",
    server="ftp://213.184.6.172/home/gnss/products/ionex/",
    earth_rot=0.5,
    ha_limit=1 * np.pi,
    radec=pointing,
    timestep=150,
    timerange=[a, b + 150],
    stat_positions=[
        statpos,
    ],
    prefix="ILTF",
)
tec = gt.getRM(
    ionexPath="/home/users/mevius/IONEXdata/",
    server="ftp://213.184.6.172/home/gnss/products/ionex/",
    earth_rot=0.5,
    ha_limit=1 * np.pi,
    radec=pointing,
    timestep=150,
    timerange=[a, b + 150],
    stat_positions=[
        statpos,
    ],
    prefix="ILTQ",
)
# tec3=gt.getRM(ionexPath='/home/mevius/IONEXdata/',earth_rot=1,ha_limit=1*np.pi,radec=pointing,timestep=100, timerange = [a, b],stat_positions=[statpos,])
times = tec["times"]
flags = tec["flags"]["st1"]
timess = [tm / (3600 * 24.0) for tm in times]
dates = tab.taql("calc ctod($timess)")

myf = open("/home/users/mevius/RMextract/examples/dataMB2.txt", "r")
dataMB2 = [[float(i) for i in j.split()] for j in myf]
dataMB2 = array(dataMB2)


rtec2 = np.average(tec2["RM"]["st1"].reshape((-1, 4)), axis=1)
rtec = np.average(tec["RM"]["st1"].reshape((-1, 4)), axis=1)
rtec3 = np.average(tec3["RM"]["st1"].reshape((-1, 4)), axis=1)

if "array" in dates.keys():
    dates = dates["array"]
else:
    dates = dates[
        dates.keys()[0]
    ]  # backward compatibility with older casacore vresions

format = "%Y/%m/%d/%H:%M:%S.%f"
mydatetimes = [datetime.datetime.strptime(mydate, format) for mydate in dates]
maskeddata = np.ma.array(tec2["RM"]["st1"], mask=np.logical_not(flags))
plot_date(mydatetimes, maskeddata, "-")
maskeddata = np.ma.array(tec["RM"]["st1"], mask=np.logical_not(flags))
plot_date(mydatetimes, maskeddata, "-")
maskeddata = np.ma.array(tec3["RM"]["st1"], mask=np.logical_not(flags))
plot_date(mydatetimes, maskeddata, "-")
plt.gcf().autofmt_xdate()
ylabel("RM (rad/m^2)")
# title("RM variation direction 3C196 2015/07/22")
title("RM variation direction 3C196 2013/03/02")
figure(2)
plot_date(mydatetimes, tec2["RM"]["st1"] - np.average(rtec2 - dataMB2[:, 1]), "b-")
plot_date(mydatetimes, tec["RM"]["st1"] - np.average(rtec - dataMB2[:, 1]), "g-")
plot_date(mydatetimes, tec3["RM"]["st1"] - np.average(rtec3 - dataMB2[:, 1]), "r-")
plot_date(mydatetimes[2::4], dataMB2[:, 1], "ko")
plot_date(mydatetimes[2::4], dataMB2[:, 1], "k-")
legend(["ILTF", "ILTQ", "IGSG", "LOFAR"])
ylabel("RM rad/m^2")
figure(3)
# plot_date(mydatetimes[2::4],dataMB2[:,1],'-')
plot_date(mydatetimes, tec2["STEC"]["st1"] / tec2["AirMass"]["st1"], "-")
plot_date(mydatetimes, tec["STEC"]["st1"] / tec3["AirMass"]["st1"], "-")
plot_date(mydatetimes, tec3["STEC"]["st1"] / tec["AirMass"]["st1"], "-")
legend(["ILTG", "ILTQ", "IGSG"])
ylabel("vTEC (TECU)")

figure(4)
# plot_date(mydatetimes[2::4],tec2['RM']['st1'][2::4]-dataMB2[:,1]-np.average(tec2['RM']['st1'][2::4]-dataMB2[:,1]),'b-')
# plot_date(mydatetimes[2::4],tec['RM']['st1'][2::4]-dataMB2[:,1]-np.average(tec3['RM']['st1'][2::4]-dataMB2[:,1]),'g-')
# plot_date(mydatetimes[2::4],tec3['RM']['st1'][2::4]-dataMB2[:,1]-np.average(tec['RM']['st1'][2::4]-dataMB2[:,1]),'r-')
plot_date(
    mydatetimes[2::4], rtec2 - dataMB2[:, 1] - np.average(rtec2 - dataMB2[:, 1]), "b-"
)
plot_date(
    mydatetimes[2::4], rtec - dataMB2[:, 1] - np.average(rtec - dataMB2[:, 1]), "g-"
)
plot_date(
    mydatetimes[2::4], rtec3 - dataMB2[:, 1] - np.average(rtec3 - dataMB2[:, 1]), "r-"
)
legend(["ILTG", "ILTQ", "IGSG"])
ylabel("RM rad/m^2")
title("differential RM GPS-LOFAR")
show()
