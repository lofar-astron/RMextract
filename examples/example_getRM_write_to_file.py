import os

import numpy as np

import RMextract.PosTools as PosTools
from RMextract.getRM import getRM

OBJECT="EoR0"
START_TIME="2013/11/28 00:00:00"
END_TIME="2013/11/28 23:59:59"
END_TIME="2013/11/28 01:59:59"
OBJECT="MWA_test"
START_TIME="2014/12/01 00:00:00"
END_TIME="2014/12/01 23:59:59"

TIME_STEP = 300.0
TIME_OFFSET=120.
out_file='test_RMextract_report_' + OBJECT
use_mean=True
use_azel = True

MWA_antennas = np.array([[-2559314.23084924,5095535.90961438,-2848889.57667157],
   [-2559293.10717106,5095498.79164383,-2848974.05801863],
   [-2559156.42269442,5095418.83230373,-2849233.34162414],
   [-2559094.63804600,5095526.84526066,-2849097.63284488],
   [-2559042.54106487,5095566.88538445,-2849072.42535023],
   [-2559068.53757350,5095654.59288871,-2848892.60844473],
   [-2559161.70851932,5095607.73286033,-2848894.91011893],
   [-2559173.78330034,5095643.10464650,-2848820.20086397],
   [-2559782.26851932,5095304.54438001,-2848875.78669410],
   [-2559644.22829760,5095369.93521424,-2848885.21756417],
   [-2559507.77695003,5095490.23883646,-2848797.39833977],
   [-2559467.43177484,5095508.01973328,-2848802.30233654],
   [-2559460.59086333,5095515.74910944,-2848794.76371318],
   [-2559491.68457220,5095527.67486954,-2848745.44773601],
   [-2559603.60609646,5095563.73884050,-2848579.72258876],
   [-2559631.28428317,5095541.41922988,-2848594.98830325],
   [-2559113.92023486,5095854.59042124,-2848492.05455485],
   [-2559133.51844911,5095831.00304170,-2848517.14718873],
   [-2559018.96896708,5095793.67783611,-2848686.69023686],
   [-2558906.48396095,5095592.28259425,-2849148.93562390],
   [-2558894.77687225,5095720.00453191,-2848930.82517056],
   [-2558880.58102582,5095762.06255238,-2848868.27661380],
   [-2558503.88881043,5095891.11710898,-2848981.31195756],
   [-2558648.85477276,5096060.47633611,-2848544.49069260],
   [-2558998.73468649,5095390.06352995,-2849423.09595365],
   [-2559238.04568324,5095263.75775157,-2849432.88470164],
   [-2558856.49159020,5095257.96516587,-2849788.57821277],
   [-2558761.92575271,5095281.91134845,-2849829.99130606],
   [-2558719.21221208,5095416.28342253,-2849628.99110746],
   [-2558836.79342206,5095555.42415917,-2849277.33903756],
   [-2558850.45931999,5095586.71918979,-2849209.71070222],
   [-2558890.31919482,5095521.92810583,-2849288.42518348]])


result = getRM(use_azel=use_azel,use_mean=use_mean,object=OBJECT,start_time=START_TIME,end_time=END_TIME, timestep=TIME_STEP,stat_positions=MWA_antennas,useEMM=True,TIME_OFFSET=TIME_OFFSET)

timerange=[result['times'][0],result['times'][-1]]
timegrid=result['times']
stat_pos=result['stat_pos']
reference_time=result['reference_time'] 
str_start_time=PosTools.obtain_observation_year_month_day_hms(reference_time)
if os.path.exists(out_file):
  os.remove(out_file)
log = open(out_file, 'a')
log.write ('Observing %s\n' % OBJECT)
if use_azel: 
  log.write('observing at a fixed azimuth and elevation\n')
if use_mean:
  log.write ('station_positions %s \n' % MWA_antennas)
  log.write ('mean of station positions %s \n' % stat_pos)
else:
  log.write ('station_positions %s \n' % MWA_antennas)
log.write ('start and end times %s %s \n' % (timerange[0], timerange[1]))
log.write ('reference time for rel_time=0: year month day hr min sec %s %s %s %s %s %s \n' % str_start_time)
log.write ('\n')
k = 0
for key in result['station_names']:
    seq_no = 0
    if use_mean:
      log.write ('data for station mean position at %s\n' % (stat_pos[k]))
    else:
      log.write ('data for station %s  at position %s\n' % (key, stat_pos[k]))
    log.write ('seq  rel_time time_width El         Az         STEC           RM (rad/m2)   VTEC factor  \n')
    for i in range (timegrid.shape[0]):
       el = result['elev'][key][i]
       if el < 0 :
         ok = 1
         stec = 0.0
         rm = 0.0
         vtec_factor = 1.0
       else:
         ok = 0
         stec =result['STEC'][key][i]
         rm = result['RM'][key][i]
         vtec_factor = 1.0 / result['AirMass'][key][i]
       az = result['azimuth'][key][i]
       rel_time = timegrid[i] - reference_time
       if i  == 0:
         time_width = reference_time - timegrid[i] 
       else:
         time_width = timegrid[i] - timegrid[i-1]
       log.write("%s : %s %s %s %s %s %s %s %s\n" % (seq_no, ok, rel_time, time_width, el, az, stec, rm, vtec_factor))
       seq_no = seq_no + 1
    k = k + 1
    try:
     if use_mean:
       break
    except: 
      pass
    log.write (' \n')
log.close()
