from RMextract.getRM import getRM
import RMextract.PosTools as PosTools
import numpy as np
import os

TIME_STEP = 300.0
TIME_OFFSET=120.0
result = getRM(MS="/home/twillis1/ASKAP_related_ionosphere/ATCA_obs.MS",use_mean=True,object= 'ATCA_dec_obs_300_mean',useEMM=True,timestep=TIME_STEP,TIME_OFFSET=TIME_OFFSET)
#result = getRM(MS="/home/twillis1/ASKAP_related_ionosphere/ATCA_obs.MS",object= 'ATCA_dec_obs_300',useEMM=True,timestep=TIME_STEP,TIME_OFFSET=TIME_OFFSET)
