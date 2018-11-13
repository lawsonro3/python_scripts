import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from glob import glob
from importlib import reload
import sowfa_precursor

sowfa_precursor = reload(sowfa_precursor)
plt.close('all')

# Creat sowfa_precursor object and enter neccisary user inputs
infPer_001m_5m =  sowfa_precursor.Sim(
    dir='/Users/mlawson/Peregrine/windsim/wake_steering/stableABLRuns/infPer_0.001m_5m',
    log='log.3.ABLSolver',
    time_dir='25000',
    avg_time=29000,
    avg_width=2000,
    z_level=90.0)

# Enter heights to do stuff... ;)
infPer_001m_5m.input['heights'] = np.array([0,1/2,1,3/2,2,3,4,5])*infPer_001m_5m.input['windHeight']

# Make plots
infPer_001m_5m.theta_w_avg_cell()
infPer_001m_5m.Umean_avg_nonnormalized()
infPer_001m_5m.Tmean_avg_nonnormalized()
infPer_001m_5m.variances_avg_cell()
