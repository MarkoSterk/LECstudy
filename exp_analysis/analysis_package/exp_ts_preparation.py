"""
Analysis of experimental signals
"""
import os
from scipy.signal import detrend
import numpy as np
import matplotlib.pyplot as plt

from exp_meta_data import select_exp_names_and_configs
from helper_funcs import smooth_ts, normalize_ts

# Matplotlib configuration
SMALL_SIZE = 6
MEDIUM_SIZE = 6
BIGGER_SIZE = 8

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

PANEL_WIDTH = 5
PANEL_HEIGHT = 5
CONVERSION = 2.54
PANEL_WIDTH = PANEL_WIDTH/CONVERSION
PANEL_HEIGHT = PANEL_HEIGHT/CONVERSION
###########################################

# choose experiment type [apiraza, cbx, control]
exp_type = input("Select experiment name [apiraza, cbx, control]: ")

##Points for sliding-window smoothing of time series
SMOOTH_POINTS = int(input("Input number of time series smoothing points: "))
##Relative time series size for cut-off at the end.
REL_CUT_OFF = float(input("Select relative time series cut-off [0.0-1.0]: "))

exp_names, _ = select_exp_names_and_configs(exp_type)

for exp_name in exp_names:
    print('Preparing TS for: ', exp_name)
    results_folder = f'Experimental_data/{exp_type}/{exp_name}'
    if not os.path.exists(f'{results_folder}/traces'):
        os.makedirs(f'{results_folder}/traces')

    data = np.loadtxt(f'{results_folder}/series.txt')

    CUT_OFF = int(REL_CUT_OFF*len(data))
    data = data[:-CUT_OFF,:]
    detrended_data = np.zeros(data.shape, float)
    smooth_data = np.zeros(data.shape, float)
    for j in range(len(data[0])):
        detrended_data[:,j] = detrend(data[:,j])
        smooth_data[:,j] = smooth_ts(detrended_data[:,j], SMOOTH_POINTS)
        detrended_data[:,j] = normalize_ts(detrended_data[:,j])
        smooth_data[:,j] = normalize_ts(smooth_data[:,j])

    time = range(len(data))
    for j in range(len(data[0])):
        fig = plt.figure(figsize=(2.0*PANEL_WIDTH, 2.0*PANEL_HEIGHT))
        ax1 = fig.add_subplot(2,1,1)
        ax1.plot(time, data[:,j], c='black', linewidth=0.5)
        ax1.set_xticks([])
        ax1.set_ylabel('Raw signal')
        ax2 = fig.add_subplot(2,1,2)
        ax2.plot(time, detrended_data[:,j], c='gray', linewidth=0.5, label=f'Detrended TS {j}')
        ax2.plot(time, smooth_data[:,j], c='lightgray', linewidth=0.5, label=f'Smoothed TS {j}')
        ax2.legend(loc='upper right')
        ax2.set_xlabel('time (s)')
        ax2.set_ylabel('Cell signal $i$')
        fig.savefig(f'{results_folder}/traces/detrended_ts_{j}.png',
                    dpi=600, bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)

    np.savetxt(f'{results_folder}/ts_cutoff.txt', [CUT_OFF], fmt='%d')
    np.savetxt(f'{results_folder}/norm_time_series.txt', detrended_data, fmt='%.4lf')
    np.savetxt(f'{results_folder}/norm_smooth_series.txt', smooth_data, fmt='%.4lf')
