"""
Analysis of time series:
- activation time/end time
- active time duration
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from exp_meta_data import select_exp_names_and_configs, get_analysis_configs
from helper_funcs import cell_activity

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

##Loads list of experiment names and corresponding configurations for analysis
exp_names, configs = select_exp_names_and_configs(exp_type)

for exp_name in exp_names:
    print('Analysing: ', exp_name)
    results_folder = f'Experimental_data/{exp_type}/{exp_name}'
    if not os.path.exists(f'{results_folder}/results'):
        os.makedirs(f'{results_folder}/results')

    fakt = np.loadtxt(f'{results_folder}/fakt.txt')
    raw_series = np.loadtxt(f'{results_folder}/series.txt')
    time_series = np.loadtxt(f'{results_folder}/norm_time_series.txt')
    smooth_series = np.loadtxt(f'{results_folder}/norm_smooth_series.txt')
    pos = np.loadtxt(f'{results_folder}/koordinate.txt')[:, :2]*fakt
    stim_pos_x, stim_pos_y, stim_time = np.loadtxt(
        f'{results_folder}/Vbod XY.txt')
    stim_pos_x = stim_pos_x * fakt
    stim_pos_y = stim_pos_y * fakt

    dist_from_stim = np.array(
        [round(np.sqrt((pos[i, 0]-stim_pos_x)**2+(pos[i, 1]-stim_pos_y)**2), 2) for i in range(len(pos))])

    configs['STIM_TIME'] = stim_time
    act_times, response_times, fraction, amplitudes, durations, deact_times, peak_times = cell_activity(
        smooth_series, raw_series, configs)

    for i in range(len(smooth_series[0])):
        fig = plt.figure(figsize=(2.0*PANEL_WIDTH, PANEL_HEIGHT))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(range(len(smooth_series)),
                smooth_series[:, i], c='gray', linewidth=0.5)
        if not np.isnan(peak_times[i]):
            ax.scatter(peak_times[i], smooth_series[int(
                peak_times[i]), i], c='red', marker='X', s=6, label='Peak')
            ax.axvline(act_times[i], c='red', linewidth=0.4, label='Act. time')
            ax.axvline(deact_times[i], c='blue',
                       linewidth=0.4, label='Deact. time')
            ax.legend(loc='upper right')
        ax.set_xlabel('time (s)')
        ax.set_ylabel('Cell signal $i$')
        fig.savefig(f'{results_folder}/results/act_time_{i}.png',
                    dpi=600, bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)

    signal_parameters = {
        'dist_from_stim': dist_from_stim,
        'act_times': act_times,
        'fraction': fraction,
        'response_times': response_times,
        'amplitudes': amplitudes,
        'durations': durations
    }
    signal_parameters = pd.DataFrame.from_dict(signal_parameters, 'columns')
    signal_parameters.to_csv(
        f'{results_folder}/results/signal_params.txt', sep=' ', index=False)
