"""
Analysis of time series:
- activation time/end time
- active time duration
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from exp_meta_data import select_exp_names_and_configs
from helper_funcs import cell_activity_manual

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

exp_type = input("Select experiment name [apiraza, cbx, control]: ")

exp_names, _ = select_exp_names_and_configs(exp_type)

for exp_name in exp_names:
    print('Analysing: ',exp_name)
    results_folder = f'Experimental_data/{exp_type}/{exp_name}'
    if not os.path.exists(f'{results_folder}/results'):
        os.makedirs(f'{results_folder}/results')

    fakt = np.loadtxt(f'{results_folder}/fakt.txt')
    raw_series = np.loadtxt(f'{results_folder}/series.txt')
    time_series = np.loadtxt(f'{results_folder}/norm_time_series.txt')
    smooth_series = np.loadtxt(f'{results_folder}/norm_smooth_series.txt')
    pos = np.loadtxt(f'{results_folder}/koordinate.txt')[:,:2]*fakt
    stim_pos_x, stim_pos_y, stim_time = np.loadtxt(f'{results_folder}/Vbod XY.txt')
    stim_pos_x = stim_pos_x * fakt
    stim_pos_y = stim_pos_y * fakt

    dist_from_stim = np.array(
        [round(np.sqrt((pos[i,0]-stim_pos_x)**2+(pos[i,1]-stim_pos_y)**2),2) for i in range(len(pos))])

    cell_data = cell_activity_manual(smooth_series)
    np.savetxt(f'{results_folder}/act_times.txt', cell_data['act_times'])
    np.savetxt(f'{results_folder}/deact_times.txt', cell_data['deact_times'])
