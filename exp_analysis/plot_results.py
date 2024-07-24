"""
Plots all results
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score 
from exp_meta_data import select_exp_names_and_configs

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

def lin_fit_func(x, k, n):
    """
    Linear fit function
    """
    return k*x+n

exp_type = input("Select experiment name [apiraza, cbx, control]: ")
analysis_type = input("Select analysis type [auto, manual]: ")
MAX_DIST = float(input("Input max. distance for binning results [example: 120.0]: "))
NUM_OF_BINS = int(input("Input number of binns for results [example: 5]: "))
MIN_DIST = 0.0
DRAW_BINS_LIMIT =  input("Input binn cut off number [example: 3 (<number of binns)]: ")
if DRAW_BINS_LIMIT in ['0', '', None, False, ' ', '\n']:
    DRAW_BINS_LIMIT = 0
else:
    DRAW_BINS_LIMIT = int(DRAW_BINS_LIMIT)

dist_from_stim_all = [[] for _ in range(NUM_OF_BINS)]
act_times_all = [[] for _ in range(NUM_OF_BINS)]
fraction_all = [[] for _ in range(NUM_OF_BINS)]
fraction_all_inc_nonactive = [[] for _ in range(NUM_OF_BINS)]
fraction_all_diffs = [[] for _ in range(NUM_OF_BINS)]
response_times_all = [[] for _ in range(NUM_OF_BINS)]
amplitudes_all = [[] for _ in range(NUM_OF_BINS)]
durations_all = [[] for _ in range(NUM_OF_BINS)]
time_to_peak_all = [[] for _ in range(NUM_OF_BINS)]

exp_names, _ = select_exp_names_and_configs(exp_type)

digitize_bins = np.linspace(MIN_DIST, MAX_DIST, NUM_OF_BINS)
for exp_name in exp_names:
    print(exp_name)
    results_folder = f'Experimental_data/{exp_type}/{exp_name}'
    if analysis_type == 'auto':
        data = pd.read_csv(f'{results_folder}/results/signal_params.txt', delimiter=' ')
        dist_from_stim = data['dist_from_stim'].to_numpy(copy=True)
        act_times = data['act_times'].to_numpy(copy=True)
        fraction = data['fraction'].to_numpy(copy=True)
        response_times = data['response_times'].to_numpy(copy=True)
        amplitudes = data['amplitudes'].to_numpy(copy=True)
        durations = data['durations'].to_numpy(copy=True)
    else:
        fakt = np.loadtxt(f'{results_folder}/fakt.txt')
        signals = np.loadtxt(f'{results_folder}/series.txt')
        pos = np.loadtxt(f'{results_folder}/koordinate.txt')[:,:2]*fakt
        stim_pos_x, stim_pos_y, stim_time = np.loadtxt(f'{results_folder}/Vbod XY.txt')
        stim_pos_x = stim_pos_x * fakt
        stim_pos_y = stim_pos_y * fakt

        act_times = np.loadtxt(f'{results_folder}/act_times.txt')
        deact_times = np.loadtxt(f'{results_folder}/deact_times.txt')
        valid_idxs = np.where((~np.isnan(act_times)) & (~np.isnan(deact_times) & (deact_times-act_times>0)))[0]
        invalid_idxs = np.where((np.isnan(act_times)) & (np.isnan(deact_times) | (deact_times-act_times==0)))[0]
        fraction = np.zeros(len(pos), float)
        fraction[:] = np.nan
        fraction[valid_idxs] = 1
        amplitudes = np.zeros(len(pos), float)
        amplitudes[:] = np.nan
        response_times = np.zeros(len(pos), float)
        response_times[:] = np.nan

        durations = deact_times - act_times
        amplitudes[valid_idxs] = np.array([np.amax(signals[int(act_times[cell]):int(deact_times[cell]),cell])-signals[int(act_times[cell]),cell] for cell in valid_idxs])
        response_times[valid_idxs] = act_times[valid_idxs] - np.nanmin(act_times)
        dist_from_stim = np.array([np.sqrt((stim_pos_x-pos[i,0])**2+(stim_pos_y-pos[i,1])**2) for i in range(len(pos))])

        time_to_peak = np.zeros(len(pos), float)
        time_to_peak[:] = np.nan
        time_to_peak[valid_idxs] = np.array([np.argmax(signals[int(act_times[cell]):int(deact_times[cell]),cell]) for cell in valid_idxs])
    
    
    #digitize_bins = np.linspace(MIN_DIST, 1.05*np.amax(dist_from_stim), NUM_OF_BINS+1)
    digitized = np.digitize(dist_from_stim, bins=digitize_bins)
    unique_bins = np.unique(digitized)
    #print(unique_bins)
    for bin_num in unique_bins:
        idxs = np.where(digitized==bin_num)[0]
        len_idxs = len(idxs)
        dist_from_stim_all[bin_num-1].append(np.nanmean(dist_from_stim[idxs]))
        if not np.isnan(np.nanmean(act_times[idxs])):
            act_times_all[bin_num-1].append(np.nanmean(act_times[idxs]))
        fraction_all[bin_num-1].append(np.nansum(fraction[idxs]))
        fraction_all_inc_nonactive[bin_num-1].append(len_idxs)
        fraction_all_diffs[bin_num-1].append((fraction_all_inc_nonactive[bin_num-1][-1]-fraction_all[bin_num-1][-1])/fraction_all_inc_nonactive[bin_num-1][-1])
        if not np.isnan(np.nanmean(response_times[idxs])):
            response_times_all[bin_num-1].append(np.nanmean(response_times[idxs]))
        if not np.isnan(np.nanmean(amplitudes[idxs])):
            amplitudes_all[bin_num-1].append(np.nanmean(amplitudes[idxs]))
        if not np.isnan(np.nanmean(durations[idxs])):
            durations_all[bin_num-1].append(np.nanmean(durations[idxs]))
        if not np.isnan(np.nanmean(time_to_peak[idxs])):
            time_to_peak_all[bin_num-1].append(np.nanmean(time_to_peak[idxs]))


dist_from_stim_avgs = [np.nanmean(elem) for elem in dist_from_stim_all]
act_times_avgs = [np.nanmean(elem) for elem in act_times_all]
response_times_avgs = [np.nanmean(elem)-np.nanmean(response_times_all[0]) for elem in response_times_all]
amplitudes_avgs = [np.nanmean(elem) for elem in amplitudes_all]
durations_avgs = [np.nanmean(elem) for elem in durations_all]
fraction_avgs = [np.nansum(elem)/np.sum(all_cells) for elem, all_cells in zip(fraction_all, fraction_all_inc_nonactive)]
time_to_peak_avgs = [np.nanmean(elem) for elem in time_to_peak_all]

#print(response_times_all[0])
#sets response times of cells in first bin to 0 relative to all other bins
response_times_all[0] = [0 for _ in response_times_all[0]]

dist_from_stim_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in dist_from_stim_all]
act_times_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in act_times_all]
response_times_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in response_times_all]
amplitudes_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in amplitudes_all]
durations_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in durations_all]
fraction_stds = [np.std(elem)/np.sqrt(len(elem)) for elem in fraction_all_diffs]
time_to_peak_stds = [np.nanstd(elem)/np.sqrt(len(elem)) for elem in time_to_peak_all]

labels = [
    'Time to peak (s)',
    'Cell response times (s)',
    'Signal amplitudes',
    'Signal durations (s)',
    'Fraction of activated cells'
]

results_folder = f'Experimental_data/{exp_type}/results'
if not os.path.exists(f'{results_folder}'):
    os.makedirs(f'{results_folder}')

plot_items = [
    (time_to_peak_avgs, time_to_peak_stds, "time_to_peak"),
    (response_times_avgs, response_times_stds, "response_times"),
    (amplitudes_avgs, amplitudes_stds, "amplitudes"),
    (durations_avgs, durations_stds, "durations"),
    (fraction_avgs, fraction_stds, "fractions")
]

fig = plt.figure(figsize=(3.0*PANEL_WIDTH, 2.0*PANEL_HEIGHT))
axes = [fig.add_subplot(2,3,i+1) for i in range(5)]
for i, (item, std_errors, filename) in enumerate(plot_items):
    item = np.array(item)
    std_errors = np.array(std_errors)
    if DRAW_BINS_LIMIT:
        axes[i].errorbar(dist_from_stim_avgs[:DRAW_BINS_LIMIT], item[:DRAW_BINS_LIMIT],
                         std_errors[:DRAW_BINS_LIMIT], c='gray', linewidth=2)
    else:
        axes[i].errorbar(dist_from_stim_avgs, item, std_errors, c='gray', linewidth=2)
    np.savetxt(f'{results_folder}/results_plot_{exp_type}_{filename}.txt', [dist_from_stim_avgs[:DRAW_BINS_LIMIT],
                                                                            item[:DRAW_BINS_LIMIT],
                                                                            std_errors[:DRAW_BINS_LIMIT]])

    x_start, x_stop = 0.94*np.nanmin(dist_from_stim_avgs), 0.95*np.nanmax(dist_from_stim_avgs)
    print(x_start, x_stop)
    x_ticks = [round(value,-1) for value in np.linspace(x_start, x_stop, NUM_OF_BINS)]
    x_ticklabels = [round(value) for value in x_ticks]
    axes[i].set_xlabel('Distance from stimulation')
    axes[i].set_ylabel(labels[i])
    axes[i].set_xlim(x_start, x_stop)
    axes[i].set_xticks(x_ticks)
    axes[i].set_xticklabels(x_ticklabels)
    axes[i].set_ylim(0.95*np.nanmin(item-std_errors), 1.05*np.nanmax(item+std_errors))

plt.subplots_adjust(wspace=0.4, hspace=0.4)
fig.savefig(f'{results_folder}/results_plot_{exp_type}.png',
            dpi=600, bbox_inches='tight', pad_inches=0.01)
plt.close(fig)
