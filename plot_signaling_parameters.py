"""
Plots cell signal parameters
"""
# pylint: disable=C0103

import numpy as np
import matplotlib.pyplot as plt
from plot_configurations import PANEL_WIDTH
from model.model_parameters import ModelParameters as MP

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
###########################################

def plot_signaling_parameters(data: dict[np.ndarray], pos: np.ndarray, stim_cell: int):
    """
    Plots cellular signaling parameters
    """
    distances = np.array([np.sqrt((pos[i,0]-pos[stim_cell,0])**2+(pos[i,1]-pos[stim_cell,1])**2) for i in range(len(pos))])
    response_times = data['response_times']
    durations = data['durations']
    fraction = data['fraction']
    amplitudes = data['amplitudes']

    NUM_OF_BINS = 5
    bins = np.linspace(0, 100, NUM_OF_BINS)
    digitized = np.digitize(distances, bins)

    bin_nums = np.unique(digitized)
    distance_bins = [distances[np.where(digitized==bin_num)] for bin_num in bin_nums]
    response_times_bins = [response_times[np.where(digitized==bin_num)] for bin_num in bin_nums]
    durations_bins = [durations[np.where(digitized==bin_num)] for bin_num in bin_nums]
    fractions_bins = [fraction[np.where(digitized==bin_num)] for bin_num in bin_nums]
    amplitudes_bins = [amplitudes[np.where(digitized==bin_num)] for bin_num in bin_nums]

    distances_avgs = np.array([np.nanmean(elem) for elem in distance_bins])
    response_times_avgs = np.array([np.nanmean(elem) for elem in response_times_bins])
    durations_avgs = np.array([np.nanmean(elem) for elem in durations_bins])
    fractions_avgs = np.array([np.sum(elem)/len(elem) for elem in fractions_bins])
    amplitudes_avgs = np.array([np.nanmean(elem) for elem in amplitudes_bins])
    
    distances_stds = np.array([np.std(elem)/np.sqrt(len(elem)) for elem in distance_bins])
    response_times_stds = np.array([np.std(elem)/np.sqrt(len(elem)) for elem in response_times_bins])
    durations_stds = np.array([np.std(elem)/np.sqrt(len(elem)) for elem in durations_bins])
    fractions_stds = np.array([0 for elem in fractions_bins])
    amplitudes_stds = np.array([np.std(elem)/np.sqrt(len(elem)) for elem in amplitudes_bins])
    
    xlim = [0.98*np.amin(distances_avgs), 1.02*np.amax(distances_avgs)]
    xlabel = 'Distance from stimulations (um)'
    plot_data = {
        'response_times': {
            'ylabel': 'Response time (s)',
            'xlabel': xlabel,
            'x': distances_avgs,
            'y': response_times_avgs,
            'yerr': response_times_stds,
            'ylim': [0.98*np.amin(response_times_avgs), 1.02*np.amax(response_times_avgs)],
            'xlim': xlim
        },
        'durations': {
            'ylabel': 'Signal duration (s)',
            'xlabel': xlabel,
            'x': distances_avgs,
            'y': durations_avgs,
            'yerr': durations_stds,
            'ylim': [0.98*np.amin(durations_avgs), 1.02*np.amax(durations_avgs)],
            'xlim': xlim
        },
        'fractions': {
            'ylabel': 'Fraction of activated cells',
            'xlabel': xlabel,
            'x': distances_avgs,
            'y': fractions_avgs,
            'yerr': fractions_stds,
            'ylim': [0, 1.05],
            'xlim': xlim
        },
        'amplitudes': {
            'ylabel': 'Signal amplitude',
            'xlabel': xlabel,
            'x': distances_avgs,
            'y': amplitudes_avgs,
            'yerr': amplitudes_stds,
            'ylim': [0.98*np.amin(amplitudes_avgs), 1.02*np.amax(amplitudes_avgs)],
            'xlim': xlim
        }
    }
    xticks = [round(value) for value in distances_avgs]
    results_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}'
    fig = plt.figure(figsize=(2.0*PANEL_WIDTH, 2.0*PANEL_WIDTH))
    axes = [fig.add_subplot(2,2,i+1) for i in range(len(plot_data.keys()))]
    for i, (_, item) in enumerate(plot_data.items()):
        #axes[i].plot(item['x'], item['y'], c='gray', linewidth=1)
        axes[i].errorbar(item['x'], item['y'], item['yerr'], c='gray', linewidth=1)
        axes[i].set_xlabel(item['xlabel'])
        axes[i].set_ylabel(item['ylabel'])
        axes[i].set_xlim(*item['xlim'])
        axes[i].set_ylim(*item['ylim'])
        axes[i].set_xticks(xticks)
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(f'{results_folder}/signal_parameters.png',
                dpi=600, bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

