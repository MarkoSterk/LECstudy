"""
Plot results: Ca time series
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from plot_configurations import PANEL_WIDTH, PANEL_HEIGHT
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


def plot_time_series():
        """
        Plots cellular time series
        """
        activation_time = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/act_times.txt')
        signal_duration = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/durations.txt')
        ca_time_series = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/ca_time_series.txt')
        ip3_time_series = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/ip3_time_series.txt')
        atp_time_series = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/atp_time_series.txt')
        cmat = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cmat.txt')
        pos = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cm_cells.txt')
        time = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/time.txt')

        cell_distances = distance.cdist(pos, pos, 'euclidean')
        cell_num = len(pos)
        sampling = 1.0/(MP.record_every*MP.dt)
        stimulated_cell = MP.stimulated_cell
        max_dist = np.amax(cell_distances[:, stimulated_cell])
        nearest_cell = np.where(cmat[:, stimulated_cell])[0][0]
        furthest_cell = np.where(cell_distances[:, stimulated_cell] == max_dist)[0][0]

        fig = plt.figure(figsize=(2.0*PANEL_WIDTH, 2.0*PANEL_HEIGHT))
        ax1 = fig.add_subplot(3, 1, 1)
        ax1.plot(time, ca_time_series[:, stimulated_cell],
                c='black', label='Stimulated cell')
        ax1.plot(time, ca_time_series[:, nearest_cell], c='blue', label='Nearest cell')
        ax1.plot(time, ca_time_series[:, furthest_cell], c='red', label='Furthest cell')
        ax1.legend(loc='upper right')
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Ca2+ signal (uM)')

        ax2 = fig.add_subplot(3, 1, 2)
        ax2.plot(time, ip3_time_series[:, stimulated_cell],
                c='black', label='Stimulated cell')
        ax2.plot(time, ip3_time_series[:, nearest_cell], c='blue', label='Nearest cell')
        ax2.plot(time, ip3_time_series[:, furthest_cell], c='red', label='Furthest cell')
        ax2.legend(loc='upper right')
        ax2.set_xlabel('time (s)')
        ax2.set_ylabel('IP3 signal (uM)')

        ax3 = fig.add_subplot(3, 1, 3)
        ax3.plot(time, atp_time_series[:, stimulated_cell],
                c='black', label='Stimulated cell')
        ax3.plot(time, atp_time_series[:, nearest_cell], c='blue', label='Nearest cell')
        ax3.plot(time, atp_time_series[:, furthest_cell], c='red', label='Furthest cell')
        ax3.legend(loc='upper right')
        ax3.set_xlabel('time (s)')
        ax3.set_ylabel('ATP signal (uM)')

        fig.savefig(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/stimulated_nearest_furthest_ca_ip3_atp_signal.png',
                dpi=600, bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)
        
        results_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/time_series'
        if not os.path.exists(results_folder):
                os.makedirs(results_folder)
        
        for i in range(cell_num):
                bin_sig = None
                if not np.isnan(signal_duration[i]):
                        bin_sig = np.full(len(ca_time_series), np.amin(ca_time_series[:,i]), float)
                        start_frame = int(activation_time[i]*sampling)
                        end_frame = start_frame + int((signal_duration[i]*sampling))
                        bin_sig[start_frame:end_frame+1] = 1.01*np.amax(ca_time_series[:,i])
                fig = plt.figure(figsize=(2.0*PANEL_WIDTH, PANEL_HEIGHT))
                ax = fig.add_subplot(1,1,1)
                ax.plot(time, ca_time_series[:,i], c='black')
                if not np.isnan(signal_duration[i]):
                        ax.plot(time, bin_sig, c='red', linewidth=0.5)
                ax.set_ylabel('Ca2+ signal / bin. signal')
                ax.set_xlabel('time (s)')
                fig.savefig(f'{results_folder}/bin_sig_{i}.png', dpi=300, bbox_inches='tight', pad_inches=0.01)
                plt.close(fig)
                
