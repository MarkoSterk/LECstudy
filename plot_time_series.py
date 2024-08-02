"""
Plot results: Ca time series
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.spatial import distance
from scipy.stats import rankdata
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

def get_stimulated_cell(pos: np.ndarray):
        """
        Gets stimulated cell
        """      
        center_x = np.amin(pos[:,0]) + (np.amax(pos[:,0])-np.amin(pos[:,0]))/2.0
        center_y = np.amin(pos[:,1]) + (np.amax(pos[:,1])-np.amin(pos[:,1]))/2.0
        stimulated_cell = MP.find_stimulated_cell(pos, (center_x, center_y))  # number of cell for stimulation
        return stimulated_cell

def plot_all_ts_and_binarized_ts():
        """
        Plots all time series on same panel and 
        """
        ca_time_series = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/ca_time_series.txt')
        act_times = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/response_times.txt')
        deact_times = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/deact_times.txt')
        pos = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cm_cells.txt')
        time = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/time.txt')
        stimulated_cell: int = get_stimulated_cell(pos)
        cell_num = len(pos)

        distance_from_origin = np.array([np.sqrt((pos[i,0]-pos[stimulated_cell,0])**2+(pos[i,1]-pos[stimulated_cell,1])**2) for i in range(len(pos))])
        distance_ranks = np.argsort(distance_from_origin)

        fig = plt.figure(figsize=(2.0*PANEL_WIDTH, 4.0*PANEL_HEIGHT))
        ax1 = fig.add_subplot(2,1,1)
        ax1.set_xlabel("time (s)")
        ax1.set_ylabel("Norm. cell signal/Distance $i$")
        ax1.set_ylim(0,cell_num+1)
        ax1.set_xlim(0, 20.0)
        count = 0
        for indx in distance_ranks:
                norm_ts = (ca_time_series[:, indx]-np.min(ca_time_series[:, indx]))/(np.max(ca_time_series[:, indx])-np.min(ca_time_series[:, indx])) + count
                ax1.plot(time, norm_ts, c='lightgray', linewidth=0.5)
                count+=1
        ax1.axvline(MP.time_of_stimulation, c='black', linewidth=0.25)

        ax2 = fig.add_subplot(2,1,2)
        ax2.set_ylabel("Distance from origin/cell $i$")
        ax2.set_xlabel("time (s)")

        bin_signal = np.zeros(ca_time_series.shape, int)
        for i in range(cell_num):
                if(np.isnan(act_times[i]) or np.isnan(deact_times[i])):
                        pass
                else:
                        start = int(act_times[i]*MP.sampling)
                        end = int(deact_times[i]*MP.sampling)
                        bin_signal[start:end, i] = 1
        
        sorted_bin_sig = np.zeros(bin_signal.shape, int)

        for count, indx in enumerate(distance_ranks):
                sorted_bin_sig[:,count] = bin_signal[:,indx]

        color_list = ['white', 'red','black']
        cmap = colors.ListedColormap(color_list)
        bounds=[0,1,2]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        ax2.pcolormesh(time,np.arange(len(sorted_bin_sig[0])),
                np.transpose(sorted_bin_sig),cmap=cmap,norm=norm)

        plt.subplots_adjust(wspace=0.2, hspace=0.2)
        fig.savefig(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/all_traces_and_bin_traces.png',
                dpi=600, bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)
