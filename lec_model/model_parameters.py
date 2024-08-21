"""
All relevant model parameters
"""
# pylint: disable=W0719
import os
import json
from enum import Enum
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from tissue_generator import load_capsule_data



class ModelParameters:
    """
    All model parameters
    """

    class ATPMode(Enum):
        """
        Enum class for ATP deffusion modes
        """
        DECOUPLED = "decoupled"
        POINT = "point"
        PARTIALLY_REGENERATIVE = "partially-regenerative"
        REGENERATIVE = "regenerative"

        def __str__(self):
            return self.value

        def __repr__(self):
            return str(self.value)

    class INITMode(Enum):
        """
        Enum class for initiator modes
        """
        POINT = "point"
        NEAREST = "nearest"

        def __str__(self):
            return self.value
        
        def __repr__(self):
            return str(self.value)
    
    class EnumEncoder(json.JSONEncoder):
        """
        JSON Encoder for Enums
        """
        def default(self, obj):
            """
            Default encoder method
            """
            if isinstance(obj, Enum):
                return obj.value
            return super().default(obj)
    
    capsule: int = 3
    t: float = 0.0 #initial time
    dt: float = 0.01 #integration step

    record_every: int = 1 ###record every 5th simulation steps
    sampling: float = 1.0/dt

    """
    Cell heterogeneity can be either True or False.
    If True, some selected parameters are varied according to a normal distribution
    """
    cell_heterogeneity: bool = True

    """
    Two modes of stimulation
    1) nearest: stimulated cell and its nearest neighbours are effected
    2) point: only stimulated cell is effected
    """
    stimulated_cell: int = None
    valid_initiator_cell_modes: list[str] = ["point", "nearest"]
    initiator_cell_mode: str

    # Ligand diffusion parameters
    """
    Diffusion modes:
    1) point
    2) partially-regenerative
    3) regenerative
    4) decoupled

    partially-regenerative:
    Each sequential cell that activates secretes less ATP according to an
    exponential decay function L0(x)=L0init*e^(x/char_dist)

    regenerative:
    Each sequential cell that activates secretes the same amount of ATP as
    the initiator cell L0(x)=L0init

    decoupled:
    No ATP is secreted. Cells are not coupled with extracellular ligands.

    point:
    Only stimulated cell secretes ATP. L0_stim = L0init
    """
    valid_diffusion_modes: list[str] = ["partially-regenerative", "regenerative", "decoupled", "point"]
    diffusion_mode: str|None = None
    Datp: float = 236.0 #difusion constant for ATP (um^2/s)
    film_thickness: float = 100.0 ##thickness of thin film for diffusion (um)
    L0init: float = 0.2 #1.0 #concentration of secreted ATP by stimulated cell (fmol)
    char_dist: float = 30.0 #characteristic distance of ATP release decrease
    apyrase_deg: bool = False #True for apiraza, False no apiraza
    apyrase_char_time: float = 0.01 #s

    ###Difusion of Ca2+ and IP3 through gap junctions
    Dca=512.7 #difusion constant of Ca2+ through GJ
    Dip3=913.9 #difusion constant for IP3 through GJ
    dif_fact = Dca/Dip3

    # Threshold for cellular activity
    Cth = 0.15 ##threshold amplitude of normalized calcium
    time_interval_for_slope = 4.0 # seconds
    slopeTh = 0.003#0.01 ## uM/s
    Cth_act = 0.125 + 0.025

    # White noise amplitude for calcium and ip3
    ca_noise_amp = 0.000
    ip3_noise_amp = 0.00

    def __init__(self, model_parameters: dict[str, str|float]):
        """
        Initialize model with desired parameters
        """
        for key, value in model_parameters.items():
            setattr(self, key, value)
        
        if "diffusion_mode" in model_parameters.keys():
            self.diffusion_mode = model_parameters["diffusion_mode"].value
        if "initiator_cell_mode" in model_parameters.keys():
            self.initiator_cell_mode = model_parameters["initiator_cell_mode"].value

        self.check_diffusion_mode()
        self.check_initiator_cells_mode()

        self.stim_dur_steps: int = int(self.stim_dur/self.dt) #number of stimulation steps

        self.Cgjca=self.dif_fact*model_parameters["Cgjip3"] #apparent constant for Ca2+ difusion

        self.results_path: str = f"results/capsule_{self.capsule}"
        if not os.path.exists(self.results_path):
            os.makedirs(self.results_path)

    def check_diffusion_mode(self):
        """
        Checks if selected diffusion mode is valid
        """
        if self.diffusion_mode not in self.valid_diffusion_modes:
            raise Exception(f"Diffusion mode must be one of {self.valid_diffusion_modes}")
    
    def check_initiator_cells_mode(self):
        """
        Checks if selected initiator/stimulated cell mode is valid.
        Either the poked cell only or their nearest neighbours are stimulated
        """
        if self.initiator_cell_mode not in self.valid_initiator_cell_modes:
            raise Exception(f"Initiator cell mode must be one of: {self.valid_initiator_cell_modes}")

    def get_initiatior_cells(self, stimulated_cell: int, conn_mat: np.ndarray, mode: str) -> list[int]:
        """
        Finds all cells that are next to the stimulated cell
        If mode == 'point' only stimulated cell is considered an initiator
        If mode == 'nearest' neighbours of stimulated cell are also considered initiators
        """
        init_cells = []
        init_cells.append(stimulated_cell)
        if mode == 'point':
            return init_cells

        for i in np.nonzero(conn_mat[:, stimulated_cell])[0]:
            init_cells.append(int(i))
        return init_cells
    
    def save_ts_data(self, ca_ts: np.ndarray, ca_bin_ts: np.ndarray, ip3_ts: np.ndarray,
                     atp_ts: np.ndarray, act_times: np.ndarray,
                     jgjca_ts: np.ndarray, jgjip3_ts: np.ndarray):
        """
        Saves all data to results/capsule_{<INT>}
        """
        np.savetxt(f"{self.results_path}/ca_ts.txt", ca_ts, fmt="%.4lf")
        np.savetxt(f"{self.results_path}/ca_bin_ts.txt", ca_bin_ts, fmt="%d")
        np.savetxt(f"{self.results_path}/ip3_ts.txt", ip3_ts, fmt="%.4lf")
        np.savetxt(f"{self.results_path}/atp_ts.txt", atp_ts, fmt="%.4lf")
        np.savetxt(f"{self.results_path}/act_times.txt", act_times, fmt="%.4lf")
        np.savetxt(f"{self.results_path}/jgjca_ts.txt", jgjca_ts, fmt="%.4lf")
        np.savetxt(f"{self.results_path}/jgjip3_ts.txt", jgjip3_ts, fmt="%.4lf")
    
    def plot_time_series(self, time_series: np.ndarray, bin_time_series, ylabel: str = "Calcium signal", xlabel: str = "time (s)"):
        """
        Plots provided time series and saves them to the results/capsule_{INT}/time_series folder
        """
        ts_path = f"{self.results_path}/time_series"
        if not os.path.exists(ts_path):
            os.makedirs(ts_path)

        time: np.ndarray = time_series[:,0]
        cell_num: int = len(time_series[0]) - 1
        for i in range(cell_num):
            max_v = np.max(time_series[:, i+1])
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(1,1,1)
            ax.plot(time, time_series[:, i+1], linewidth=1, c='black')
            ax.plot(time, bin_time_series[:, i+1]*max_v, linewidth=1, c='red')
            ax.set_xlim(0, np.amax(time))
            ax.set_ylim(0.99*np.amin(time_series[:,i+1]), 1.01*np.amax(time_series[:,i+1]))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            fig.savefig(f"{ts_path}/cell_{i}.png",
                        dpi=600, bbox_inches="tight", pad_inches=0.01)
            plt.close(fig)
    
    def run_simulation(self, cells: list, cell_distances: np.ndarray):
        """
        Runs actual simulation
        """
        simulation_time = 0.0
        simulation_step = 0
        while simulation_time < self.tend:
            if simulation_step%100 == 0:
                print(simulation_time)
            calcium_coupling, ip3_coupling = self.calculate_ca_ip3_coupling(cells)
            atp_coupling = self.calculate_atp_coupling(
                                                        cells,
                                                        cell_distances,
                                                        self.stimulated_cell,
                                                        simulation_time,
                                                        mode=self.diffusion_mode
                                                    )
            #atp_coupling = [0 for i in range(len(cells))]
            for i, cell in enumerate(cells):
                cell.run_model_step(simulation_step,
                                    simulation_time,
                                    atp_coupling[i],
                                    calcium_coupling[i],
                                    ip3_coupling[i])
            simulation_step += 1
            simulation_time += self.dt

    @staticmethod
    def find_stimulated_cell(cell_pos: np.ndarray, target: tuple) -> int:
        """
        Finds cell closest to provided coordinates
        """
        dist = np.array([np.sqrt((target[0]-cell_pos[i,0])**2+(target[1]-cell_pos[i,1])**2) for i in range(len(cell_pos))])
        closest_cell = int(np.argmin(dist))
        return closest_cell

    @staticmethod
    def calculate_ca_ip3_coupling(cells: list) -> tuple[list]:
        """
        Calculates coupling currents for Ca and IP3
        """

        ca_currents = []
        ip3_currents = []
        for cell in cells:
            ca_current = 0.0
            ip3_current = 0.0
            for neighbour, weight in cell.neighbours:
                ip3_current += (cell.model.Cgjip3*weight)*(cells[neighbour].P-cell.P)
                ca_current += (cell.model.Cgjca*weight)*(cells[neighbour].C-cell.C)
            ca_currents.append(float(ca_current))
            ip3_currents.append(float(ip3_current))

        #print("Ca cell coupling: ", ca_currents)
        return ca_currents, ip3_currents
    
    @staticmethod
    def calculate_atp_coupling(cells: list,
                           cell_distances: np.ndarray,
                           stimulated_cell: int,
                           time: float,
                           mode: str = 'partially-regenerative') -> list[float]:
        """
        Calculates the diffusion ligand coupling between cells
        Four modes are available (point, partially-regenerative, regenerative, decoupled):
        1) point source mode: only stimulated cell secretes the atp
        2) partially-regenerative mode: activated cells secrete sequentially less atp
        3) regenerative mode: activated cells secrete the same amount of ligand
        4) decoupled mode: There is no ligand diffusion
        """
        if mode == "decoupled":
            return [0.0 for _ in cells]
        all_cells_atp = []
        for i in range(len(cells)):
            L_i = 0.0
            for j, cell_j in enumerate(cells):
                if cell_j.time_of_activation:
                    L0_j = 0.0
                    if mode == 'point' and j == stimulated_cell:
                        L0_j = cell_j.model.L0init
                    elif mode == 'partially-regenerative':
                        #Activated cells secrete sequentially less ATP
                        L0_j = cell_j.model.L0init * \
                            np.e**(-cell_distances[stimulated_cell, j]/cell_j.model.char_dist)
                    elif mode == 'regenerative':
                        #Activated cells secrete the same amount of ATP
                        L0_j = cell_j.model.L0init

                    L0_j = ((L0_j*10**6)/(cell_j.model.film_thickness*4.0*np.pi*cell_j.model.Datp*(time-cell_j.time_of_activation))) * \
                        np.e**(-(cell_distances[i, j]**2)/(4.0 * cell_j.model.Datp*(time-cell_j.time_of_activation)))
                    
                    if(cell_j.model.apyrase_deg):
                        #char_time_atp_deg = (cell_j.apyrase_Km*np.log(3)+(2/3)*L0_j)/(cell_j.apyrase_Vmax)*60.0 #characteristic time of atp degradation in s
                        char_time_atp_deg = cell_j.model.apyrase_char_time
                        L0_j = L0_j * np.e**(-((time-cell_j.time_of_activation)/(char_time_atp_deg)))
                    
                    L_i += L0_j
            all_cells_atp.append(L_i)
        return all_cells_atp
    

    @staticmethod
    def extract_time_series_data(cells: list) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Extracts and returns all time series data as numpy arrays
        """

        ts_length: int = len(cells[0].calcium_time_series)
        cell_num: int = len(cells)
        time: list[float] = cells[0].simulation_time
        ca_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)
        ip3_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)
        atp_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)
        jgjca_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)
        jgjip3_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)


        ca_ts[:,0] = np.array(time)
        ip3_ts[:,0] = np.array(time)
        atp_ts[:,0] = np.array(time)
        jgjca_ts[:,0] = np.array(time)
        jgjip3_ts[:,0] = np.array(time)

        for i, cell in enumerate(cells):
            ca_ts[:, i+1] = np.array(cell.calcium_time_series)
            ip3_ts[:, i+1] = np.array(cell.ip3_time_series)
            atp_ts[:, i+1] = np.array(cell.atp_time_series)
            jgjca_ts[:, i+1] = np.array(cell.jgjca_time_series)
            jgjip3_ts[:, i+1] = np.array(cell.jgjip3_time_series)

        return ca_ts, ip3_ts, atp_ts, jgjca_ts, jgjip3_ts
    
    @staticmethod
    def extract_bin_signals(ca_ts: np.ndarray, act_frames: np.ndarray, deact_frames: np.ndarray):
        """
        Extracts binarized signals based on ca ts, activation time and activation amplitude
        """
        ca_bin_ts = np.zeros(ca_ts.shape, float)
        ts_length: int = len(ca_ts)
        ca_bin_ts[:,0] = ca_ts[:,0]
        cell_num: int = len(ca_ts[0])-1
        for i in range(cell_num):
            if not np.isnan(act_frames[i]):
                start: int = int(act_frames[i])
                end: int = int(deact_frames[i])
                ca_bin_ts[start:end, i+1] = 1
        return ca_bin_ts
        

    @staticmethod
    def extracts_activation_times(cells: list):
        """
        Extracts activation times of cells
        """
        act_times: np.ndarray = np.zeros(len(cells), float)
        act_amps: np.ndarray = np.zeros(len(cells), float)
        act_frames: np.ndarray = np.zeros(len(cells), float)
        for i, cell in enumerate(cells):
            act_times[i] = np.nan
            act_amps[i] = np.nan
            act_frames[i] = np.nan
            if(cell.time_of_activation):
                act_times[i] = cell.time_of_activation
                act_amps[i] = cell.activation_amp
                act_frames[i] = cell.index_of_activation

        #act_times = act_times - np.nanmin(act_times)
        return act_times, act_amps, act_frames
    
    def plot_activation_sequence(self, pos: np.ndarray, act_times: np.ndarray,
                                 capsule_data):
        """
        Plots activation sequence of cells
        """
        cell_num = len(pos)
        vmin=np.nanmin(act_times, axis=None)
        vmax=np.nanmax(act_times, axis=None)

        tru_col_map=[]
        for i in range(cell_num):
            if np.isnan(act_times[i]):
                tru_col_map.append(plt.cm.gray(0.4))
            else:
                tru_col_map.append(plt.cm.jet_r(float(act_times[i]-vmin)/float(vmax-vmin)))

        dx=10.0

        minx = np.amin(pos[:,0]) - 6
        miny = np.amin(pos[:,1]) - 6

        fig=plt.figure(figsize=(2.0*5, 2.0*5))
        ax=fig.add_subplot(111)
        for i, vor_region in capsule_data["cells"].items():
            polygon = Polygon(vor_region["vertices_coord"], closed=True,
                              fc=tru_col_map[i], ec='black', zorder=0)
            ax.add_patch(polygon)
        ax.scatter(pos[:,0], pos[:,1], s=12, c="black")
        ax.scatter(pos[:,0], pos[:,1], s=12, c='black', zorder=10)
        ax.scatter(pos[self.stimulated_cell,0], pos[self.stimulated_cell,1], s=30, color='black', marker='x')

        ax.plot([minx,minx+dx], [miny, miny], linewidth=3, color='black')
        ax.text(minx+1, miny-4, r'10 $\mathrm{\mu}$m', fontsize=14)

        ax.set_axis_off()
        fig.savefig(
            f"results/capsule_{self.capsule}/activation_sequence.png",
            dpi=600, bbox_inches='tight')
        plt.close(fig)
    
    def plot_gj_currents(self, stimulated_cell: int, 
                         stim_cell_neighbours: list[int],
                         jgjca_ts: np.ndarray, jgjip3_ts: np.ndarray):
        """
        Draws and saved GJ currents for stimulated cell and their neighbours.
        """
        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        for cell in [stimulated_cell, *stim_cell_neighbours]:
            ax1.plot(jgjca_ts[:,0], jgjca_ts[:,cell], linewidth=1.5, label=f"{cell=}")
            ax2.plot(jgjip3_ts[:,0], jgjip3_ts[:,cell], linewidth=1.5, label=f"{cell=}")
        ax1.legend(loc="upper right")
        ax2.legend(loc="upper right")
        ax1.set_xlabel("time (s)")
        ax2.set_xlabel("time (s)")

        ax1.set_ylabel("GJ Ca2+ current")
        ax2.set_ylabel("GJ IP3 current")

        ax1.set_xlim(0, 50.0)
        ax2.set_xlim(0, 50.0)

        plt.subplots_adjust(wspace=0.2, hspace=0.2)
        fig.savefig(f"results/capsule_{self.capsule}/GJ_currents.png", dpi=600,
                                        bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)
    
    def save_model_parameters(self, model_parameters: dict):
        """
        Saves provided model parameters as JSON file
        to results folder
        """
        with open(f"results/capsule_{self.capsule}/model_parameters.json",
                  'w', encoding="utf-8") as json_file:
            json.dump(model_parameters, json_file, indent=4, cls=self.EnumEncoder)
    
    def calculate_activity_params(self, ca_ts: np.ndarray, act_times: np.ndarray,
                                  start_amps: np.ndarray, act_frames: np.ndarray, time: np.ndarray) -> tuple[np.ndarray]:
        """
        Calculates cell parameters: signal duration, signal amplitude and response time
        """
        max_amps_indx: list[int] = [np.argmax(ca_ts[int(act_frame):, i]) + int(act_frame) if not np.isnan(act_frame) else np.nan for i, act_frame in enumerate(act_frames)]
        max_amps: list[float] = [ca_ts[max_indx, i] if not np.isnan(max_indx) else np.nan for i, max_indx in enumerate(max_amps_indx)]
        deact_frames: list[np.ndarray] = [np.where(ca_ts[max_amp_indx:, i] <= (max_amps[i]+start_amps[i])/2.0)[0] if not np.isnan(max_amp_indx) else np.nan for i, max_amp_indx in enumerate(max_amps_indx)]
        deact_frames: list[int] = [deact[0] if isinstance(deact, np.ndarray) else np.nan for deact in deact_frames]
        deact_frames: list[int] = [deact+max_amp_indx if not np.isnan(deact) and not np.isnan(max_amp_indx) else np.nan for deact, max_amp_indx in zip(deact_frames, max_amps_indx)]
        deact_frames: list[int] = [len(time)-1 if not np.isnan(act_frame) and np.isnan(deact) else deact for act_frame, deact in zip(act_frames, deact_frames)]

        sig_durs: list[float] = [(time[int(end_frame)] - time[int(start_frame)]) if not np.isnan(start_frame) and not np.isnan(end_frame) else np.nan
                                 for start_frame, end_frame in zip(act_frames, deact_frames)]
        response_times: np.ndarray = act_times - np.nanmin(act_times)

        return np.array(sig_durs), response_times, np.array(max_amps), np.array(deact_frames)

    def save_activity_params(self, durations: np.ndarray,
                             resp_times: np.ndarray, max_amps: np.ndarray,
                             fractions_act_cells: np.ndarray):
        """
        Saves activity params
        """
        np.savetxt(f"results/capsule_{self.capsule}/signal_durations.txt", durations, fmt="%.4lf")
        np.savetxt(f"results/capsule_{self.capsule}/response_times.txt", resp_times, fmt="%.4lf")
        np.savetxt(f"results/capsule_{self.capsule}/max_amplitudes.txt", max_amps, fmt="%.4lf")
        np.savetxt(f"results/capsule_{self.capsule}/fraction_of_active_cells.txt",
                   fractions_act_cells, fmt="%.2lf")

    def cluster_cells(self, pos: np.ndarray, bins: int) -> tuple[dict[list[int]], list[float]]:
        """
        Clusters cells into groups based on distance from stimulated cell
        """
        distances = np.sqrt(np.sum((pos - pos[self.stimulated_cell])**2, axis=1))
        sorted_dist: np.ndarray = np.argsort(distances)
        elements_in_bin: int = int(len(pos)/bins)
        groups: dict[list[int]] = {}
        for current_bin in range(bins):
            group_start: int = int(current_bin*elements_in_bin)
            group_end: int = group_start + elements_in_bin
            groups[current_bin] = list(sorted_dist[group_start:group_end+1])

        group_distances: list[float] = [np.average(distances[members])
                                        for members in groups.values()]
        return groups, group_distances

    def calculate_bin_avg_and_stderr(self, data: np.ndarray, groups: dict[list[int]]) -> np.ndarray:
        """
        Calculates average value and stderr for each group of cells
        """
        avg_stderr: np.ndarray = np.zeros((len(groups), 2), float)
        for group, members in groups.items():
            valid_cells: int = sum([1 if not np.isnan(data[member]) else 0 for member in members])
            avg_value: float = np.nanmean(data[members]) if valid_cells>0 else np.nan
            std_err: float = np.nanstd(data[members])/np.sqrt(len(members)) if valid_cells>0 else np.nan
            avg_stderr[group, 0] = avg_value
            avg_stderr[group, 1] = std_err
        return avg_stderr

    def get_fraction_of_act_cells(self, cell_groups: dict[list[int]],
                                  resp_times: np.ndarray) -> list[float]:
        """
        Calculates number of active cells for each group
        """
        fractions: list[float] = []
        for _, members in cell_groups.items():
            fraction: float = sum([1 if not np.isnan(resp_time) else 0 for resp_time in resp_times[members]])/len(members)
            fractions.append(fraction)
        return fractions

    def plot_errorbars(self, ax, x: list[int], data: np.ndarray,
                       ylabel, xlabel=r"Distance from stimulation ($\mathrm{\mu}$m)"):
        """
        Adds errorbar plot to axes object
        """
        ax.errorbar(x, data[:,0], yerr=data[:,1], fmt='-o',
                    capsize=5, capthick=2, elinewidth=2, c='black')
        ax.set_ylabel(ylabel)
        ax.set_xlim(0.95*np.amin(x), 1.02*np.amax(x))
        ax.set_xlabel(xlabel)

        return ax

    def plot_activity_params(self, durations: np.ndarray, resp_times: np.ndarray,
                             amps: np.ndarray, pos: np.ndarray, bins: int = 5):
        """
        Plots all cell activity params
        """
        cell_groups: dict[list[int]]
        group_distances: list[float]
        cell_groups, group_distances = self.cluster_cells(pos, bins)
        clustered_durations: np.ndarray = self.calculate_bin_avg_and_stderr(durations, cell_groups)
        clustered_resp_times: np.ndarray = self.calculate_bin_avg_and_stderr(resp_times,
                                                                             cell_groups)
        clustered_amps: np.ndarray = self.calculate_bin_avg_and_stderr(amps, cell_groups)
        fraction_act_cells: list[float] = self.get_fraction_of_act_cells(cell_groups, resp_times)
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        axes[0,1] = self.plot_errorbars(axes[0,1], group_distances, clustered_durations,
                                        "Signal duration (s)")
        axes[0,0] = self.plot_errorbars(axes[0,0], group_distances, clustered_resp_times,
                                        "Response time (s)")
        axes[1,0].plot(group_distances, fraction_act_cells, c='black', linewidth=2.0)
        axes[1,0].set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
        axes[1,0].set_ylabel("Fraction of active cells")
        axes[1,0].set_xlim(0.95*np.amin(group_distances), 1.02*np.amax(group_distances))
        axes[1,0].set_ylim(0,1.02)
        axes[1,1] = self.plot_errorbars(axes[1,1], group_distances, clustered_amps,
                                        "Signal amplitude")
        fig.subplots_adjust(wspace=0.2, hspace=0.2)
        fig.savefig(f"results/capsule_{self.capsule}/cell_activity_params.png",
                    dpi=600, bbox_inches='tight', pad_inches=0.01)
        plt.close(fig)
        return fraction_act_cells
