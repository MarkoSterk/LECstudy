"""
All relevant model parameters
"""
# pylint: disable=W0719
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

class ModelParameters:
    """
    All model parameters
    """
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
    
    def save_ts_data(self, ca_ts: np.ndarray, ca_bin_ts: np.ndarray, ip3_ts: np.ndarray, atp_ts: np.ndarray, act_times: np.ndarray):
        """
        Saves all data to results/capsule_{<INT>}
        """
        np.savetxt(f"{self.results_path}/ca_ts.txt", ca_ts)
        np.savetxt(f"{self.results_path}/ca_bin_ts.txt", ca_bin_ts, fmt="%d")
        np.savetxt(f"{self.results_path}/ip3_ts.txt", ip3_ts)
        np.savetxt(f"{self.results_path}/atp_ts.txt", atp_ts)
        np.savetxt(f"{self.results_path}/act_times.txt", act_times)
    
    def plot_time_series(self, time_series: np.ndarray, ylabel: str = "Calcium signal", xlabel: str = "time (s)"):
        """
        Plots provided time series and saves them to the results/capsule_{INT}/time_series folder
        """
        ts_path = f"{self.results_path}/time_series"
        if not os.path.exists(ts_path):
            os.makedirs(ts_path)

        time: np.ndarray = time_series[:,0]
        cell_num: int = len(time_series[0]) - 1
        for i in range(cell_num):
            fig = plt.figure(figsize=(6, 4))
            ax = fig.add_subplot(1,1,1)
            ax.plot(time, time_series[:, i+1], linewidth=1, c='black')
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
                    L0_j = 0.0 ##decoupled
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
                        np.e**(-(cell_distances[i, j]**2)/(4.0 *
                            cell_j.model.Datp*(time-cell_j.time_of_activation)))
                    
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
        ca_bin_ts: np.ndarray = np.zeros((ts_length, cell_num+1), int)
        ip3_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)
        atp_ts: np.ndarray = np.zeros((ts_length, cell_num+1), float)


        ca_ts[:,0] = np.array(time)
        ip3_ts[:,0] = np.array(time)
        atp_ts[:,0] = np.array(time)

        for i, cell in enumerate(cells):
            ca_ts[:, i+1] = np.array(cell.calcium_time_series)
            ip3_ts[:, i+1] = np.array(cell.ip3_time_series)
            atp_ts[:, i+1] = np.array(cell.atp_time_series)
        
        ca_bin_ts = ModelParameters.extract_bin_signals(cells, ca_ts)

        return ca_ts, ca_bin_ts, ip3_ts, atp_ts
    
    @staticmethod
    def extract_bin_signals(cells: list, ca_ts: np.ndarray):
        """
        Extracts binarized signals based on ca ts, activation time and activation amplitude
        """
        ca_bin_ts = np.zeros(ca_ts.shape, float)
        ts_length: int = len(ca_ts)
        ca_bin_ts[:,0] = ca_ts[:,0]
        for i, cell in enumerate(cells):
            if cell.index_of_activation is not None:
                start: int = cell.index_of_activation
                max_amp = np.amax(ca_ts[start:, i+1])
                max_time = np.where(ca_ts[start:, i+1] == max_amp)[0]
                end = None
                for j in range(start+1, ts_length, 1):
                    if ca_ts[j, i+1] < (cell.activation_amp+max_amp)/2.0 and j > max_time:
                        end = j
                        break
                ca_bin_ts[start:, i+1] = 1
                if end is not None:
                    ca_bin_ts[end:, i+1] = 0
        return ca_bin_ts
        

    @staticmethod
    def extracts_activation_times(cells: list):
        """
        Extracts activation times of cells
        """
        act_times: np.ndarray = np.zeros(len(cells), float)
        for i, cell in enumerate(cells):
            act_times[i] = np.nan
            if(cell.time_of_activation):
                act_times[i] = cell.time_of_activation

        return act_times
    
    @staticmethod
    def plot_activation_sequence(model):
        """
        Plots activation sequence of cells
        """
        stimulated_cell = model.stimulated_cell
        pos = np.loadtxt(f'cell_data/capsule_{model.capsule}/cm_cells.txt')
        edgepoints = np.loadtxt(f'cell_data/capsule_{model.capsule}/cell_edge_points_x_y_cellnum.txt')
        act_times = np.loadtxt(f'results/capsule_{model.capsule}/act_times.txt')
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
        dy=10.0
        x0=10.0
        y0=10.0

        fig=plt.figure(figsize=(2.0*5, 2.0*5))
        ax=fig.add_subplot(111)
        for i in range(cell_num):
            data=edgepoints[np.where(edgepoints[:,2]==i)][:,:2]
            polygon=Polygon(list(data), fc=tru_col_map[i], ec='black', zorder=0)
            ax.add_artist(polygon)
        ax.scatter(pos[:,0], pos[:,1], s=12, c='black', zorder=10)
        ax.scatter(pos[stimulated_cell,0], pos[stimulated_cell,1], s=30, color='black', marker='x')

        ax.plot([x0,x0+dx], [y0, y0], linewidth=2, color='black')
        ax.text(x0, y0-10, '10 um', fontsize=8)
        ax.plot([x0, x0], [y0,y0+dy], linewidth=2, color='black')
        ax.text(x0-10, y0, '10 um', fontsize=8, rotation='vertical')

        ax.set_axis_off()
        fig.savefig(
            f"results/capsule_{model.capsule}/activation_sequence.png",
            dpi=600, bbox_inches='tight')
        plt.close(fig)
