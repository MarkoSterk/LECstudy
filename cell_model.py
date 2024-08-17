"""
Model for single cell
"""
# pylint: disable=R0913, R0902, C0103, E1101
import random
import numpy as np

from cell_parameters import CellParameters

class CellModel(CellParameters):
    """
    CellModel class holds the entire model for each cell
    """

    def __init__(self,
                 model_parameters,
                 custom_parameters: dict[str, float],
                 cell_num: int,
                 volume: float,
                 dist_from_stimulation: float,
                 init_cells: list[int],
                 neighbours: list[int] = None
                 ):
        """
        Initializer for single cell
        """
        self.model = model_parameters

        for key, value in custom_parameters.items():
            setattr(self, key, value)

        self.cell_num: int = cell_num
        self.volume: int = volume
        self.dist_from_stimulation: int = dist_from_stimulation
        self.init_cells: list[int] = init_cells

        self.neighbours: list[int] = []
        if neighbours:
            self.neighbours = neighbours

        self.C: float = 0.125
        self.P: float = 0.01
        self.Osscc: float = 0.0
        self.stretch: float = 0.0
        self.Rs: float = 17000.0
        self.Rsp: float = 0.0
        self.G: float = 14.0
        self.PIP: float = 49997.0
        self.k_ip3_in: float = 0.0
        self.k_ca_out: float = 0.0

        self.CA_OUT: float = 0.0 #amount of added Ca2+ into stimulated cell(s)
        self.IP3_IN: float = 0.5771 #amount of added IP3 into stimulated cells

        #data storage
        self.time_of_activation = None
        self.activation_amp = None
        self.index_of_activation = None
        self.simulation_time: list[float] = []
        self.calcium_time_series: list[float] = []
        self.ca_bin_time_series: list[float] = []
        self.ip3_time_series: list[float] = []
        self.atp_time_series: list[float] = []
        self.jgjca_time_series: list[float] = []
        self.jgjip3_time_series: list[float] = []

    def __repr__(self) -> str:
        """
        String representation for printing of class instance
        """
        return f"Cell(index={self.cell_num})"

    def run_model_step(self, step: int, time: float, L: float, jgjca: float, jgjip3: float):
        """
        Runs next model step
        """
        self.store_data(step, time, L, jgjca, jgjip3)
        self.calculate_activation_time(time)

        dstretch: float = self.calculate_strech(time)
        currents: dict[str, float] = {**self.calculate_currents(L), "jgjca": jgjca, "jgjip3": jgjip3, "L": L}
        differentials: dict[str, float] = {**self.calculate_differentials(currents), "dstretch": dstretch}
        self.calculate_new_state(differentials)

    def calculate_strech(self, time: float) -> float:
        """
        Calculates current strech of cell due to stimulation
        """
        dstretch: float = 0.0
        # Stimulation
        if self.model.time_of_stimulation <= time < (self.model.time_of_stimulation + self.model.stim_dur):
            self.stretch = self.Astretch * \
                np.e**(-self.stretch_exp*self.dist_from_stimulation)
            if self.cell_num in self.init_cells:
                self.k_ip3_in = self.IP3_IN/self.model.stim_dur
        else:
            dstretch = (-self.k_stretch*self.stretch)*self.model.dt
            self.k_ip3_in = 0.0
        return dstretch

    def calculate_currents(self, L: float) -> dict[str, float]:
        """
        Calculates all currents
        """
        jpump: float = (self.Vpump*self.C**2) / (self.Kpump**2+self.C**2)

        w_inf = ((1.0+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)) /
            (1.0+(1.0/self.Kc)+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)))
        pryr = ((w_inf*(1.0+((self.C/self.Kb)**3))) / (1.0+((self.Ka/self.C)**4)+((self.C/self.Kb)**3)))
        jryr: float = self.kryr*pryr

        alpha4 = self.Aalpha4 * np.e**(self.exp_alpha4*self.dist_from_stimulation)
        k1 = (self.alpha1*(self.C**3))/((self.beta1**3)+(self.C**3))
        k4 = (alpha4*self.P)/(self.beta4+self.P)
        phi = (1.0)/(1.0+((self.k2)/(self.k3+k4))*(1.0+k4/self.k5))
        O = (phi*self.P)/(((self.k1m+self.k2)/(k1))*phi+self.P)
        jip3r3: float = self.kip3r3*O**4

        jsscc: float = self.ksscc*self.Osscc
        pr: float = (L*self.Rs)/(self.epsilon*self.Rtot*(self.K1+L))
        rh: float = self.alpha*((self.C)/(self.K3+self.C))*self.G

        return {
            "jpump": jpump,
            "jryr": jryr,
            "jip3r3": jip3r3,
            "jsscc": jsscc,
            "pr": pr,
            "rh": rh
        }

    def calculate_activation_time(self, time: float):
        """
        Calculates cells activation time from TS slope
        """
        time_step = self.model.record_every * self.model.dt ## time between two data points
        points_for_slope = int(self.model.time_interval_for_slope/time_step)
        if len(self.calcium_time_series)>points_for_slope+5 and not self.time_of_activation:
            slope = (self.calcium_time_series[-1] - self.calcium_time_series[-points_for_slope])/(points_for_slope*time_step)
            if slope>self.model.slopeTh and self.calcium_time_series[-1]>self.model.Cth_act:
                self.time_of_activation = time
                self.index_of_activation = len(self.calcium_time_series)
                self.activation_amp = self.calcium_time_series[-1]
                #print(f"Cell {self.cell_num} activated with amplitude {self.activation_amp} ({self.model.Cth_act})")


    def calculate_differentials(self, currents: dict[str, float]):
        """
        Calculates all differentials based on time step and currents
        """
        ca_noise: float = self.model.ca_noise_amp*np.random.uniform(-1, 1)
        p_noise = self.model.ip3_noise_amp*np.random.uniform(-1, 1)

        dC = (currents["jgjca"]+self.Jleak-currents["jpump"]+currents["jryr"] +
                   currents["jip3r3"]+currents["jsscc"]+self.k_ca_out)*self.model.dt + ca_noise
        dG = (self.ka*(self.sigma+currents["pr"]) * (self.Gtot-self.G)-self.kd*self.G)*self.model.dt
        dRsp = (currents["L"]*(((self.kp*self.Rs)/(self.K1+currents["L"])) -((self.ke*self.Rsp)/(self.K2+currents["L"]))))*self.model.dt
        dRs = (self.kr*self.Rtot-(self.kr+((self.kp*currents["L"]) /
                            (self.K1+currents["L"])))*self.Rs-self.kr*self.Rsp)*self.model.dt
        dP = (currents["rh"]*((1.0)/(self.Na*self.volume))*(10**(6)) *
                        self.PIP-self.kdeg*self.P+currents["jgjip3"]+self.k_ip3_in)*self.model.dt + p_noise
        dPIP = (self.rr*self.PIP2tot-(currents["rh"]+self.rr)*self.PIP -
                            self.rr*self.Na*self.volume*self.P*0.000001)*self.model.dt
        dOsscc = (self.stretch*self.kf-(self.stretch * self.kf+self.kb)*self.Osscc)*self.model.dt

        return {
            "dC": dC,
            "dG": dG,
            "dRsp": dRsp,
            "dRs": dRs,
            "dP": dP,
            "dPIP": dPIP,
            "dOsscc": dOsscc
        }

    def calculate_new_state(self, differentials: dict[str, float]) -> None:
        """
        Calculates new state values for cell
        """
        if differentials["dstretch"]:
            self.stretch+=differentials["dstretch"]
        self.Osscc += differentials["dOsscc"]
        self.C += differentials["dC"]
        self.G += differentials["dG"]
        self.Rsp += differentials["dRsp"]
        self.Rs += differentials["dRs"]
        self.P += differentials["dP"]
        self.PIP += differentials["dPIP"]

    def store_data(self, step: int, time: float, L: int, jgjca: float, jgjip3: float):
        """
        Stores current data to designated lists
        """
        if step%self.model.record_every == 0:
            self.simulation_time.append(time)
            self.calcium_time_series.append(self.C)
            self.ip3_time_series.append(self.P)
            self.atp_time_series.append(L)
            self.jgjca_time_series.append(jgjca)
            self.jgjip3_time_series.append(jgjip3)

    @classmethod
    def generate_cells(cls, model, stimulated_cell: int, init_cells: list[int],
                       pos: np.ndarray, volumes: np.ndarray,
                       dist_from_stimulation: np.ndarray, cmat: np.ndarray,
                       weights: np.ndarray, noise: bool = False) -> list:
        """
        Generates all cells for model
        """
        cells = []
        for i in range(len(pos)):
            volume = volumes[i]*10**(-18) #conversions
            neighbours = CellModel.get_neighbours(i, cmat, weights)
            cells.append(cls(
                model, CellModel.generate_cell_parameters(noise),
                i, volume, dist_from_stimulation[i, stimulated_cell],
                init_cells, neighbours
            ))
        return cells
    
    @staticmethod
    def generate_cell_parameters(noise: bool) -> dict[str, float]:
        """
        Generates missing cell parameters with or without noise
        """
        return {
            "ksscc": 1.20*CellModel.random_normal_number(noise),
            "Kpump": 0.5030,#heterogenost povzroča nestabilnost
            "kryr": 16.04*CellModel.random_normal_number(noise),
            "kip3r3": 155.8*CellModel.random_normal_number(noise),
            "kdeg": 1.25*CellModel.random_normal_number(noise)
        }
    
    @staticmethod
    def random_normal_number(noise=True, mean=1.0, std_dev=0.0):
        """
        Generates a random number
        """
        if not noise:
            return 1.0
        while True:
            num = random.normalvariate(mean, std_dev)
            if 0.9 <= num <= 1.1:
                return num
    
    @staticmethod
    def get_neighbours(cell: int, conn_mat: np.ndarray, weights: np.ndarray) -> list[int]:
        """
        Finds all neighbours of cell
        Returns list of neighbours indices
        """
        neighbours = []
        for i in np.where(conn_mat[:, cell] == 1.0)[0]:
            if i != cell:
                neighbours.append((int(i), weights[i, cell]))
        return neighbours
