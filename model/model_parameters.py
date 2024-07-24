"""
Numerical model parameters
"""
# pylint: disable=R0903
import numpy as np

class ModelParameters:
    """
    All model parameters
    """
    capsule: int = 0
    t: float = 0.0 #initial time
    dt: float = 0.01 #integration step
    tend: float = 150.0 #end time
    time_of_stimulation: float = 5.0 #time of stimulation

    stim_dur: float = 3.0 #stimulation duration
    stim_dur_steps: int = int(stim_dur/dt) #number of stimulation steps

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
    initiator_cell_mode: str = 'nearest'

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
    diffusion_mode: str = 'point'
    Datp: float = 236.0 #difusion constant for ATP (um^2/s)
    film_thickness: float = 100.0 ##thickness of thin film for diffusion (um)
    L0init: float = 1.0 #concentration of secreted ATP by stimulated cell (fmol)
    char_dist: float = 30.0 #characteristic distance of ATP release decay

    # Threshold for cellular activity
    Cth = 0.15 ##threshold amplitude of normalized calcium
    time_interval_for_slope = 2.0 # seconds
    slopeTh = 0.01 ## uM/s

    # White noise amplitude for calcium and ip3
    ca_noise_amp = 0.001
    ip3_noise_amp = 0.001

    @staticmethod
    def find_stimulated_cell(cell_pos: np.ndarray, target: tuple) -> int:
        """
        Finds cell closest to provided coordinates
        """
        dist = np.array([np.sqrt((target[0]-cell_pos[i,0])**2+(target[1]-cell_pos[i,1])**2) for i in range(len(cell_pos))])
        closest_cell = np.argmin(dist)
        return closest_cell
