"""
All coupling functions
"""
# pylint: disable=C0103

from typing import List, Tuple
import numpy as np

from .cell_model import CellModel

def get_neighbours(cell: int, conn_mat: np.ndarray, weights: np.ndarray) -> List[int]:
    """
    Finds all neighbours of cell
    Returns list of neighbours indices
    """
    neighbours = []
    for i in np.where(conn_mat[:, cell] == 1.0)[0]:
        if i != cell:
            neighbours.append((int(i), weights[i, cell]))
    return neighbours


def calculate_ca_ip3_coupling(cells: List[CellModel]) -> Tuple[list]:
    """
    Calculates coupling currents for Ca and IP3
    """

    ca_currents = []
    ip3_currents = []
    for cell in cells:
        ca_current = 0.0
        ip3_current = 0.0
        for neighbour, weight in cell.neighbours:
            ip3_current += (cell.Cgjip3*weight)*(cells[neighbour].P-cell.P)
            ca_current += (cell.Cgjca*weight)*(cells[neighbour].C-cell.C)
        ca_currents.append(ca_current)
        ip3_currents.append(ip3_current)

    return ca_currents, ip3_currents


def calculate_atp_coupling(cells: List[CellModel],
                           cell_distances: np.ndarray,
                           stimulated_cell: int,
                           time: float,
                           mode: str = 'partially-regenerative') -> List[float]:
    """
    Calculates the diffusion ligand coupling between cells
    Four modes are available (point, partially-regenerative, regenerative, decoupled):
    1) point source mode: only stimulated cell secretes the atp
    2) partially-regenerative mode: activated cells secrete sequentially less atp
    3) regenerative mode: activated cells secrete the same amount of ligand
    4) decoupled mode: There is no ligand diffusion
    """
    all_cells_atp = []
    for i in range(len(cells)):
        L_i = 0.0
        for j, cell_j in enumerate(cells):
            if cell_j.time_of_activation:
                L0_j = 0.0 ##decoupled
                if mode == 'point' and j == stimulated_cell:
                    L0_j = cell_j.L0init
                elif mode == 'partially-regenerative':
                    #Activated cells secrete sequentially less ATP
                    L0_j = cell_j.L0init * \
                        np.e**(-cell_distances[stimulated_cell, j]/cell_j.char_dist)
                elif mode == 'regenerative':
                    #Activated cells secrete the same amount of ATP
                    L0_j = cell_j.L0init

                L_i += ((L0_j*10**6)/(cell_j.film_thickness*4.0*np.pi*cell_j.Datp*(time-cell_j.time_of_activation))) * \
                    np.e**(-(cell_distances[i, j]**2)/(4.0 *
                           cell_j.Datp*(time-cell_j.time_of_activation)))
        all_cells_atp.append(L_i)
    return all_cells_atp


def get_initiatior_cells(stimulated_cell: int, conn_mat: np.ndarray, mode: str) -> List[int]:
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
