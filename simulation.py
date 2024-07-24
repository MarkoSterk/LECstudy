"""
Main simulation file
"""
# pylint: disable=C0103
# pylint: disable=E0401
import os
import random
import numpy as np
from scipy.spatial import distance

from model.cell_model import CellModel
from model.cell_coupling import (get_neighbours,
                                 calculate_ca_ip3_coupling,
                                 calculate_atp_coupling,
                                 get_initiatior_cells)
from model.model_parameters import ModelParameters as MP
from model.cell_heterogeneity import cell_heterogeneity
from helper_functions.get_time_series import get_time_series

from plot_activation_sequence import plot_activation_sequence
from plot_time_series import plot_time_series
from plot_movie_frames import plot_movie_frames
from plot_signaling_parameters import plot_signaling_parameters

def random_normal_number(mean=1.0, std_dev=0.0):
    while True:
        num = random.normalvariate(mean, std_dev)
        if 0.9 <= num <= 1.1:
            return num

# Importing cell data
# coordinates of cms
pos = np.loadtxt(f"cell_data/capsule_{MP.capsule}/cm_cells.txt")
# weights for connections
weights = np.loadtxt(f"cell_data/capsule_{MP.capsule}/cmat_weight.txt")*CellModel.cell_height
# perimeter, area and volume of cells
per_area_vol = np.loadtxt(f"cell_data/capsule_{MP.capsule}/perimeters_area_volume.txt")
# binarized connectivity matrix of cells
cmat = np.loadtxt(f"cell_data/capsule_{MP.capsule}/cmat.txt")
cell_distances = distance.cdist(pos, pos, 'euclidean')


cell_num = len(pos)  # number of all cells
center_x = np.amin(pos[:,0]) + (np.amax(pos[:,0])-np.amin(pos[:,0]))/2.0
center_y = np.amin(pos[:,1]) + (np.amax(pos[:,1])-np.amin(pos[:,1]))/2.0
stimulated_cell = MP.find_stimulated_cell(pos, (center_x, center_y))  # number of cell for stimulation
MP.stimulated_cell = stimulated_cell
init_cells = get_initiatior_cells(stimulated_cell,
                                  cmat,
                                  mode=MP.initiator_cell_mode
                                  )


# Creates all cells and adds them to a list
cells = []
for i in range(cell_num):
    volume = per_area_vol[i, 2]*10**(-18)
    neighbours = get_neighbours(i, cmat, weights)
    cells.append(CellModel(i, volume,
                           cell_distances[i, stimulated_cell],
                           neighbours=neighbours,
                           parameters={
                                "ksscc": 1.20*random_normal_number(),
                                "Kpump": 0.5030,#heterogenost povzroča nestabilnost
                                "kryr": 16.04*random_normal_number(),
                                "kip3r3": 155.8*random_normal_number(),
                                "kdeg": 1.25*random_normal_number()
                            }
                           ))

# for cell in cells:
#     print(getattr(cell, "ksscc"))


# Actual simulation
simulation_time = 0.0
simulation_step = 0
while simulation_time < MP.tend:
    if simulation_step%100 == 0:
        print(simulation_time)
    calcium_coupling, ip3_coupling = calculate_ca_ip3_coupling(cells)
    atp_coupling = calculate_atp_coupling(
        cells,
        cell_distances,
        MP.stimulated_cell,
        simulation_time,
        mode=MP.diffusion_mode
    )
    #atp_coupling = [0 for i in range(len(cells))]
    for i, cell in enumerate(cells):
        cell.run_model_step(simulation_time,
                            simulation_step,
                            calcium_coupling[i],
                            ip3_coupling[i],
                            atp_coupling[i],
                            init_cells)
    simulation_step += 1
    simulation_time += MP.dt


# Collects all time series
results_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}'
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
results = get_time_series(cells, MP.sampling)
for key, item in results.items():
    np.savetxt(f'{results_folder}/{key}.txt', item, fmt='%.3lf')


#plot_movie_frames()
plot_activation_sequence()
plot_signaling_parameters(results, pos, stimulated_cell)
plot_time_series()
