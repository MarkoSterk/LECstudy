"""
Entry point for simulation
"""
# pylint: disable=E1101
import numpy as np
from scipy.spatial import distance
from model_parameters import ModelParameters as MP
from cell_model import CellModel

#CHANGE MODEL PARAMETERS HERE
model_parameters: dict = dict(
    capsule=2,
    tend = 150.0, #end time
    time_of_stimulation = 15.0, #time of stimulation
    stim_dur = 3.0, #stimulation duration
    initiator_cell_mode = "point",
    cell_height=8.0,
    diffusion_mode = "decoupled",
    Cgjip3 = 0.05,
    L0init = 0.5, #0.2 concentration of secreted ATP by stimulated cell (fmol)
    char_dist= 30.0, #characteristic distance of ATP release decrease
    apyrase_deg = True, #True for apiraza, False (default) no apiraza
    apyrase_char_time = 1.0, #s
    Cth = 0.15, ##threshold amplitude of normalized calcium
    time_interval_for_slope = 4.0, # seconds
    slopeTh = 0.003,#0.01 ## uM/s
    Cth_act = 0.125 + 0.025,
    cell_heterogeneity = False,
    ca_noise_amp = 0.00,
    ip3_noise_amp = 0.00
)
#DON'T CHANGE AFTER THIS POINT

model: MP = MP(model_parameters)

pos = np.loadtxt(f"cell_data/capsule_{model.capsule}/cm_cells.txt")
weights = np.loadtxt(f"cell_data/capsule_{model.capsule}/cmat_weight.txt")*model.cell_height
per_area_vol = np.loadtxt(f"cell_data/capsule_{model.capsule}/perimeters_area_volume.txt")
cmat = np.loadtxt(f"cell_data/capsule_{model.capsule}/cmat.txt")
cell_distances = distance.cdist(pos, pos, 'euclidean')

cell_num = len(pos)  # number of all cells
center_x = np.amin(pos[:,0]) + (np.amax(pos[:,0])-np.amin(pos[:,0]))/2.0
center_y = np.amin(pos[:,1]) + (np.amax(pos[:,1])-np.amin(pos[:,1]))/2.0
stimulated_cell = MP.find_stimulated_cell(pos, (center_x, center_y))
model.stimulated_cell = stimulated_cell
init_cells = model.get_initiatior_cells(stimulated_cell,
                                  cmat,
                                  mode=model.initiator_cell_mode
                                  )
print(f"Initiator cells: {init_cells}")
cells: list[CellModel] = CellModel.generate_cells(model,
                                                  stimulated_cell,
                                                  init_cells,
                                                  pos,
                                                  per_area_vol,
                                                  cell_distances,
                                                  cmat,
                                                  weights,
                                                  model.cell_heterogeneity)

time: list[int] = model.run_simulation(cells, cell_distances)

ca_ts, ca_bin_ts, ip3_ts, atp_ts = MP.extract_time_series_data(cells)
act_times = MP.extracts_activation_times(cells)

model.save_ts_data(ca_ts, ca_bin_ts, ip3_ts, atp_ts, act_times)
model.plot_activation_sequence(model)
model.plot_time_series(ca_ts)
