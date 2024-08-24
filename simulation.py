"""
Entry point for simulation
"""
# pylint: disable=E1101
import numpy as np
from scipy.spatial import distance
from lec_model import ModelParameters as MP
from lec_model import CellModel
from tissue_generator import load_capsule_data

#CHANGE MODEL PARAMETERS HERE
model_parameters: dict = dict(
    capsule=0,
    tend = 150.0, #end time
    time_of_stimulation = 15.0, #time of stimulation
    stim_dur = 3.0, #stimulation duration
    initiator_cell_mode = MP.INITMode.POINT,
    cell_height=8.0,
    diffusion_mode = MP.ATPMode.DECOUPLED,
    Cgjip3 = 0.05,
    L0init = 0.2, #0.2 concentration of secreted ATP by stimulated cell (fmol)
    char_dist= 30.0, #characteristic distance of ATP release decrease
    apyrase_deg = False, #True for apiraza, False (default) no apiraza
    apyrase_char_time = 1.0, #s
    Cth = 0.15, ##threshold amplitude of normalized calcium
    time_interval_for_slope = 4.0, # seconds
    slopeTh = 0.002,#0.01 ## uM/s
    Cth_act = 0.125 + 0.025,
    cell_heterogeneity = True,
    ca_noise_amp = 0.00,
    ip3_noise_amp = 0.00
)
#DON'T CHANGE AFTER THIS POINT

model: MP = MP(model_parameters)

capsule_data = load_capsule_data(model.capsule)

pos = np.array([v["cm_xy"] for v in capsule_data["cells"].values()])
weights = capsule_data["weights"]*capsule_data["cell_height"]
volumes = np.array([v["volume"] for v in capsule_data["cells"].values()])
cmat = capsule_data["bin_conn_mat"]
cell_distances = distance.cdist(pos, pos, 'euclidean')

center_x = np.amin(pos[:,0]) + (np.amax(pos[:,0])-np.amin(pos[:,0]))/2.0
center_y = np.amin(pos[:,1]) + (np.amax(pos[:,1])-np.amin(pos[:,1]))/2.0
stimulated_cell = MP.find_stimulated_cell(pos, (center_x, center_y))
model.stimulated_cell = stimulated_cell
init_cells = model.get_initiatior_cells(stimulated_cell,
                                  cmat,
                                  mode=model.initiator_cell_mode
                                  )
print(f"Stimulated cells: {init_cells}")
cells: list[CellModel] = CellModel.generate_cells(model,
                                                  stimulated_cell,
                                                  init_cells,
                                                  pos,
                                                  volumes,
                                                  cell_distances,
                                                  cmat,
                                                  weights,
                                                  model.cell_heterogeneity)

model.run_simulation(cells, cell_distances)

ca_ts, ip3_ts, atp_ts, jgjca_ts, jgjip3_ts = MP.extract_time_series_data(cells)
act_times, act_amps, act_frames = MP.extracts_activation_times(cells)


model.plot_activation_sequence(pos, act_times, capsule_data)

durations, resp_times, amps, deact_frames = model.calculate_activity_params(ca_ts[:,1:],
                                                              act_times,
                                                              act_amps,
                                                              act_frames,
                                                              ca_ts[:,0].flatten())
ca_bin_ts = model.extract_bin_signals(ca_ts, act_frames, deact_frames)
model.save_ts_data(ca_ts, ca_bin_ts, ip3_ts, atp_ts, act_times, jgjca_ts, jgjip3_ts)
model.plot_time_series(ca_ts, ca_bin_ts)
group_distances, fractions_act_cells, clustered_durations, clustered_amps, clustered_resp_times = model.plot_activity_params(durations, resp_times, amps, pos, 5)
model.save_activity_params(durations, resp_times, amps, fractions_act_cells)
model.save_clustered_activity_params(group_distances, fractions_act_cells, clustered_durations, clustered_amps, clustered_resp_times)

stim_cell_neighbours = [int(neigh) for neigh in np.nonzero(cmat[:, stimulated_cell])[0]]

model_parameters["stimulated_cell"] = stimulated_cell
model_parameters["stim_cell_neighbours"] = stim_cell_neighbours
model.save_model_parameters(model_parameters)
model.plot_gj_currents(stimulated_cell, stim_cell_neighbours, jgjca_ts, jgjip3_ts)
