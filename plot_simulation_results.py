"""
Plots all results
"""
import os
import numpy as np
import pandas as pd
from model.model_parameters import ModelParameters as MP
from plot_signaling_parameters import plot_signaling_parameters

results_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}'
pos = np.loadtxt(f"cell_data/capsule_{MP.capsule}/cm_cells.txt")
data = {}
data['response_times'] = np.loadtxt(f'{results_folder}/response_times.txt')
data['durations'] = np.loadtxt(f'{results_folder}/durations.txt')
data['fraction'] = np.loadtxt(f'{results_folder}/fraction.txt')
data['amplitudes'] = np.loadtxt(f'{results_folder}/amplitudes.txt')

plot_signaling_parameters(data, pos, np.argmin(data['response_times']))
