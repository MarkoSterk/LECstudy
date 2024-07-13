"""
Calculates cellular activation times
"""
import numpy as np
from model.model_parameters import ModelParameters as MP


def calculate_activation_times():
    """
    Calculates activation times (slope) of cells
    """

    def cell_activation_time(series: np.array):
        """
        Calculates act. time of cell
        """
        time_step = MP.record_every * MP.dt ## time between two data points
        points_for_slope = int(MP.time_interval_for_slope/time_step)
        left_points = int(points_for_slope/2)
        right_points = points_for_slope-left_points+1
        for t in range(left_points, len(series)-right_points, 1):
            #slope = (series[t+right_points]-series[t-left_points])/points_for_slope
            avg_sig = np.average(series[t-left_points:t+right_points])
            #if slope >= MP.slopeTh and series[t]>=MP.Cth:
            if avg_sig>MP.Cth:
                return t*time_step
        return np.nan

    time_series = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/ca_time_series_i.txt')
    cell_num = len(time_series[0])

    slope_act_times = np.full(cell_num, np.nan)
    norm_time_series = np.zeros(time_series.shape, float)
    for i in range(cell_num):
        norm_time_series[:,i] = (time_series[:,i] - np.min(time_series[:,i]))/(np.max(time_series[:,i]) - np.min(time_series[:,i]))
        slope_act_times[i] = cell_activation_time(norm_time_series[:,i])

    np.savetxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/slope_activation_times.txt', slope_act_times, fmt="%.2lf")
    np.savetxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/norm_time_series.txt', norm_time_series, fmt="%.2lf")

if __name__ == "__main__":
    calculate_activation_times()
