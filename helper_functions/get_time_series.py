"""
Collects all time series
"""

import numpy as np
from typing import List
from model.cell_model import CellModel

def get_time_series(cells: List[CellModel], sampling: float) -> dict[np.ndarray]:
    """
    Gets all time series from cell data
    """
    def calculate_signal_duration(act_time: float, series: np.array, sampling: float):
        """
        Calculates total signal duration
        """
        if np.isnan(act_time):
            return np.nan
        for t in reversed(range(int(act_time*sampling), len(series), 1)):
            if series[t]>=(((np.amax(series)-series[int(act_time*sampling)])/2) + series[int(act_time*sampling)]):
                osc_dur = (t - int(act_time*sampling))/sampling
                return osc_dur
        return np.nan
    
    def calculate_amp_increase(act_time: float, series: np.ndarray, sampling: float):
        """
        Calculates the amplitude increase of cell
        """
        if np.isnan(act_time):
            return np.nan
        return np.amax(series) - series[int(act_time*sampling)]

    activation_times = np.array([np.nan for _ in range(len(cells))])
    response_times = np.array([np.nan for _ in range(len(cells))])
    deact_times = np.array([np.nan for _ in range(len(cells))])
    durations = np.array([np.nan for _ in range(len(cells))])
    amplitudes = np.array([np.nan for _ in range(len(cells))])
    fraction = np.array([0 for _ in range(len(cells))])
    ca_time_series = np.zeros(
                (len(cells[0].calcium_time_series), len(cells)), float)
    ip3_time_series = np.zeros(
                (len(cells[0].calcium_time_series), len(cells)), float)
    atp_time_series = np.zeros(
                (len(cells[0].calcium_time_series), len(cells)), float)
    time = cells[0].simulation_time

    for i, cell in enumerate(cells):
        ca_time_series[:,i] = cell.calcium_time_series[:]
        ip3_time_series[:,i] = cell.ip3_time_series[:]
        atp_time_series[:,i] = cell.atp_time_series[:]
        if cell.time_of_activation:
            activation_times[i] = cell.time_of_activation
            fraction[i] = 1
            durations[i] = calculate_signal_duration(activation_times[i],
                                                     ca_time_series[:,i],
                                                     sampling)
            deact_times[i] = activation_times[i] + durations[i]
            amplitudes[i] = calculate_amp_increase(activation_times[i],
                                                   ca_time_series[:,i],
                                                   sampling)
            

    for i in range(len(cells)):
        if not np.isnan(activation_times[i]):
            response_times[i] = activation_times[i] - np.nanmin(activation_times)

    return {
            'time': time,
            'act_times': activation_times,
            'deact_times': deact_times,
            'response_times': response_times,
            'fraction': fraction,
            'durations': durations,
            'amplitudes': amplitudes,
            'ca_time_series': ca_time_series,
            'ip3_time_series': ip3_time_series,
            'atp_time_series': atp_time_series
            }
