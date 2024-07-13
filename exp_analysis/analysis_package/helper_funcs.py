"""
Helper functions
"""
import sys
from typing import List, Tuple
import numpy as np
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

def smooth_ts(time_series: np.ndarray, N: int) -> np.ndarray:
    """
    Smooths time series with N points
    """
    right = int(N/2)
    left = N - right
    smooth_time_series = np.zeros(time_series.shape, float)
    for i in range(len(time_series)):
        if left<=i<=len(time_series)-right:
            smooth_time_series[i] = np.average(time_series[i-left:i+right])
        elif i<left:
            smooth_time_series[i] = np.average(time_series[:i+right])
        elif i>len(time_series)-right:
            smooth_time_series[i] = np.average(time_series[i-left:])
    return smooth_time_series

def normalize_ts(time_series, amp=1, offset=0) -> np.ndarray:
    """
    Normalizes time series between offset and amp+offset
    :param: time_series - time series for normalization
    :param: amp - amplitude of normalization default=1
    :param: offset - offset of normalized ts default=0
    """
    vmin_data = np.min(time_series[:])
    vmax_data = np.max(time_series[:])
    new_ts = amp*(time_series[:] - vmin_data)/(vmax_data-vmin_data)+offset
    return new_ts

def calculate_slope(time_series: np.ndarray, step: int, points: int) -> float:
    """
    Calculates current slope
    """
    slope = (time_series[step+points]-time_series[step])/points
    return slope
    

def calculate_act_times(data: np.ndarray, configs: dict) -> np.ndarray:
    """
    Calculates activation times of cells
    """
    points_for_slope = configs['points_for_slope']
    slope_th = configs['slope_th']
    amp_th = configs['amp_th']
    stim_time = configs['stim_time']

    act_times = np.zeros(len(data[0]), float)
    act_times[:] = np.nan

    for i in range(len(data[0])):
        for j in range(int(stim_time), len(data)-points_for_slope, 1):
            current_slope = calculate_slope(data[:,i], j, points_for_slope)
            if np.isnan(act_times[i]) and current_slope>slope_th and data[j,i]>amp_th:
                act_times[i] = int(j+1)
    return act_times

def calculate_peak_times(data: np.ndarray, stim_time: int) -> np.ndarray:
    """
    Finds peak times
    """         
    peak_times = np.zeros(len(data[0]))
    peak_amps = np.zeros(len(data[0]))
    for i in range(len(data[0])):
        peak_times[i] = np.argmax(data[stim_time:,i]) + stim_time
        peak_amps[i] = data[int(peak_times[i]),i]
    return peak_times, peak_amps


def calculate_deactivation_times(data: np.ndarray,
                                 act_times: np.ndarray,
                                 peak_times: np.ndarray) -> List[np.ndarray]:
    """
    Calculates deactivation times
    """
    deact_times = np.zeros(len(data[0]), float)
    deact_amps = np.zeros(len(data[0]), float)
    for i in range(len(data[0])):
        deact_times[i] = len(data)
        deact_amps[i] = data[-1,i]
        if not np.isnan(act_times[i]) and not np.isnan(peak_times[i]):
            deact_amp = data[int(peak_times[i]),i] - (data[int(peak_times[i]),i] - data[int(act_times[i]),i])/2
            for j in range(int(peak_times[i]), len(data), 1):
                deact_times[i] = j
                deact_amps[i] = deact_amp
                if data[j,i]<=deact_amp:
                    break
    return deact_times, deact_amps

def calculate_durations(act_times: np.ndarray, deact_times: np.ndarray) -> np.ndarray:
    """
    Calculates durations of cellular activities
    """
    durations = np.zeros(len(act_times), float)
    durations[:] = deact_times - act_times
    return durations

def cell_activity(data: np.ndarray, raw_series: np.ndarray, configs: dict) -> Tuple[np.ndarray]:
    """
    Determines time interval of cellular activation
    """
    amp_fact = configs['AMP_FACT']
    distance=configs['DISTANCE']
    width = configs['WIDTH']
    prominence=configs['PROMINENCE']
    rel_height = configs['REL_HEIGHT']
    stim_time = configs['STIM_TIME']

    act_times = np.zeros(len(data[0]))
    response_times = np.zeros(len(data[0]))
    fraction = np.zeros(len(data[0]))
    amplitudes = np.zeros(len(data[0]))
    durations = np.zeros(len(data[0]))

    deact_times = np.zeros(len(data[0]))
    peak_times = np.zeros(len(data[0]))

    act_times[:] = np.nan
    response_times[:] = np.nan
    amplitudes[:] = np.nan
    durations[:] = np.nan

    peak_times[:] = np.nan
    deact_times[:] = np.nan

    for i in range(len(data[0])):
        peaks, _ = find_peaks(data[:, i],
                            height=amp_fact,
                            distance=distance,
                            width=width,
                            prominence=prominence)

        _, _, left_ips, right_ips = peak_widths(data[:, i], peaks, rel_height=rel_height)
        if len(peaks)==1:
            act_times[i] = int(left_ips[0])
            durations[i] = int(right_ips[0]) - act_times[i]
            fraction[i] = 1
            peak_times[i] = int(peaks[0])
            amplitudes[i] = raw_series[int(peaks[0]), i]-raw_series[int(left_ips[0]),i]
            deact_times[i] = int(right_ips[0])

    for i in range(len(data[0])):
        if not np.isnan(act_times[i]):
            response_times[i] = act_times[i] - np.nanmin(act_times)

    return act_times, response_times, fraction, amplitudes, durations, deact_times, peak_times


def cell_activity_manual(data: np.ndarray):
    """
    Manual setting of activation and deactivation times for cells
    """
    cell_data = {
        'act_times': [np.nan for _ in range(len(data[0]))],
        'deact_times': [np.nan for _ in range(len(data[0]))]
    }
    
    click_params = {
        'next_cell': None,
        'key_d': False,
        'key_a': False
    }

    current_cell = 0
    cell_num = len(data[0])
    time = range(len(data))
    def on_click(event, cell):
        """
        on-click event handler
        """
        if event.inaxes == ax:
            if event.button==plt.MouseButton.RIGHT:
                if click_params['key_a']:
                    cell_data['act_times'][cell]=event.xdata
                    print(f'Activation time of cell {cell}/{cell_num}:', cell_data['act_times'][cell])
                if click_params['key_d']:
                    cell_data['deact_times'][cell]=event.xdata
                    print(f'Deactivation time of cell {cell}/{cell_num}:', cell_data['deact_times'][cell])

    def on_key_press(event, cell):
        """
        On-press event handler
        """
        if str(event.key) == 'escape':
            ###Exits the whole thing
            sys.exit()
        if str(event.key) == 'left':
            ###If the current cell is not the first (0)
            ###the cell counter is turned back. It gets turned back by 2
            ###because it gets increased again by one in the main while
            ### loop after the figure closes
            if cell > 0:
                click_params['next_cell'] = cell-2
                plt.close()
        if str(event.key) == 'right':
            ###If the response time of this cell was not already set it sets the time to Nan
            ### and continous to the next cell
            ### else it just closes the figure and continous on.
            click_params['next_cell'] = cell
            plt.close()
        if str(event.key) in ['d', 'D']:
            click_params['key_d'] = True
        if str(event.key) in ['a', 'A']:
            click_params['key_a'] = True
    def on_key_release(event):
        if str(event.key) in ['d', 'D']:
            click_params['key_d'] = False
        if str(event.key) in ['a', 'A']:
            click_params['key_a'] = False

    plt.rcParams['backend'] = 'TkAgg'
    plt.rcParams["figure.figsize"] = [8, 5]
    plt.rcParams["figure.autolayout"] = True

    while current_cell < cell_num:
        fig = plt.figure()
        fig.set_tight_layout(False)
        fig.canvas.mpl_connect('button_release_event', lambda event: on_click(event, current_cell))
        fig.canvas.mpl_connect('key_press_event', lambda event: on_key_press(event, current_cell))
        fig.canvas.mpl_connect('key_release_event', on_key_release)
        ax = fig.add_subplot(1,1,1)
        # pylint: disable-next=W0612
        cursor = Cursor(ax, horizOn=True, vertOn=True, color='green', linewidth=1.0, useblit=True)
        fig.suptitle(f'Cell {current_cell}')
        ax.plot(time, data[:,current_cell], linewidth=0.4, c='gray')
        ax.set_xlabel('time (s)')
        ax.set_ylabel('Cell signal (a.u.)')
        plt.get_current_fig_manager().window.wm_geometry("+10+10") # move the window
        plt.show()
        current_cell = click_params['next_cell'] + 1
    plt.close()

    return cell_data
