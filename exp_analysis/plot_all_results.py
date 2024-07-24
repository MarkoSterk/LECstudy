"""
Plots results (panels) with all three experiment types
"""
#pylint: disable=E0401,W0611
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore

# Matplotlib configuration
SMALL_SIZE = 6
MEDIUM_SIZE = 6
BIGGER_SIZE = 8

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

PANEL_WIDTH = 6.5
PANEL_HEIGHT = 6.5
CONVERSION = 2.54
PANEL_WIDTH = PANEL_WIDTH/CONVERSION
PANEL_HEIGHT = PANEL_HEIGHT/CONVERSION
###########################################

exp_types = ["control", "apiraza", "cbx"]

files = [("time_to_peak", "Time to peak (s)", "upper right"),
        ("response_times", "Response time (s)", "upper left"),
        ("amplitudes", r'Amplitude $([\mathrm{Ca}^{2+}])$', "upper right"),
        ("durations", "Signal duration (s)", "upper right"),
        ("fractions", "Fraction of activated cells", "lower left")]

colors = {
    exp_types[0]: "black",
    exp_types[1]: "red",
    exp_types[2]: "green"
}

for filename, label, leg_loc in files:
    fig = plt.figure(figsize=(PANEL_WIDTH, PANEL_HEIGHT))
    ax = fig.add_subplot(1,1,1)
    for exp_type in exp_types:
        results_folder = f'Experimental_data/{exp_type}/results'
        path = f'{results_folder}/results_plot_{exp_type}_{filename}.txt'
        data = np.loadtxt(path)
        ax.errorbar(data[0,:], data[1,:], data[2,:], c=colors[exp_type], linewidth=2, label=exp_type)

    ax.legend(loc=leg_loc)
    ax.set_xlabel("Distance from stimulation (um)")
    ax.set_ylabel(label)
    ax.set_xlim(10, 100)
    fig.savefig(f"Experimental_data/{filename}.png", pad_inches=0.01, bbox_inches='tight', dpi=600)
    plt.close(fig)
