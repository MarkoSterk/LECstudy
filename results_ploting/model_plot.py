""" Plots for model data """

import numpy as np
import matplotlib.pyplot as plt

from configs import get_dimensions

data_folder: str = "data/model"
results_folder: str ="results/model"

amps: np.ndarray = np.loadtxt(f"{data_folder}/amplitude.txt")
durs: np.ndarray = np.loadtxt(f"{data_folder}/duration.txt")
resp_times: np.ndarray = np.loadtxt(f"{data_folder}/responsetime.txt")
fracs: np.ndarray = np.loadtxt(f"{data_folder}/fraction.txt")


fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.plot(amps[:,0], amps[:,1], c='black', linewidth=1.0, marker="o", label="control")
ax.plot(amps[:,0], amps[:,2], c='red', linewidth=1.0, marker="o", label="apyrase")
ax.plot(amps[:,0], amps[:,3], c='green', linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Signal amplitude ($\mathrm{\mu}$M)")
ax.set_xlim(15, 80)
ax.set_ylim(0.15, 0.32)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/amplitudes.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.plot(durs[:,0], durs[:,1], c='black', linewidth=1.0, marker="o", label="control")
ax.plot(durs[:,0], durs[:,2], c='red', linewidth=1.0, marker="o", label="apyrase")
ax.plot(durs[:,0], durs[:,3], c='green', linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Signal duration (s)")
ax.set_xlim(15, 80)
ax.set_ylim(24, 36)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/durations.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)

fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.plot(resp_times[:,0], resp_times[:,1], c='black', linewidth=1.0, marker="o", label="control")
ax.plot(resp_times[:,0], resp_times[:,2], c='red', linewidth=1.0, marker="o", label="apyrase")
ax.plot(resp_times[:,0], resp_times[:,3], c='green', linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Response time (s)")
ax.set_xlim(15, 80)
#ax.set_ylim(24, 36)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/response_times.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)

fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.plot(fracs[:,0], fracs[:,1], c='black', linewidth=1.0, marker="o", label="control")
ax.plot(fracs[:,0], fracs[:,2], c='red', linewidth=1.0, marker="o", label="apyrase")
ax.plot(fracs[:,0], fracs[:,3], c='green', linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Fraction of active cells")
ax.set_xlim(15, 80)
ax.set_ylim(0, 1.05)
ax.legend(loc="lower left")
fig.savefig(f"{results_folder}/fractions.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)
