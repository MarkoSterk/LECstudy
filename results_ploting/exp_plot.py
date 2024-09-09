""" Plots for experimental data """

import numpy as np
import matplotlib.pyplot as plt

from configs import get_dimensions

data_folder: str = "data/experiment"
results_folder: str ="results/experiment"

control: np.ndarray = np.loadtxt(f"{data_folder}/all_values_dist_frac_resptime_dur_amp_value_std_control.txt")
apyrase: np.ndarray = np.loadtxt(f"{data_folder}/all_values_dist_frac_resptime_dur_amp_value_std_apiraza.txt")
cbx: np.ndarray = np.loadtxt(f"{data_folder}/all_values_dist_frac_resptime_dur_amp_value_std_cbx.txt")

##Fractions
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,1], control[:,2], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,1], apyrase[:,2], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:,0], cbx[:,1], cbx[:,2], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Fraction of active cells")
#ax.set_xlim(15, 80)
ax.set_ylim(0, 1.05)
ax.legend(loc="lower left")
fig.savefig(f"{results_folder}/fractions.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)

##Response times
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,3], control[:,4], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,3], apyrase[:,4], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,3], cbx[:2,4], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Response time (s)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/response_times.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


##Signal durations
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,5], control[:,6], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,5], apyrase[:,6], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,5], cbx[:2,6], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Signal duration (s)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/durations.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


##Signal amplitudes
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,7], control[:,8], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,7], apyrase[:,8], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,7], cbx[:2,8], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Relative amplitude ($\mathrm{\mu}$M)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/amplitudes.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


##POOLED DATA PLOT
control: np.ndarray = np.loadtxt(f"{data_folder}/pooled_data_dist_fracs_resptimes_durs_amps_avg_stderr_control.txt")
apyrase: np.ndarray = np.loadtxt(f"{data_folder}/pooled_data_dist_fracs_resptimes_durs_amps_avg_stderr_apiraza.txt")
cbx: np.ndarray = np.loadtxt(f"{data_folder}/pooled_data_dist_fracs_resptimes_durs_amps_avg_stderr_cbx.txt")

##Fractions
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.plot(control[:,0], control[:,1], c='black',
            linewidth=1.0, marker="o", label="control")
ax.plot(apyrase[:,0], apyrase[:,1], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.plot(cbx[:,0], cbx[:,1], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Fraction of active cells")
#ax.set_xlim(15, 80)
ax.set_ylim(0, 1.05)
ax.legend(loc="lower left")
fig.savefig(f"{results_folder}/fractions_pooled.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)

##Response times
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,2]-np.nanmin(control[:2,2]), control[:,3], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,2]-np.nanmin(apyrase[:2,2]), apyrase[:,3], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,2]-np.nanmin(cbx[:2,2]), cbx[:2,3], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Response time (s)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/response_times_pooled.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


##Signal durations
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,4], control[:,5], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,4], apyrase[:,5], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,4], cbx[:2,5], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Signal duration (s)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/durations_pooled.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)


##Signal amplitudes
fig = plt.figure(figsize=get_dimensions())
ax = fig.add_subplot(1,1,1)
ax.errorbar(control[:,0], control[:,6], control[:,7], c='black',
            linewidth=1.0, marker="o", label="control")
ax.errorbar(apyrase[:,0], apyrase[:,6], apyrase[:,7], c='red',
            linewidth=1.0, marker="o", label="apyrase")
ax.errorbar(cbx[:2,0], cbx[:2,6], cbx[:2,7], c='green',
            linewidth=1.0, marker="o", label="cbx")
ax.set_xlabel(r"Distance from stimulation ($\mathrm{\mu}$m)")
ax.set_ylabel("Relative amplitude ($\mathrm{\mu}$M)")
#ax.set_xlim(15, 80)
#ax.set_ylim(0, 1.05)
ax.legend(loc="upper right")
fig.savefig(f"{results_folder}/amplitudes_pooled.png", dpi=600,
            bbox_inches="tight", pad_inches=0.01)
plt.close(fig)
