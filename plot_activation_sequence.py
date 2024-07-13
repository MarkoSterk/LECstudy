"""
Plots cell activation sequence
"""
# pylint: disable=C0103

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from plot_configurations import PANEL_WIDTH
from model.model_parameters import ModelParameters as MP

def plot_activation_sequence():
    """
    Plots activation sequence of cells
    """
    stimulated_cell = MP.stimulated_cell
    pos = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cm_cells.txt')
    edgepoints = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cell_edge_points_x_y_cellnum.txt')
    act_times = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/act_times.txt')
    cell_num = len(pos)

    vmin=np.nanmin(act_times, axis=None)
    vmax=np.nanmax(act_times, axis=None)

    tru_col_map=[]
    for i in range(cell_num):
        if np.isnan(act_times[i]):
            tru_col_map.append(plt.cm.gray(0.4))
        else:
            tru_col_map.append(plt.cm.jet_r(float(act_times[i]-vmin)/float(vmax-vmin)))

    dx=10.0
    dy=10.0
    x0=10.0
    y0=10.0

    fig=plt.figure(figsize=(2.0*PANEL_WIDTH, 2.0*PANEL_WIDTH))
    ax=fig.add_subplot(111)
    for i in range(cell_num):
        data=edgepoints[np.where(edgepoints[:,2]==i)][:,:2]
        polygon=Polygon(list(data), fc=tru_col_map[i], ec='black', zorder=0)
        ax.add_artist(polygon)
    ax.scatter(pos[:,0], pos[:,1], s=12, c='black', zorder=10)
    ax.scatter(pos[stimulated_cell,0], pos[stimulated_cell,1], s=30, color='black', marker='x')

    ax.plot([x0,x0+dx], [y0, y0], linewidth=2, color='black')
    ax.text(x0, y0-10, '10 um', fontsize=8)
    ax.plot([x0, x0], [y0,y0+dy], linewidth=2, color='black')
    ax.text(x0-10, y0, '10 um', fontsize=8, rotation='vertical')

    ax.set_axis_off()
    fig.savefig(
        f"results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/activation_sequence.png",
        dpi=600, bbox_inches='tight')
    plt.close(fig)
