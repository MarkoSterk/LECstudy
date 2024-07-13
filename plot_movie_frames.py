"""
Plots movie frames
"""
# pylint: disable=C0103
import os
import cv2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from plot_configurations import PANEL_WIDTH
from model.model_parameters import ModelParameters as MP

def plot_movie_frames(): 
    """
    Plots movie frames and constructs the movie
    """
    stimulated_cell = MP.stimulated_cell
    pos = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cm_cells.txt')
    edgepoints = np.loadtxt(f'cell_data/capsule_{MP.capsule}/cell_edge_points_x_y_cellnum.txt')
    act_times = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/activation_times.txt')
    time = np.loadtxt(f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/time.txt')
    cell_num = len(pos)

    results_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/frames'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    dx=10.0
    dy=10.0
    x0=10.0
    y0=10.0

    activated = []
    frames_num = 0

    for i, t in enumerate(time):
        print(f'Frame: {i}')
        tru_col_map=[]
        for j in range(cell_num):
            if not np.isnan(act_times[j]):
                if t<act_times[j]:
                    tru_col_map.append('lightgray')
                else:
                    tru_col_map.append('red')
                    if j not in activated:
                        activated.append(j)
            else:
                tru_col_map.append('lightgray')

        fig=plt.figure(figsize=(2.0*PANEL_WIDTH, 2.0*PANEL_WIDTH))
        ax=fig.add_subplot(111)
        for j in range(cell_num):
            data=edgepoints[np.where(edgepoints[:,2]==j)][:,:2]
            polygon=Polygon(list(data), fc=tru_col_map[j], ec='black', zorder=0)
            ax.add_artist(polygon)
        ax.scatter(pos[:,0], pos[:,1], s=12, c='black', zorder=10)
        ax.scatter(pos[stimulated_cell,0], pos[stimulated_cell,1], s=30, color='black', marker='x')

        ax.plot([x0,x0+dx], [y0, y0], linewidth=2, color='black')
        ax.text(x0, y0-10, '10 um', fontsize=8)
        ax.plot([x0, x0], [y0,y0+dy], linewidth=2, color='black')
        ax.text(x0-10, y0, '10 um', fontsize=8, rotation='vertical')

        ax.set_axis_off()
        fig.savefig(
            f"results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/frames/frame_{i}.png",
            dpi=300, bbox_inches='tight')
        plt.close(fig)

        frames_num+=1
        if len(activated)==len(act_times[~np.isnan(act_times)]):
            break


    image_folder = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/frames'
    video_name = f'results/capsule_{MP.capsule}/init-cell_{MP.initiator_cell_mode}_diff-mode_{MP.diffusion_mode}/activation_sequence.avi'

    images = [f'frame_{i}.png' for i in range(frames_num)]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, _ = frame.shape

    video = cv2.VideoWriter(video_name, 0, 3, (width,height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()
