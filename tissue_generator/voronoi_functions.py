"""
Voronoi methods for constructing the 
cellular tissue
"""
#pylint: disable=E0611
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import Voronoi, voronoi_plot_2d
from tissue_generator.helper_functions import (generate_random_coordinates,
                              distances_between_neighbours,
                              create_bin_conn_matrix)

def shoelace_formula(polygon_points: np.ndarray):
    """
    Calculates surface area of polygon
    with polygon vertices coordinates
    """
    x = polygon_points[:, 0]
    y = polygon_points[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

def generate_voronoi_network(pos: np.ndarray) -> tuple[Voronoi,
                                                       np.ndarray,
                                                       np.ndarray]:
    """
    Generates voronoi netowkr from provided seed coordinates (pos)
    Returns voronoi points, regions and vertices
    """
    voronoi = Voronoi(pos)
    vor_points: np.ndarray = np.array(voronoi.points)
    vor_vertices: np.ndarray = np.array(voronoi.vertices)

    return voronoi, vor_points, vor_vertices

def find_common_edge_vertices(polygon1, polygon2)-> np.ndarray:
    """
    Finds common edge vertices of two polygons
    """
    # Step 1: View the arrays as structured arrays
    polygon1_view = polygon1.view([('', polygon1.dtype)] * polygon1.shape[1])
    polygon2_view = polygon2.view([('', polygon2.dtype)] * polygon2.shape[1])

    # Step 2: Use numpy.intersect1d to find shared vertices
    shared_vertices = np.intersect1d(polygon1_view, polygon2_view)

    # Step 3: Convert the result back to the original 2D coordinates
    shared_vertices = shared_vertices.view(polygon1.dtype).reshape(-1, 2)
    return shared_vertices

def find_connected_cells(cell_data: dict[str, list[int]|float|np.ndarray]
                        ) -> dict[str, list[int]|float|np.ndarray]:
    """
    Finds all neighbouring cells and calculates their shared border size
    """
    for cell_i, data_i in cell_data.items():
        verts_i = data_i["vertices_coord"]
        for cell_j, data_j in cell_data.items():
            verts_j = data_j["vertices_coord"]
            if cell_i != cell_j:
                shared_edge_vertices: np.ndarray = find_common_edge_vertices(verts_i, verts_j)
                if len(shared_edge_vertices) != 0:
                    x1, y1 = shared_edge_vertices[0]
                    x2, y2 = shared_edge_vertices[1]
                    length: float = np.sqrt((x1-x2)**2+(y1-y2)**2)
                    cell_data[cell_i]["neighbours"][cell_j] = length
    return cell_data

#pylint: disable-next=R0913,R0914
def crop_voronoi_network(vor_regions: list[list[int]], vor_vertices: np.ndarray,
                         min_coord: float, max_coord: float, cell_height: float,
                         crop: float = 20) -> dict[str, list[int]|float|np.ndarray]:
    """
    Filters voronoi network to remove all regions with vertices outside the
    threshold coordinates. The default crop value is 20 (um).

    Returns a dictionary of dictionaries with keys of all valid regions,
    vertices numbers for valid regions, vertices coordinates for valid regions,
    and center of mass, surface area and volume for each valid region.
    """
    #minxy: float = min_coord+crop
    #maxxy: float = max_coord-crop
    ok_regions: dict[str, list[int]|np.ndarray] = {}
    cell_count = 0
    for _, region in enumerate(vor_regions):
        check: int = 0
        verts = []
        for j in region:
            x, y = vor_vertices[j, 0], vor_vertices[j, 1]
            verts.append((x, y))
            if((min_coord+crop <= x <=max_coord-crop) and (min_coord+crop <= y <= max_coord-crop)):
                check+=1
        if check == len(region) and len(region)>2:
            verts = np.array(verts)
            area: float = shoelace_formula(verts)
            volume: float = area * cell_height
            ok_regions[cell_count] = {
                "index": cell_count,
                "vertices_coord": verts,
                "cm_xy": (np.average(verts[:,0]), np.average(verts[:,1])),
                "surface_area": area,
                "volume": volume,
                "neighbours": {}
            }
            cell_count+=1
    ok_regions = find_connected_cells(ok_regions)
    return {"cells": ok_regions}


def plot_voronoi_network(voronoi: Voronoi, show_vertices: bool = False):
    """
    Plots the voronoi network and returns matplotlib figure object
    """
    figure = plt.figure(figsize=(5,5))
    axes = figure.add_subplot(1,1,1)
    voronoi_plot_2d(voronoi, axes, show_vertices=show_vertices)
    return figure

def plot_polygon_to_figure(axes, polygon_vertices: np.ndarray, cm_xy: tuple[float], **kwargs):
    """
    Plots a polygon to the provided figure axes
    """
    polygon = Polygon(polygon_vertices, closed=True, edgecolor="black",
                      facecolor="lightgray", linewidth=2, **kwargs)
    axes.add_patch(polygon)
    axes.scatter(cm_xy[0], cm_xy[1], s=4, c="black")


def plot_area_and_volume(cell_data: dict[str, list[int]|float|np.ndarray]):
    """
    Creates a figures with two subplots. One for areas and
    for volumes of cells.
    """
    volumes = [v["volume"] for v in cell_data.values()]
    areas = [v["surface_area"] for v in cell_data.values()]
    distances = distances_between_neighbours(cell_data).flatten()
    degree_dist = np.array([len(v["neighbours"].keys()) for v in cell_data.values()])

    figure, axes = plt.subplots(2,2, figsize=(10, 10))
    #Histogram chart for surface area
    axes[0,0].hist(areas, bins=20, color='blue',
                 edgecolor='black', linewidth=1.2,
                 rwidth=0.9, density=True)
    axes[0,0].set_title('Surface Areas Distribution')
    axes[0,0].set_ylabel('Frequency')
    axes[0,0].set_xlabel('Surface Area')
    # Histogram chart for volume
    axes[0,1].hist(volumes, bins=20, color='green',
                 edgecolor='black', linewidth=1.2,
                 rwidth=0.9, density=True)
    axes[0,1].set_title('Volumes Distribution')
    axes[0,1].set_ylabel('Frequency')
    axes[0,1].set_xlabel('Volume')
    # Histogram chart of intercellular distances
    axes[1,0].hist(distances, bins=20, color='gray',
                 edgecolor='black', linewidth=1.2,
                 rwidth=0.9, density=True)
    axes[1,0].set_title('Distance Distribution')
    axes[1,0].set_ylabel('Frequency')
    axes[1,0].set_xlabel('Distance')
    # Histogram chart of intercellular distances
    axes[1,1].hist(degree_dist, bins=8, color='cyan',
                 edgecolor='black', linewidth=1.2,
                 rwidth=0.9, density=True)
    axes[1,1].set_title('Degree Distribution')
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].set_xlabel('Degree')

    return figure

#pylint: disable-next= R0913,R0914
def create_cell_network(capsule: int, threshold_dist: float, real_intercellular_dist: float,
                        crop: float, number_of_cells: int = 700,
                        min_xy: float = 0.0, max_xy: float = 200):
    """
    Creates network and saves it to results folder
    returns model
    """
    results_path: str = f"cell_data/capsule_{capsule}"
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    positions: np.ndarray = generate_random_coordinates(number_of_cells,
                                                  threshold_dist,
                                                  min_xy,
                                                  max_xy)
    vor, _, vertices = generate_voronoi_network(positions)
    fig = plot_voronoi_network(vor)

    fig.savefig(f"{results_path}/voronoi_network.png",
                dpi=600, bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

    model = crop_voronoi_network(vor.regions, vertices, min_xy, max_xy,
                                           cell_height=real_intercellular_dist, crop=crop)
    model["bin_conn_mat"] = create_bin_conn_matrix(model["cells"])
    with open(f"{results_path}/model_data.pkl", "wb") as file:
        pickle.dump(model, file)

    cropped_regions = model["cells"]

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(1,1,1)
    for _, vor_region in cropped_regions.items():
        plot_polygon_to_figure(ax, vor_region["vertices_coord"], vor_region["cm_xy"])
    ax.set_axis_off()
    fig.savefig(f"{results_path}/cropped_voronoi.png", dpi=600,
                bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)

    fig = plot_area_and_volume(cropped_regions)
    fig.subplots_adjust(wspace=0.3, hspace=0.2)
    fig.savefig(f"{results_path}/cell_area_volume_intercelldist_degree_distribution.png",
                dpi=600, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

    print(f"SECCUESS: successfully created tissue network with {len(cropped_regions)} cells.")
    return model

if __name__ == "__main__":

    CAPSULE: int = 1
    NUMBER_OF_CELLS: int = 700
    MIN_XY: float = 0.0
    MAX_XY: float = 200
    THRESHOLD_DIST: float = 6.0
    REAL_INTERCELL_DIST: float = 8.0 #fine-tuning of cell valumes
    CROP: float = 60.0

    create_cell_network(CAPSULE, THRESHOLD_DIST, REAL_INTERCELL_DIST, CROP,
                        number_of_cells=NUMBER_OF_CELLS, min_xy=MIN_XY, max_xy=MAX_XY)
