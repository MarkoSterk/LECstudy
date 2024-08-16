"""
Constructs cellular networks with realistic 
cellular properties and intercellular connections
"""

from tissue_generator import create_cell_network

#All results are saved to the "cell_data/capsule_<CAPSULE>" folder
#Because of the randomness of the cell generator it is sometimes not possible to generate
#all cells (NUMBER_OF_CELLS). In this case the generator exits with a warning but
#continues with the number of cells it managed to generate. The warning is printed
#to the console. Often only a few cells are missing and the results are still ok
#
#
#CROP: the crop parameter influences the final number of cells. The higher this parameter
#the fewer cells are in the final model. This is because only cells whos edge vertices
#fall into the interval (MIN_XY+CROP) and (MAX_XY-CROP) are kept. This cropping is neccessary
#because the Voronoi network adds edge vertices to all seed points on the edge of the system.
#These voronoi cells usually have a very large surface area because the added edge vertices are far
#from the system edge (look at voronoi_network.png in the cell_data/capsule folder for clarification)
#
#The final cell network has parameters that are ploted in the
#cell_area_volume_intercelldist_degree_distribution.png image.
#A general rule is that cells have a typical area of around 80 um2
#and a volume of around 520 um3. This corresponds to a
#cell radius of 5 um (diameter = 10 um)
#

CAPSULE: int = 0 #capsule number
THRESHOLD_DIST: float = 6.0 #minimal distance between two points in the voronoi network
REAL_INTERCELL_DIST: float = 8.0 #fine-tuning of cell valumes
CROP: float = 40.0 #no default - number of final cells depends on this
NUMBER_OF_CELLS: int = 700 #default 700
MIN_XY: float = 0.0 #default 0.0
MAX_XY: float = 200 #default 200.0

model = create_cell_network(CAPSULE,
                            THRESHOLD_DIST,
                            REAL_INTERCELL_DIST,
                            CROP,
                            NUMBER_OF_CELLS,
                            MIN_XY,
                            MAX_XY)

#model data structure
#{
#   "cells": [
#               {
#                   "index": cell_count,
#                   "vertices_coord": np.array([(x1, y1), (x2, y2), (x3, y3)...]),
#                   "cm_xy": (cm_x, cm_y),
#                   "surface_area": float,
#                   "volume": float,
#                   "neighbours": {
#                                   1: float (length of connected polygon side)
#                                   2: ------------
#                                   n: ------------
#                               }
#               }
#           ],
#   "bin_conn_mat": np.ndarray (2D, int)
#   "weights": np.ndarray (2D, float)
#   "cell_height": float
#   "crop": float
#   "threshold_dist": float
#}
#
