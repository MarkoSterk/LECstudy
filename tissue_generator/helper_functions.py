"""
Helper functions for tissue constructor
"""
import numpy as np
from scipy.spatial.distance import cdist

class TriesOverflowException(Exception):
    """
    Tries overflow exception.
    Probable cause: to high threshold distance
    """
    def __init__(self, append_msg: str = ""):
        super().__init__(f"""Could not construct positions array.
                    Try decreasing the threshold distance or lower cell number or try again.
                    Because of the randomness of the process it sometimes can't complete.
                    {append_msg}""")

def random_xy_coord(min_coord: float, max_coord: float):
    """
    Generates a random pair of x,y coordinates
    """
    x = np.random.uniform(min_coord, max_coord)
    y = np.random.uniform(min_coord, max_coord)
    return x, y

def generate_random_coordinates(n: int, r_th: float,
                                max_coord: float,
                                min_coord: float = 200.0) -> np.ndarray:
    """
    Generates random coordinates for points in voronoi network
    with a set minimal threshold distance between points.
    """
    TRIES_LIMIT: int = 3000
    pos: list[tuple[float]] = []
    x0, y0 = random_xy_coord(min_coord, max_coord)
    pos.append((x0, y0))
    counter = 1
    exit_builder: bool = False
    for _ in range(1, n):
        checked_num = 0
        tries: int = 0
        while checked_num < counter:
            x, y = random_xy_coord(min_coord, max_coord)
            for j in range(counter):
                dist: float = np.sqrt((x-pos[j][0])**2+(y-pos[j][1])**2)
                if dist < r_th:
                    checked_num = 0
                    tries+=1
                    break
                checked_num+=1
            if tries > TRIES_LIMIT:
                exit_builder = True
        if exit_builder:
            print(f"""WARNING: exited random coordinates generator prematurely because of tries ({TRIES_LIMIT}) overflow.
                  Finished with {counter}/{n} cells.""")
            break
        pos.append((x, y))
        counter+=1
    return np.array(pos)

def distances_between_cells(pos: np.ndarray) -> np.ndarray:
    """
    Calculates distances between all cell pairs in pos array
    Returns matrix 2D with distances for all pairs
    """
    return cdist(pos, pos, 'euclidean')

def distances_between_neighbours(cell_data: dict[str, list[int]|float|np.ndarray]) -> np.ndarray:
    """
    Calculates distances between all connected cell pairs
    """
    intercell_distances = []
    for cell_i, data_i in cell_data.items():
        xi, yi = data_i["cm_xy"]
        for cell_j, data_j in cell_data.items():
            xj, yj = data_j["cm_xy"]
            if cell_i != cell_j and cell_j in data_i["neighbours"]:
                dist: float = np.sqrt((xi-xj)**2+(yi-yj)**2)
                intercell_distances.append(dist)
    return np.array(intercell_distances)

def create_bin_conn_matrix(cell_data: dict[str, list[int]|float|np.ndarray]) -> np.ndarray:
    """
    Creates 2D matrix of binarized intercellular connections
    """
    cell_num: int = len(cell_data)
    matrix: np.ndarray = np.zeros((cell_num, cell_num), int)
    for cell_i, data_i in cell_data.items():
        for cell_j in cell_data:
            if cell_j in data_i["neighbours"]:
                matrix[cell_i, cell_j] = 1
    return matrix

if __name__ == "__main__":
    NUMBER_OF_CELLS: int = 400
    MIN_XY: float = 0.0
    MAX_XY: float = 1.0
    THRESHOLD_DIST: float = 0.04
    positions: np.ndarray = generate_random_coordinates(NUMBER_OF_CELLS,
                                                  THRESHOLD_DIST,
                                                  MIN_XY,
                                                  MAX_XY)
    distances: np.ndarray = distances_between_cells(positions)
    check_dist: np.ndarray = distances[np.where((distances<THRESHOLD_DIST) & (distances>0.0))]
    if len(check_dist) > 0:
        print(f"WARNING: minimal distance between cells is smaller then the set {THRESHOLD_DIST=}")
    else:
        print("SUCCESS: cell coordinates generated successfully with given parameters")
