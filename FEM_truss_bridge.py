"""
    This script is an example for the utilization of the FEM_scratch program, and the subsequent post_processing using FEM_post_processing.
    The case analyzed in this example case is that of a simple truss bridge.
"""

import numpy as np
from FEM_scratch import FEM
from FEM_post_processing import FEM_visualization


# Parameters
E = 250e9
A = 0.00589
I = 0.00004545

# Node coordinates
node_positions = np.array([
    [0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0], [6, 0], [7, 0], [8, 0],
    [0.5, 0.5], [3.5, 0.5], [4, 0.5], [4.5, 0.5], [7.5, 0.5], [1, 1], [3, 1],
    [4, 1], [5, 1], [7, 1], [1.5, 1.5], [2.5, 1.5], [4, 1.5], [5.5, 1.5],
    [6.5, 1.5], [2, 2], [3, 2], [4, 2], [5, 2], [6, 2]
])

# Connectivity matrix
connectivity = np.array([
    [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9],
    [1, 10], [11, 5], [12, 5], [5, 13], [14, 9], [10, 15], [16, 11],
    [17, 12], [13, 18], [19, 14], [15, 20], [21, 16], [22, 17], [18, 23],
    [24, 19], [20, 25], [25, 21], [27, 22], [23, 29], [29, 24], [25, 26],
    [26, 27], [27, 28], [28, 29]
])

# Fixed nodes (1 means that that node is fixed)
fixed_nodes = np.array([1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# Create element positions
elements = []
for i in range(len(connectivity)):
    a, b = connectivity[i] - 1
    elements.append([*node_positions[a], *node_positions[b]])
elements = np.array(elements)

# Global matrix parameters
node_dof = 3
total_dof = node_dof * len(node_positions)

# Force vector
global_forces = np.zeros(total_dof)
loaded_node_1 = 3
loaded_node_2 = 7
global_forces[(loaded_node_1 - 1) * node_dof + 1] = -500000
global_forces[(loaded_node_2 - 1) * node_dof + 1] = -500000

# Execution of the FEM program
node_displacements, node_forces = FEM(elements, connectivity, fixed_nodes, node_dof, total_dof, global_forces,  E, A, I)

# Post processing to visualize FEM analysis
FEM_visualization(elements, node_displacements, node_dof, connectivity)

