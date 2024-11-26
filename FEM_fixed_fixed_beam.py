"""
    This script is an example for the utilization of the FEM_scratch program, and the subsequent post_processing using FEM_post_processing.
    The case analyzed in this example case is that of a fixed-fixed beam
"""

import numpy as np
from FEM_scratch import FEM
from FEM_post_processing import FEM_visualization


# Parameters
E = 200e9       # Young's Modulus
A = 0.00589     # Cross-sectional area
I = 0.00001     # Second area moment of the cross-section

# Node coordinates
node_positions = np.array([
    [0, 0], [1, 0], [2, 0], [3, 0], [4, 0],
    [5, 0], [6, 0], [7, 0], [8, 0], [9, 0],
    [10, 0], [11, 0], [12, 0], [13, 0], [14, 0],
    [15, 0], [16, 0], [17, 0], [18, 0], [19, 0],
    [20, 0]
])

# Connectivity matrix
connectivity = np.array([
    [1, 2], [2, 3], [3, 4], [4, 5], [5, 6],
    [6, 7], [7, 8], [8, 9], [9, 10], [10, 11],
    [11, 12], [12, 13], [13, 14], [14, 15], [15, 16],
    [16, 17], [17, 18], [18, 19], [19, 20], [20, 21]
])

# Fixed nodes (1 means that that node is fixed)
fixed_nodes = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])

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
loaded_node = 11
global_forces[(loaded_node - 1) * node_dof + 1] = -10000

# Execution of the FEM program
node_displacements, node_forces = FEM(elements, connectivity, fixed_nodes, node_dof, total_dof, global_forces,  E, A, I)

# Post processing to visualize FEM analysis
FEM_visualization(elements, node_displacements, node_dof, connectivity)

