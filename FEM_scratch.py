import numpy as np
import sympy as sp


def FEM(elements, connectivity, fixed_nodes, node_dof, total_dof, global_forces,  E, A, I):

    # Symbolic variables
    global chi
    chi = sp.Symbol('chi')

    # Local matrix assembly
    local_matrices = [] # Stores all rotated local matrices
    for k in range(len(elements)):
        v1 = v1Func(k, elements)
        v2 = v2Func(k, elements)

        element_angle = np.arctan2((elements[k, 3] - elements[k, 1]), (elements[k, 2] - elements[k, 0]))
        he = np.sqrt((elements[k, 0] - elements[k, 2])**2 + (elements[k, 1] - elements[k, 3])**2) # Element length
        
        # Matrix K1
        K1 = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                d1 = sp.diff(v1[i], chi)
                d2 = sp.diff(v1[j], chi)
                K1[i, j] = E * A * sp.integrate(d1 * d2, (chi, 0, 1))
        
        # Matrix K2
        K2 = np.zeros((4, 4))
        for i in range(4):
            for j in range(4):
                dd1 = sp.diff(sp.diff(v2[i], chi), chi)
                dd2 = sp.diff(sp.diff(v2[j], chi), chi)
                K2[i, j] = E * I * sp.integrate(dd1 * dd2, (chi, 0, 1))

        # Combine matrices and rearrange matrix rows/columns 
        Kinter = np.block([
            [K1, np.zeros((2, 4))],
            [np.zeros((4, 2)), K2]
        ])
        Kinter = np.block([
            [Kinter[0]],
            [Kinter[4]],
            [Kinter[3]],
            [Kinter[1]],
            [Kinter[2]],
            [Kinter[5]]
        ])

        Kinter = np.block([
            [
                Kinter[:, 0:1],
                Kinter[:, 4:5],
                Kinter[:, 3:4],
                Kinter[:, 1:2],
                Kinter[:, 2:3],
                Kinter[:, 5:6]
            ]
        ])

        # Rotate the matrix
        Rot = np.array([
            [np.cos(element_angle), np.sin(element_angle), 0, 0, 0, 0],
            [-np.sin(element_angle), np.cos(element_angle), 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, np.cos(element_angle), np.sin(element_angle), 0],
            [0, 0, 0, -np.sin(element_angle), np.cos(element_angle), 0],
            [0, 0, 0, 0, 0, 1]
        ])
        KTotal = Rot.T @ Kinter @ Rot
        local_matrices.append(KTotal)

    # Construction of the global stiffness matrix
    global_matrix = np.zeros((total_dof, total_dof))
    for e in range(len(connectivity)):
        i = connectivity[e, 0]  # First node of the element
        j = connectivity[e, 1]  # Second node of the element

        # Positioning at the top-left corner of the submatrix
        iG = (i - 1) * node_dof  # Global index for node i
        jG = (j - 1) * node_dof  # Global index for node j

        # Insert values from submatrices into the global_matrix matrix
        global_matrix[iG:iG + 3, iG:iG + 3] += local_matrices[e][0:3, 0:3]
        global_matrix[iG:iG + 3, jG:jG + 3] += local_matrices[e][0:3, 3:6]
        global_matrix[jG:jG + 3, iG:iG + 3] += local_matrices[e][3:6, 0:3]
        global_matrix[jG:jG + 3, jG:jG + 3] += local_matrices[e][3:6, 3:6]


    # Creating new matrix variable, to accomodate for the boundary conditions
    BC_global_matrix = global_matrix
    BC_global_forces = global_forces

    # Iterate over fixed_nodes in reverse order (Evaluating the boundary conditions)
    for i in range(len(fixed_nodes)):
        # Start from the end of fixed_nodes
        k = len(fixed_nodes) - i - 1

        if fixed_nodes[k] == 1:
            M = k * node_dof  # Starting index for the fixed node in Python
            BC_global_matrix = np.delete(BC_global_matrix, slice(M, M + 3), axis=0)  # Remove rows
            BC_global_matrix = np.delete(BC_global_matrix, slice(M, M + 3), axis=1)  # Remove columns
            BC_global_forces = np.delete(BC_global_forces, slice(M, M + 3))          # Remove entries

    # Calculate the free degrees of freedom
    free_node_displacements = np.linalg.solve(BC_global_matrix, BC_global_forces)

    # Rebuild the vector of degrees of freedom, including fixed nodes
    node_displacements = np.zeros(total_dof)
    j = 0  # Index for `free_node_displacements`

    for i in range(len(fixed_nodes)):
        k = i * node_dof  # Compute the global degree of freedom index
        if fixed_nodes[i] == 1:
            node_displacements[k] = 0
            node_displacements[k + 1] = 0
            node_displacements[k + 2] = 0
        else:
            node_displacements[k] = free_node_displacements[j]
            node_displacements[k + 1] = free_node_displacements[j + 1]
            node_displacements[k + 2] = free_node_displacements[j + 2]
            j += 3

    # Calculating the force vector, including the fixed nodes 
    node_forces = global_matrix @ node_displacements

    return node_displacements, node_forces


# Define v1 function (Lagrangian interpolation functions)
def v1Func(num, elem):
    he = np.sqrt((elem[num, 0] - elem[num, 2])**2 + (elem[num, 1] - elem[num, 3])**2)
    v1 = [1 - chi, chi]
    return np.array(v1)

# Define v2 function (Hermitic functions)
def v2Func(num, elem):
    he = np.sqrt((elem[num, 0] - elem[num, 2])**2 + (elem[num, 1] - elem[num, 3])**2)
    v2 = [1 - 3 * chi**2 + 2 * chi**3,
          -he * chi * (1 - chi)**2,
          3 * chi**2 - 2 * chi**3,
          -he * chi * (chi**2 - chi)]
    return np.array(v2)

