
import matplotlib.pyplot as plt


def FEM_visualization(elements, node_displacements, node_dof, connectivity):
    """
    Plot the initial and final positions of elements in a structure.

    Parameters:
        elements (ndarray): Array of element coordinates.
        node_displacements (ndarray): Displacement vector for all nodes.
        node_dof (int): Degrees of freedom per node.
        connectivity (ndarray): Element connectivity array.
    """
    plt.figure()
    plt.title("Initial and Deformed Structure")
    size = 30
    
    for i in range(len(elements)):
        # Initial positions of nodes
        plt.plot([elements[i, 0], elements[i, 2]], [elements[i, 1], elements[i, 3]], 'g-', linewidth=1)
        plt.scatter([elements[i, 0], elements[i, 2]], [elements[i, 1], elements[i, 3]], s=size, c='blue')
        
        # Final positions of nodes
        c = connectivity[i]
        ka = (c[0] - 1) * node_dof
        kb = (c[1] - 1) * node_dof
        PosNueva = [
            node_displacements[ka] + elements[i, 0], node_displacements[ka + 1] + elements[i, 1],
            node_displacements[kb] + elements[i, 2], node_displacements[kb + 1] + elements[i, 3]
        ]
        plt.plot([PosNueva[0], PosNueva[2]], [PosNueva[1], PosNueva[3]], 'r-', linewidth=1)
        plt.scatter([PosNueva[0], PosNueva[2]], [PosNueva[1], PosNueva[3]], s=size, c='red')
    
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid(True)
    plt.show()

def PostAnalisis(elements, node_displacements, node_dof, E, connectivity, he):
    """
    Perform post-analysis to calculate maximum deformation and stress.

    Parameters:
        elements (ndarray): Array of element coordinates.
        node_displacements (ndarray): Displacement vector for all nodes.
        node_dof (int): Degrees of freedom per node.
        E (float): Young's modulus.
        connectivity (ndarray): Element connectivity array.
        he (float): Element length.

    Returns:
        MaxDeformacion (ndarray): Maximum deformation for each element.
        MaxEsfuerzo (ndarray): Maximum stress for each element.
    """
    chi = sp.symbols('chi')
    MaxDeformacion = []
    MaxEsfuerzo = []
    
    for k in range(len(elements)):
        # Define v1 and v2 functions
        v1 = v1Func(k, elements) / he**3
        v2 = v2Func(k, elements) / he**3
        
        a = connectivity[k, 0]
        b = connectivity[k, 1]
        
        # Calculate deformation
        eps0 = node_displacements[(a - 1) * node_dof] * sp.diff(v1[0], chi) + \
               node_displacements[(b - 1) * node_dof] * sp.diff(v1[1], chi)
        eps1 = (node_displacements[(a - 1) * node_dof + 1] * sp.diff(sp.diff(v2[0], chi), chi) + 
                node_displacements[(a - 1) * node_dof + 2] * sp.diff(sp.diff(v2[1], chi), chi) + 
                node_displacements[(b - 1) * node_dof + 1] * sp.diff(sp.diff(v2[2], chi), chi) + 
                node_displacements[(b - 1) * node_dof + 2] * sp.diff(sp.diff(v2[3], chi), chi))
        deformacion = eps0 + eps1
        
        # Calculate stress
        esfuerzo = E * eps0 - E * (0.203 / 2) * eps1
        
        # Evaluate deformation and stress at points
        gradDeformacion = [deformacion.subs(chi, i / 10 - 1 / 10) for i in range(1, 12)]
        gradEsfuerzo = [esfuerzo.subs(chi, i / 10) for i in range(1, 12)]
        
        # Store maximum deformation and stress
        MaxDeformacion.append(float(max(gradDeformacion)))
        MaxEsfuerzo.append(float(max(gradEsfuerzo)))
    
    return np.array(MaxDeformacion), np.array(MaxEsfuerzo)