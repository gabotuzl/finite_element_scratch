# Finite Element Scratch

**Finite Element Scratch** is a functioning finite element method (FEM) analysis tool built entirely from scratch, using only `numpy` and `sympy` for computation and `matplotlib` for visualization. This project demonstrates the fundamental steps in FEM, including the creation of stiffness matrices, solving global equations, and visualizing results for different cases.
**Note: All units used in this project are SI units**

## Features

- **Stiffness Matrix Construction**:
  - Local matrices are derived and integrated using interpolation functions.
  - Lagrange functions are used for axial displacement.
  - Hermitic functions are used for bending displacement.
- **Global Matrix Assembly**:
  - Connectivity of elements is taken into account to construct the global matrix.
- **Solution for Displacements and Forces**:
  - Compute the displacements of all nodes.
  - Calculate the forces exerted on each node.
- **Visualization**:
  - Use `matplotlib` to visualize the results.

## Example Use Cases

This repository includes two example cases:

1. **Truss Bridge Analysis**: `FEM_truss_bridge.py`
2. **Fixed-Fixed Beam Analysis**: `FEM_fixed_fixed_beam.py`

## Getting Started

### Prerequisites

Ensure you have Python installed along with the following libraries:

- `numpy`
- `sympy`
- `matplotlib`

You can install the required libraries using pip:

```bash
pip install numpy sympy matplotlib
```

### Usage

1. Clone the repository:

```bash
git clone https://github.com/gabotuzl/finite_element_scratch.git
cd finite_element_scratch
```

2. Run one of the example cases:

```bash
python FEM_truss_bridge.py
```

or

```bash
python FEM_fixed_fixed_beam.py
```

3. Customize your analysis by modifying specific fields as described below.

### Customizing Cases

To define your custom FEM analysis, update the following fields in your script:

- **`node_positions` (np.array)**: An array containing the global coordinates of each node. The position in the array corresponds to the node number.
  
  Example:
  ```python
  node_positions = np.array([[0, 0], [1, 0], [2, 0]])
  ```

- **`connectivity` (np.array)**: An array defining the connections between nodes. For example, if node 1 connects to node 2, the element is `[1, 2]`. Nodes are numbered from 1 to N.

  Example:
  ```python
  connectivity = np.array([[1, 2], [2, 3]])
  ```

- **`fixed_nodes` (np.array)**: An array specifying boundary conditions. A `1` at a position indicates the corresponding node is fixed.

  Example:
  ```python
  fixed_nodes = np.array([1, 0, 1])  # Node 1 and Node 3 are fixed
  ```

- **`global_forces` (np.array)**: An array containing applied forces. Use the same format as shown in the examples to ensure correct application of point loads.

  Example:
  ```python
  global_forces = np.zeros(total_dof)
  loaded_node_1 = 3
  loaded_node_2 = 7
  global_forces[(loaded_node_1 - 1) * node_dof + 1] = -500000
  global_forces[(loaded_node_2 - 1) * node_dof + 1] = -500000
  ```

## Contributing

Contributions are welcome! If you have suggestions, improvements, or new features to add, feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.


