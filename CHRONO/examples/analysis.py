import numpy as np
from scipy.linalg import solve, eigh
from beam_properties import dof_per_node, total_dof, n_nodes


def solve_static_analysis(K_global, force_vector):
    
    # Apply cantilever boundary conditions (fix first node)
    # We'll remove the first 6 DOFs from the system
    fixed_dofs = list(range(dof_per_node))  # First node's DOFs
    free_dofs = list(range(dof_per_node, total_dof))  # All other DOFs
    
    # Partition the stiffness matrix and force vector
    K_ff = K_global[np.ix_(free_dofs, free_dofs)]
    F_f = force_vector[free_dofs]
    
    # Solve for displacements
    u_f = solve(K_ff, F_f)
    
    # Create full displacement vector (with zeros at constrained DOFs)
    u_full = np.zeros(total_dof)
    u_full[free_dofs] = u_f
    
    # Calculate reaction forces
    reaction_forces = K_global[np.ix_(fixed_dofs, range(total_dof))] @ u_full
    
    return u_full, reaction_forces


def modal_analysis(K_global, M_global, num_modes=10):

    # Apply cantilever boundary conditions
    fixed_dofs = list(range(dof_per_node))  # First node's DOFs
    free_dofs = list(range(dof_per_node, total_dof))  # All other DOFs
    
    # Partition the stiffness and mass matrices
    K_ff = K_global[np.ix_(free_dofs, free_dofs)]
    M_ff = M_global[np.ix_(free_dofs, free_dofs)]
    
    # Solve the generalized eigenvalue problem
    # The eigenvalues are ω², where ω is the natural frequency in rad/s
    eigvals, eigvecs = eigh(K_ff, M_ff, eigvals=(0, min(num_modes-1, len(free_dofs)-1)))
    
    # Calculate natural frequencies in Hz
    frequencies = np.sqrt(eigvals) / (2 * np.pi)
    
    # Create full mode shapes (with zeros at constrained DOFs)
    mode_shapes = np.zeros((total_dof, len(eigvals)))
    mode_shapes[free_dofs, :] = eigvecs
    
    return frequencies, mode_shapes


def extract_displacements(u_full):
    
    # Extract nodal displacements and rotations from the solution vector
    displacements = []
    for i in range(n_nodes):
        start_dof = i * dof_per_node
        node_disp = {
            'node': i,
            'x': u_full[start_dof],     # x displacement
            'y': u_full[start_dof + 1], # y displacement
            'z': u_full[start_dof + 2], # z displacement
            'rx': u_full[start_dof + 3], # rotation around x
            'ry': u_full[start_dof + 4], # rotation around y
            'rz': u_full[start_dof + 5]  # rotation around z
        }
        displacements.append(node_disp)
    return displacements


def create_force_vector(load_magnitude=-100.0):

    # Create a force vector with a vertical load at the last node
    force_vector = np.zeros(total_dof)
    force_vector[n_nodes * dof_per_node - 5] = load_magnitude  # Y-direction force at last node
    return force_vector