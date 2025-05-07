import numpy as np

import beam_properties as bp

def calculate_element_stiffness_matrix(nodes_stiffness, L):
    """
    Calculate the stiffness matrix for a beam element using the complete sectional stiffness matrix.
    
    Parameters:
    -----------
    L : float
        Length of the beam element
    sectional_stiffness_matrix : numpy.ndarray
        6x6 sectional stiffness matrix containing all stiffness coefficients

    [K] = ∫[B(x)]ᵀ[D][B(x)]dx

    Returns:
    --------
    K : numpy.ndarray
        12x12 element stiffness matrix (including all DOFs)
    """
    
    # Initialize stiffness matrix
    K_element = np.zeros((12, 12))
    
    # Define matrices of shape functions derivatives
    # We need to integrate B^T * C * B over the element length
    
    # Use Gauss-Legendre quadrature for integration
    gauss_points = np.array([-np.sqrt(1/3), np.sqrt(1/3)])
    gauss_weights = np.array([1.0, 1.0])
    
    # Map Gauss points from [-1, 1] to [0, L]
    x_points = L/2 * (gauss_points + 1)
    
    # Perform numerical integration
    for i, (x, w) in enumerate(zip(x_points, gauss_weights)):
        # Calculate the B matrix at this Gauss point
        B = calculate_B_matrix(x, L)
        
        # Compute contribution to stiffness matrix: B^T * C * B * w * det(J)
        # The determinant of the Jacobian for the mapping from [0,L] to [-1,1] is L/2
        K_element += np.dot(np.dot(B.T, nodes_stiffness), B) * w * (L/2)
    
    return K_element

def calculate_B_matrix(x, L):
    """Calculates the strain-displacement matrix B for a Timoshenko beam element
    using cubic Hermite functions for transverse displacements and linear functions for rotations.
    Args:
        x: position along the element (local coordinate, 0 <= x <= L)
        L: length of the element
    Returns:
        B: 6x12 strain-displacement matrix
    """
    B = np.zeros((6, 12))

    # Normalized coordinate
    xi = x / L

    # --- Linear Shape Functions (for u, θx, θy, θz) ---
    # N1(xi) = 1 - xi, N2(xi) = xi
    # dN1/dx = -1/L, dN2/dx = 1/L
    N1_val = 1 - xi
    N2_val = xi
    dN1dx = -1/L
    dN2dx =  1/L

    # --- Cubic Hermite Shape Functions (for v, w) ---
    H1 = 1 - 3 * xi**2 + 2 * xi**3
    H2 = L * (xi - 2 * xi**2 + xi**3)
    H3 = 3 * xi**2 - 2 * xi**3
    H4 = L * (-xi**2 + xi**3)

    # First derivatives w.r.t. x
    dH1dx = (-6 * xi + 6 * xi**2) / L
    dH2dx = (1 - 4 * xi + 3 * xi**2)
    dH3dx = (6 * xi - 6 * xi**2) / L
    dH4dx = (-2 * xi + 3 * xi**2)

    # Second derivatives w.r.t. x (Not used for Timoshenko B matrix directly)
    # d2H1dx2 = (-6 + 12 * xi) / L**2
    # d2H2dx2 = (-4 + 6 * xi) / L
    # d2H3dx2 = (6 - 12 * xi) / L**2
    # d2H4dx2 = (-2 + 6 * xi) / L

    # DOF ordering:
    # Node 1: [u1, v1, w1, θx1, θy1, θz1]
    # Node 2: [u2, v2, w2, θx2, θy2, θz2]

    # Strain/curvature vector: [εx, κz, κy, γxy, γxz, κx]^T

    # 0. Axial strain: εx = du/dx (Linear N)
    B[0, 0] = dN1dx  # u1
    B[0, 6] = dN2dx  # u2

    # 1. Bending curvature about z: κz = dθz/dx (Linear N for θz)
    B[1, 5] = dN1dx  # θz1
    B[1, 11] = dN2dx # θz2

    # 2. Bending curvature about y: κy = dθy/dx (Linear N for θy)
    B[2, 4] = dN1dx  # θy1
    B[2, 10] = dN2dx # θy2

    # 3. Shear strain γxy = dv/dx - θz (Cubic H for dv/dx, Linear N for θz)
    B[3, 1] = dH1dx   # v1 contribution to dv/dx
    B[3, 5] = dH2dx   # θz1 contribution to dv/dx (This is the H term for v related to θz1)
    B[3, 7] = dH3dx   # v2 contribution to dv/dx
    B[3, 11] = dH4dx  # θz2 contribution to dv/dx (This is the H term for v related to θz2)
    # Subtract θz interpolated linearly
    B[3, 5] -= N1_val # θz1 (contribution from -θz term)
    B[3, 11] -= N2_val # θz2 (contribution from -θz term)

    # 4. Shear strain γxz = dw/dx - θy (Cubic H for dw/dx, Linear N for θy)
    B[4, 2] = dH1dx   # w1 contribution to dw/dx
    B[4, 4] = dH2dx   # θy1 contribution to dw/dx (This is the H term for w related to θy1)
    B[4, 8] = dH3dx   # w2 contribution to dw/dx
    B[4, 10] = dH4dx  # θy2 contribution to dw/dx (This is the H term for w related to θy2)
    # Subtract θy interpolated linearly
    B[4, 4] -= N1_val # θy1 (contribution from -θy term)
    B[4, 10] -= N2_val # θy2 (contribution from -θy term)

    # 5. Torsion κx = dθx/dx (Linear N for θx)
    B[5, 3] = dN1dx  # θx1
    B[5, 9] = dN2dx  # θx2

    return B


def assemble_global_stiffness_matrix():
    """
    Calculate the global stiffness matrix for the entire beam model
    """
    # Get beam model from properties
    beam_model = bp.create_beam_model()
    
    # Initialize global stiffness matrix
    K_global = np.zeros((bp.total_dof, bp.total_dof))
    
    # Loop through all elements
    for i, element in enumerate(beam_model["elements"]):
        # Get nodes of this element
        start_node = element["nodes"][0]
        end_node = element["nodes"][1]
        
        # Get sectional stiffness matrix
        sectional_stiffness = element["stiffness"]
        
        # Calculate element stiffness matrix
        K_element = calculate_element_stiffness_matrix(
            sectional_stiffness, 
            bp.element_length
        )
        
        # Assemble into global stiffness matrix
        for local_i in range(bp.dof_per_element):
            # Map local DOF to global DOF
            if local_i < bp.dof_per_node:
                global_i = start_node * bp.dof_per_node + local_i
            else:
                global_i = end_node * bp.dof_per_node + (local_i - bp.dof_per_node)
            
            for local_j in range(bp.dof_per_element):
                # Map local DOF to global DOF
                if local_j < bp.dof_per_node:
                    global_j = start_node * bp.dof_per_node + local_j
                else:
                    global_j = end_node * bp.dof_per_node + (local_j - bp.dof_per_node)
                
                # Add contribution to global stiffness matrix
                K_global[global_i, global_j] += K_element[local_i, local_j]
    
    return K_global

if __name__ == "__main__":
    # Print beam parameters
    print(f"Beam length: {bp.beam_length} m")
    print(f"Number of elements: {bp.n_el}")
    print(f"Element length: {bp.element_length} m")
    print(f"Total DOFs: {bp.total_dof}")
    
    # Calculate global stiffness matrix
    K_global = assemble_global_stiffness_matrix()
    
    print(f"\nGlobal stiffness matrix size: {K_global.shape}")
    
    # Check conditioning of the stiffness matrix
    try:
        cond_num = np.linalg.cond(K_global)
        print(f"Condition number of global stiffness matrix: {cond_num:.2e}")
    except:
        print("Could not calculate condition number - matrix may be singular")
    
    # Show a small portion of the stiffness matrix for verification
    print("\nSample of global stiffness matrix (10x10 corner):")
    np.set_printoptions(precision=2, suppress=True)
    print(K_global[:10, :10])
    
    # Save matrix to file for further analysis
    np.save("K_global.npy", K_global)
    
    print("\nStiffness matrix saved to K_global.npy")