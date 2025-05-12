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

def extract_coupling_factors(sonata_matrix):
    """Extract coupling factors from Sonata matrix"""
    coupling_factors = {}
    
    # Axial-shear coupling (terms [0,4] and [4,0])
    coupling_factors['axial_shear'] = sonata_matrix[0,4] / sonata_matrix[0,0]
    
    # Bending mode coupling (terms [1,2] and [2,1])
    coupling_factors['bending_coupling'] = sonata_matrix[1,2] / sonata_matrix[1,1]
    
    # Bending-shear coupling (terms [1,3] and [3,1])
    coupling_factors['bending_shear'] = sonata_matrix[1,3] / sonata_matrix[1,1]
    
    # Shear-torsion coupling (terms [4,5] and [5,4])
    coupling_factors['shear_torsion'] = sonata_matrix[4,5] / sonata_matrix[4,4]
    
    return coupling_factors

def calculate_B_matrix(x, L):
    """Calculates the strain-displacement matrix B for a beam element
    using the strain ordering expected by Sonata
    
    Args:
        x: position along the element (local coordinate, 0 <= x <= L)
        L: length of the element
        
    Returns:
        B: 6x12 strain-displacement matrix with Sonata's strain ordering
    """
    B = np.zeros((6, 12))
    # Normalized coordinate
    xi = x / L
    
    # Basic shape functions
    N1_val = 1 - xi
    N2_val = xi
    dN1dx = -1/L
    dN2dx = 1/L
    
    # Cubic Hermite shape functions
    H1 = 1 - 3 * xi**2 + 2 * xi**3
    H2 = L * (xi - 2 * xi**2 + xi**3)
    H3 = 3 * xi**2 - 2 * xi**3
    H4 = L * (-xi**2 + xi**3)
    
    # Derivatives
    dH1dx = (-6 * xi + 6 * xi**2) / L
    dH2dx = (1 - 4 * xi + 3 * xi**2)
    dH3dx = (6 * xi - 6 * xi**2) / L
    dH4dx = (-2 * xi + 3 * xi**2)
    
    # Based on the reordering observed in Sonata matrix, I'll reorder the strains
    # Assume Sonata strain order: [εx, γxy, γxz, κx, κy, κz]
    
    # 0. Axial strain: εx = du/dx (same as original)
    B[0, 0] = dN1dx  # u1
    B[0, 6] = dN2dx  # u2
    
    # 1. Shear strain γxy = dv/dx - θz (was index 3 in original)
    B[1, 1] = dH1dx   # v1
    B[1, 5] = dH2dx - N1_val  # θz1
    B[1, 7] = dH3dx   # v2
    B[1, 11] = dH4dx - N2_val  # θz2
    
    # 2. Shear strain γxz = dw/dx - θy (was index 4 in original)
    B[2, 2] = dH1dx   # w1
    B[2, 4] = dH2dx - N1_val  # θy1
    B[2, 8] = dH3dx   # w2
    B[2, 10] = dH4dx - N2_val  # θy2
    
    # 3. Torsion κx = dθx/dx (was index 5 in original)
    B[3, 3] = dN1dx  # θx1
    B[3, 9] = dN2dx  # θx2
    
    # 4. Bending curvature about y: κy = dθy/dx (was index 2 in original)
    B[4, 4] = dN1dx  # θy1
    B[4, 10] = dN2dx # θy2
    
    # 5. Bending curvature about z: κz = dθz/dx (was index 1 in original)
    B[5, 5] = dN1dx  # θz1
    B[5, 11] = dN2dx # θz2
    
    # Add coupling terms based on Sonata matrix structure
    
    # 1. Axial-shear coupling (terms [0,4] and [4,0])
    # This represents warping effects
    coupling_factors = extract_coupling_factors(bp.K1)
    cf_axial_shear = coupling_factors.get('axial_shear', 0)
    if cf_axial_shear != 0:
        # Axial strain affected by θy rotation
        B[0, 4] += cf_axial_shear * N1_val  # θy1
        B[0, 10] += cf_axial_shear * N2_val  # θy2
        
        # Shear strain γxz affected by axial displacement
        B[4, 0] += cf_axial_shear * N1_val  # u1
        B[4, 6] += cf_axial_shear * N2_val  # u2
    
    # 2. Bending mode coupling (terms [1,2] and [2,1])
    # This represents asymmetry in the cross-section
    cf_bending = coupling_factors.get('bending_coupling', 0)
    if cf_bending != 0:
        # Bending about z affected by θy rotation
        B[1, 4] += cf_bending * N1_val  # θy1
        B[1, 10] += cf_bending * N2_val  # θy2
        
        # Bending about y affected by θz rotation
        B[2, 5] += cf_bending * N1_val  # θz1
        B[2, 11] += cf_bending * N2_val  # θz2
    
    # 3. Bending-shear coupling (terms [1,3], [3,1], etc.)
    cf_bend_shear = coupling_factors.get('bending_shear', 0)
    if cf_bend_shear != 0:
        # Bending curvature affected by w displacement
        B[1, 2] += cf_bend_shear * dH1dx  # w1
        B[1, 8] += cf_bend_shear * dH3dx  # w2
        
        # Shear strain affected by bending rotation
        B[3, 4] += cf_bend_shear * dN1dx  # θy1
        B[3, 10] += cf_bend_shear * dN2dx  # θy2
    
    # 4. Shear-torsion coupling (terms [4,5] and [5,4])
    # This represents non-coincident shear center
    cf_shear_torsion = coupling_factors.get('shear_torsion', 0)
    if cf_shear_torsion != 0:
        # Shear γxz affected by θx rotation
        B[4, 3] += cf_shear_torsion * N1_val  # θx1
        B[4, 9] += cf_shear_torsion * N2_val  # θx2
        
        # Torsion κx affected by θy rotation
        B[5, 4] += cf_shear_torsion * N1_val  # θy1
        B[5, 10] += cf_shear_torsion * N2_val  # θy2
    
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
