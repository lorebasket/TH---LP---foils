import numpy as np
import beam_properties as bp

def calculate_element_stiffness_matrix(sectional_stiffness, L):
    """
    Compute the 12x12 element stiffness matrix using Bauchau's linear formulation:
    Ke = ∫ B(x)^T * C * B(x) dx

    Parameters:
    - sectional_stiffness: 6x6 matrix (from SONATA)
    - L: length of the beam element

    Returns:
    - 12x12 element stiffness matrix
    """
    Ke = np.zeros((12, 12))

    # Two-point Gauss-Legendre quadrature
    gauss_pts = np.array([-np.sqrt(1/3), np.sqrt(1/3)])
    weights = np.array([1.0, 1.0])

    for xi_hat, w in zip(gauss_pts, weights):
        # Map xi_hat from [-1, 1] to physical x in [0, L]
        x = (L / 2) * (xi_hat + 1)
        B = calculate_B_matrix(x, L, sectional_stiffness)
        Ke += B.T @ sectional_stiffness @ B * w * (L / 2)

    return Ke

def calculate_B_matrix(x, L, sectional_stiffness):
    """
    Construct the 6x12 strain-displacement matrix B(x)
    consistent with Bauchau's linear theory (warping via beta)
    """
    B = np.zeros((6, 12))
    xi = x / L

    # Linear shape functions
    N1 = 1 - xi
    N2 = xi
    dN1 = -1 / L
    dN2 = 1 / L

    # Cubic Hermite shape functions (bending)
    H1 = 1 - 3*xi**2 + 2*xi**3
    H2 = L * (xi - 2*xi**2 + xi**3)
    H3 = 3*xi**2 - 2*xi**3
    H4 = L * (-xi**2 + xi**3)

    dH1 = (-6*xi + 6*xi**2) / L
    dH2 = (1 - 4*xi + 3*xi**2)
    dH3 = (6*xi - 6*xi**2) / L
    dH4 = (-2*xi + 3*xi**2)

    beta = calculate_warping_parameter(sectional_stiffness)

    # εx = du/dx - beta*θx
    B[0, 0] = dN1
    B[0, 6] = dN2
    B[0, 3] = -beta * N1
    B[0, 9] = -beta * N2

    # γxy = dv/dx - θz
    B[1, 1] = dH1
    B[1, 7] = dH3
    B[1, 5] = -N1
    B[1, 11] = -N2

    # γxz = dw/dx - θy
    B[2, 2] = dH1
    B[2, 8] = dH3
    B[2, 4] = -N1
    B[2, 10] = -N2

    # κx = dθx/dx - beta*du/dx
    B[3, 3] = dN1
    B[3, 9] = dN2
    B[3, 0] = -beta * N1
    B[3, 6] = -beta * N2

    # κy = dθy/dx
    B[4, 4] = dN1
    B[4, 10] = dN2

    # κz = dθz/dx
    B[5, 5] = dN1
    B[5, 11] = dN2

    return B

def calculate_warping_parameter(C):
    """
    Compute the warping parameter beta using stiffness matrix C:
    beta = (1 - alpha) / (1 + alpha)
    where alpha = (c/d) * (t_v/t_h) * (Gv/Gh)
    and Gv and Gh are estimated from the shear terms in the stiffness matrix.
    """
    c = bp.b
    d = bp.h
    t_v = bp.s
    t_h = bp.s

    # Estimate shear modulus in vertical and horizontal walls
    Gv = C[4, 4] / bp.A
    Gh = C[5, 5] / bp.A

    alpha = (c/d) * (t_v/t_h) * (Gv/Gh)
    beta = (1 - alpha) / (1 + alpha)
    return beta

def assemble_global_stiffness_matrix():
    """
    Assemble the global stiffness matrix from all elements using Bauchau's procedure.
    """
    beam_model = bp.create_beam_model()
    K_global = np.zeros((bp.total_dof, bp.total_dof))

    for element in beam_model["elements"]:
        n1, n2 = element["nodes"]
        C = element["stiffness"]
        Ke = calculate_element_stiffness_matrix(C, bp.element_length)

        for i in range(12):
            gi = n1 * 6 + i if i < 6 else n2 * 6 + (i - 6)
            for j in range(12):
                gj = n1 * 6 + j if j < 6 else n2 * 6 + (j - 6)
                K_global[gi, gj] += Ke[i, j]

    return K_global

if __name__ == "__main__":
    K_global = assemble_global_stiffness_matrix()
    print("Global stiffness matrix shape:", K_global.shape)
    np.save("K_global.npy", K_global)
    
    # Save matrix to file for further analysis
    np.save("K_global.npy", K_global)
    
    print("\nStiffness matrix saved to K_global.npy")
