import numpy as np

# Beam parameters
n_el = 100
n_nodes = n_el + 1
beam_length = 2  # m
b = 100e-3            # m
h = 100e-3 
s = 40e-3            # m
b_int = b-2*s
h_int = h-2*s
A = b * h -((b-2*s)*(h-2*s)) # Area
Iyy = (b * h**3 - (b_int)*(h_int)**3) / 12 # Second moment of area 
Izz = (b**3 * h - (b_int)**3*(h_int)) / 12 # Second moment of area 

rho = 7850           # density [kg/m^3]
sectional_mass = rho * A
E = 210e9             # Young Modulus [Pa]
nu = 0.3            # Poisson's ratio
G = E /(2*(1+nu))    # Shear ratio
Ixx = Iyy + Izz
k_shear = 1      #motivation # Shear correction factor (APPROXIMATION for hollow square, adjust if needed)

element_length = beam_length / n_el
phi = (12 * E * Iyy) / (G * A * element_length**2) # Note: This phi might not be directly used in Timoshenko element formulation
dof_per_node = 6  # 3 displacement + 3 rotation
dof_per_element = 12  # 6 DOF per node * 2 nodes per element
total_dof = (n_nodes) * dof_per_node

# Input stiffness matrices (6x6 sectional stiffness for different composite sections)
# Order of terms corresponds to B-matrix strains: [εx, κz, κy, γxy, γxz, κx]^T
K1 = np.array([
    [E*A, 0, 0, 0, 0, 0],          # Axial stiffness (εx)
    [0, E*Izz, 0, 0, 0, 0],        # Bending stiffness Mz (κz) - related to v displacement
    [0, 0, E*Iyy, 0, 0, 0],        # Bending stiffness My (κy) - related to w displacement
    [0, 0, 0, G*A*k_shear, 0, 0],  # Shear stiffness (γxy)
    [0, 0, 0, 0, G*A*k_shear, 0],  # Shear stiffness (γxz)
    [0, 0, 0, 0, 0, G*Ixx]         # Torsional stiffness G*J (κx) - using Ixx as approx for J
])

K2 = np.array([
    [E*A, 0, 0, 0, 0, 0],
    [0, E*Izz, 0, 0, 0, 0],
    [0, 0, E*Iyy, 0, 0, 0],
    [0, 0, 0, G*A*k_shear, 0, 0],
    [0, 0, 0, 0, G*A*k_shear, 0],
    [0, 0, 0, 0, 0, G*Ixx]
])

K3 = np.array([
    [E*A, 0, 0, 0, 0, 0],
    [0, E*Izz, 0, 0, 0, 0],
    [0, 0, E*Iyy, 0, 0, 0],
    [0, 0, 0, G*A*k_shear, 0, 0],
    [0, 0, 0, 0, G*A*k_shear, 0],
    [0, 0, 0, 0, 0, G*Ixx]
])

## Input stiffness matrices (6x6 sectional stiffness for different composite sections)
#K1 = np.array([
#    [238781114.19621196,2.605195384951224e-15,-1.5895551477921026e-18,-2.886296777791002e-10,-7383.2802423778785,-10.541715031703083],
#    [2.605195384066693e-15,58044624.04659345,-382.75020421727334,19.81182234254027,1.5498648569585233e-10,-1.0309163156718082e-12],
#    [-1.5895548352787762e-18,-382.7502031456955,20686102.974347055,123796.030915383,5.000632457103843e-13,1.671669294580594e-10],
#    [-2.88629677779102e-10,19.81179151763727,123796.03091211985,91836348042.67496,-2.6806653469861883e-12,1.674238127515022e-13],
#    [-7383.28021992029,1.54986485695838e-10,5.00063245717894e-13,-2.6806653461113677e-12,99588187701.10786,-973381.7038727745],
#    [-10.541746261310754,-1.030916315648805e-12,1.6716692945804869e-10,1.6742381286401418e-13,-973381.70364537,302990024131.3176]
#])
#
#K2 = np.array([
#    [238781114.19621196,2.605195384951224e-15,-1.5895551477921026e-18,-2.886296777791002e-10,-7383.2802423778785,-10.541715031703083],
#    [2.605195384066693e-15,58044624.04659345,-382.75020421727334,19.81182234254027,1.5498648569585233e-10,-1.0309163156718082e-12],
#    [-1.5895548352787762e-18,-382.7502031456955,20686102.974347055,123796.030915383,5.000632457103843e-13,1.671669294580594e-10],
#    [-2.88629677779102e-10,19.81179151763727,123796.03091211985,91836348042.67496,-2.6806653469861883e-12,1.674238127515022e-13],
#    [-7383.28021992029,1.54986485695838e-10,5.00063245717894e-13,-2.6806653461113677e-12,99588187701.10786,-973381.7038727745],
#    [-10.541746261310754,-1.030916315648805e-12,1.6716692945804869e-10,1.6742381286401418e-13,-973381.70364537,302990024131.3176]         # Bending z (Ey)
#])
#
#K3 = np.array([
#    [238781114.19621196,2.605195384951224e-15,-1.5895551477921026e-18,-2.886296777791002e-10,-7383.2802423778785,-10.541715031703083],
#    [2.605195384066693e-15,58044624.04659345,-382.75020421727334,19.81182234254027,1.5498648569585233e-10,-1.0309163156718082e-12],
#    [-1.5895548352787762e-18,-382.7502031456955,20686102.974347055,123796.030915383,5.000632457103843e-13,1.671669294580594e-10],
#    [-2.88629677779102e-10,19.81179151763727,123796.03091211985,91836348042.67496,-2.6806653469861883e-12,1.674238127515022e-13],
#    [-7383.28021992029,1.54986485695838e-10,5.00063245717894e-13,-2.6806653461113677e-12,99588187701.10786,-973381.7038727745],
#    [-10.541746261310754,-1.030916315648805e-12,1.6716692945804869e-10,1.6742381286401418e-13,-973381.70364537,302990024131.3176]         # Bending z (Ey)
#])

# old and asymmetric
#K2 = np.array([
#    [181439704539787.12,0.0,0.0,0.0,-7382811147.34076,-44.64380383642636],
#    [0.0,43956959328385.055,-1438047210.7943156,22664527253.159916,0.0,0.0],
#    [0.0,-1438047210.7563467,15869126778442.68,39180528954.10046,0.0,0.0],
#    [0.0,22664527292.75935,39180528965.639404,7.240174421959046e+16,0.0,0.0],
#    [-7382811099.6284895,0.0,0.0,0.0,7.867994941373424e+16,-4090130555585.707],
#    [-34.56007546292466,0.0,0.0,0.0,-4090130556046.7227,2.3550723413355952e+17]         # Bending z (Ey)
#])



# Input stiffness matrices (6x6 sectional stiffness for different composite sections)
M1 = np.diag([
    rho*A,            # Axial mass
    rho*A,            # Transverse mass y
    rho*A,            # Transverse mass Z
    rho*A*Ixx,        # Rotational mass around x
    rho*A*Izz,        # Rotational mass around y
    rho*A*Iyy         # Rotational mass around z
])

M2 = np.diag([
    rho*A,            # Axial mass
    rho*A,            # Transverse mass y
    rho*A,            # Transverse mass Z
    rho*A*Ixx,        # Rotational mass around x
    rho*A*Izz,        # Rotational mass around y
    rho*A*Iyy         # Rotational mass around z
])

M3 = np.diag([
    rho*A,            # Axial mass
    rho*A,            # Transverse mass y
    rho*A,            # Transverse mass Z
    rho*A*Ixx,        # Rotational mass around x
    rho*A*Izz,        # Rotational mass around y
    rho*A*Iyy         # Rotational mass around z
])

## Input stiffness matrices (6x6 sectional stiffness for different composite sections)
#M1 = np.array([
#    [8.860789025763818,0.0,0.0,0.0,-0.00027421869856204383,1.1601830607332886e-14],        # Axial mass
#    [0.0,8.860789025763818,0.0,0.00027421869856204383,0.0,0.0],        # Transverse mass y
#    [0.0,0.0,8.860789025763818,-1.1601830607332886e-14,0.0,0.0],        # Transverse mass z
#    [0.0,0.00027421869856204383,-1.1601830607332886e-14,14936.320926087838,0.0,0.0],        # Rotational mass around x
#    [-0.00027421869856204383,0.0,0.0,0.0,3694.5307479602657,1.2370104940373494e-12],        # Rotational mass around y
#    [1.1601830607332886e-14,0.0,0.0,0.0,1.2370104940373494e-12,11241.790178127612]         # Rotational mass around z
#])
#
#M2 = np.array([
#    [8.860789025763818,0.0,0.0,0.0,-0.00027421869856204383,1.1601830607332886e-14],        # Axial mass
#    [0.0,8.860789025763818,0.0,0.00027421869856204383,0.0,0.0],        # Transverse mass y
#    [0.0,0.0,8.860789025763818,-1.1601830607332886e-14,0.0,0.0],        # Transverse mass z
#    [0.0,0.00027421869856204383,-1.1601830607332886e-14,14936.320926087838,0.0,0.0],        # Rotational mass around x
#    [-0.00027421869856204383,0.0,0.0,0.0,3694.5307479602657,1.2370104940373494e-12],        # Rotational mass around y
#    [1.1601830607332886e-14,0.0,0.0,0.0,1.2370104940373494e-12,11241.790178127612]         # Rotational mass around z
#])
#
#M3 = np.array([
#    [8.860789025763818,0.0,0.0,0.0,-0.00027421869856204383,1.1601830607332886e-14],        # Axial mass
#    [0.0,8.860789025763818,0.0,0.00027421869856204383,0.0,0.0],        # Transverse mass y
#    [0.0,0.0,8.860789025763818,-1.1601830607332886e-14,0.0,0.0],        # Transverse mass z
#    [0.0,0.00027421869856204383,-1.1601830607332886e-14,14936.320926087838,0.0,0.0],        # Rotational mass around x
#    [-0.00027421869856204383,0.0,0.0,0.0,3694.5307479602657,1.2370104940373494e-12],        # Rotational mass around y
#    [1.1601830607332886e-14,0.0,0.0,0.0,1.2370104940373494e-12,11241.790178127612]         # Rotational mass around z
#])

def create_beam_model():
    """
    Create the beam model with nodes, elements and stiffness assignments
    """
    # Define the beam portions for each stiffness matrix
    first_portion_end = int(n_nodes * 0.33)
    second_portion_end = int(n_nodes * 0.66)

    first_portion = list(range(first_portion_end))
    second_portion = list(range(first_portion_end, second_portion_end))
    third_portion = list(range(second_portion_end, n_nodes))

    print(f"Beam portions: {len(first_portion)} nodes in first, {len(second_portion)} in second, {len(third_portion)} in third")

    # Map nodes to stiffness matrices
    nodes_stiffness = []
    for i in range(n_nodes):
        if i < first_portion_end - 1:
            nodes_stiffness.append(K1)
        elif i < second_portion_end - 1:
            nodes_stiffness.append(K2)
        else:
            nodes_stiffness.append(K3)

    # Map nodes to mass matrices
    nodes_mass = []
    for i in range(n_nodes):
        if i < first_portion_end - 1:
            nodes_mass.append(M1)
        elif i < second_portion_end - 1:
            nodes_mass.append(M2)
        else:
            nodes_mass.append(M3)


    # Generate node positions
    nodes = []
    for i in range(n_nodes):
        x = i * element_length
        nodes.append({"position": [x, 0, 0], "index": i})

    # Define element connectivity
    elements = []
    for i in range(n_el):
        elements.append({"nodes": [i, i+1], "stiffness":nodes_stiffness[i], "mass": nodes_mass[i]})
    
    # Return model details
    return {
        "nodes": nodes,
        "elements": elements,
        "first_portion_end": first_portion_end,
        "second_portion_end": second_portion_end
    }

beam_model = create_beam_model()
print(f'Stiffness matrix for comparison : {beam_model["elements"][0]["stiffness"]}')