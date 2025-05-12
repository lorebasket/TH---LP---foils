import numpy as np

# Beam parameters
n_el = 100
n_nodes = n_el + 1
beam_length = 2000   # mm
b = 100              # mm
h = 100              # mm
s = 10               # mm
b_int = b-2*s        # mm
h_int = h-2*s        # mm
A = b * h -((b-2*s)*(h-2*s)) # mm^2
Iyy = (b * h**3 - (b_int)*(h_int)**3) / 12 # mm^4 
Izz = (b**3 * h - (b_int)**3*(h_int)) / 12 # mm^4 

rho = 7800e-6             # density [kg/mm^3]
sectional_mass = rho * A
E = 0.210e6               # Young Modulus SI kg,mm
nu = 0.3                  # Poisson's ratio
G = E /(2*(1+nu))         # Shear ratio
Ixx = Iyy + Izz

element_length = beam_length / n_el
phi = (12 * E * Iyy) / (G * A * element_length**2) # Note: This phi might not be directly used in Timoshenko element formulation
dof_per_node = 6  # 3 displacement + 3 rotation
dof_per_element = 12  # 6 DOF per node * 2 nodes per element
total_dof = (n_nodes) * dof_per_node

## Theoretical Timoshenko stiffness matrix to test the script validity
#K1 = np.array([
#    [E*A, 0, 0, 0, 0, 0],        # Axial stiffness (εx)
#    [0, G*A, 0, 0, 0, 0],        # Shear stiffness (γxy) 
#    [0, 0, G*A, 0, 0, 0],        # Shear stiffness (γxz)
#    [0, 0, 0, G*Ixx, 0, 0],      # Torsional stiffness G*J (κx) - using Ixx as approx for J
#    [0, 0, 0, 0, E*Izz, 0],      # Bending stiffness Mz (κz) - related to v displacement
#    [0, 0, 0, 0, 0, E*Iyy]       # Bending stiffness My (κy) - related to w displacement 
#])
#
#K2= np.array([
#    [E*A, 0, 0, 0, 0, 0],  
#    [0, G*A, 0, 0, 0, 0],  
#    [0, 0, G*A, 0, 0, 0],  
#    [0, 0, 0, G*Ixx, 0, 0],
#    [0, 0, 0, 0, E*Izz, 0],
#    [0, 0, 0, 0, 0, E*Iyy] 
#])
#
#K3 = np.array([
#    [E*A, 0, 0, 0, 0, 0],          
#    [0, G*A, 0, 0, 0, 0],    
#    [0, 0, G*A, 0, 0, 0],    
#    [0, 0, 0, G*Ixx, 0, 0],  
#    [0, 0, 0, 0, E*Izz, 0],  
#    [0, 0, 0, 0, 0, E*Iyy]    
#])

# Input stiffness matrices from SONATA
K1 = np.array([
    [756962212.3838665,0.0,0.0,0.0,28409144.149016462,-9999473.55783619],
    [0.0,132674369.51376443,15590.244307241333,-12599172.913559467,0.0,0.0],
    [0.0,15590.244308950809,132587213.76207954,4227114.887939889,0.0,0.0],
    [0.0,-12599172.913557041,4227114.887936292,638797135601.2262,0.0,0.0],
    [28409144.148775954,0.0,0.0,0.0,1035861823584.9604,-48678841.60016193],
    [-9999473.557804847,0.0,0.0,0.0,-48678841.6030593,1035711351955.2004]
])
K2 = np.array([
    [756962212.3838665,0.0,0.0,0.0,28409144.149016462,-9999473.55783619],
    [0.0,132674369.51376443,15590.244307241333,-12599172.913559467,0.0,0.0],
    [0.0,15590.244308950809,132587213.76207954,4227114.887939889,0.0,0.0],
    [0.0,-12599172.913557041,4227114.887936292,638797135601.2262,0.0,0.0],
    [28409144.148775954,0.0,0.0,0.0,1035861823584.9604,-48678841.60016193],
    [-9999473.557804847,0.0,0.0,0.0,-48678841.6030593,1035711351955.2004]
])
K3 = np.array([
    [756962212.3838665,0.0,0.0,0.0,28409144.149016462,-9999473.55783619],
    [0.0,132674369.51376443,15590.244307241333,-12599172.913559467,0.0,0.0],
    [0.0,15590.244308950809,132587213.76207954,4227114.887939889,0.0,0.0],
    [0.0,-12599172.913557041,4227114.887936292,638797135601.2262,0.0,0.0],
    [28409144.148775954,0.0,0.0,0.0,1035861823584.9604,-48678841.60016193],
    [-9999473.557804847,0.0,0.0,0.0,-48678841.6030593,1035711351955.2004]
])


## Input stiffness matrices (6x6 sectional stiffness for different composite sections)
#M1 = np.diag([
#    rho*A,            # Axial mass
#    rho*A,            # Transverse mass y
#    rho*A,            # Transverse mass Z
#    rho*Ixx,        # Rotational mass around x
#    rho*Izz,        # Rotational mass around y
#    rho*Iyy         # Rotational mass around z
#])
#
#M2 = np.diag([
#    rho*A,            # Axial mass
#    rho*A,            # Transverse mass y
#    rho*A,            # Transverse mass Z
#    rho*Ixx,        # Rotational mass around x
#    rho*Izz,        # Rotational mass around y
#    rho*Iyy         # Rotational mass around z
#])
#
#M3 = np.diag([
#    rho*A,            # Axial mass
#    rho*A,            # Transverse mass y
#    rho*A,            # Transverse mass Z
#    rho*Ixx,        # Rotational mass around x
#    rho*Izz,        # Rotational mass around y
#    rho*Iyy         # Rotational mass around z
#])

# Input stiffness matrices from SONATA
M1 = np.array([
    [0.028115739317115272,0.0,0.0,0.0,0.0010551967826830364,-0.0003714090178660503],        # Axial mass
    [0.0,0.028115739317115272,0.0,-0.0010551967826830364,0.0,0.0],        # Transverse mass y
    [0.0,0.0,0.028115739317115272,0.0003714090178660503,0.0,0.0],        # Transverse mass z
    [0.0,-0.0010551967826830364,0.0003714090178660503,76.8243449915628,0.0,0.0],        # Rotational mass around x
    [0.0010551967826830364,0.0,0.0,0.0,38.415258246604175,-0.0008645649916275387],        # Rotational mass around y
    [-0.0003714090178660503,0.0,0.0,0.0,-0.0008645649916275387,38.409086744958636]         # Rotational mass around z
])

M2 = np.array([
    [0.028115739317115272,0.0,0.0,0.0,0.0010551967826830364,-0.0003714090178660503],        # Axial mass
    [0.0,0.028115739317115272,0.0,-0.0010551967826830364,0.0,0.0],        # Transverse mass y
    [0.0,0.0,0.028115739317115272,0.0003714090178660503,0.0,0.0],        # Transverse mass z
    [0.0,-0.0010551967826830364,0.0003714090178660503,76.8243449915628,0.0,0.0],        # Rotational mass around x
    [0.0010551967826830364,0.0,0.0,0.0,38.415258246604175,-0.0008645649916275387],        # Rotational mass around y
    [-0.0003714090178660503,0.0,0.0,0.0,-0.0008645649916275387,38.409086744958636]         # Rotational mass around z
])

M3 = np.array([
    [0.028115739317115272,0.0,0.0,0.0,0.0010551967826830364,-0.0003714090178660503],        # Axial mass
    [0.0,0.028115739317115272,0.0,-0.0010551967826830364,0.0,0.0],        # Transverse mass y
    [0.0,0.0,0.028115739317115272,0.0003714090178660503,0.0,0.0],        # Transverse mass z
    [0.0,-0.0010551967826830364,0.0003714090178660503,76.8243449915628,0.0,0.0],        # Rotational mass around x
    [0.0010551967826830364,0.0,0.0,0.0,38.415258246604175,-0.0008645649916275387],        # Rotational mass around y
    [-0.0003714090178660503,0.0,0.0,0.0,-0.0008645649916275387,38.409086744958636]         # Rotational mass around z
])


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
