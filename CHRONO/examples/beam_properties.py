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

element_length = beamimport numpy as np

# Beam parameters
n_el = 100
n_nodes = n_el + 1
from_inches_to_mm = 25.4

steel_box_beam = False

if steel_box_beam == True:
    beam_length = 2000
    b = 100
    h = 100
    n_plyes = 1
    ply_thickness = 10
    rho = 7800e-9
    E = 210000000000e-7
    G = 80769230769.23077e-7
    nu = 0.3                  # Poisson's ratio

else:
    beam_length = 30*from_inches_to_mm   # mm
    b = 0.953*from_inches_to_mm              # mm
    h = 0.537*from_inches_to_mm              # mm
    n_plyes = 6
    ply_thickness = 0.005*from_inches_to_mm
    rho = 1 
    E = 144762.02644794726            # Young Modulus SI kg,mm
    G = 6116.705342871011
    nu = 0.42                  # Poisson's ratio

s = ply_thickness*n_plyes               # mm
b_int = b-2*s        # mm
h_int = h-2*s        # mm
A = b * h -((b-2*s)*(h-2*s)) # mm^2
Iyy = (b * h**3 - (b_int)*(h_int)**3) / 12 # mm^4 
Izz = (b**3 * h - (b_int)**3*(h_int)) / 12 # mm^4 
Ixx = Iyy + Izz

sectional_mass = rho * A

element_length = beam_length / n_el
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

# example case COBRA with 6plies (15)6
K1 = np.array([
    [5647609.943325554,    0.7344734729634163,-0.013125052835588003,-5996718.293929776,72.63992870172287,267.7663925945215],
    [0.7344734641663502,   445159.81875698175,-20.553886404203837,-1.6027123957783094,3038463.4295200594,-1593.194778146545],
    [-0.013125057727511767,-20.553886403233136,192355.74932248812,-3.4295366441119155,329.94794714224264,3203794.916084049],
    [-5996718.293930015,   -1.6027124691283574,-3.429536616475426,51810554.03693702,-21.67275858924122,-80.2091344840405],
    [72.6399288658993,     3038463.4295201497,329.9479470686943,-21.672759088486035,179634073.53118718,-4144.554411196352],
    [267.766392657707,     -1593.1947780183302,3203794.916083931,-80.20913561431072,-4144.554411406311,444285081.35465825]
])

K2 = K1
K3 = K1

# ======= MASS MATRICES ======== #

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

# Input stiffness matrices (6x6 sectional stiffness for different composite sections)
M1 = np.array([
    [55.35444882732644,0.0,0.0,0.0,0.0006865039404114752,0.0016752440667648291],        # Axial mass
    [0.0,55.35444882732644,0.0,-0.0006865039404114752,0.0,0.0],        # Transverse mass y
    [0.0,0.0,55.35444882732644,-0.0016752440667648291,0.0,0.0],        # Transverse mass z
    [0.0,-0.0006865039404114752,-0.0016752440667648291,6096.399021178555,0.0,0.0],        # Rotational mass around x
    [0.0006865039404114752,0.0,0.0,0.0,1757.0920258485057,-0.011118257653655535],        # Rotational mass around y
    [0.0016752440667648291,0.0,0.0,0.0,-0.011118257653655535,4339.306995330096]         # Rotational mass around z
])

M2 = np.array([
    [55.35444882732644,0.0,0.0,0.0,0.0006865039404114752,0.0016752440667648291],        # Axial mass
    [0.0,55.35444882732644,0.0,-0.0006865039404114752,0.0,0.0],        # Transverse mass y
    [0.0,0.0,55.35444882732644,-0.0016752440667648291,0.0,0.0],        # Transverse mass z
    [0.0,-0.0006865039404114752,-0.0016752440667648291,6096.399021178555,0.0,0.0],        # Rotational mass around x
    [0.0006865039404114752,0.0,0.0,0.0,1757.0920258485057,-0.011118257653655535],        # Rotational mass around y
    [0.0016752440667648291,0.0,0.0,0.0,-0.011118257653655535,4339.306995330096]         # Rotational mass around z
])

M3 = np.array([
    [55.35444882732644,0.0,0.0,0.0,0.0006865039404114752,0.0016752440667648291],        # Axial mass
    [0.0,55.35444882732644,0.0,-0.0006865039404114752,0.0,0.0],        # Transverse mass y
    [0.0,0.0,55.35444882732644,-0.0016752440667648291,0.0,0.0],        # Transverse mass z
    [0.0,-0.0006865039404114752,-0.0016752440667648291,6096.399021178555,0.0,0.0],        # Rotational mass around x
    [0.0006865039404114752,0.0,0.0,0.0,1757.0920258485057,-0.011118257653655535],        # Rotational mass around y
    [0.0016752440667648291,0.0,0.0,0.0,-0.011118257653655535,4339.306995330096]         # Rotational mass around z
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
