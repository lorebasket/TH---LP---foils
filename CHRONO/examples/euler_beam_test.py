import numpy as np
import pychrono as chrono
import pychrono.fea as fea

# Initialize Chrono system
system = chrono.ChSystemNSC()
system.Set_G_acc(chrono.ChVectorD(0, 0, -9.81))

# Beam parameters
length = 1.0  # Beam length [m]
n_elements = 10  # Number of finite elements
section_area = 0.01  # Cross-section area [m^2]
density = 7800  # Density [kg/m^3]

# Sectional stiffness matrix (6x6) for Euler-Bernoulli beam
# Format: [EA, GAy, GAz, GJ, EIy, EIz]
E = 2.1e11  # Young's modulus [Pa]
G = 7.9e10  # Shear modulus [Pa]
Iy = 8.33e-7  # Area moment of inertia y [m^4]
Iz = 8.33e-7  # Area moment of inertia z [m^4]
J = 1.67e-6   # Torsional constant [m^4]

K = np.diag([
    E * section_area,  # EA
    G * section_area,  # GAy
    G * section_area,  # GAz
    G * J,            # GJ
    E * Iy,           # EIy
    E * Iz             # EIz
])

# Create mesh
mesh = fea.ChMesh()
system.Add(mesh)

# Create beam section with stiffness matrix
beam_section = fea.ChBeamSectionEuler()
beam_section.SetAxialRigidity(K[0,0])
beam_section.SetTorsionRigidity(K[3,3])
beam_section.SetBendingRigidityY(K[4,4])
beam_section.SetBendingRigidityZ(K[5,5])
beam_section.SetMassPerUnitLength(density * section_area)

# Create nodes
nodes = []
for i in range(n_elements + 1):
    node = fea.ChNodeFEAxyzrot(chrono.ChFrameD(chrono.ChVectorD(i * length/n_elements, 0, 0)))
    nodes.append(node)
    mesh.AddNode(node)
    \
    # Fix first node (cantilever constraint)
    if i == 0:
        node.SetFixed(True)

# Create elements
for i in range(n_elements):
    element = fea.ChElementBeamEuler()
    element.SetNodes(nodes[i], nodes[i+1])
    element.SetSection(beam_section)
    mesh.AddElement(element)

# Apply homogeneous distributed force
force_per_unit_length = chrono.ChVectorD(0, 0, -100)  # [N/m]
for element in mesh.GetElements():
    if isinstance(element, fea.ChElementBeamEuler):
        element.SetDistributedLoad(force_per_unit_length)

# Visualization (optional)
vis = chrono.ChVisualSystemIrrlicht()
vis.AttachSystem(system)
vis.SetWindowSize(1024, 768)
vis.SetWindowTitle('Euler Beam with Distributed Load')
vis.Initialize()
vis.AddCamera(chrono.ChVectorD(0.5, -1, 0.2), chrono.ChVectorD(0.5, 0, 0))
vis.AddTypicalLights()

# Simulation loop
time_step = 1e-4
while vis.Run():
    vis.BeginScene()
    vis.Render()
    vis.EndScene()
    system.DoStepDynamics(time_step)

# Print node displacements
for i, node in enumerate(nodes):
    print(f"Node {i}: Displacement = {node.GetPos().z} m")

#Add data logging to save results
with open("beam_displacements.csv", "w") as f:
    f.write("Node,Displacement\n")
    for i, node in enumerate(nodes):
        f.write(f"{i},{node.GetPos().z}\n")
