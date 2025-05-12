import numpy as np
import time

# Import from our modules
from beam_properties import create_beam_model, n_nodes, E, Iyy, beam_length, b, h, s
from stiffness_matrix_assembly import assemble_global_stiffness_matrix 
from mass_matrix_assembly import assemble_global_mass_matrix
from analysis import solve_static_analysis, modal_analysis, extract_displacements, create_force_vector
from visualization import plot_deformed_shape, plot_mode_shapes

def main():
    print("=== Composite Cantilever Beam Analysis ===")
    
    # Create the beam model
    print("\nCreating beam model...")
    start_time = time.time()
    beam_model = create_beam_model()
    nodes = beam_model["nodes"]
    elements = beam_model["elements"]
    section_breaks = [
        nodes[beam_model["first_portion_end"]]["position"][0],
        nodes[beam_model["second_portion_end"]]["position"][0]
    ]
    print(f"Model created in {time.time() - start_time:.3f} seconds")
    
    # Assemble global stiffness matrix
    print("\nAssembling global stiffness matrix...")
    start_time = time.time()
    K_global = assemble_global_stiffness_matrix()
    print(f"Stiffness matrix assembled in {time.time() - start_time:.3f} seconds")
    
    # Assemble global mass matrix for modal analysis
    print("\nAssembling global mass matrix...")
    start_time = time.time()
    M_global = assemble_global_mass_matrix()
    print(f"Mass matrix assembled in {time.time() - start_time:.3f} seconds")
    
    # Create force vector (100N downward at the last node)
    force_vector = create_force_vector(-100.0)
    
    # Solve static analysis
    print("\nSolving static analysis...")
    start_time = time.time()
    u_full, reaction_forces = solve_static_analysis(K_global, force_vector)
    print(f"Static analysis completed in {time.time() - start_time:.3f} seconds")
    
    # Extract displacements
    print("\nExtracting displacements...")
    displacements = extract_displacements(u_full)
    
    # Print maximum displacement
    tip_disp = displacements[-1]["y"]
    print(f"Tip vertical displacement: {tip_disp:.6f} m")
    
    # Print reaction forces
    print("\nReaction forces at fixed end:")
    for i, force in enumerate(reaction_forces):
        dof_name = ["Fx", "Fy", "Fz", "Mx", "My", "Mz"][i]
        print(f"  {dof_name}: {force:.2f}")
    
    # Plot deformed shape
    print("\nPlotting deformed shape...")
    plot_deformed_shape(nodes, displacements, section_breaks)
    
    # Calculate theoretical tip displacement (Euler-Bernoulli)
    tip_force = 100.0  # Magnitude of force
    L = beam_length
    I = Iyy
    tip_theoretical_EB = tip_force * L**3 / (3 * E * I)

    # Calculate theoretical tip displacement (Timoshenko)
    # Add shear correction factor (approximate for hollow square section)
    G = E / (2 * (1 + 0.3))  # Shear modulus
    A = b * h - ((b - 2 * s) * (h - 2 * s))  # Area
    k = 0.5  # Approximate shear correction factor for hollow square
    tip_theoretical_T = tip_theoretical_EB * (1 + (E * I) / (k * G * A * L**2))

    # Print comparison
    print(f"Theoretical tip displacement (Euler-Bernoulli): {tip_theoretical_EB*1000:.4f} mm")
    print(f"Theoretical tip displacement (Timoshenko): {tip_theoretical_T*1000:.4f} mm")
    print(f"FEA tip displacement: {abs(tip_disp)*1000:.4f} mm")

    # Perform modal analysis
    print("\nPerforming modal analysis...")
    start_time = time.time()
    frequencies, mode_shapes = modal_analysis(K_global, M_global, num_modes=10)
    print(f"Modal analysis completed in {time.time() - start_time:.3f} seconds")
    
    # Print natural frequencies
    print("\nNatural frequencies:")
    for i, freq in enumerate(frequencies):
        print(f"  Mode {i+1}: {freq:.4f} Hz")
    
    # Plot mode shapes
    print("\nPlotting mode shapes...")
    plot_mode_shapes(nodes, mode_shapes, frequencies)

    
    print("\n=== Analysis completed successfully ===")

if __name__ == "__main__":
    main()
