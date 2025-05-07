import numpy as np
import matplotlib.pyplot as plt
from beam_properties import dof_per_node


def plot_deformed_shape(nodes, displacements, section_breaks, scale_factor=1):
    """
    Plot the deformed shape of the beam
    
    Args:
        nodes: List of node dictionaries with positions
        displacements: List of displacement dictionaries
        section_breaks: List of x-positions where beam sections change
        scale_factor: Scale factor for displaying deformations
    """
    # Extract original positions
    x_orig = np.array([node["position"][0] for node in nodes])
    y_orig = np.array([node["position"][1] for node in nodes])
    z_orig = np.array([node["position"][2] for node in nodes])
    
    # Extract displacements
    x_disp = np.array([disp["x"] for disp in displacements])
    y_disp = np.array([disp["y"] for disp in displacements])
    z_disp = np.array([disp["z"] for disp in displacements])
    
    # Calculate deformed positions
    x_def = x_orig + scale_factor * x_disp
    y_def = y_orig + scale_factor * y_disp
    z_def = z_orig + scale_factor * z_disp
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Plot original and deformed shape
    plt.plot(x_orig, y_orig, 'k-', label='Original')
    plt.plot(x_def, y_def, 'r-', label=f'Deformed (scale={scale_factor})')
    
    # Mark different sections
    for sb in section_breaks:
        plt.axvline(x=sb, color='g', linestyle='--', alpha=0.5)
    
    plt.grid(True)
    plt.xlabel('X position (m)')
    plt.ylabel('Y position (m)')
    plt.title('Beam Deformation under Vertical Load')
    plt.legend()
    
    plt.figure(figsize=(12, 6))
    plt.plot(x_orig, x_disp, label='X displacement')
    plt.plot(x_orig, y_disp, label='Y displacement')
    plt.plot(x_orig, z_disp, label='Z displacement')
    plt.grid(True)
    plt.xlabel('X position (m)')
    plt.ylabel('Displacement (m)')
    plt.title('Displacements along the beam')
    plt.legend()
    
    plt.show()


def plot_mode_shapes(nodes, mode_shapes, frequencies, num_modes_to_plot=4):
    """
    Plot the first few mode shapes of the beam
    
    Args:
        nodes: List of node dictionaries with positions
        mode_shapes: Mode shapes from modal analysis
        frequencies: Natural frequencies in Hz
        num_modes_to_plot: Number of modes to plot
    """
    x_orig = np.array([node["position"][0] for node in nodes])
    
    plt.figure(figsize=(12, 8))
    for i in range(min(num_modes_to_plot, len(frequencies))):
        # Extract y-direction displacements for this mode
        mode_y = np.array([mode_shapes[n*dof_per_node + 1, i] for n in range(len(nodes))])
        
        # Normalize the mode shape
        mode_y = mode_y / np.max(np.abs(mode_y))
        
        plt.subplot(num_modes_to_plot, 1, i+1)
        plt.plot(x_orig, mode_y)
        plt.grid(True)
        plt.ylabel(f'Mode {i+1}')
        plt.title(f'f = {frequencies[i]:.2f} Hz')
        
        if i < num_modes_to_plot - 1:
            plt.xticks([])
    
    plt.xlabel('X position (m)')
    plt.tight_layout()
    plt.show()