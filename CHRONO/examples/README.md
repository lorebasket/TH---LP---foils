# Composite Beam Analysis

This project implements a custom finite element analysis for composite beams with orthotropic material properties. It allows specifying different sectional stiffness matrices for different portions of the beam.

## Project Structure

The code is organized into several Python modules:

- `beam_properties.py`: Defines the beam geometry, stiffness matrices, and creates the basic model
- `matrix_assembly.py`: Handles assembly of global stiffness and mass matrices
- `analysis.py`: Contains functions for static and modal analysis
- `visualization.py`: Provides functions for plotting results
- `main.py`: Main script that runs the full analysis

## How to Run

Simply run the main script:

```bash
python main.py
```

This will:
1. Create the beam model with three different sections
2. Assemble the stiffness and mass matrices
3. Apply a 100N vertical load at the free end
4. Solve for static displacements and reaction forces
5. Calculate natural frequencies and mode shapes
6. Display plots of the deformed shape and mode shapes

## Customization

To modify the beam properties, edit `beam_properties.py`:
- Change stiffness matrices to match your composite material
- Adjust beam length and number of elements
- Modify the section distribution along the beam

For different loading conditions, modify the `create_force_vector` function in `analysis.py`.

## Dependencies

- NO pychrono dependency for the moment
- NumPy
- SciPy
- Matplotlib

## Analysis Capabilities

- Static analysis with custom boundary conditions
- Modal analysis for natural frequencies and mode shapes
- Visualization of deformed shapes and mode shapes
- Support for orthotropic material properties

## N.B. Stiffness and mass matrix evaluation
- The stiffness and mass matrices are evaluated using Gaussian quadrature
- The integration is performed over the beam element
- [K] = ∫[B(x)]ᵀ[D][B(x)]dx
- [M] = ∫[N(x)]ᵀ[M][N(x)]dx
- [B(x)] is the strain-displacement matrix
- [N(x)] is the shape function matrix
- [D] is the sectional stiffness matrix
- [M] is the sectional mass matrix
