import os

import numpy as np
from Bio.PDB import PDBParser
import pyvista as pv


# Load PDB file
# pdb_file = "./data/pdb/example.pdb"  # Replace with your PDB file path
pdb_file = "./data/pdb/1EN2.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("structure", pdb_file)


# Parameters for visualization
cutoff = 30.0  # Cutoff distance in Angstroms
resolution = 100  # Resolution for 3D grid

scale_factor = 0.8  # Scale factor for van der Waals radii

atom_vdw_radii = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "H": 1.2,  # Add H for completeness
}
# Toggle for heavy atoms only
include_only_heavy_atoms = True  # Set to True to exclude H, False to include all atoms

# Extract atomic positions and van der Waals radii
positions = []
vdw_radii = []

# Example van der Waals radii values (can be replaced with a proper mapping)
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if include_only_heavy_atoms and atom.element == "H":
                    continue  # Skip hydrogen atoms if toggle is on
                positions.append(atom.coord)
                vdw_radii.append(
                    atom_vdw_radii.get(atom.element, 1.5) * scale_factor
                )  # Apply scale factor
positions = np.array(positions)
vdw_radii = np.array(vdw_radii)

# Center positions
center = np.mean(positions, axis=0)
positions -= center


def generate_vdw_surface(positions, vdw_radii, cutoff, resolution):
    """Generate a van der Waals surface representation."""
    # Create a 3D grid
    grid = pv.ImageData()  # Changed from UniformGrid to ImageData
    grid.dimensions = [resolution, resolution, resolution]
    grid.origin = [-cutoff, -cutoff, -cutoff]
    grid.spacing = [2 * cutoff / resolution] * 3  # Equal spacing in all dimensions

    # Initialize the density grid
    density = np.zeros(grid.dimensions)

    # Create grid points
    x = np.linspace(-cutoff, cutoff, resolution)
    y = np.linspace(-cutoff, cutoff, resolution)
    z = np.linspace(-cutoff, cutoff, resolution)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    # Calculate density
    for pos, radius in zip(positions, vdw_radii):
        distances = np.linalg.norm(grid_points - pos, axis=1)
        density += np.exp(-((distances / radius) ** 2)).reshape(
            resolution, resolution, resolution
        )

    # Add the density data to the grid
    grid.point_data["values"] = density.flatten(
        order="F"
    )  # PyVista expects Fortran order
    return grid


# Generate the van der Waals surface
grid = generate_vdw_surface(positions, vdw_radii, cutoff, resolution=resolution)

# Create isosurface directly using PyVista
surface = grid.contour([np.max(grid.point_data["values"]) * 0.1])

# Optional: Smooth and simplify the surface
surface.smooth(n_iter=100, relaxation_factor=0.1)
surface.decimate(0.5)

# Create a plotter with optimized settings
plotter = pv.Plotter(
    lighting="three lights", off_screen=True
)  # off_screen=True for screenshot
plotter.add_mesh(
    surface,
    color="gray",
    smooth_shading=True,
    specular=0.5,
    ambient=0.3,
    diffuse=0.8,
    show_edges=False,
)

# Optimize rendering performance
plotter.enable_anti_aliasing()
plotter.enable_eye_dome_lighting()
plotter.camera.parallel_projection = True
plotter.camera_position = "xy"
plotter.camera.zoom(1.5)

# Add axes
plotter.add_axes()

# # Show the plot
# plotter.show()

# Save the plot
os.makedirs("./results/images/pdb", exist_ok=True)
plotter.screenshot("./results/images/pdb/1EN2_vdw.png")
