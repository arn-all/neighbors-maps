# coding: utf-8
import ase.lattice
import ase.lattice.tetragonal
import ase.lattice.monoclinic
import ase.lattice.triclinic
import ase.lattice.hexagonal
import ase.lattice.orthorhombic
import matplotlib.pyplot as plt
import math

import ase
from neighbors_map import Atoms

# List of Bravais lattices with lattice constants
bravais_lattices = [
    ("cubic", "SimpleCubic", 3.0),
    ("cubic", "FaceCenteredCubic", 3.0),
    ("cubic", "BodyCenteredCubic", 3.0),
    ("cubic", "Diamond", 3.0),
    ("tetragonal", "SimpleTetragonal", {'a': 3.0, 'c/a': 1.5}),
    ("tetragonal", "CenteredTetragonal", {'a': 3.0, 'c/a': 1.5}),
    ("orthorhombic", "SimpleOrthorhombic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5}),
    ("orthorhombic", "BaseCenteredOrthorhombic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5}),
    ("orthorhombic", "FaceCenteredOrthorhombic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5}),
    ("orthorhombic", "BodyCenteredOrthorhombic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5}),
    ("monoclinic", "SimpleMonoclinic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5, 'alpha': 90}),
    ("monoclinic", "BaseCenteredMonoclinic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5, 'alpha': 90}),
    ("triclinic", "Triclinic", {'a': 3.0, 'b/a': 1.2, 'c/a': 1.5, 'alpha': 90, 'beta': 90, 'gamma': 90}),
    ("hexagonal", "Hexagonal", {'a': 3.0, 'c/a': 1.5}),
    ("hexagonal", "HexagonalClosedPacked", {'a': 3.0, 'c/a': 1.5}),
    ("hexagonal", "Graphite", {'a': 3.0, 'c/a': 1.5})
]

# Function to generate the descriptor for a given lattice type
def generate_descriptor(lattice_module, lattice_class, element, size, img_target_size, lattice_constants):
    lattice = getattr(getattr(ase.lattice, lattice_module), lattice_class)(element, size=size, latticeconstant=lattice_constants)
    atoms = Atoms(lattice)
    central_atom = atoms.tree.query(atoms.positions.mean(axis=0))[1]
    d_max = atoms.tree.query(atoms.positions[0], k=128)[0][-1]
    print(lattice_class, f'{d_max:.2f}')

    descriptor = atoms.get_neighbor_map(atom_id=central_atom, r_cut=d_max, img_target_size=img_target_size)
    return descriptor

# Parameters
element = "Fe"  # Use an element compatible with the lattice constants
size = (12, 12, 12)
img_target_size = 32

# Calculate the number of subplots needed
num_subplots = len(bravais_lattices)
num_rows = int(math.ceil(num_subplots / 4))
fig, axes = plt.subplots(num_rows, 4, figsize=(14, 4 * num_rows))

# Generate and plot descriptors for all Bravais lattices
for i, (lattice_type, lattice_class, lattice_constants) in enumerate(bravais_lattices):
    descriptor = generate_descriptor(lattice_type, lattice_class, element, size, img_target_size, lattice_constants)
    
    row = i // 4
    col = i % 4
    
    # If there's only one row, axes is a 1D array
    if num_rows == 1:
        axes[col].imshow(descriptor)
        axes[col].set_title(f"{lattice_class}", size=16)
        axes[col].axis("off")
    else:
        axes[row, col].imshow(descriptor)
        axes[row, col].set_title(f"{lattice_class}", size=16)
        axes[row, col].axis("off")

# Hide any remaining empty subplots
for i in range(num_subplots, num_rows * 4):
    row = i // 4
    col = i % 4
    
    if num_rows == 1:
        axes[col].axis("off")
    else:
        axes[row, col].axis("off")

plt.tight_layout()
plt.savefig("bravais_lattices.pdf")
