import ase
from scipy.spatial import KDTree
import numpy as np
import warnings

class Atoms(ase.Atoms):
    """
    Implements the get_descriptor() attribute to ase.Atoms. 
    It returns the Neighbor Map descriptor of an atom, as described in Allera et al. (2023).
    Notes on current implementation: 
        * Periodic boundary conditions are applied by default for cubic cells. 
          Non-cubic cells will fallback to non-periodic implementation. 
          If the atom is close to the cell boundary it will affect the result. 
        * Limited to single-element systems.
    """

    def __init__(self, *args, **kwargs):
        """
        Inheritance from ase.Atoms
        """
        super().__init__(*args, **kwargs)
        self.wrap()
        self.tree = self.build_periodic_kdtree(self.positions, boxsize=self.cell.lengths())

    @staticmethod
    def build_periodic_kdtree(positions, boxsize):
        try:
          return KDTree(positions, boxsize=boxsize)
        except ValueError:
          warnings.warn("Neighbor Map: Non cubic cell. Falling back to non-periodic implementation.")
          return KDTree(positions)
        
    @staticmethod
    def pad(x: np.ndarray, image_size: int, target_size: int) -> np.ndarray:
        """
        Pad the input numpy array with zeroes on the right and bottom to match the target size.
        x: Numpy array to be padded.
        image_size: Current size of the image.
        target_size: Desired size of the image.
        Returns: Padded numpy array.
        """
        return np.pad(x, ((0, target_size-image_size), (0, target_size-image_size)))

    def get_descriptor(self, atom_id: int, r_cut: float, img_target_size: int = 32) -> np.ndarray:
        """
        Compute the descriptor for the atom at the given ID.
        atom_id: ID of the atom for which to compute the descriptor.
        r_cut: Radius for limiting neighbors lookup.
        img_target_size: Target size of the image.
        Returns: Descriptor of the atom.
        Example:
        >>> import ase.lattice
        >>> import matplotlib.pyplot as plt
        >>> p = Atoms(ase.lattice.cubic.BodyCenteredCubic("Fe", size=(4,4,4)))
        >>> plt.imshow(p.get_descriptor(atom_id=0, r_cut=6.0, img_target_size=32))
        """
        atom = self.positions[atom_id]
        
        # Limit neighbors lookup to a ball of radius r_cut around the central atom 
        ids_in_rcut = self.tree.query_ball_point(atom, r_cut)
        atoms_in_rcut = self.positions[ids_in_rcut]
        number_atoms_in_rcut = atoms_in_rcut.shape[0]
        tree_in_ball = self.build_periodic_kdtree(atoms_in_rcut,  boxsize=self.cell.lengths())

        # to avoid returning the atom itself and avoid returning more atoms than there are in the ball
        k_n = np.arange(2, min(img_target_size, len(ids_in_rcut))+2) 

        # for first row of matrix:
        d_to_central, id_neigh_central = tree_in_ball.query(atom, k=k_n)

        # for the next n_neigh-1 rows
        d_to_neigh_j, id_neigh_of_neigh_j = tree_in_ball.query(atoms_in_rcut[id_neigh_central[:-1]], k=k_n)
        v = np.vstack((d_to_central,d_to_neigh_j))

        # Weights: [1, d_0j**m, ... ]
        img_size = min(number_atoms_in_rcut, img_target_size)
        m = -1/3 * np.ones(img_size)
        A = np.concatenate(([1], v[0, :-1]))**m
        
        # transform into a column vector so weights are applied row-wise
        A = A.reshape((-1,1))

        # pad with zeroes on the right and bottom
        return self.pad(np.sqrt(A/v**2), img_size, img_target_size)
