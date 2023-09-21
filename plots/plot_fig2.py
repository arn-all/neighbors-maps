import numpy as np
import ase.lattice
import ase.lattice.hexagonal
import matplotlib.pyplot as plt
from scipy.spatial import KDTree, distance
import scipy.spatial
import neighbors_map

from scipy.stats import qmc
radius = 0.08
engine = qmc.PoissonDisk(d=3, radius=radius, seed=1)
poisson_data = engine.random(1000)*10

# We try two different levels of noise
noise_amplitude = [0.0, 0.15]

normal_cell = (8, 8, 8)
small_cell = (6, 6, 6)

# We define several typical strucutres
systems = [ase.lattice.cubic.BodyCenteredCubic("Fe", size=normal_cell),
          ase.lattice.cubic.FaceCenteredCubic("Al", size=small_cell),
          ase.lattice.hexagonal.HexagonalClosedPacked(size=normal_cell,
                                 symbol="Ti", pbc=(1,1,1)),
          ase.lattice.cubic.Diamond(directions=[[1,0,0], [0,1,0], [0,0,1]],
                                 size=small_cell, symbol="Si", pbc=(1,1,1)),
          ase.lattice.cubic.SimpleCubic(directions=[[1,0,0], [0,1,0], [0,0,1]],
                                 size=normal_cell, symbol="Au",
                                 pbc=(1,1,1), latticeconstant=2.0)]

atoms_data =  [s.positions for s in systems] + [poisson_data]
phases = ["bcc", "fcc", "hcp", "diamond", "simple cubic", "non-crystalline"]

### FIGURE
# We try two different levels of noise
noise_amplitude = [0.0, 0.15]

def get_centre(sample):
  t = scipy.spatial.KDTree(sample)
  _, i = t.query(sample.mean(axis=0))
  return i

def plot_neighbor_maps(systems, noise_amplitude, filename):
    fig = plt.figure(figsize=(11,5))

    # create 2x1 subfigs
    subfigs = fig.subfigures(nrows=2, ncols=1)
    row_titles = ["Perfect structures (0 Kelvin)", "Structures with Gaussian noise ($\sigma = 0.15 \mathrm{\AA}$)"]

    for row_i, (row_title, subfig) in enumerate(zip(row_titles, subfigs)):
        subfig.suptitle(row_title, y=1.05, size=13)

        # create 1x6 subplots per subfig
        axs = subfig.subplots(nrows=1, ncols=6)
        for col, ax in enumerate(axs):

            for idx, s in enumerate(systems + [neighbors_map.Atoms(positions=poisson_data)]):

              np.random.seed(0)
              noise = np.random.normal(loc=0, scale=noise_amplitude[row_i], size=(s.positions.shape[0], 3))
              saved_pos = s.positions.copy()
              s.positions = s.positions+noise

              centre = get_centre(s.positions)

              d = neighbors_map.Atoms(s).get_neighbor_map(centre, r_cut=6.0)
              #d = Atoms(s).get_descriptor(0, r_cut=6.0)
              #_, _, d  = compute_descriptor_test(centre, s, 6.0)
              axs[idx].imshow(d, cmap="viridis")
              axs[idx].axis("off")
              axs[idx].set_title(phases[idx])

              s.positions = saved_pos

    fig.tight_layout() 

    plt.savefig(f"plots/pdf/{filename}.pdf", bbox_inches='tight')
    plt.savefig(f"plots/png/{filename}.png", dpi=300, bbox_inches='tight')




if __name__ == "__main__":
    plot_neighbor_maps(systems, noise_amplitude, "all_structures")