import ase.lattice                                                       
import matplotlib.pyplot as plt                                          

from neighbors_map import Atoms

def neighbor_map(beta, gamma, r_cut=6.0, img_target_size=32):
  """Calculates the neighbor map of a given atom.

  Args:
    beta: The beta parameter.
    gamma: The gamma parameter.
    r_cut: The cutoff radius.
    img_target_size: The target size of the image.

  Returns:
    A numpy array representing the neighbor map.
  """

  p = Atoms(ase.lattice.cubic.BodyCenteredCubic("Fe", size=(4,4,4)))
  p.neighbor_map.beta = beta
  p.neighbor_map.gamma = gamma
  return p.get_neighbor_map(atom_id=0, r_cut=r_cut, img_target_size=img_target_size)

# Generate a 5x5 grid of plots.
fig, axes = plt.subplots(5, 5, figsize=(8, 8), sharex=True, sharey=True)

hyperparams = [1, 2, 3, 4, 5]

for i, beta in enumerate(hyperparams):
  for j, gamma in enumerate(hyperparams):
    axes[i, j].imshow(neighbor_map(beta, gamma))
    axes[i, j].tick_params(
                axis='both',
                which='both',
                bottom=False,
                top=False,
                labelbottom=False)
    # axes[i, j].set_title(f"beta={beta:.1f}, gamma={gamma:.1f}")

for i, h in enumerate(hyperparams):
    axes[i, 0].set_ylabel(f"beta = {h}")
    axes[-1, i].set_xlabel(f"gamma = {h}")

plt.tight_layout()
plt.savefig("plots/pdf/gamma_beta.pdf")
plt.savefig("plots/png/gamma_beta.png", dpi=300)
