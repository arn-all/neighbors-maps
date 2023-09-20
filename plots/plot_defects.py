from ase.build import bulk
from neighbors_map import Atoms
from ase import Atom
import matplotlib.pyplot as plt

def get_cutoff(atoms, n_neightbours, cutoff=10.):
    """
    Function to get the distance of n_th neighbour using matscipy neighbour list.
    """
    from matscipy.neighbours import neighbour_list

    i, d = neighbour_list("id", atoms, cutoff=cutoff)
    first_atom_d = d[i == 0]
    first_atom_d.sort()
    return first_atom_d[n_neightbours]

symbols = {"Simple Cubic": "Po", 
           "BCC": "Fe", 
           "FCC": "Ni", 
           "Diamond": "C"}

unit_cells = {}

for structure_type, symbol in symbols.items():
    unit_cells[structure_type] = bulk(symbol, cubic=True)    

bulks = {}

for structure_type, unit_cell in unit_cells.items():
    bulks[structure_type] = Atoms(unit_cell * [4, 4, 4])

for structure_type, symbol in symbols.items():
    unit_cells[structure_type] = bulk(symbol, cubic=True)    

vacancies = {}
interstitials = {}

for structure_type, bulk in bulks.items():
    vacancy = bulk.copy()
    # delete atom with index 1
    del vacancy[1]
    vacancies[structure_type] = Atoms(vacancy)

    # add an atom
    interstitial = bulk.copy()
    interstitial.extend(Atom(symbol, (2.5, 2.5, 1.5)))
    interstitials[structure_type] = Atoms(interstitial)

fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(10, 8))

for (structure_type, structure), ax in zip(bulks.items(), axes[0]):

    ax.imshow(structure.get_neighbor_map(atom_id=0, r_cut=1.2 * get_cutoff(structure, 32), img_target_size=32))
    ax.set_title("Bulk " + structure_type)
    ax.axis("off")

for (structure_type, structure), ax in zip(vacancies.items(), axes[1]):

    ax.imshow(structure.get_neighbor_map(atom_id=0, r_cut=1.2 * get_cutoff(structure, 32), img_target_size=32))
    ax.set_title("Vacancy in " + structure_type)
    ax.axis("off")


for (structure_type, structure), ax in zip(interstitials.items(), axes[2]):

    ax.imshow(structure.get_neighbor_map(atom_id=0, r_cut=1.2 * get_cutoff(structure, 32), img_target_size=32))
    ax.set_title("Interstitial in " + structure_type)
    ax.axis("off")

    
fig.tight_layout()
plt.savefig("plots/pdf/plot_defects.pdf")
plt.savefig("plots/pdf/plot_defects.png", dpi=300)
plt.show()