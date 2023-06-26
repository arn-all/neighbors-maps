# Neighbors maps

Companion repository for Allera et al. (2023).


**Example**

```py
import ase.lattice
import matplotlib.pyplot as plt

from neighbor_maps import Atoms

p = Atoms(ase.lattice.cubic.BodyCenteredCubic("Fe", size=(4,4,4)))
plt.imshow(p.get_descriptor(atom_id=0, r_cut=6.0, img_target_size=32))
plt.axis("off")
```
![gist_example](https://user-images.githubusercontent.com/45487966/243297715-8303bd6d-6199-40ef-b3bd-89984103183d.png)
