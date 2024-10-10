import MDAnalysis as mda
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import math as m

u = mda.Universe("4pti.psf", "4pti_eq.dcd")
print(u)
print(len(u.trajectory))

previous_positions = None
r = []
velocities = []
for ts in u.trajectory[:50]:
    current_positions = u.atoms.positions

    if previous_positions is not None:
        l = current_positions - previous_positions
        r.append(l)
    previous_positions = current_positions.copy()
r_np = np.array(r)
distances = np.linalg.norm(r_np, axis=2)
v = distances/2 #dist in Angstroem/fs
velocities.append(v)
print(velocities)