import MDAnalysis as mda 
import pandas as pd
from matplotlib import pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np

# Loading a structure or trajectory

u = mda.Universe("4pti.psf", "4pti_eq.dcd", "4pti_prod.dcd")

# Hbonds
print("hbonds selection")

hbonds = HydrogenBondAnalysis(universe=u)

protein_hydrogens_sel = hbonds.guess_hydrogens("protein")
protein_acceptors_sel = hbonds.guess_acceptors("protein")

water_hydrogens_sel = "resname TIP3 and name H1 H2"
water_acceptors_sel = "resname TIP3 and name OH2"

hbonds.hydrogens_sel = f"({protein_hydrogens_sel}) or ({water_hydrogens_sel} and around 10 not resname TIP3)"
hbonds.acceptors_sel = f"({protein_acceptors_sel}) or ({water_acceptors_sel} and around 10 not resname TIP3)"

print(f"hydrogen_sel = {protein_hydrogens_sel}")
print(f"acceptors_sel = {protein_acceptors_sel}")

hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel=None,
    hydrogens_sel=hbonds.hydrogens_sel,
    acceptors_sel=hbonds.acceptors_sel,
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=150,
    update_selections=False
)

hbonds.run(start=15,verbose=True)

print(hbonds.results.hbonds.shape)
hbonds.results.hbonds.dtype
n_hbond = hbonds.results.hbonds

print("end of hbonds calculation")

n_hbonds = hbonds.count_by_time()

np.save("4pti_hbonds.npy", n_hbond)

# Counts the number of hydrogen bonds for each frame
plt.plot(hbonds.times, n_hbonds, lw=2)

plt.title("Number of hydrogen bonds over time", weight="bold")
plt.xlabel("Time (ps)")
plt.ylabel(r"$N_{HB}$")
plt.grid(True)
plt.savefig("hbonds.jpg")

plt.show()