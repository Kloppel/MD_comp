import MDAnalysis as mda 
import pandas as pd
from matplotlib import pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.dihedrals import Dihedral


universe = mda.Universe("1eru.psf", "1eru_eq.dcd", "1eru_prod.dcd")

# Loading a structure or trajectory


u = mda.Universe("1eru.psf", "1eru_eq.dcd", "1eru_prod.dcd")
print(u)
print(len(u.trajectory))

# Working with trajectories
print(len(u.trajectory))

# Radius of gyration

for ts in u.trajectory[:1]:
    time = u.trajectory.time
    rgyr = u.atoms.radius_of_gyration()
    print(f"Frame: {ts.frame:3d}, Time: {time:4.0f} ps, Rgyr: {rgyr:.4f} A")

print(u.trajectory.frame)
print(u.trajectory[1].frame)
frame = u.trajectory.frame
time = u.trajectory.time
rgyr = u.atoms.radius_of_gyration()
print("Frame: {:3d}, Time: {:4.0f} ps, Rgyr: {:.4f} A".format(frame, time, rgyr))

rgyr = []
time = []
protein = u.select_atoms("protein")
for ts in u.trajectory:
    time.append(u.trajectory.time)
    rgyr.append(protein.radius_of_gyration())

rgyr_df = pd.DataFrame(rgyr, columns=['Radius of gyration (A)'], index=time)
rgyr_df.index.name = 'Time (ps)'

rgyr_df.head()

rgyr_df.plot(title='Radius of gyration')
plt.savefig("1Eru_RoG.jpg")
plt.show()

# Writing out coordinates
# RMSD

ca = u.select_atoms('protein')
with mda.Writer('4pti.xtc', ca.n_atoms) as w:
    for ts in u.trajectory:
        w.write(ca)

bb = u.select_atoms('protein')

u.trajectory[0] # first frame
first = bb.positions
u.trajectory[-1] #last frame
last = bb.positions
rms.rmsd(first, last)
u.trajectory[0] # set to first frame

rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['protein'])
rmsd_analysis.run()
print(rmsd_analysis.results.rmsd.shape)

rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                       columns=['Backbone', 'Protein'],
                       index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df.index.name = 'Time (ps)'
rmsd_df.head()
rmsd_df.plot(title='RMSD')
plt.savefig("1Eru_RMSD.jpg")
plt.show()
plt.savefig("1Eru_RMSD.jpg")
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

hbonds.run(start=0, verbose=True)

print(hbonds.results.hbonds.shape)
hbonds.results.hbonds.dtype
n_hbond = hbonds.results.hbonds

print("end of hbonds calculation")

n_hbonds = hbonds.count_by_time()

np.save("1Eru_hbonds.npy", n_hbond)

# Counts the number of hydrogen bonds for each frame
plt.plot(hbonds.times, n_hbonds, lw=2)

plt.title("Number of hydrogon bonds over time", weight="bold")
plt.xlabel("Time (ps)")
plt.ylabel(r"$N_{HB}$")
plt.grid(True)
plt.savefig("1eru_hbonds.jpg")
plt.show()