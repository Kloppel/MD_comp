import MDAnalysis as mda 
import pandas as pd
from matplotlib import pyplot as plt
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import numpy as np
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.dihedrals import Dihedral

used_program = "namd/"

universe = mda.Universe(f"{used_program}1eru.psf", f"{used_program}1eru_eq1.dcd", f"{used_program}1eru_prod1.dcd")

# Loading a structure or trajectory


u = mda.Universe(f"{used_program}1eru.psf", f"{used_program}1eru_eq1.dcd", f"{used_program}1eru_prod1.dcd")
print(u)
print(len(u.trajectory))

# Working with trajectories
print(len(u.trajectory))

# Output Values
output_folder = f"{used_program}output/"
data = []
with open(f"{used_program}1eru_sim.out", "r") as file:
    for line in file:
        #if "ENERGY:" in line:
        if line.startswith("ENERGY: "):
            values = line.split()[1:]
#            values = [float(i) for i in values]
            data.append(values)
columns = ["TS", "BOND", "ANGLE", "DIHED", "IMPRP", "ELECT", "VDW", "BOUNDARY", "MISC", "KINETIC", "TOTAL", "TEMP", "POTENTIAL", "TOTAL3", "TEMPAVG", "PRESSURE", "GPRESSURE", "VOLUME", "PRESSAVG", "GPRESSAVG"]
data, new_list = data[4:], []
for line in data:
    dataline = [float(i) for i in line]
    new_list.append(dataline)
data = new_list
new_list = new_list[10:]
df = pd.DataFrame(new_list, columns=columns)
df = df[4:]

plt.plot(df['TS'], df['TEMP'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Temperature (K)")
plt.title("Temperature over Time") 
plt.savefig(f"{output_folder}temperature.jpg")
plt.clf()
plt.plot(df['TS'], df['TOTAL'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Total Energy")
plt.title("Total Energy over Time") 
plt.savefig(f"{output_folder}total_energy.jpg")
plt.clf()
plt.plot(df['TS'], df['KINETIC'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Kinetic Energy")
plt.title("Kinetic Energy over Time") 
plt.savefig(f"{output_folder}kinetic_energy.jpg")
plt.clf()
plt.plot(df['TS'], df['POTENTIAL'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Potentielle Energy")
plt.title("Potential Energy over Time") 
plt.savefig(f"{output_folder}potential_energy.jpg")
plt.clf()
plt.plot(df['TS'], df['VOLUME'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Volumen")
plt.title("Volume over Time") 
plt.savefig(f"{output_folder}volume.jpg")
plt.clf()
plt.plot(df['TS'], df['PRESSURE'])
plt.xlabel("Time Step (2fs/step)")
plt.ylabel("Pressure (bar)")
plt.title("Pressure over Time") 
plt.savefig(f"{output_folder}pressure.jpg")
plt.clf()

# Radius of gyration

for ts in u.trajectory[:1]:
    time = ts.time
    rgyr = u.atoms.radius_of_gyration()
    print(f"Frame: {ts.frame:3d}, Time: {time:4.0f} ps, Rgyr: {rgyr:.4f} A")

print(u.trajectory.frame)
print(u.trajectory[1].frame)
frame = u.trajectory.frame
time = ts.time
rgyr = u.atoms.radius_of_gyration()
print("Frame: {:3d}, Time: {:4.0f} ps, Rgyr: {:.4f} A".format(frame, time, rgyr))

rgyr = []
time = []
protein = u.select_atoms("protein")
for ts in u.trajectory:
    time.append(ts.time)
    rgyr.append(protein.radius_of_gyration())

data = [{"Time (ps)": ts.time, "Radius of gyration (A)": protein.radius_of_gyration()} for ts in u.trajectory]
data = data[10:]
rgyr_df = pd.DataFrame(data).set_index("Time (ps)")
rgyr_df.index.name = 'Time (ps)'

rgyr_df.head()

rgyr_df.plot(title="Radius of Gyration", ylabel="Radius of gyration (A)", xlabel="Time (ps)", grid=True)
plt.savefig(f"{output_folder}1Eru_RoG.jpg", dpi=300)
plt.show()

# Writing out coordinates
# RMSD

ca = u.select_atoms('protein')
with mda.Writer('1eru.xtc', ca.n_atoms) as w:
    for ts in u.trajectory:
        w.write(ca)

bb = u.select_atoms('backbone')

u.trajectory[0] # first frame
first = ca.positions
u.trajectory[-1] #last frame
last = ca.positions
rms.rmsd(first, last)
u.trajectory[0] # set to first frame

rmsd_analysis = rms.RMSD(u, select='backbone', groupselections=['protein'])
rmsd_analysis.run()
print(rmsd_analysis.results.rmsd.shape)

rmsd_df = pd.DataFrame(rmsd_analysis.results.rmsd[:, 2:],
                       columns=['backbone', 'protein'],
                       index=rmsd_analysis.results.rmsd[:, 1])
rmsd_df.index.name = 'Time (ps)'
rmsd_df.head()
rmsd_df.plot(title='RMSD')
plt.savefig(f"{output_folder}1Eru_RMSD.jpg")
plt.show()

# RMSF

R = rms.RMSF(ca).run()

rmsf_df = pd.DataFrame({'Residue ID': ca.resids, 'RMSF (A)': R.results.rmsf})
rmsf_df.to_csv("1Eru_RMSF.csv", index=False)

plt.plot(ca.resids, R.results.rmsf, label="RMSF")
plt.title("RMSF of 1Eru")
plt.xlabel('Residue number')
plt.ylabel('RMSD ($\AA$)')
plt.legend()
plt.grid(True)
plt.savefig(f"{output_folder}1Eru_RMSF.jpg")
plt.show()

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

hbonds.run(start=12,stop=540, verbose=True)
print(hbonds.results.hbonds.shape)
n_hbonds = hbonds.results.hbonds
np.save(f"{output_folder}1Eru_hbonds.npy", n_hbonds)
print("end of hbonds calculation")

n_hbonds_time = hbonds.count_by_time()
np.save(f"{output_folder}1Eru_hbonds_count.npy", n_hbonds_time)

# Counts the number of hydrogen bonds for each frame
plt.plot(hbonds.times, n_hbonds_time)
plt.title("Number of hydrogon bonds over time", weight="bold")
plt.xlabel("Time (ps)")
plt.ylabel(r"$N_{HB}$")
plt.grid(True)
plt.savefig(f"{output_folder}1eru_hbonds.jpg", dpi=300)
plt.show()