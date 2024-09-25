from openmm import *
from openmm.app import *
from openmm.unit import *
from sys import stdout
from simtk.openmm.app import Modeller
import os

# Load PDB file
pdb = PDBFile('1eru.pdb')

# Load CHARMM force field and parameters
param_list=['toppar_water_ions.str', 'toppar_ions_won.str', 'toppar_dum_noble_gases.str', 'toppar_all36_synthetic_polymer_patch.str', 'toppar_all36_synthetic_polymer.str', 'toppar_all36_prot_retinol.str', 'toppar_all36_prot_na_combined.str', 'toppar_all36_prot_modify_res.str', 'toppar_all36_prot_model.str', 'toppar_all36_prot_heme.str', 'toppar_all36_prot_fluoro_alkanes.str', 'toppar_all36_prot_c36m_d_aminoacids.str', 'toppar_all36_prot_arg0.str', 'toppar_all36_polymer_solvent.str', 'toppar_all36_na_rna_modified.str', 'toppar_all36_nano_lig_patch.str', 'toppar_all36_nano_lig.str', 'toppar_all36_na_nad_ppi.str', 'toppar_all36_moreions.str', 'toppar_all36_lipid_yeast.str', 'toppar_all36_lipid_tag.str', 'toppar_all36_lipid_sphingo.str', 'toppar_all36_lipid_prot.str', 'toppar_all36_lipid_oxidized.str', 'toppar_all36_lipid_mycobacterial.str', 'toppar_all36_lipid_model.str', 'toppar_all36_lipid_miscellaneous.str', 'toppar_all36_lipid_lps.str', 'toppar_all36_lipid_lnp.str', 'toppar_all36_lipid_inositol.str', 'toppar_all36_lipid_hmmm.str', 'toppar_all36_lipid_ether.str', 'toppar_all36_lipid_detergent.str', 'toppar_all36_lipid_dag.str', 'toppar_all36_lipid_cholesterol.str', 'toppar_all36_lipid_cardiolipin.str', 'toppar_all36_lipid_bacterial.str', 'toppar_all36_lipid_archaeal.str', 'toppar_all36_label_spin.str', 'toppar_all36_label_fluorophore.str', 'toppar_all36_carb_imlab.str', 'toppar_all36_carb_glycopeptide.str', 'toppar_all36_carb_glycolipid.str', 'top_interface.rtf', 'top_all36_prot.rtf', 'top_all36_na.rtf', 'top_all36_lipid.rtf', 'top_all36_cgenff.rtf', 'top_all36_carb.rtf', 'par_interface.prm', 'par_all36_na.prm', 'par_all36m_prot.prm', 'par_all36_lipid.prm', 'par_all36_cgenff.prm', 'par_all36_carb.prm', 'cam.str']
param_files= ["../params/"+ filename for filename in param_list]
params = CharmmParameterSet(*param_files)

print("PDB and Forcefield loaded")

# Finds possible sulfidbonds between cys with a maximum distance of 2.4 Angström
def find_disulfide_bonds(pdb_structure):
    disulfide_bonds = []
    atom_pairs = []

    for chain in pdb_structure.topology.chains():
        cysteines = [res for res in chain.residues() if res.name == 'CYS']
        for i, res1 in enumerate(cysteines):
            for res2 in cysteines[i+1:]:
                # Find the SG atoms of the cysteines
                sg1 = next(atom for atom in res1.atoms() if atom.name == 'SG')
                sg2 = next(atom for atom in res2.atoms() if atom.name == 'SG')

                # Calculate the distance between the two sulfur atoms (SG)
                distance = norm(pdb_structure.positions[sg1.index] - pdb_structure.positions[sg2.index])
                if distance < 0.24*nanometers:  # 2.4 Å is typical for disulfide bonds
                    disulfide_bonds.append((res1, res2))
                    atom_pairs.append((sg1.index, sg2.index))
    
    return disulfide_bonds, atom_pairs

print("Possible sulfidbonds detected")

# Detect disulfide bonds in the PDB structure
disulfides, atom_pairs = find_disulfide_bonds(pdb)
print(f"Found {len(disulfides)} disulfide bonds.")

print("disulfidbons detected in pdb")

# Apply disulfide bond patches
modeller = Modeller(pdb.topology, pdb.positions)

disulfides = [(res1, res2)]  # Add your cysteine residue pairs

# Loop through the list of disulfide pairs and add bonds
for sg1, sg2 in disulfides:
    modeller.addBond(res1, res2)
    print(f"Disulfide bond added between residue {res1} and residue {res2}")

# Create the system using the patched force field
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=NoCutoff, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

print("Creating System")
print("applying disulfidbond patches")

# Create the system using the patched force field
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=NoCutoff, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

print("Creating System")

# Write out PSF and PDB files
psf = CharmmPsfFile()
psf.writeFile(pdb.topology, pdb.positions, '1eru_sul.psf')
PDBFile.writeFile(pdb.topology, pdb.positions, open('1eru_sul.pdb', 'w'))

print("writing pdb and psf file")

# Set up the integrator and simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

print("setting up the integrator and the simulation")

# Minimize energy
simulation.minimizeEnergy()

print("minimization")

# Add reporters
simulation.reporters.append(PDBReporter('1eru_sul_output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

print("wrinting output")

# Run the simulation for 10,000 steps
simulation.step(10000)

print("simulating")