from openmm import *
from openmm.app import *
from openmm.unit import *
from simtk.openmm.app import Modeller
from simtk.openmm.app import CharmmPsfFile, CharmmParameterSet
import os
from numpy.linalg import norm  # For distance calculation

# Load PDB file
pdb = PDBFile('1eru.pdb')

# Load CHARMM force field and parameters
param_list = [
    'toppar_water_ions.str', 'toppar_ions_won.str', 'toppar_dum_noble_gases.str',
    'toppar_all36_synthetic_polymer_patch.str', 'toppar_all36_synthetic_polymer.str',
    'toppar_all36_prot_retinol.str', 'toppar_all36_prot_na_combined.str',
    'toppar_all36_prot_modify_res.str', 'toppar_all36_prot_model.str',
    'toppar_all36_prot_heme.str', 'toppar_all36_prot_fluoro_alkanes.str',
    'toppar_all36_prot_c36m_d_aminoacids.str', 'toppar_all36_prot_arg0.str',
    'toppar_all36_polymer_solvent.str', 'toppar_all36_na_rna_modified.str',
    'toppar_all36_nano_lig_patch.str', 'toppar_all36_nano_lig.str',
    'toppar_all36_na_nad_ppi.str', 'toppar_all36_moreions.str', 'toppar_all36_lipid_yeast.str',
    'toppar_all36_lipid_tag.str', 'toppar_all36_lipid_sphingo.str',
    'toppar_all36_lipid_prot.str', 'toppar_all36_lipid_oxidized.str',
    'toppar_all36_lipid_mycobacterial.str', 'toppar_all36_lipid_model.str',
    'toppar_all36_lipid_miscellaneous.str', 'toppar_all36_lipid_lps.str',
    'toppar_all36_lipid_lnp.str', 'toppar_all36_lipid_inositol.str',
    'toppar_all36_lipid_hmmm.str', 'toppar_all36_lipid_ether.str',
    'toppar_all36_lipid_detergent.str', 'toppar_all36_lipid_dag.str',
    'toppar_all36_lipid_cholesterol.str', 'toppar_all36_lipid_cardiolipin.str',
    'toppar_all36_lipid_bacterial.str', 'toppar_all36_lipid_archaeal.str',
    'toppar_all36_label_spin.str', 'toppar_all36_label_fluorophore.str',
    'toppar_all36_carb_imlab.str', 'toppar_all36_carb_glycopeptide.str',
    'toppar_all36_carb_glycolipid.str', 'top_interface.rtf', 'top_all36_prot.rtf',
    'top_all36_na.rtf', 'top_all36_lipid.rtf', 'top_all36_cgenff.rtf',
    'top_all36_carb.rtf', 'par_interface.prm', 'par_all36_na.prm',
    'par_all36m_prot.prm', 'par_all36_lipid.prm', 'par_all36_cgenff.prm',
    'par_all36_carb.prm', 'cam.str'
]
param_files = ["../params/" + filename for filename in param_list]
params = CharmmParameterSet(*param_files)

print("PDB and Forcefield loaded")

# Finds possible disulfide bonds between CYS residues with a max distance of 2.4 Å
def find_disulfide_bonds(pdb_structure):
    disulfide_bonds = []
    for chain in pdb_structure.topology.chains():
        cysteines = [res for res in chain.residues() if res.name == 'CYS']
        for i, res1 in enumerate(cysteines):
            for res2 in cysteines[i+1:]:
                # Find the SG atoms of the cysteines
                sg1 = next(atom for atom in res1.atoms() if atom.name == 'SG')
                sg2 = next(atom for atom in res2.atoms() if atom.name == 'SG')

                # Calculate the distance between the two sulfur atoms (SG)
                distance = norm(pdb_structure.positions[sg1.index] - pdb_structure.positions[sg2.index])
                if distance < 0.24 * nanometers:  # 2.4 Å typical for disulfide bonds
                    disulfide_bonds.append((res1, res2))
    
    return disulfide_bonds

# Detect disulfide bonds in the PDB structure
disulfides = find_disulfide_bonds(pdb)
print(f"Found {len(disulfides)} disulfide bonds.")

# Apply disulfide bond patches
modeller = Modeller(pdb.topology, pdb.positions)

# Add bonds for detected disulfide bonds
for res1, res2 in disulfides:
    sg1 = next(atom for atom in res1.atoms() if atom.name == 'SG')
    sg2 = next(atom for atom in res2.atoms() if atom.name == 'SG')
    modeller.topology.addBond(sg1, sg2)
    print(f"Disulfide bond added between residue {res1.index} and residue {res2.index}")

forcefield = ForceField('params')
# Create the system using the patched force field
system = forcefield.createSystem(params,modeller.topology, 
                             nonbondedMethod=NoCutoff, 
                             nonbondedCutoff=1*nanometer, 
                             constraints=HBonds)

print("Creating System")

# Write out PSF and PDB files
psf = CharmmPsfFile()  # You need to initialize with a proper structure
psf.writeFile(modeller.topology, modeller.positions, '1eru_sul.psf')
PDBFile.writeFile(modeller.topology, modeller.positions, open('1eru_sul.pdb', 'w'))

print("Writing PDB and PSF files")

# Set up the integrator and simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print("Setting up the integrator and the simulation")

# Run the simulation for 10,000 steps
simulation.step(10000)
print("Simulating")

