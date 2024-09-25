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
param_list = ['charmm36.xml'
]
param_files = ["../charmm36/" + filename for filename in param_list]
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


# Create the system using the patched force field
system = psf.createSystem(params,modeller.topology, 
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

