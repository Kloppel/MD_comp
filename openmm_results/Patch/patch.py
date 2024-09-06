from openmm import *
from openmm.app import *
from openmm.unit import *
from sys import stdout
import os

# Load PDB file
pdb = PDBFile('1eru.pdb')

# Load CHARMM force field and parameters
forcefield = ForceField('../charmm.xml/charmm36.xml')

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

# Detect disulfide bonds in the PDB structure
disulfides, atom_pairs = find_disulfide_bonds(pdb)
print(f"Found {len(disulfides)} disulfide bonds.")

# Apply disulfide bond patches
for res1, res2 in disulfides:
    forcefield.applyPatch('disulfide', res1, res2)

# Create the system using the patched force field
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=NoCutoff, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

# Write out PSF and PDB files
psf = CharmmPsfFile()
psf.writeFile(pdb.topology, pdb.positions, '1eru_sul.psf')
PDBFile.writeFile(pdb.topology, pdb.positions, open('1eru_sul.pdb', 'w'))

# Set up the integrator and simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Minimize energy
simulation.minimizeEnergy()

# Add reporters
simulation.reporters.append(PDBReporter('1eru_sul_output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# Run the simulation for 10,000 steps
simulation.step(10000)
