"""Regularize 1LYZ geometry: rebuild missing atoms, add H, minimize."""
from pdbfixer import PDBFixer
from openmm.app import PDBFile, ForceField, NoCutoff
from openmm import Platform, LangevinMiddleIntegrator
from openmm.app import Simulation, NoCutoff
import openmm.unit as unit

fixer = PDBFixer(filename='/tmp/rec-eval/eval-3/without_skill/1LYZ_protein.pdb')
fixer.findMissingResidues()
fixer.missingResidues = {}          # don't build terminal extensions
fixer.findMissingAtoms()
print('missing heavy atoms:', {r: [a.name for a in v] for r, v in fixer.missingAtoms.items()} or 'none')
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

ff = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
system = ff.createSystem(fixer.topology, nonbondedMethod=NoCutoff)
integrator = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
sim = Simulation(fixer.topology, system, integrator, Platform.getPlatformByName('CPU'))
sim.context.setPositions(fixer.positions)
print('initial energy:', sim.context.getState(getEnergy=True).getPotentialEnergy())
sim.minimizeEnergy(maxIterations=2000, tolerance=10*unit.kilojoule_per_mole/unit.nanometer)
print('final energy:  ', sim.context.getState(getEnergy=True).getPotentialEnergy())

with open('/tmp/rec-eval/eval-3/without_skill/1LYZ_fixed_H.pdb', 'w') as f:
    PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), f, keepIds=True)

# heavy-atom-only version for meeko (it adds its own hydrogens from templates)
with open('/tmp/rec-eval/eval-3/without_skill/1LYZ_fixed.pdb', 'w') as out:
    for line in open('/tmp/rec-eval/eval-3/without_skill/1LYZ_fixed_H.pdb'):
        if line.startswith(('ATOM', 'HETATM')) and line[76:78].strip() == 'H':
            continue
        out.write(line)
print('done')
