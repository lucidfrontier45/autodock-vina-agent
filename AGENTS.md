# AGENTS.md

## Python env

- use uv as manager

## Skills

- `.agents/skills/ligand-prep` - ligand preparation for AutoDock/Vina with Meeko
  (SMILES/SDF -> 3D, H-explicit SDF -> PDBQT; atom typing, rotatable bonds, charges,
  macrocycles, reactive and tethered docking). Read its `SKILL.md` before preparing
  ligands.
- `.agents/skills/receptor-prep` - receptor preparation for AutoDock/Vina with Meeko
  (PDB/mmCIF/PQR -> PDBQT; residue templates, protonation states, altlocs, waters/HETATM
  cleanup, flexible and reactive (covalent) sidechains, grid box/GPF). Read its `SKILL.md`
  before preparing receptors.
