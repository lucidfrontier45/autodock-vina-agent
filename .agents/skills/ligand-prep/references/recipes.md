# Recipes

Concrete, copy-pasteable workflows. All paths assume the project root.

## 1. Prepare a single ligand

```bash
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --smiles "CC(=O)Nc1ccc(O)cc1" --name paracetamol \
    --out ligands/paracetamol.sdf \
    --pdbqt ligands/paracetamol.pdbqt

uv run python .agents/skills/ligand-prep/scripts/check_pdbqt.py ligands/paracetamol.pdbqt
```

## 2. Prepare a library from a multi-molecule SD file

```bash
# normalize everything to 3D + explicit H first
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --sdf library.sdf --out library_3d.sdf

# then one PDBQT per molecule
uv run mk_prepare_ligand.py -i library_3d.sdf --multimol_outdir pdbqt/
```

or in one step:

```bash
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --sdf library.sdf --out library_3d.sdf --pdbqt-dir pdbqt/
```

Large libraries: add `-z` to `mk_prepare_ligand.py` to get `.tar.gz` archives
(`--multimol_targz_size` controls how many PDBQT files per archive, default 10000).

## 3. Keep polar hydrogens (don't merge them)

```bash
uv run mk_prepare_ligand.py -i mol.sdf -o mol.pdbqt --merge_these_atom_types
```

## 4. Reusable configuration file

```json
{
  "merge_these_atom_types": ["H"],
  "charge_model": "gasteiger",
  "rigid_macrocycles": false,
  "min_ring_size": 7
}
```

```bash
uv run mk_prepare_ligand.py -i mol.sdf -c prep.json -o mol.pdbqt
```

CLI flags override the file; `--add_atom_types` extends it.

## 5. Custom atom type (SMARTS)

Type aromatic C-H hydrogens as `HX` instead of `H`:

```bash
uv run mk_prepare_ligand.py -i mol.sdf -o mol.pdbqt \
    --add_atom_types '[{"smarts": "[H]c", "atype": "HX"}]'
```

Python equivalent (note the 0-based `IDX` and the swapped atom order):

```python
mk_prep = MoleculePreparation(add_atom_types=[{"smarts": "c[H]", "atype": "HX", "IDX": [1]}])
```

## 6. Rigidify a conjugated bond / an amide

```bash
uv run mk_prepare_ligand.py -i mol.sdf -o mol.pdbqt \
    --rigidify_bonds_smarts "C=CC=C" --rigidify_bonds_indices 1 2 \
    --rigidify_bonds_smarts "[CX3](=O)[NX3]" --rigidify_bonds_indices 1 3
```

Python (0-based pairs):

```python
mk_prep = MoleculePreparation(
    rigidify_bonds_smarts=["C=CC=C", "[CX3](=O)[NX3]"],
    rigidify_bonds_indices=[(0, 1), (0, 2)],
)
```

## 7. Vina-style parameters / OpenFF

```python
mk_prep = MoleculePreparation(load_atom_params=["vina_params", "openff-2.3.0"])
```

`openff-forcefields` must be installed. Local JSON files need the `.json` suffix.

## 8. Macrocycles

Default: macrocycles (7-33 atoms) are broken open and become flexible branches.

```bash
uv run mk_prepare_ligand.py -i macro.sdf -o macro.pdbqt --rigid_macrocycles
uv run mk_prepare_ligand.py -i macro.sdf -o macro.pdbqt --min_ring_size 10 --max_ring_size 30
```

## 9. Reactive docking

```bash
uv run mk_prepare_ligand.py -i warhead.sdf -o warhead.pdbqt \
    --reactive_smarts "[CX3]=[SX2]" --reactive_smarts_idx 1
```

The reactive atom is retyped; Meeko returns one setup per reactive atom.

## 10. Tethered (covalent) docking - needs ProDy

```bash
uv run mk_prepare_ligand.py -i mol.sdf -o mol.pdbqt \
    --receptor receptor.pdb --rec_residue "A:LYS:204" \
    --tether_smarts "[CX4][NX3]" --tether_smarts_indices 1 2
```

## 11. Batch preparation + validation loop in Python

```python
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

mk_prep = MoleculePreparation()
outdir = Path("pdbqt"); outdir.mkdir(exist_ok=True)

for i, smi in enumerate(smiles_list):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    setups = mk_prep(mol)
    pdbqt, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        print(f"FAIL {smi}: {err}")
        continue
    (outdir / f"lig_{i:04d}.pdbqt").write_text(pdbqt)
```

## 12. Prototyping with Vina without touching the filesystem

```python
from vina import Vina
v = Vina(sf_name="vina")
v.set_receptor("rec.pdbqt")
v.set_ligand_from_string(pdbqt_string)
v.compute_vina_maps(center=[x, y, z], box_size=[20, 20, 20])
v.optimize()
print(v.score())
```
