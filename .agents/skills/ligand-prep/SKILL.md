---
name: ligand-prep
description: Prepare small-molecule ligands for AutoDock/Vina docking with Meeko - convert SMILES/SDF/MOL2 into PDBQT, control atom typing, hydrogen merging, rotatable bonds, charges, macrocycles, reactive and tethered (covalent) docking. Use when the task mentions ligand preparation, PDBQT generation, mk_prepare_ligand.py, MoleculePreparation, atom typing, rotatable bonds, charges, or preparing a ligand library for docking.
---

# Ligand Preparation (Meeko)

Produce AutoDock-ready **PDBQT** files from ligands. Meeko is the parameterizer:
it does **not** add hydrogens and does **not** generate 3D coordinates - RDKit must do
that first (use `scripts/ligprep.py` in this skill).

## Hard requirements for the input molecule

The molecule handed to Meeko MUST have:

1. **All hydrogens as real atoms** (no implicit H).
2. **3D coordinates** on every atom.
3. **Valid valences / formal charges** (RDKit sanitization).

Prefer **SDF** over MOL2. Format: SD files (.sdf) with multiple molecules are supported.

## Standard workflow

```bash
# 1. build a 3D, H-explicit SDF from SMILES (this skill's helper)
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --smiles "CC(=O)Nc1ccc(O)cc1" --name paracetamol --out ligands/paracetamol.sdf

# 2. parameterize to PDBQT
uv run mk_prepare_ligand.py -i ligands/paracetamol.sdf -o ligands/paracetamol.pdbqt
```

Batch, one PDBQT per molecule in an SD file:

```bash
uv run mk_prepare_ligand.py -i ligands/library.sdf --multimol_outdir ligands/pdbqt
```

Print to STDOUT instead of a file (single molecule only): use `-` instead of `-o`.

## Python API (preferred when a script is needed)

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
AllChem.EmbedMolecule(mol, randomSeed=42)

mk_prep = MoleculePreparation()          # defaults: merge H, gasteiger charges, AD4 types
molsetup = mk_prep(mol)[0]                # returns a LIST (multi-entry only for reactive docking)
pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(molsetup)
```

Notes:
- `mk_prep(mol)` returns a list; take `[0]` unless you configured reactive docking.
- `PDBQTWriterLegacy.write_string` returns `(pdbqt_string, is_ok, error_message)` - always check `is_ok`.
- The PDBQT string can be handed straight to Vina's Python API - no file needed.
- Inspect results: `[a.atom_type for a in molsetup.atoms]`, `molsetup.atom_params`
  (dict of per-atom parameter lists, e.g. `vina_ri`, `rmin_half`), `molsetup.rings`,
  `molsetup.flexibility_model`.

## Configuration

Everything in `MoleculePreparation` is settable from a JSON config file:

```json
{"merge_these_atom_types": ["H"], "charge_model": "gasteiger"}
```

```bash
uv run mk_prepare_ligand.py -i mol.sdf -c config.json
```

Command-line args **override** the config file, except `--add_atom_types`, which
**extends** it. In Python: `MoleculePreparation.from_config(config_dict)`.

## The knobs that matter most

| Goal | CLI | Python |
|---|---|---|
| Keep polar hydrogens | `--merge_these_atom_types` (empty) | `merge_these_atom_types=()` |
| Merge a different type | `--merge_these_atom_types HD` | `merge_these_atom_types=("HD",)` |
| Custom atom type | `--add_atom_types '[{"smarts":"[H]c","atype":"HX"}]'` | `add_atom_types=[{"smarts":"c[H]","atype":"HX","IDX":[1]}]` |
| Full atom-type table | `-p ad4_types.json `/ `vina_params` | `load_atom_params=["vina_params"]` |
| Rigidify a bond | `-r "C=CC=C" -b 1 2` | `rigidify_bonds_smarts=["C=CC=C"], rigidify_bonds_indices=[(0,1)]` |
| Non-rotatable amide | default (leave off `-a`) | `flexible_amides=False` |
| Zero charges | `--charge_model zero` | `charge_model="zero"` |
| OpenFF charges | `--charge_model espaloma` (+ `openff-forcefields`) | `charge_model="espaloma"` |
| Rigid macrocycle | `--rigid_macrocycles` | `rigid_macrocycles=True` |
| Keep polar H in output | default merged; see merging above | - |

**Index base differs**: SMARTS bond indices are **1-based on the CLI** and **0-based in
the Python API**. `--add_atom_types` `"IDX"` defaults to `[0]` and is 0-based in both.

## Macrocycles / rings

- Default: macrocycles (rings of `min_ring_size`..`max_ring_size`, 7-33 atoms) are
  opened into flexible branches.
- `--rigid_macrocycles` keeps them in the input conformation.
- `--min_ring_size N` / `--max_ring_size N` control what counts as a macrocycle.
- `--keep_chorded_rings`, `--keep_equivalent_rings` change ring perception.

## Reactive docking (covalent bond formation)

```bash
uv run mk_prepare_ligand.py -i mol.sdf --reactive_smarts "[C]=[C]" --reactive_smarts_idx 1
```

The reactive atom gets a special type; `mk_prep(mol)` then returns **multiple**
molecule setups (one per reactive atom). Index is 1-based.

## Tethered / covalent docking (attach to a receptor residue)

Requires **ProDy**.

```bash
uv run mk_prepare_ligand.py -i mol.sdf \
    --receptor receptor.pdb --rec_residue "A:LYS:204" \
    --tether_smarts "[CX4][NX3]" --tether_smarts_indices 1 2
```

## After preparation - validate before docking

1. `TORSDOF` in the PDBQT matches the expected number of rotatable bonds.
2. Total charge is sensible (sum of partial charges).
3. No `nan`/`inf` charges (use `--bad_charge_ok` only if you know why).
4. Atom types are the expected AD4 types (`C`, `A`, `N`, `NA`, `OA`, `SA`, `HD` ...).
5. `--add_index_map` writes `REMARK SMILES IDX` so PDBQT atoms can be mapped back to input atoms.
6. The heavy-atom coordinates are the same as the input (Meeko does not move atoms).

## Common problems

- **"molecule has no 3D coordinates" / missing H** -> run RDKit `AddHs` + `EmbedMolecule`
  first (see `scripts/ligprep.py`).
- **Wrong atom count after merging** -> merged H are added to the parent and not written;
  check `molsetup.atom_params` / `is_ignore`.
- **Amide bond rotates and geometry is bad** -> do not pass `-a`; keep amides rigid.
- **Vina complains about charges** -> try `--charge_model espaloma` or `zero`.
- **Macrocycle opens unexpectedly** -> use `--rigid_macrocycles`.

## Reference files

- `references/cli-options.md` - every option of `mk_prepare_ligand.py`
- `references/python-api.md` - `MoleculePreparation` / `MoleculeSetup` API details
- `references/recipes.md` - copy-pasteable workflows (library prep, custom typing, macrocycles, reactive/tethered)
- `scripts/ligprep.py` - SMILES/SDF -> H-explicit 3D SDF -> PDBQT helper
- `scripts/check_pdbqt.py` - validate a PDBQT (TORSDOF, atom types, charge sanity)
- Upstream docs: <https://meeko.readthedocs.io/en/develop/lig_overview.html>,
  `lig_prep_basic.html`, `lig_prep_advanced.html`, `lig_cli_options.html`
