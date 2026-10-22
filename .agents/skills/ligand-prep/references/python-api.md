# Meeko Python API (ligand preparation)

## Pipeline

```
RDKit mol (H-explicit, 3D)
  -> MoleculePreparation.prepare(mol)  -> list[MoleculeSetup]
  -> PDBQTWriterLegacy.write_string(setup) -> (pdbqt_string, is_ok, error)
```

## MoleculePreparation

```python
MoleculePreparation(
    merge_these_atom_types=("H",),   # atom types folded into their parent
    merge_these_atom_params=(),
    merge_rmin_half=False,
    hydrate=False,
    flexible_amides=False,           # amides rigid by default
    rigid_macrocycles=False,
    untyped_macrocycles=False,
    min_ring_size=7,
    max_ring_size=33,
    keep_chorded_rings=False,
    keep_equivalent_rings=False,
    double_bond_penalty=50,
    macrocycle_allow_A=False,
    rigidify_bonds_smarts=[],
    rigidify_bonds_indices=[],
    input_atom_params=None,
    load_atom_params="ad4_types",
    add_atom_types=None,
    input_offatom_params=None,
    load_offatom_params=None,
    charge_model="gasteiger",         # gasteiger | nagl | espaloma | zero | read
    charge_atom_prop=None,
    dihedral_model=None,
    reactive_smarts=None,
    reactive_smarts_idx=None,
    add_index_map=False,
    remove_smiles=False,
    compute_charges=False,
    crippen=False,
    crippen_as_solpar=False,
    override_ad4sol_par_including_q=False,
    override_ad4sol_par_including_q_qasp=0.0,
)
```

Other class methods:
- `MoleculePreparation.from_config(dict)` - build from a config dict.
- `MoleculePreparation.from_json_file(path)`
- `MoleculePreparation.get_defaults_dict()`
- `MoleculePreparation.packaged_params` - names of the packaged parameter files.
- `mk_prep(mol)` (i.e. `__call__`) - returns **a list** of `MoleculeSetup`.
  More than one entry only with reactive docking.

## MoleculeSetup (the parameterized molecule)

Useful attributes / methods:

| Member | Meaning |
|---|---|
| `molsetup.atoms` | list of `Atom` objects (`.atom_type`, `.charge`, `.is_ignore`, `.element`, `.position`) |
| `molsetup.atom_params` | dict: parameter name -> list of per-atom values (`atype`, `vina_ri`, `rmin_half`, `epsilon`, ...) |
| `molsetup.rings` | perceived rings |
| `molsetup.flexibility_model` | rotatable-bond model |
| `molsetup.rotamers` | dihedral rotamers |
| `molsetup.name` | molecule name |
| `molsetup.get_smiles_and_order()` | SMILES + atom order |
| `molsetup.get_coord(atom_idx)` | coordinates |
| `molsetup.get_charge(idx)` | partial charge |
| `molsetup.get_is_ignore(idx)` | whether the atom is merged away |
| `molsetup.to_json()` / `from_json()` | serialise the setup |

## PDBQTWriterLegacy

```python
pdbqt_string, is_ok, error_message = PDBQTWriterLegacy.write_string(
    molsetup, add_index_map=False, remove_smiles=False, bad_charge_ok=False
)
```

Other methods: `write_from_polymer`, `write_string_from_polymer`,
`adapt_pdbqt_for_autodock4_flexres`, `remark_index_map`.

## Atom typing

Default table: `meeko/data/params/ad4_types.json` (AutoDock4 types).
Other packaged files: `ad4_hb.json`, `ad4_vdw.json`, `ad4_desolv_volume.json`,
`ad4_desolv_param.json`, `vina_params.json`.

```python
mk_prep = MoleculePreparation(load_atom_params=["vina_params", "openff-2.3.0"])
```

Local files need the `.json` suffix. Entries are matched in order - the last matching
SMARTS wins. Each entry: `{"smarts": ..., "atype": ..., "IDX": [...]}` (`IDX` defaults
to `[0]`, 0-based, and selects which atoms of the SMARTS get the type).

OpenFF vdW parameters require the `openff-forcefields` package.

## Charges

- `gasteiger` (default) - fast, from RDKit.
- `nagl` - NAGL model (needs the `nagl` package).
- `espaloma` - ESPALOMA charges (needs `espaloma_charge`).
- `zero` - all charges zero.
- `read` - take charges from the input file property (`PartialCharge` for SDF,
  `_TriposPartialCharge` for MOL2; override with `charge_atom_prop`).

## Hydrated docking

```python
from meeko import Hydrate
```

CLI: `-w/--hydrate`.

## Other useful imports

```python
from meeko import (
    MoleculePreparation, MoleculeSetup, PDBQTWriterLegacy,
    PDBQTMolecule, PDBQTReceptor, RDKitMolCreate, RDKitMoleculeSetup,
    Polymer, Monomer, ResidueTemplate, ResidueChemTemplates, ResiduePadder,
    Hydrate, AtomTyper, UniqAtomParams,
)
```
