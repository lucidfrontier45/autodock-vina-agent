# mk_prepare_ligand.py options (Meeko 0.8.x)

## Input / output

| Option | Meaning |
|---|---|
| `-i, --mol FILE` | **Required.** Input molecule file (SDF preferred, MOL2 supported). |
| `-o, --out FILE` | Output PDBQT filename. Single-molecule input only. |
| `-,` | Do not write a file; write PDBQT to STDOUT. Ignores `-o`. Single molecule only. |
| `--multimol_outdir DIR` | Directory for PDBQT files from multi-molecule input. Incompatible with `-o` and `-`. |
| `--multimol_prefix PREFIX` | Replace internal molecule name with prefix (multi-mol only). |
| `-z, --multimol_targz` | Compress the PDBQT outputs into `.tar.gz` archives. |
| `--multimol_targz_size N` | Number of PDBQT files per archive. Default 10000. |
| `--name_from_prop PROP` | Set the molecule name from an RDKit/SDF property. |
| `-v, --verbose` | Print information about the molecule setup. |

## Molecule preparation

| Option | Meaning |
|---|---|
| `-c, --config_file FILE` | JSON config for `MoleculePreparation`. CLI args override it. |
| `--merge_these_atom_types [T ...]` | Atom types to merge into their parent (polar H merged by default). Pass no value to merge nothing. |
| `-r, --rigidify_bonds_smarts SMARTS` | SMARTS pattern(s) whose bonds become non-rotatable. |
| `-b, --rigidify_bonds_indices i j [...]` | 1-based indices of the two SMARTS atoms forming the bond. |
| `-a, --flexible_amides` | Let amide bonds rotate (non-planar). Not recommended. |
| `-p, --load_atom_params NAME [NAME ...]` | Atom typing SMARTS tables. Packaged: `ad4_types`, `ad4_hb`, `ad4_vdw`, `ad4_desolv_volume`, `ad4_desolv_param`, `vina_params`. Local files need the `.json` suffix. |
| `-aa, --add_atom_types JSON` | JSON list of `{"smarts": ..., "atype": ..., "IDX": [...]}` added on top of the config's types. |
| `--charge_model MODEL` | `gasteiger` (default), `nagl`, `espaloma`, `zero`, `read`. |
| `--charge_atom_prop PROP` | Read partial charges from an input atom property. Default `PartialCharge` (SDF), `_TriposPartialCharge` (MOL2). |
| `--bad_charge_ok` | Allow NaN/Inf charges in the PDBQT. |
| `--double_bond_penalty N` | Penalty > 100 prevents breaking double bonds (default 50). |
| `--rigid_macrocycles` | Keep macrocycles rigid in the input conformation. |
| `--macrocycle_allow_A` | Allow bond break with atom type A; the atom is retyped as C. |
| `--min_ring_size N` | Minimum ring size considered for opening (default 7). |
| `--keep_chorded_rings` | Keep all rings from exhaustive ring perception. |
| `--keep_equivalent_rings` | Keep rings with equivalent size and neighbours. |
| `-w, --hydrate` | Add water molecules for hydrated docking. |
| `--add_index_map` | Write a `REMARK SMILES IDX` map from input atoms to PDBQT atoms. |
| `--remove_smiles` | Do not write the SMILES remark into the PDBQT. |
| `--rename_atoms` | Rename atoms to `<original name><1-based index>`. |

## Reactive docking

| Option | Meaning |
|---|---|
| `--reactive_smarts SMARTS` | SMARTS defining the reactive group. |
| `--reactive_smarts_idx N` | 1-based index of the reactive atom inside the SMARTS. |

## Covalent docking (tethered) - needs ProDy

| Option | Meaning |
|---|---|
| `--receptor FILE` | Receptor file (`.pdb`, `.mmcif`), requires ProDy. |
| `--rec_residue RES` | Attachment residue, e.g. `"A:LYS:204"`, `"A:HIS:"`, `":LYS:"`. |
| `--tether_smarts SMARTS` | SMARTS defining the ligand atoms used for attachment. |
| `--tether_smarts_indices I J` | 1-based indices of the two attached SMARTS atoms. Default `1 2`. |

## Behaviour notes

- Command-line arguments override `-c/--config_file`; `--add_atom_types` **extends** it.
- Merging an atom type adds its partial charge to the parent and marks it `is_ignore`,
  so it is not written to the PDBQT. Changing atom typing changes what gets merged.
- SMARTS bond indices: 1-based on the CLI, 0-based in the Python API.
- Atom typing SMARTS entries are matched in file order; the last match wins.
