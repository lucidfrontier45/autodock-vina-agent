# mk_prepare_receptor.py options (Meeko 0.8.0)

Verified against the installed version's `--help`. The develop-branch docs rename
`-x/--delete_bad_res` to `-d/--delete_bad_res`; check `--help` if flags clash.

## Input

| Option | Meaning |
|---|---|
| `--read_pdb FILE` | Read PDB with RDKit's parser (no ProDy needed). |
| `--read_pqr FILE` | Read PQR (no ProDy); pair with `--charge_model read`. |
| `-i, --read_with_prody FILE` | Read PDB **or mmCIF** via ProDy (must be installed: `uv add prody`). |
| `--read_json FILE` | Read a Meeko-parameterized JSON receptor. Parameters and flexres state are lost. |

## Output (write flags)

Flags starting with `--write`/`-p` accept either an explicit filename or, with no
argument, derive it from `-o/--output_basename`:

```bash
mk_prepare_receptor.py -i sys.pdb -o myenzyme -p -j            # myenzyme.pdbqt + myenzyme.json
mk_prepare_receptor.py -i sys.pdb --write_pdbqt myenzyme.pdbqt --write_json myenzyme.json
```

| Option | Meaning |
|---|---|
| `-o, --output_basename BASE` | Default basename for write flags without filename. |
| `-p, --write_pdbqt [FILE]` | PDBQT output. `_rigid`/`_flex` suffixes appear when flexres exist. |
| `-j, --write_json [FILE]` | Parameterized receptor as JSON (round-trip via `--read_json`). |
| `--write_pdb FILE` | Prepared receptor as PDB (filename required). |
| `-g, --write_gpf [FILE]` | AutoGrid input file. **Requires `-p`** (GPF points at the PDBQT) and a box. |
| `-v, --write_vina_box [FILE]` | Vina box config (`.box.txt`) + `.box.pdb` for visualization. Requires a box. |
| `--debug_fn FILE` | Log debug level to file. |
| `--ignore_https_cert` | Skip cert check when fetching CCD definitions from PDB (last resort). |

## Receptor perception

| Option | Meaning |
|---|---|
| `-n, --set_template SEL=TMPL,...` | Force templates, e.g. `A:5,7=CYX,B:17=HID`. |
| `-d, --delete_residues SEL` | Delete residues, e.g. `A:350,B:15,16,17`. |
| `-b, --blunt_ends POS` | Blunt ends: `A:123,200=2,A:1=0` (atom index after `=`). |
| `--add_templates FILE` | Extra residue templates: JSON, or `resname:file.sdf`. Repeatable. |
| `--cache_templates [FILE]` | Cache templates to `$HOME/.meeko_residue_chem_templates_cached.json` (cumulative updates). Default OFF. |
| `--config_preset NAME` | Packaged config (only choice: `scofu1`). |
| `--config_file FILE` | Local JSON config; overrides preset option-wise, CLI args override both. |
| `-x, --delete_bad_res` | Delete (instead of erroring on) residues that fail **template matching**. Does NOT rescue RDKit mol-build failures. |
| `--delete_bad_res_from_box_radius R` | Delete unmatched residues >R A outside the box; near/inside still errors. Needs a box. Incompatible with `-x`. |
| `--default_altloc LOC` | Default alternate location for all residues. |
| `--wanted_altloc SEL=LOC` | Per-residue altloc, e.g. `:5=B,B:17=A`; overrides `--default_altloc`. |
| `--forgive_extra_bonds` | Tolerate excess inter-residue bonds (clashed structures) - may match wrong templates. |
| `--adj_padding` | Pad each residue with adjacent residues instead of template reactions. **Requires `--compute_charges`.** |

## Charges

| Option | Meaning |
|---|---|
| `--compute_charges` | Compute from scratch instead of reading pre-calculated template charges (slower). |
| `--charge_model M` | `gasteiger` (default), `nagl`, `espaloma`, `zero`, `read` (needs `--read_pqr`). nagl/espaloma are slow. |

## Grid box (exactly one combo required for `-v`/`-g`)

| Option | Meaning |
|---|---|
| `--box_size X Y Z` | Box size in Angstrom. |
| `--box_center X Y Z` | Box center in Angstrom. |
| `--box_center_off_reactive_res` | Center 5 A along CA-CB beyond CB; exactly one reactive flexres. Needs `--box_size`. |
| `--box_enveloping FILE` | Box envelops atoms of this file (`.sdf .mol .mol2 .pdb .pdbqt`). |
| `--padding P` | Padding around `--box_enveloping` atoms (A). |

## Flexible / reactive residues

| Option | Meaning |
|---|---|
| `-f, --flexres SEL` | Flexible sidechains, e.g. `-f ":42,B:23"`. Repeatable. |
| `-t, --rot_terminal_group SEL` | Make terminal functional group rotatable (e.g. terminal Glu/Asp carboxylate). |
| `-r, --reactive_flexres SEL` | Reactive flexible residues. **Max 8.** |
| `--reactive_name RES:ATOM` | Reactive atom per residue TYPE, e.g. `TRP:NE1`. Repeatable. Overridden by `-s`. |
| `-s, --reactive_name_specific SEL=ATOM` | Reactive atom per residue ID, e.g. `A:42=NE2`; residue becomes reactive. |
| `--r_eq_12 V` | Equilibrium distance for reactive 1-2 interactions (default 1.8 A). |
| `--eps_12 V` | Epsilon for reactive 1-2 interactions (default 2.5). |
| `--r_eq_13_scaling F` | r_eq scaling for 1-3 (default 0.5). |
| `--r_eq_14_scaling F` | r_eq scaling for 1-4 (default 0.5). |

## Default reactive atoms (from mk_prepare_receptor source)

```text
SER:OG   LYS:NZ   TYR:OH   CYS:SG   HIE:NE2   HID:ND1   GLU:OE2   THR:OG1   MET:SD
```

## Built-in template sets

- `residue_chem_templates.json`: standard protein/DNA/RNA residues, Amber24 modified
  residues (amino19, mod_amino19, RNA/DNA libs, ...), CCD entries with conflicting names.
- `NAKB_templates.json` (add via `--add_templates`): modified nucleotides from NAKB.
- Runtime CCD fetch: unknown residue names are resolved by downloading a definition CIF
  from rcsb.org (internet required) - supports noncovalent ligands and unmodified
  backbones only; use `chemtempgen` for the rest (see `python-api.md`).

