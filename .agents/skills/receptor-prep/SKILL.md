---
name: receptor-prep
description: Prepare protein/receptor structures for AutoDock/Vina docking with Meeko - PDB/mmCIF/PQR to docking-ready PDBQT, residue template matching, protonation-state templates, alternate locations (altloc), waters/HETATM cleanup, flexible sidechains (flexres), reactive/covalent residues, grid box and GPF files. Use when the task mentions receptor preparation, mk_prepare_receptor.py, Polymer, flexres, flexible residues, reactive docking receptor, GPF, grid box, altloc, blunt ends, or residue templates.
---

# Receptor Preparation (Meeko)

Turn a PDB/mmCIF/PQR structure into docking-ready **PDBQT** with `mk_prepare_receptor.py`.
Tested against the Meeko version installed in this repo (check `uv run python -c "import meeko; print(meeko.__version__)"`; current: 0.8.0).

## Mental model

A receptor is a `Polymer` of `Monomer`s. **Chemistry comes from residue templates; only
positions and inter-residue connectivity come from the input file.** Each residue is
matched against templates (813 built-in: standard protein/DNA/RNA, Amber24 modified
residues, conflicting CCD entries) using **connectivity and elements only** - not bond
orders, not atom names. Heavy atoms must match exactly (missing or excess = failure);
hydrogens are optional and get added from template SMILES without force-field
minimization.

Consequences:

- Meeko does **not** predict protonation states - it adds H matching whichever template
  matched. His tautomer/protonation, neutral LYS (LYN), neutral ASP/GLU (ASH/GLH) and
  disulfide CYX templates are auto-tried; force one with `-n`.
- Input coordinates must be sane. Bad geometry breaks template matching or, worse,
  makes RDKit fail to even build the residue.

## Standard workflow (validated)

```bash
# 0. (recommended) inspect the file: chains, altlocs, waters, HETATM
uv run python .agents/skills/receptor-prep/scripts/recprep.py structure.pdb

# 1. clean: pick altloc, drop waters + HETATM, delete hopeless residues
uv run python .agents/skills/receptor-prep/scripts/recprep.py structure.pdb \
    --remove-waters --remove-het --altloc A --delete A:119 --out rec_clean.pdb

# 2. prepare: PDBQT + Vina box config
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    --box_center 15 19 47 --box_size 30 30 30

# 3. validate
uv run python .agents/skills/receptor-prep/scripts/check_receptor.py rec.pdbqt rec.box.txt
```

`-o rec -p` writes `rec.pdbqt`; with flexible residues it writes `rec_rigid.pdbqt` +
`rec_flex.pdbqt`. `-v` writes `rec.box.txt` (center/size for Vina) + `rec.box.pdb` (visualize).
Writing a GPF (`-g`) requires `-p` and a box (see combos below).

## The two failure stages - different remedies  **(critical)**

| Stage | Error signature | Bypass with `-x`? | Fix |
|---|---|---|---|
| 1. RDKit mol build (per residue) | `unable to build rdkit mol for residue X corresponding to key A:N`, `AtomValenceException` | **NO** | fix the input file: `recprep.py --delete A:N`, correct altloc, or fix coordinates |
| 2. Template matching | `Template matching failed for: ['A:119', ...]`, `No template matched` | **YES** (`-x/--delete_bad_res`) | usually waters/HETATM near protein (remove them), altlocs (`--wanted_altloc`), missing sidechain atoms, steric clashes; or force a template with `-n` |

Stage-1 failures come from bad geometry (e.g. 1990s structures with 0.9 A C-O bonds ->
RDKit perceives extra bonds -> valence error). Stage-2 failures name the residues and
list what each candidate template expected (`heavy_miss`, `bond_excess`) - read the log.
Old/deposited structures with clash clusters (REMARK 500) often fail stage 2 at exactly
those residues; prefer re-refined PDB entries.

`--delete_bad_res_from_box_radius R` is a targeted alternative to `-x`: unmatched residues
more than R outside the box are deleted, anything near/inside the box still errors.
`--forgive_extra_bonds` tolerates excess inter-residue bonds (risky).

## Residue selection & assignment language

```text
A:5,7,BB:12C          chain A residues 5 and 7 + chain BB residue 12C (insertion code)
:42                   residue 42 with EMPTY chain ID
```

Selection options: `-d/--delete_residues`, `-f/--flexres`, `-r/--reactive_flexres`,
`-t/--rot_terminal_group`.

Assignment options take `selection=value`:

```text
-n "A:5,7=CYX,A:19,B:17=HID"     force templates (protonation/state)
-b "A:123,200=2,A:1=0"           blunt ends (atom index; artificial chain breaks)
--wanted_altloc ":5=B,B:17=A"    altloc per residue (overrides --default_altloc)
-s "A:42=NE2"                    reactive atom for one residue (makes it reactive)
```

## Protonation / templates

- Auto: all candidate templates for a residue name are tried (HIS -> HID/HIE/HIP ...).
- Force: `-n "A:57=HIP,A:15=HID,A:6=LYN,A:220=ASH"`.
- Unknown residues (cofactors, modified AAs): Meeko auto-fetches a definition CIF from
  the PDB by CCD code (needs internet; `--ignore_https_cert` last resort), or supply
  templates with `--add_templates file.json` (repeatable; also `resname:file.sdf`).
  `--cache_templates` caches to `~/.meeko_residue_chem_templates_cached.json`.
- Custom/exotic residues (modified backbones): build templates with `meeko.chemtempgen`
  - see `references/python-api.md`.

## Flexible sidechains (induced fit, Vina `--flex`)

```bash
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -x \
    -f "A:20,A:103" --box_center 15 19 47 --box_size 30 30 30
# -> rec_rigid.pdbqt + rec_flex.pdbqt
uv run python .agents/skills/receptor-prep/scripts/check_receptor.py rec_flex.pdbqt
```

Dock: `vina --receptor rec_rigid.pdbqt --flex rec_flex.pdbqt --ligand lig.pdbqt ...`.
Keep it to few solvent-exposed residues (2-5 sidechains; never buried/core residues).
`-t "A:231"` also rotates terminal functional groups (e.g. carboxylate of a terminal Glu/Asp).

## Reactive docking (covalent)

```bash
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -x \
    -r "A:23" --box_center_off_reactive_res --box_size 30 30 30
```

- Max **8** reactive residues. Default reactive atoms: `SER:OG, LYS:NZ, TYR:OH, CYS:SG,
  HIE:NE2, HID:ND1, GLU:OE2, THR:OG1, MET:SD`; override per type
  `--reactive_name "TRP:NE1"` (repeatable) or per residue `-s "A:42=NE2"` (that residue
  becomes reactive).
- Writes `rec_rigid.pdbqt`, `rec_flex.pdbqt` (reactive atoms get types like `1O4`, `1H5`)
  and `rec.reactive_config` (for AutoDock-GPU; the ligand side is prepared with
  `mk_prepare_ligand.py --reactive_smarts ...`, see the ligand-prep skill).
- `--r_eq_12 1.8 --eps_12 2.5 --r_eq_13_scaling 0.5 --r_eq_14_scaling 0.5` tune the
  covalent restraints. `--box_center_off_reactive_res` needs exactly one reactive residue.

## Charges

Default: **Gasteiger charges read from the template file** (fast). `--compute_charges`
recomputes with `--charge_model {gasteiger, nagl, espaloma, zero, read}` (`read` needs
`--read_pqr`; nagl/espaloma are slow). `--adj_padding` (pad residues with neighbors
instead of template reactions) requires `--compute_charges`. `--config_preset scofu1` /
`--config_file x.json` override templates+charges wholesale (CLI args win).

## Grid box: the three valid combos

```text
--box_center X Y Z --box_size X Y Z
--box_center_off_reactive_res --box_size X Y Z     (one reactive residue only)
--box_enveloping FILE --padding N                  (FILE: .sdf/.mol/.mol2/.pdb/.pdbqt)
```

## Input formats

| Flag | Format | Needs |
|---|---|---|
| `--read_pdb` | PDB | RDKit only |
| `-i/--read_with_prody` | PDB **or mmCIF** | **ProDy (not installed here)** |
| `--read_pqr` | PQR | charges `read` with `--charge_model read` |
| `--read_json` | Meeko JSON | round-trip; flexres state is lost |

For mmCIF: `uv add prody` (or `uv pip install prody`) then `-i`. Altlocs and multi-MODEL
files need handling regardless of format.

## Validation checklist

1. `check_receptor.py` on all outputs: sane type histogram (C/A/N/NA/OA/SA/HD ...),
   no NaN charges, total charge near an integer (net protein charge).
2. Atom count: rigid PDBQT atoms = input heavy atoms + added polar H - deleted residues.
3. Read the run log: every `tried N templates` line that did NOT end in success.
4. Flex blocks: one `BEGIN_RES` per requested residue; reactive types (`1O4` etc.) only
   where intended.
5. Box: `rec.box.pdb` overlaps the binding site in a viewer before docking.

## Common problems

- `unable to build rdkit mol` -> stage-1: delete the residue via `recprep.py --delete` or
  fix its coordinates; `-x` will NOT save you.
- `Template matching failed ... HOH` -> remove waters (`recprep.py --remove-waters`).
- `matched with excess inter-residue bond(s)` -> steric clash / close contact; clean the
  structure or `-x` if far from the site.
- Altloc crash -> `recprep.py --altloc A` or `--default_altloc A` / `--wanted_altloc ":5=B"`.
- `missing box center to write .gpf` -> use one of the three combos above.
- Missing sidechain atoms -> model them in (e.g. PDBFixer) or delete the residue.
- Very old entry with dozens of REMARK 500 clashes -> find a re-refined structure.

## Reference files

- `references/cli-options.md` - every `mk_prepare_receptor.py` option (verified on 0.8.0)
- `references/python-api.md` - `Polymer` API, `PDBQTWriterLegacy`, `chemtempgen` templates
- `references/recipes.md` - copy-pasteable end-to-end workflows
- `scripts/recprep.py` - inspect/clean PDB before Meeko (altloc, waters, HETATM, deletions)
- `scripts/check_receptor.py` - validate PDBQT/GPF/box/reactive_config outputs
- Sibling skill: `.agents/skills/ligand-prep` for the ligand side.
- Upstream docs: <https://meeko.readthedocs.io/en/develop/rec_overview.html>,
  `cli_rec_prep.html`, `py_rec_prep.html`, `rec_cli_options.html`, `py_build_temp.html`
