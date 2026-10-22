# Receptor preparation recipes

Copy-pasteable workflows, all validated against Meeko 0.8.0. Paths assume project root;
`recprep.py` / `check_receptor.py` live in this skill's `scripts/` directory.

## 1. Clean pipeline: PDB -> PDBQT + box (no flexres)

```bash
P=structure.pdb   # e.g. downloaded from rcsb.org

# inspect first - shows chains, altlocs, waters, HETATM, hydrogens
uv run python .agents/skills/receptor-prep/scripts/recprep.py $P

# clean: one altloc, no waters, no HETATM (add --delete A:119 for hopeless residues)
uv run python .agents/skills/receptor-prep/scripts/recprep.py $P \
    --altloc A --remove-waters --remove-het --out rec_clean.pdb

uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    --box_center 15 19 47 --box_size 30 30 30

uv run python .agents/skills/receptor-prep/scripts/check_receptor.py rec.pdbqt rec.box.txt
```

Real-world example (lysozyme 4LYZ): removing waters cleared all template failures
except a 3-residue clash cluster (ASP119/ALA122/ARG125, CB-NH1 = 1.36 A in the crystal);
those are template-MATCH failures, so `-x` bypasses them cleanly:

```bash
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -x -v \
    --box_center 15 19 47 --box_size 30 30 30
```

## 2. Decide: bypass or fix? (the two failure stages)

```text
"unable to build rdkit mol for residue X"   -> stage 1, -x will NOT work.
   fix: recprep.py --delete CHAIN:RESNUM   (or fix coords / pick altloc)
"Template matching failed for: [...]"       -> stage 2, -x works.
   better: read the log (heavy_miss/bond_excess per template) and fix the cause:
     - HOH in the list            -> --remove-waters
     - altloc mentioned           -> --altloc / --wanted_altloc / --default_altloc
     - heavy_miss on one residue  -> missing sidechain atoms: model in or delete
     - bond_excess / clash cluster-> delete residues far from site, or remodel
```

Targeted bypass (keep residues near the binding site strict):

```bash
mk_prepare_receptor.py ... --box_center 15 19 47 --box_size 30 30 30 \
    --delete_bad_res_from_box_radius 10
```

## 3. Box around a ligand (co-crystal) instead of manual center

```bash
# box envelops the co-crystallized ligand, 6 A padding
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    --box_enveloping cocrystal_ligand.pdb --padding 6

# then dock with Vina using rec.box.txt values (center_x/size_x ...)
```

## 4. Flexible sidechain docking (Vina --flex)

```bash
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -x -v \
    -f "A:20,A:103" \
    --box_enveloping cocrystal_ligand.pdb --padding 6

uv run python .agents/skills/receptor-prep/scripts/check_receptor.py \
    rec_rigid.pdbqt rec_flex.pdbqt

uv run vina --receptor rec_rigid.pdbqt --flex rec_flex.pdbqt \
    --ligand ligand.pdbqt --config rec.box.txt --out out.pdbqt
```

Guidelines: 2-5 solvent-exposed residues near the site; never buried/core residues;
avoid ARG/LYS unless the charge interacts with the ligand (big, floppy).

## 5. Reactive (covalent) docking setup

```bash
# receptor side: TYR A:23 becomes reactive (default reactive atom OH)
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -x \
    -r "A:23" --box_center_off_reactive_res --box_size 30 30 30

uv run python .agents/skills/receptor-prep/scripts/check_receptor.py \
    rec_rigid.pdbqt rec_flex.pdbqt rec.reactive_config
```

- `rec.reactive_config` + `rec_rigid.pdbqt` go to **AutoDock-GPU**; the map grid GPF-like
  file lists ligand types and `intnbp` covalent pair parameters.
- Ligand side: `mk_prepare_ligand.py --reactive_smarts ...` (see ligand-prep skill).
- Custom reactive atom for a residue type / instance:

```bash
... --reactive_name "TRP:NE1" ...      # any TRP listed in -r reacts via NE1
... -s "A:42=NE2" ...                  # exactly HIS A:42 reacts via NE2
```

Defaults: SER:OG, LYS:NZ, TYR:OH, CYS:SG, HIE:NE2, HID:ND1, GLU:OE2, THR:OG1, MET:SD.

## 6. Protonation states and special residues

```bash
# force: ASP A:6 neutral (ASH), HIS A:57 doubly protonated (HIP),
#        HIS A:15 delta tautomer (HID), CYS A:220 free thiol (CYS not CYX)
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    -n "A:6=ASH,A:57=HIP,A:15=HID,A:220=CYS" \
    --box_center 15 19 47 --box_size 30 30 30
```

Template names for common states: `HID` (ND1-H), `HIE` (NE2-H), `HIP` (+1, both),
`LYN` (neutral LYS), `ASH`/`GLH` (neutral ASP/GLU), `CYX` (disulfide), `CYS` (thiol).
Without `-n`, Meeko tries all candidates and picks by connectivity (positions are not
energy-minimized, so H placement follows the template, not H-bond optimization).

## 7. Chain breaks / blunt ends

Artificial breaks (missing loops, ligand-like polymer fragments):

```bash
# blunt end at residue A:123 atom index 2, and A:200 atom 2; A:1 atom 0
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    -b "A:123,200=2,A:1=0" --box_center 15 19 47 --box_size 30 30 30
```

## 8. mmCIF input

```bash
uv add prody    # or: uv pip install prody
uv run mk_prepare_receptor.py -i structure.cif -o rec -p -v \
    --box_center 15 19 47 --box_size 30 30 30
```

## 9. AutoGrid (AutoDock4) GPF

```bash
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -g \
    --box_center 15 19 47 --box_size 30 30 30
# -> rec_rigid.gpf (or rec.gpf), rec.pdbqt must exist: -g requires -p

autogrid4 -p rec.gpf -l rec.glg
```

## 10. All-Python receptor prep

```python
from meeko import Polymer, MoleculePreparation, PDBQTWriterLegacy

mk_prep = MoleculePreparation()
polymer = Polymer.from_pdb_string(open("rec_clean.pdb").read(), allow_bad_res=True)
polymer.flexibilize_sidechain("A:20", mk_prep)
polymer.flexibilize_sidechain("A:103", mk_prep)

rigid, flex = PDBQTWriterLegacy.write_from_polymer(polymer)
open("rec_rigid.pdbqt", "w").write(rigid)
open("rec_flex.pdbqt", "w").write("".join(flex.values()))
```

Box/center still comes from the CLI or your own geometry math; `polymer.to_json()`
lets you persist the parameterized receptor for later runs.

## 11. Adding a cofactor/ligand template (CCD)

If a cofactor (e.g. NAP, FAD variant) fails to match:

```bash
# internet: Meeko fetches the CCD definition and builds templates automatically.
# If it cannot (modified backbone, covalent links): build templates, then:
uv run mk_prepare_receptor.py --read_pdb rec_clean.pdb -o rec -p -v \
    --add_templates my_cofactor.json \
    --box_center 15 19 47 --box_size 30 30 30
```

Template JSON building is in `references/python-api.md` (chemtempgen).

