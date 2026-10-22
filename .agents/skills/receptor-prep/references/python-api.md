# Receptor preparation from Python (Meeko 0.8.0)

## Build a Polymer

```python
from meeko import Polymer, MoleculePreparation

polymer = Polymer.from_pdb_file("rec.pdb")
# or:
polymer = Polymer.from_pdb_string(pdb_string)
# or from ProDy object:
# polymer = Polymer.from_prody(mol)
# or from PQR string:
# polymer = Polymer.from_pqr_string(pqr_string)
```

`from_pdb_string` signature (0.8.0) - the same knobs as the CLI:

```python
Polymer.from_pdb_string(
    pdb_string,
    chem_templates=None,        # ResidueChemTemplates, or loaded via add_templates JSON
    mk_prep=None,               # MoleculePreparation instance (atom typing / charges)
    set_template=None,          # e.g. {"A:57": "HIP"}
    residues_to_delete=None,    # e.g. {"A:119"}
    ignore_https_cert=False,
    allow_bad_res=False,        # CLI -x equivalent (template-match failures only!)
    bonds_to_delete=None,
    blunt_ends=None,
    wanted_altloc=None,         # e.g. {"A:5": "B"}
    default_altloc=None,        # e.g. "A"
    forgive_extra_bonds=False,
    delete_bad_res_from_box_radius=None,
    box_size=None, box_center=None,
    adj_padding=False,
)
```

The parameterization is residue-by-residue, using the same `MoleculePreparation` as
ligands:

```python
mk_prep = MoleculePreparation(
    load_atom_params=["openff-2.3.0"],
    charge_model="nagl",          # needs recent openff-toolkit
    merge_these_atom_types=[],    # keep polar H for receptor H-bonds
    rigid_macrocycles=True,
)
polymer = Polymer.from_pdb_string(pdb_string, mk_prep=mk_prep)
```

Note: CLI default charges come from pre-calculated template values; passing a plain
`MoleculePreparation()` computes Gasteiger charges instead.

## Flexibilize sidechains

```python
polymer.flexibilize_sidechain("A:20", mk_prep)   # residue ID "chain:resnum"
polymer.rigidify_sidechain("A:20", mk_prep)
polymer.rigidify_all(mk_prep)
polymer.parameterize(mk_prep)                    # re-parameterize everything
polymer.monomers["A:3"].parameterize(mk_prep, "A:3")  # one residue
```

The `mk_prep` must be passed: flexibilization sets bonds rotatable/non-rotatable
according to its config. Use the same config as the polymer was built with.

## Write PDBQT

```python
from meeko import Polymer, MoleculePreparation, PDBQTWriterLegacy

mk_prep = MoleculePreparation()
polymer = Polymer.from_pdb_string(open("rec.pdb").read(), allow_bad_res=True)
polymer.flexibilize_sidechain("A:20", mk_prep)

rigid_pdbqt, flex_pdbqt = PDBQTWriterLegacy.write_from_polymer(polymer)
# rigid_pdbqt: str
# flex_pdbqt:  {"A:20": "<flex pdbqt string>"}  (dict residue-id -> string)
```

Concatenate all flex strings into one `rec_flex.pdbqt` for Vina `--flex`.

## Inspect / export

```python
polymer.monomers                  # dict "A:20" -> Monomer
polymer.get_valid_monomers()
polymer.get_ignored_monomers()
polymer.to_json()                  # round-trip: Polymer.from_json(...)
polymer.to_pdb(new_positions=None) # PDB string of prepared model (H added)
polymer.to_rdkit_mol()             # RDKit molecule of the whole polymer
```

## Residue templates

```python
from meeko import ResidueChemTemplates

rct = ResidueChemTemplates.create_from_defaults()
rct.residue_templates   # dict name -> {"smiles", "atom_name", "link_labels"}  (813 entries)
rct.ambiguous           # resname -> list of candidate template names (e.g. HIS: HID/HIE/HIP)
```

Pass `chem_templates=rct` to `Polymer.from_pdb_string`.

## Building custom templates (chemtempgen)

For residues Meeko cannot auto-derive (modified backbones, cross-links). Internet needed
for the CCD fetch. Sketch for `CRO` (GFP fluorophore, Ser-Tyr-Gly condensation):

```python
from meeko.chemtempgen import (ChemicalComponent, fetch_from_pdb,
                               export_chem_templates_to_json, AA_recipe)
import copy

cc0 = ChemicalComponent.from_cif(fetch_from_pdb("CRO"), "CRO")
embed_smarts = "[NX2][CX4][CX3][NX3][CX4][CX3](=O)[OX2]"   # the (modified) backbone

# mid-chain variant: drop leaving atoms (their H become implicit at the blunt ends)
cc = copy.deepcopy(cc0)
cc = cc.make_embedded(allowed_smarts=embed_smarts, leaving_names={"H2", "OXT"})
cc = (cc.make_pretty_smiles()
        .make_link_labels_from_patterns(pattern_to_label_mapping=
                                        AA_recipe.pattern_to_label_mapping_standard))
cc.ResidueTemplate_check()

# N-terminal variant: remove OXT, protonate (cap) the N
cc_N = copy.deepcopy(cc0)
cc_N = (cc_N.make_embedded(allowed_smarts=embed_smarts, leaving_names={"OXT"})
          .make_capped(allowed_smarts=embed_smarts, capping_names={"N1"}, protonate=True)
          .make_pretty_smiles()
          .make_link_labels_from_patterns(pattern_to_label_mapping=
                                          AA_recipe.pattern_to_label_mapping_standard))
cc_N.resname += "_N"

# C-terminal variant: deprotonate acid, remove H2
acidic = {"[H][O][C](=O)": 0}   # SMARTS -> proton index; see acidic_proton_loc_canonical
cc_C = copy.deepcopy(cc0)
cc_C = (cc_C.make_canonical(acidic_proton_loc=acidic)
          .make_embedded(allowed_smarts=embed_smarts, leaving_names={"H2"})
          .make_pretty_smiles()
          .make_link_labels_from_patterns(pattern_to_label_mapping=
                                          AA_recipe.pattern_to_label_mapping_standard))
cc_C.resname += "_C"

json_str = export_chem_templates_to_json([cc, cc_N, cc_C])   # or json_fname="CRO.json"
```

Resulting JSON shape (what `--add_templates` consumes):

```json
{
  "ambiguous": {"CRO": ["CRO", "CRO_N", "CRO_C"]},
  "residue_templates": {
    "CRO": {
      "smiles": "[H]NC([H])(...)[H]",
      "atom_name": ["H", "N1", "CA1", "..."],
      "link_labels": {"1": "N-term", "27": "C-term"}
    }
  }
}
```

Key rules:

- Template SMILES carry **real hydrogens** in square brackets `[H]`; atoms bonded to
  neighbors (link labels) have implicit H instead.
- `link_labels`: atom index (in SMILES order) -> `"N-term"` / `"C-term"` (or custom).
- `ambiguous` maps input resname -> candidate templates; matching uses connectivity +
  elements only, heavy atoms must match exactly.
- Load with `mk_prepare_receptor.py --add_templates CRO.json` or
  `chem_templates=ResidueChemTemplates.from_json_file("CRO.json") (see ResidueChemTemplates.add_json_file / add_dict to extend defaults).

