#!/usr/bin/env python
"""Ligand preparation helper: SMILES/SDF -> H-explicit 3D SDF -> PDBQT.

Meeko does NOT add hydrogens and does NOT generate 3D coordinates; RDKit must do
that first. This script does both, then parameterizes with Meeko.

Examples
--------
# one SMILES -> SDF + PDBQT
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --smiles "CC(=O)Nc1ccc(O)cc1" --name paracetamol \
    --out ligands/paracetamol.sdf --pdbqt ligands/paracetamol.pdbqt

# a multi-molecule SDF -> 3D SDF + one PDBQT per molecule
uv run python .agents/skills/ligand-prep/scripts/ligprep.py \
    --sdf in.sdf --out out.sdf --pdbqt-dir out_pdbqt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.warning")


# --------------------------------------------------------------------------------------
# 3D preparation
# --------------------------------------------------------------------------------------
def build_3d(mol, seed: int = 0xF00D, optimize: bool = True) -> Chem.Mol:
    """Add explicit hydrogens and generate a single 3D conformer."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        # fall back to random coordinates + distance geometry
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError("RDKit could not embed a 3D conformer")
    if optimize:
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:  # pragma: no cover - MMFF may not cover the element set
            pass
    return mol


def read_molecules(path: Path):
    """Yield (name, mol) from an SD file (multi-molecule supported)."""
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    for idx, mol in enumerate(supplier):
        if mol is None:
            print(f"[warn] skipping unreadable molecule #{idx + 1} in {path}", file=sys.stderr)
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"{path.stem}_{idx + 1}"
        yield name, mol


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--smiles", action="append", default=[],
                     help="SMILES string (repeatable). Use 'name:SMILES' or --name to name it.")
    src.add_argument("--sdf", type=Path, help="input SD file (may contain several molecules)")
    ap.add_argument("--name", action="append", default=[],
                    help="molecule name; matches the order of --smiles")
    ap.add_argument("--out", type=Path, help="output SDF (H-explicit, 3D)")
    ap.add_argument("--pdbqt", type=Path, help="output PDBQT (single molecule only)")
    ap.add_argument("--pdbqt-dir", type=Path, help="directory for PDBQT files (one per molecule)")
    ap.add_argument("--config", type=Path, help="JSON config for MoleculePreparation")
    ap.add_argument("--charge-model", default=None,
                    choices=["gasteiger", "nagl", "espaloma", "zero", "read"],
                    help="override the charge model from the config")
    ap.add_argument("--merge-these-atom-types", nargs="*", default=None,
                    help="atom types to merge into their parent (default: H). "
                         "Pass no value to merge nothing.")
    ap.add_argument("--rigid-macrocycles", action="store_true")
    ap.add_argument("--flexible-amides", action="store_true")
    ap.add_argument("--add-index-map", action="store_true")
    ap.add_argument("--bad-charge-ok", action="store_true")
    ap.add_argument("--seed", type=int, default=0xF00D)
    ap.add_argument("--no-optimize", action="store_true", help="skip MMFF optimization")
    args = ap.parse_args(argv)

    # ---- collect input molecules -----------------------------------------------------
    entries = []  # (name, mol)
    if args.sdf is not None:
        entries = list(read_molecules(args.sdf))
    else:
        names = args.name or []
        for i, smi in enumerate(args.smiles):
            name = names[i] if i < len(names) else f"lig_{i + 1}"
            if ":" in name:  # "name:SMILES" form
                name, smi = name.split(":", 1)
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                print(f"[error] could not parse SMILES: {smi}", file=sys.stderr)
                return 2
            mol.SetProp("_Name", name)
            entries.append((name, mol))

    if not entries:
        print("[error] no molecules found", file=sys.stderr)
        return 2

    # ---- 3D -------------------------------------------------------------------------
    mols = []
    for name, mol in entries:
        try:
            mol3d = build_3d(mol, seed=args.seed, optimize=not args.no_optimize)
        except Exception as exc:
            print(f"[error] 3D embedding failed for {name}: {exc}", file=sys.stderr)
            return 2
        mol3d.SetProp("_Name", name)
        mols.append((name, mol3d))

    # ---- write SDF ------------------------------------------------------------------
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(args.out))
        for _, mol in mols:
            writer.write(mol)
        writer.close()
        print(f"[ok] wrote {args.out} ({len(mols)} molecule(s))")

    # ---- PDBQT ----------------------------------------------------------------------
    if args.pdbqt is None and args.pdbqt_dir is None:
        return 0
    if args.pdbqt is not None and len(mols) > 1:
        print("[error] --pdbqt takes a single molecule; use --pdbqt-dir", file=sys.stderr)
        return 2

    from meeko import MoleculePreparation, PDBQTWriterLegacy

    config = {}
    if args.config is not None:
        config = json.loads(args.config.read_text())
    if args.charge_model is not None:
        config["charge_model"] = args.charge_model
    if args.merge_these_atom_types is not None:
        config["merge_these_atom_types"] = tuple(args.merge_these_atom_types)
    if args.rigid_macrocycles:
        config["rigid_macrocycles"] = True
    if args.flexible_amides:
        config["flexible_amides"] = True

    mk_prep = MoleculePreparation.from_config(config) if config else MoleculePreparation()

    if args.pdbqt_dir is not None:
        args.pdbqt_dir.mkdir(parents=True, exist_ok=True)

    failures = 0
    for name, mol in mols:
        try:
            setups = mk_prep(mol)
        except Exception as exc:
            print(f"[error] Meeko preparation failed for {name}: {exc}", file=sys.stderr)
            failures += 1
            continue
        for i, setup in enumerate(setups):
            pdbqt, is_ok, err = PDBQTWriterLegacy.write_string(
                setup,
                add_index_map=args.add_index_map,
                bad_charge_ok=args.bad_charge_ok,
            )
            if not is_ok:
                print(f"[error] PDBQT generation failed for {name}: {err}", file=sys.stderr)
                failures += 1
                continue
            if args.pdbqt is not None:
                path = args.pdbqt
            else:
                suffix = "" if len(setups) == 1 else f"_{i + 1}"
                path = args.pdbqt_dir / f"{name}{suffix}.pdbqt"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(pdbqt)
            torsdof = next((ln for ln in pdbqt.splitlines() if ln.startswith("TORSDOF")), "TORSDOF ?")
            print(f"[ok] {path}  ({torsdof.strip()})")

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
