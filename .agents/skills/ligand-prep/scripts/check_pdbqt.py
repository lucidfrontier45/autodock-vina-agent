#!/usr/bin/env python
"""Sanity-check Meeko-generated PDBQT files before docking.

Reports: TORSDOF, atom-type histogram, total partial charge, NaN/Inf charges,
unexpected atom types, and (with --sdf) the number of heavy atoms in the input.

Usage:
  uv run python .agents/skills/ligand-prep/scripts/check_pdbqt.py ligands/*.pdbqt
  uv run python .agents/skills/ligand-prep/scripts/check_pdbqt.py --sdf mol.sdf ligands/mol.pdbqt
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

# AutoDock4 atom types
AD4_TYPES = {
    "A", "NA", "OA", "SA", "HD", "H", "C", "N", "F", "Mg", "P", "S", "Cl",
    "Ca", "Mn", "Fe", "Zn", "Br", "I", "Si", "B", "W", "Ni", "Cu", "Fe",
    "G0", "CG0", "CG1", "W",
}


def parse_pdbqt(path: Path):
    atoms = []
    torsdof = None
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(("ATOM", "HETATM")):
                fields = line.split()
                serial = fields[1]
                name = fields[2]
                x, y, z = float(fields[6]), float(fields[7]), float(fields[8])
                charge = float(fields[10]) if fields[10] not in ("",) else float("nan")
                atype = fields[-1]
                atoms.append({"serial": serial, "name": name, "xyz": (x, y, z),
                              "charge": charge, "type": atype})
            elif line.startswith("TORSDOF"):
                torsdof = int(line.split()[1])
    return atoms, torsdof


def check(path: Path, sdf: Path | None) -> int:
    atoms, torsdof = parse_pdbqt(path)
    problems = []
    print(f"== {path}")
    if not atoms:
        print("   no ATOM records found")
        return 1
    types = Counter(a["type"] for a in atoms)
    total_charge = sum(a["charge"] for a in atoms)
    print(f"   atoms      : {len(atoms)}")
    print(f"   atom types : {dict(sorted(types.items()))}")
    print(f"   TORSDOF    : {torsdof}")
    print(f"   total charge: {total_charge:+.3f}")

    bad = [a for a in atoms if not math.isfinite(a["charge"])]
    if bad:
        problems.append(f"{len(bad)} atom(s) with NaN/Inf charge")

    unknown = sorted(set(types) - AD4_TYPES)
    if unknown:
        problems.append(f"non-AD4 atom types present: {unknown}")

    if sdf is not None:
        try:
            from rdkit import Chem
            supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
            mol = next((m for m in supplier if m is not None), None)
            if mol is not None:
                heavy = mol.GetNumHeavyAtoms()
                print(f"   input heavy atoms: {heavy} (from {sdf})")
                if heavy > len(atoms):
                    problems.append("PDBQT has fewer atoms than the input has heavy atoms")
        except Exception as exc:  # pragma: no cover
            print(f"   (could not read {sdf}: {exc})")

    if problems:
        for p in problems:
            print(f"   [warn] {p}")
        return 1
    print("   [ok] no problems detected")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbqt", nargs="+", type=Path, help="PDBQT file(s) to check")
    ap.add_argument("--sdf", type=Path, help="optional input SDF, for heavy-atom cross-check")
    args = ap.parse_args(argv)

    rc = 0
    for path in args.pdbqt:
        if not path.exists():
            print(f"[error] {path} not found", file=sys.stderr)
            rc = 2
            continue
        rc |= check(path, args.sdf)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
