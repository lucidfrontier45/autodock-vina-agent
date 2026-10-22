#!/usr/bin/env python
"""Sanity-check Meeko receptor-preparation outputs before docking.

Handles rigid PDBQT (_rigid / plain), flexible PDBQT (_flex, BEGIN_RES blocks),
Vina box config files, AutoGrid GPF files, and reactive_config files.

Usage:
  uv run python .agents/skills/receptor-prep/scripts/check_receptor.py receptor.pdbqt
  uv run python .agents/skills/receptor-prep/scripts/check_receptor.py *_rigid.pdbqt *_flex.pdbqt
  uv run python .agents/skills/receptor-prep/scripts/check_receptor.py --box rec.box.txt --gpf rec.gpf
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path

# Common AutoDock4 atom types (receptors also carry metal/halogen types).
AD4_TYPES = {
    "H", "HD", "HS", "C", "A", "N", "NA", "NS", "OA", "OS", "F", "Mg", "MG",
    "P", "SA", "S", "Cl", "CL", "Ca", "CA", "Mn", "MN", "Fe", "FE", "Zn", "ZN",
    "Br", "BR", "I", "Se", "SE", "Si", "SI", "B", "W",
}
# Reactive-docking receptor types look like 1O4, 1H5, N3, C7, S2, ... (numbered variants)
NUMBERED = re.compile(r"^[ACNOSIBP]\d?$")


def check_rigid(path: Path) -> int:
    atoms, problems = [], []
    with path.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                f = line.split()
                try:
                    q = float(f[10])
                except (ValueError, IndexError):
                    problems.append(f"unparsable charge field in: {line.rstrip()}")
                    continue
                atoms.append({"name": f[2], "res": f"{f[3]} {f[4]}{f[5]}",
                              "q": q, "t": f[-1]})
    if not atoms:
        print(f"ERROR {path.name}: no atoms found")
        return 1

    types = Counter(a["t"] for a in atoms)
    total_q = sum(a["q"] for a in atoms)
    unknown = {t: n for t, n in types.items()
               if t not in AD4_TYPES and not NUMBERED.match(t)}
    for a in atoms:
        if not math.isfinite(a["q"]):
            problems.append(f"non-finite charge on {a['name']} {a['res']}")
            break

    print(f"== {path.name} (rigid receptor) ==")
    print(f"atoms: {len(atoms)}")
    print("types: " + "  ".join(f"{t}:{n}" for t, n in types.most_common()))
    print(f"total partial charge: {total_q:+.3f}"
          + ("" if abs(total_q) < 1.0 else "   <- off by more than 1 e, investigate"))
    if unknown:
        print(f"WARNING: non-standard atom types {sorted(unknown)} "
              "(fine for reactive docking types you expect; suspicious otherwise)")
    for p in problems:
        print(f"ERROR: {p}")
    return 1 if problems else 0


def check_flex(path: Path) -> int:
    text = path.read_text()
    blocks = re.split(r"^BEGIN_RES ", text, flags=re.M)[1:]
    problems = []
    if not blocks:
        print(f"ERROR {path.name}: no BEGIN_RES blocks")
        return 1
    print(f"== {path.name} (flexible sidechains) ==")
    print(f"flexible residues: {len(blocks)}")
    for b in blocks:
        res = b.splitlines()[0].strip()
        n_atoms = sum(1 for l in b.splitlines() if l.startswith(("ATOM", "HETATM")))
        ok_end = "END_RES" in b
        n_br = sum(1 for l in b.splitlines() if l.startswith("BRANCH"))
        n_ebr = sum(1 for l in b.splitlines() if l.startswith("ENDBRANCH"))
        ok_close = n_br == n_ebr
        print(f"  {res}: {n_atoms} atoms"
              + ("" if ok_end else "  <- missing END_RES")
              + ("" if ok_close else f"  <- unbalanced BRANCH ({n_br}) / ENDBRANCH ({n_ebr})"))
        if not (ok_end and ok_close):
            problems.append(res)
    return 1 if problems else 0


def check_box(path: Path) -> int:
    center = size = None
    for line in path.read_text().splitlines():
        k, _, v = line.partition("=")
        k, v = k.strip(), v.strip()
        if k == "center_x": center = v
        elif k == "size_x": size = v
    print(f"== {path.name} (Vina box) ==")
    print(path.read_text().strip())
    return 0


def check_gpf(path: Path) -> int:
    types, npts, spacing = None, None, None
    for line in path.read_text().splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "npts":
            npts = f[1:4]
        elif f[0] == "spacing":
            spacing = f[1]
        elif f[0] == "map":
            types = types or []
            types.append(Path(f[1]).stem.split(".")[-1])
    print(f"== {path.name} (AutoGrid GPF) ==")
    print(f"grid points: {npts}, spacing: {spacing} A")
    print(f"maps ({len(types or [])}): {' '.join(sorted(set(types or [])))}")
    return 0


def check_reactive(path: Path) -> int:
    print(f"== {path.name} (AutoDock-GPU reactive config) ==")
    txt = path.read_text()
    first = txt.splitlines()[0] if txt else ""
    print(first[:120] + ("..." if len(first) > 120 else ""))
    n = sum(1 for l in txt.splitlines() if l.startswith("intnbp"))
    print(f"{n} intnbp reactive-pair lines")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", type=Path)
    args = ap.parse_args()
    if not args.files:
        ap.error("give me _rigid/_flex PDBQT, .box.txt, .gpf or .reactive_config files")

    rc = 0
    for p in args.files:
        name = p.name
        try:
            if name.endswith("_flex.pdbqt"):
                rc |= check_flex(p)
            elif name.endswith((".pdbqt",)):
                rc |= check_rigid(p)
            elif name.endswith(".gpf"):
                rc |= check_gpf(p)
            elif name.endswith("reactive_config"):
                rc |= check_reactive(p)
            elif name.endswith((".txt", ".box.txt")):
                rc |= check_box(p)
            else:
                print(f"skipping {p} (unknown type)")
        except FileNotFoundError:
            print(f"ERROR: {p} not found")
            rc |= 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
