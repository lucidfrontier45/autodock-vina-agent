#!/usr/bin/env python
"""Inspect and clean a PDB file BEFORE running mk_prepare_receptor.py.

mk_prepare_receptor.py builds one RDKit molecule per residue first, then matches
templates. Residues that cannot even be built into an RDKit mol (valence errors
from bad geometry, mixed alternate locations, ...) are FATAL and cannot be
bypassed with -x/--delete_bad_res: they must be fixed or removed in the input
file itself. This script does exactly that.

Examples
--------
# just inspect (chains, altlocs, waters, HETATM, hydrogens)
uv run python .agents/skills/receptor-prep/scripts/recprep.py structure.pdb

# delete residues, drop waters and all other HETATM, pick altloc A, write out
uv run python .agents/skills/receptor-prep/scripts/recprep.py structure.pdb \
    --delete A:119,B:15,16 --remove-waters --remove-het --altloc A \
    --out structure_clean.pdb

Selection language (same as Meeko): "A:5,7,BB:12C" = residues 5 and 7 of chain A
plus residue 12C of chain BB. Empty chain = residues without a chain ID: ":42".
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------- PDB columns
def rec(line: str) -> str:
    return line[0:6].strip()


def altloc(line: str) -> str:
    return line[16] if len(line) > 16 else " "


def resname(line: str) -> str:
    return line[17:20].strip() if len(line) >= 20 else ""


def chain(line: str) -> str:
    return line[21] if len(line) > 21 else " "


def resseq(line: str) -> str:
    return line[22:26].strip() if len(line) >= 26 else ""


def icode(line: str) -> str:
    return line[26] if len(line) > 26 else " "


def element(line: str) -> str:
    e = line[76:78].strip() if len(line) >= 78 else ""
    if not e:
        e = line[12:16].strip()[:1] if len(line) >= 13 else ""
    return e


def res_key(line: str) -> tuple[str, str, str]:
    return (chain(line), resseq(line), icode(line))


# ---------------------------------------------------------- selection language
def parse_residues(spec: str) -> set[tuple[str, str]]:
    """'A:5,7,BB:12C' -> {('A','5'), ('A','7'), ('BB','12C')} (icode appended)."""
    out: set[tuple[str, str]] = set()
    cur = ""
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            cur, _, nums = tok.partition(":")
        else:
            nums = tok
        cur = cur.strip()
        n = nums.strip()
        if n:
            out.add((cur, n))
    return out


# ---------------------------------------------------------------------- inspect
def inspect(lines: list[str]) -> None:
    models = sum(1 for l in lines if rec(l) == "MODEL")
    atom_lines = [l for l in lines if rec(l) in ("ATOM", "HETATM")]

    chains: dict[str, set[tuple[str, str]]] = defaultdict(set)
    altlocs: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    waters = 0
    het: dict[str, int] = defaultdict(int)
    hydrogens = 0
    elements: dict[str, int] = defaultdict(int)

    for l in atom_lines:
        k = res_key(l)
        chains[chain(l)].add((resseq(l), icode(l).strip()))
        if altloc(l) != " ":
            altlocs[k].add(altloc(l))
        if resname(l) in ("HOH", "DOD"):
            waters += 1
        if rec(l) == "HETATM":
            het[resname(l)] += 1
        e = element(l)
        elements[e] += 1
        if e == "H":
            hydrogens += 1

    print(f"records: {len(lines)} lines, {len(atom_lines)} ATOM/HETATM, "
          f"{models} MODEL(s)" + ("  <- multi-model file!" if models > 1 else ""))
    print(f"chains: {sorted(chains)}")
    for c in sorted(chains):
        rs = sorted(chains[c])
        print(f"  chain {c or '(blank)'}: {len(rs)} residues (incl. waters/ions)")
    print(f"waters (HOH/DOD): {waters}")
    print(f"hydrogens: {hydrogens}")
    if het:
        print("HETATM residues: " + ", ".join(f"{k}({v})" for k, v in sorted(het.items())))
    else:
        print("HETATM residues: none")
    if altlocs:
        print(f"ALTERNATE LOCATIONS in {len(altlocs)} residues "
              "(mk_prepare_receptor refuses to guess; pass --altloc to pick one):")
        for k in sorted(altlocs)[:10]:
            print(f"  {k[0]}:{k[1]}{k[2].strip()} -> {''.join(sorted(altlocs[k]))}")
        if len(altlocs) > 10:
            print(f"  ... and {len(altlocs) - 10} more")
    problems = []
    if models > 1:
        problems.append("multiple MODEL records - split into one file per model")
    if altlocs:
        problems.append("alternate locations present - pick one with --altloc")
    if waters:
        problems.append("waters present - usually removed (--remove-waters); Meeko tries to "
                        "template-match them and fails if they sit close to the protein")
    if problems:
        print("\nSuggested cleanup:")
        for p in problems:
            print(f"  - {p}")


# ------------------------------------------------------------------------ main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb", type=Path)
    ap.add_argument("--out", type=Path, help="write cleaned PDB here (omit to inspect only)")
    ap.add_argument("--delete", metavar="A:119,B:15,16",
                    help="delete these residues from the file")
    ap.add_argument("--keep-chains", metavar="A,B", help="keep only these chains")
    ap.add_argument("--remove-waters", action="store_true", help="drop HOH/DOD residues")
    ap.add_argument("--remove-het", action="store_true",
                    help="drop all HETATM records (use with --keep-het to spare some)")
    ap.add_argument("--keep-het", metavar="NAG,HEM", help="HETATM resnames to keep")
    ap.add_argument("--altloc", metavar="A",
                    help="keep only this alternate location (blank ones always kept)")
    ap.add_argument("--drop-h", action="store_true", help="drop hydrogen atoms")
    args = ap.parse_args()

    lines = args.pdb.read_text().splitlines()
    print(f"== {args.pdb} ==")
    inspect(lines)
    if args.out is None and not any([args.delete, args.keep_chains, args.remove_waters,
                                     args.remove_het, args.altloc, args.drop_h]):
        return 0

    del_res = parse_residues(args.delete) if args.delete else set()
    keep_chains = {c.strip() for c in args.keep_chains.split(",")} if args.keep_chains else None
    keep_het = {r.strip() for r in args.keep_het.split(",")} if args.keep_het else set()

    out_lines: list[str] = []
    dropped = defaultdict(int)
    TER_is_kept_context = False
    for l in lines:
        r = rec(l)
        if r in ("ATOM", "HETATM"):
            c, s, i = res_key(l)
            if (c, s + i.strip()) in del_res:
                dropped["deleted residues"] += 1
                continue
            if keep_chains is not None and c not in keep_chains:
                dropped["other chains"] += 1
                continue
            if args.remove_waters and resname(l) in ("HOH", "DOD"):
                dropped["waters"] += 1
                continue
            if args.remove_het and r == "HETATM" and resname(l) not in keep_het:
                dropped["HETATM"] += 1
                continue
            a = altloc(l)
            if args.altloc and a != " " and a != args.altloc:
                dropped["altloc"] += 1
                continue
            if args.drop_h and element(l) == "H":
                dropped["hydrogens"] += 1
                continue
            out_lines.append(l)
        elif r in ("TER",) :
            # keep TER only if the previous atom line survived
            if out_lines and rec(out_lines[-1]) in ("ATOM", "HETATM"):
                out_lines.append(l)
            else:
                continue
        else:
            out_lines.append(l)

    if args.out is None:
        print("\nNo edits requested (--out not given): wrote nothing.")
        return 0

    args.out.write_text("\n".join(out_lines) + "\n")
    print(f"\nwrote {args.out}: kept {sum(1 for l in out_lines if rec(l) in ('ATOM','HETATM'))} atoms")
    for k, v in sorted(dropped.items()):
        print(f"  dropped {v:6d} {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
