# autodock-vina-agent

A workspace for **structure-based molecular docking with AutoDock Vina**, set up to be
driven end-to-end by an AI coding agent. The agent (pi / Claude-style harness) executes
the full docking pipeline — receptor preparation, ligand preparation, docking, and
result parsing — using Vina CLI binaries bundled in the repo and
[Meeko](https://github.com/forlilab/Meeko) for PDBQT preparation. Domain expertise is
encoded as agent **skills** under `.agents/skills/`, so the agent prepares structures
and tunes dockings correctly without step-by-step instruction.

## Repository layout

| path | contents |
|---|---|
| `bin/` | AutoDock Vina binaries (`vina`, `vina_split`) — **not included; download them yourself** (see Setup) |
| `.agents/skills/` | agent skills: `ligand-prep`, `receptor-prep`, `simple-docking` (+ helper scripts) |
| `data/` | worked example: imatinib docked into c-Abl kinase (PDB [1IEP](https://www.rcsb.org/structure/1iep)) |
| `pyproject.toml`, `uv.lock` | Python environment (managed with [uv](https://docs.astral.sh/uv/)) |

## Setup

Requires Python ≥ 3.14 and [uv](https://docs.astral.sh/uv/). Works on
**Linux** (x86_64/arm64), **macOS** (Intel/Apple Silicon), and **Windows**:

    uv sync

This installs Meeko, RDKit, gemmi, scipy, and tqdm from the lockfile. Run all Python
tools with `uv run`.

### Vina binaries

The AutoDock Vina binaries are **not distributed with this repo** — you must download
them yourself and place them in `bin/`:

**Linux / macOS** — replace `<PLATFORM>` with `linux_x86_64`, `linux_aarch64`,
`mac_x86_64`, or `mac_aarch64`:

    mkdir -p bin
    curl -L -o bin/vina https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_<PLATFORM>
    curl -L -o bin/vina_split https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_split_1.2.7_<PLATFORM>
    chmod +x bin/vina bin/vina_split

**Windows** — download `vina_1.2.7_win.exe` and `vina_split_1.2.7_win.exe` from the
[releases page](https://github.com/ccsb-scripps/AutoDock-Vina/releases) into `bin\`
(e.g. as `bin\vina.exe`) and point the tools at them via the `VINA_BIN` environment
variable or the `--vina-bin` flag.

All docking steps call the binaries directly as `bin/vina` and `bin/vina_split`. The
skills and examples assume **Vina 1.2.7**; newer versions should work as well.

## Workflow

1. **Ligand preparation** — SMILES/SDF/MOL2 → 3D, H-explicit SDF → PDBQT
   (Meeko: atom typing, charges, rotatable bonds, macrocycles).
2. **Receptor preparation** — PDB/mmCIF/PQR → protonated, docking-ready PDBQT
   (residue templates, altloc handling, water/HETATM cleanup, optional flexible residues).
3. **Docking** — `bin/vina` rigid-receptor docking with a defined search box
   (config files, exhaustiveness/seed control, batch docking via `vina_split`).
4. **Results** — score tables, pose parsing, PDBQT → SDF export for viewing.

Each step has a skill with detailed docs and helper scripts:

| skill | covers |
|---|---|
| [`.agents/skills/ligand-prep`](.agents/skills/ligand-prep/SKILL.md) | `mk_prepare_ligand.py`, atom typing, charges, reactive/tethered ligands, PDBQT validation |
| [`.agents/skills/receptor-prep`](.agents/skills/receptor-prep/SKILL.md) | `mk_prepare_receptor.py`, protonation templates, flexres, grid box / GPF |
| [`.agents/skills/simple-docking`](.agents/skills/simple-docking/SKILL.md) | `bin/vina` CLI, search boxes, exhaustiveness/seeds, score-only/local-only, pose export |

## Example results

**`data/` — imatinib into c-Abl kinase (1IEP):** best mode **−13.25 kcal/mol**
(exhaustiveness 32, seed 42), reproducing the crystallographic binding mode.

## Acknowledgments

If you use Vina in published work, please cite:

- J. Eberhardt, D. Santos-Martins, A. F. Tillack, S. Forli,
  *AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings*,
  J. Chem. Inf. Model. (2021), DOI [10.1021/acs.jcim.1c00203](https://doi.org/10.1021/acs.jcim.1c00203)
- O. Trott, A. J. Olson,
  *AutoDock Vina: improving the speed and accuracy of docking with a new scoring function,
  efficient optimization, and multithreading*, J. Comput. Chem. (2010),
  DOI [10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

## License

[MIT](LICENSE) © 2026 Du Shiqiao
