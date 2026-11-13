# SEAS² Dry Lab: MD and QM/MM study of T7 endonuclease I with bulged DNA

Research record for the computational part of the BIOMOD 2026 project **SEAS²** (Team TUAT).
This repository collects the inputs, analysis definitions, result tables, figure scripts and source data
behind the "MD / QM/MM" page of the team wiki. MD and QM/MM are kept together because the QM/MM model
was built from the MD structures.

> 日本語の概要は [docs/README_ja.md](docs/README_ja.md) にあります。

## What was done

| Step | Method | System | Length / size |
|---|---|---|---|
| 1 | Classical MD (GROMACS 2025.1, CHARMM36, TIP3P) | WT T7E1 homodimer + 16-bp DNA with 0–5 inserted adenines (0A–5A) | 50 ns each |
| 2 | Classical MD | 3A complex | 500 ns |
| 3 | Active-site rebuild | 3A complex, metal–protein distances restrained toward crystal structures (PDB 1M0D, 2PFJ) | – |
| 4 | QM/MM (CP2K 2024.1, UKS PBE-D3, DZVP-MOLOPT-SR-GTH, electrostatic embedding) | 3A complex, 95,456 atoms, 52 QM atoms | reactant, transition state, cleavage-side structure |

## Main results (and what they do not show)

- A reaction path for cleavage of the A7–C8 phosphodiester bond was obtained: attack by a Mn-bound water,
  transfer of three protons, and breaking of the P–O3′ bond.
- Transition state: ΔE‡ = +57.6 kcal/mol, ΔG‡(310.15 K) = +52.7 kcal/mol (partial Hessian, 20 atoms),
  one imaginary frequency (65.5i cm⁻¹), connected to the reactant side and the cleavage side.
- **This barrier is far higher than expected for an enzyme.** No stable product was found (the structure after
  cleavage is only 0.05 kcal/mol below the transition state), and the largest force at the transition state
  (1.14 × 10⁻³ Ha/Bohr) is slightly above the target (1.0 × 10⁻³). The calculation therefore gives a candidate
  mechanism and shows the problems of the present model; it does not reproduce the real reaction quantitatively.
- MD: contacts with the inserted adenines and the geometry of the attacking water change strongly with time in
  the 500 ns run of 3A. The 50 ns comparison of 0A–5A cannot rank the bulge lengths.
- The differences in Mn–phosphate distance between conditions were already present in the starting models.

See [docs/limitations.md](docs/limitations.md) before using any number from this repository.

## Layout

```
01_model_building/   how the T7E1–DNA complexes were assembled (known and unknown steps)
02_md/
  inputs/            GROMACS inputs                          <-- to be added, see docs/todo.md
  analysis_definitions/   definitions, atom mapping, validation report of the MD analyses
  analysis_csv/      per-frame and summary tables (0A–5A, 50 ns; 3A, 500 ns blocks)
03_qmmm/
  inputs/            CP2K inputs                             <-- to be added, see docs/todo.md
  tools/             finite-difference partial Hessian and thermochemistry driver for CP2K
  results/           energies, distances along the path, transition-state validation
  structures/        reactant / transition state / after cleavage (atoms within 12 Å of the scissile P)
04_figures/
  scripts/           Python scripts that drew the figures (the figures themselves are on Google Drive)
  chimerax/          ChimeraX command files for the 3D pictures
  vmd/               VMD script for the reaction movie
docs/                methods summary, limitations, figure index, to-do list
```

## Figures and tables

The figures, their editable PowerPoint files and the source-data workbooks are kept on the team's Google Drive
(link to be added here). This repository holds the scripts and command files that produced them (`04_figures/`).

## Software

| | version |
|---|---|
| GROMACS | 2025.1 |
| CP2K | 2024.1 |
| CHARMM-GUI | v3.7 (3A system) |
| Python | 3, with numpy, scipy, matplotlib, Pillow, python-pptx, openpyxl |
| ChimeraX, VMD, PyMOL | model building and pictures |

Calculations were run on TSUBAME4.0.

## Partial Hessian tool

`03_qmmm/tools/fd52v31.py` builds a finite-difference partial Hessian with CP2K single points and gives harmonic
thermochemistry. It needs a JSON config through the environment variable `FD_CFG`
(`fd_config_example.json` shows the keys; the atom numbers in it are from an earlier model).

```
export FD_CFG=fd_config.json
python3 fd52v31.py validate TEMPLATE.inp
python3 fd52v31.py prep SET --template TEMPLATE.inp --restart R.restart --wfn R.wfn --mode full --natoms 20
python3 fd52v31.py work SET          # inside a batch job
python3 fd52v31.py hess SET
python3 fd52v31.py thermo fd_R/hess.json fd_TS/hess.json --T 310.15
```

## Not in this repository

Trajectories, wavefunction files and full-resolution cube files are too large for GitHub.
See [docs/todo.md](docs/todo.md) for where they are to be deposited.

## License

Not yet chosen. See [docs/todo.md](docs/todo.md).
