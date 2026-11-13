CP2K inputs are to be added here (see `docs/todo.md`).

Settings used: CP2K 2024.1; UKS PBE-D3; DZVP-MOLOPT-SR-GTH; GTH-PBE pseudopotentials; plane-wave cut-off 400 Ry;
QM region 52 atoms + 7 hydrogen link atoms, charge +1, multiplicity 11; QM cell 24 Å;
electrostatic embedding (GAUSS); MM: CHARMM36 / TIP3P; 95,456 atoms in total.

Reaction coordinates: RC_P = d(P–Ow) − d(P–O3′); proton transfer from the nucleophilic water to O2P;
proton transfer from the relay water to O3′.
Transition state: dimer method, then single-point forces and a finite-difference partial Hessian (20 atoms).
For finite differences and single points the QM cell was not re-centred (`CENTER NEVER`), with the mid-point of
the QM atoms and their MM link partners at the cell centre.
