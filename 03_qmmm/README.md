# QM/MM

- `results/energy_table.csv`: relative energies (kcal/mol, reference: reactant) of all structures on the path
- `results/path_distances.csv`: key distances for the same structures
- `results/ts_validation.csv`: checks of the transition state
- `results/metal_ligand_distances.csv`: metal–ligand distances
- `structures/`: reactant (`R_env.pdb`), transition state (`TS_env.pdb`), structure after cleavage (`P_env.pdb`);
  residues and waters within 12 Å of the scissile P. B-factor column: 0 surroundings (MM), 1 QM, 2 reaction core, 3 Mn.
  Occupancy column marks key atoms: 1 P, 2 O3′, 3 Ow, 4 H27260, 5 O2P, 6 O14871, 7 H14873, 8 O69876, 9 H69877, 10 Mn3, 11 Mn4, 12 O1P.
- `tools/fd52v31.py`: driver for finite-difference partial Hessians with CP2K and harmonic thermochemistry

Order of events on the path (from `results/path_distances.csv`):
1. Before the transition state: the water oxygen binds P (3.46 → 1.72 Å) and its proton (H27260) moves to O2P.
2. At the transition state: P–O3′ breaks (2.30 Å) while H69877 (Mn-bound water → relay water) and
   H14873 (relay water → O3′) are in flight. These two protons and the P–O3′ stretch make up the imaginary mode.
No intermediate was found as a minimum.
