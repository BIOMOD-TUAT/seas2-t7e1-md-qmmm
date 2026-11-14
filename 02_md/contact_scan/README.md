# Contact scan: all T7E1 residues vs. inserted adenines

Contact occupancy (%) between every protein residue (chains A and B, residues 17-145, 258 residues) and the
inserted adenines, for 1A-5A. 0A has no inserted adenine.

- contact: shortest heavy-atom distance <= 0.45 nm (heavy atom = atom name not starting with H)
- occupancy: frames in contact / 501 x 100
- frames: 501 per condition; 0-50,000 ps (3A: 400-50,400 ps, the production run starts at 400 ps)
- inserted adenines: residues C9 to C(8+n)

`contact_top10_by_range.csv`: the ten residues with the largest difference between conditions.
236 of 258 residues have 0 % in all five conditions (chain A 122 / 129, chain B 114 / 129).

Check against the earlier 50 ns values for 3A, chain B (this run / earlier): Lys67 10.38 / 10.38, Gly68 99.80 / 99.8,
Leu69 100.00 / 100.0, Glu71 50.30 / 49.9, Lys103 99.20 / 99.0.

Notes
- Ser17 of chain A is the N-terminus of the model (residues 1-16 are missing in the crystal structure) and Gly18 is next
  to it. Their contacts may be an artefact of the truncation.
- One trajectory per condition. The starting structures differ between conditions.
- Still to be added from the cluster: `contact_all_residues.csv`, `contact_matrix_wide.csv`,
  `contact_top20_per_condition.csv`, `contact_five_residues_check.csv`, `contact_all_residues_3A_500ns.csv`
  and the scripts (`cs_prep.py`, `cs_job.sh`, `cs_aggregate.py`, `cs_agg_job.sh`).
