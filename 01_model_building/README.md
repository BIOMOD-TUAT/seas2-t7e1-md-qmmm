# Model building

Sources: files, scripts and logs on the cluster ("record"), and the recollection of the person who built
the models ("builder"). Items marked *not recorded* must be confirmed before they are cited.

| Item | What is known | Source |
|---|---|---|
| Protein | WT T7E1 homodimer, chains A and B, residues 17–145 | record |
| Protein coordinates | 1M0D and 2PFJ fetched in PyMOL and superposed at the active site | builder |
| Which PDB supplied the final protein coordinates, repair of missing residues | *not recorded* | – |
| DNA | Built in ChimeraX (version not recorded); not taken from a crystal structure or a prediction | builder |
| DNA sequences | C strand 5′-TGTGGCAC-(A)n-CTGCACGC-3′ (n = 0–5), D strand 5′-GCGTGCAGGTGCCACA-3′ | record |
| Placement of the DNA | Superposed (three points, every third nucleotide) on a DNA arm of the 2PFJ junction that touches the active site; 2PFJ DNA, extra ions and sulfate then deleted; no manual adjustment | builder |
| Bulge | Adenines inserted directly after C8 of the C strand; exact procedure *not recorded* | record / – |
| System set-up | 3A: CHARMM-GUI v3.7 (CHARMM36m). 0A: ISOLDE-relaxed structure, then `gmx pdb2gmx`. 1A, 2A, 4A, 5A: `gmx pdb2gmx` | record |
| Mn²⁺ | Four ions (two per active site). 0A, 1A, 2A, 4A, 5A share one set of four positions added after `pdb2gmx`; 3A has them in the CHARMM-GUI input. Origin of the coordinates *not recorded* | record |
| Scissile bond | C strand A7–C8. Chosen because the C8 phosphate was the DNA phosphate closest to the Mn²⁺ of active site 1 in every starting model; the DNA was not placed to target it; not an experimentally determined cleavage site | record / builder |

Mn–P(C8) distances in the starting models (nm, Mn3 / Mn4):
0A 0.720 / 0.873, 1A 0.714 / 0.879, 2A 0.358 / 0.512, 3A 0.342 / 0.377, 4A 0.332 / 0.507, 5A 0.518 / 0.618.
