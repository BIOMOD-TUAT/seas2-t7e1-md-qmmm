# Figure scripts

Figures, editable PowerPoint files and source-data workbooks are on the team Google Drive (link in the top README).
This folder has what was used to make them.

## scripts/

Python 3 with numpy, scipy, matplotlib, Pillow, python-pptx, openpyxl. Run inside `scripts/`.

| script | makes | reads |
|---|---|---|
| `md_metal_site_and_stability.py` | Mn–P distances, Mn coordination table, RMSD / interface plots (MD) | `02_md/analysis_csv/Fig07_frames_*.csv`, `Fig03_timeseries_*.csv` |
| `plot_contacts.py` | contact plot (inserted A vs. T7E1 residues) | values typed in the script, `img/contacts_500ns.png` |
| `plot_md_qmmm_figures.py` | stability, contacts, metal site, attacking water, QM/MM model, structures, spin density + ELF panels | `02_md/analysis_csv/`, `elf.pkl`, pictures in `img/` |
| `elf_plane_section.py` | ELF on the plane through O3′, P and the water oxygen → `elf.pkl` | `data/*_elf_web.cube`, `data/*_env.pdb` |
| `plot_energy_profile.py` | energy and distances along the path, transition-state table | `03_qmmm/results/energy_table.csv`, `path_distances.csv` |
| `draw_reaction_scheme.py`, `scheme_base.py` | reaction scheme (PowerPoint shapes) | – |

`img/` is not in the repository: these are the pictures saved from ChimeraX (`complex_wide.png`, `qm_region.png`,
`structure_R/TS/P.png`, `spin_crop.png`) and the 500 ns contact plot (`contacts_500ns.png`).
Scripts that write PowerPoint files convert them to PNG with LibreOffice (`soffice.py` is a small wrapper around
`soffice --headless --convert-to pdf`; replace it with your own call).

## chimerax/

Command files for the 3D pictures. Start ChimeraX in a folder that contains `struct/` (the PDB files in
`03_qmmm/structures/` plus `complex_R.pdb`) and `cube/` (spin density and ELF cube files), then `open <file>.cxc`.

- `qm_region_overview.cxc`, `qm_region_only.cxc`, `qm_region_closeup.cxc`: QM region in the complex
- `structure_R.cxc`, `structure_TS.cxc`, `structure_P.cxc`: reactant, transition state, after cleavage
- `spin_density_*.cxc`, `elf_isosurface_*.cxc`: spin density and ELF isosurfaces

Colours: QM atoms ball-and-stick by element, Mn purple, nucleophilic water oxygen green, transferring protons magenta,
MM region dim ribbons and thin grey sticks. PDB B-factor column: 0 MM, 1 QM, 2 reaction core, 3 Mn.

## vmd/

`movie_t7e1.tcl`: reaction movie from a multi-frame PDB (`rqB_movie_traj.pdb`, not in the repository).
