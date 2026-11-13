# To do before this repository is complete

Files that exist only on the computing cluster and still have to be exported
(remove user names, group names and absolute paths before adding them):

- [ ] `02_md/inputs/`: `.mdp` files, topology and index files for 0A–5A and the 500 ns run of 3A; equilibration protocol
- [ ] `02_md/`: analysis scripts that produced `analysis_csv/`
- [ ] `03_qmmm/inputs/`: CP2K inputs for the reactant optimisation, restrained scans, dimer search,
      single points, finite-difference Hessian; list of the 52 QM atoms and the 7 link atoms
- [ ] `03_qmmm/results/`: Mulliken charges and spins (`charges_atoms.csv`, `charges_groups.csv`), ELF mid-point check
- [ ] `01_model_building/`: starting PDB files of the six complexes; record of how the bulged DNA was built
- [ ] Large data (trajectories, restart and cube files): deposit elsewhere (e.g. Zenodo or a shared drive) and link here
- [ ] Choose a license (for example MIT for code and CC BY 4.0 for data and figures)
- [ ] Add the Google Drive link for figures and source data to the top-level README
- [ ] Check the Fig. 4a values (read from a printed heat map) against the original contact table
