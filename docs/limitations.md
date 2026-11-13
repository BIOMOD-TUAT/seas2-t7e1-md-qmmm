# Limitations

## MD
- One trajectory per condition; no independent replicas, no statistical test.
- The models were not built in the same way: 3A through CHARMM-GUI, the others through `gmx pdb2gmx`;
  equilibration histories differ.
- The DNA was built in ChimeraX and placed by superposition on a DNA arm of the Holliday-junction complex 2PFJ.
  It is a model, not an experimental structure of T7E1 with bulged DNA.
- The Mn–phosphate distances differ between conditions already in the starting models
  (0A/1A about 0.72 nm, 2A–4A 0.33–0.36 nm, 5A 0.52 nm) and largely stay that way during 50 ns.
- Mn²⁺ is a non-bonded model. The Mn–Mn distance is about 0.50–0.52 nm (crystal: 0.352 nm), and one Mn²⁺
  binds Glu20 and Glu37, which is not seen in the crystal structure.
- Thresholds (contact 0.45 nm; water–P 0.24–0.40 nm and angle ≥ 140°; Mn–O 0.26 nm) were set in this work.
  The Mn–Mn criterion (level 2c) is strongly threshold dependent and is not used for conclusions.
- In the 500 ns run of 3A the attack-competent, Mn-bound water is absent between 30 and 120 ns and present
  54–94 % of the time after 150 ns. The 50 ns value of 3A is not representative.
- Fig. 4a values were read from a printed heat map of the original analysis and have not been checked against
  the original table.

## QM/MM
- Only the 3A model was studied.
- Barrier ΔE‡ = +57.6 kcal/mol is far above what is expected for an enzyme.
- No stable product: the structure after cleavage is 0.05 kcal/mol below the transition state;
  the energy is flat for P–O3′ = 2.3–2.5 Å.
- Largest force at the transition state 1.14 × 10⁻³ Ha/Bohr (target 1.0 × 10⁻³).
- ΔG‡ uses a harmonic partial Hessian of 20 atoms, not free-energy sampling.
- Grid ("egg-box") noise: up to about 1 kcal/mol per structure, up to 0.2 kcal/mol on the barrier.
- Mn–Mn distance in the QM/MM model is 4.06 Å (crystal 3.52 Å); the phosphate binds one Mn²⁺ only.
- Both Mn²⁺ are high spin with parallel spins (multiplicity 11); the antiparallel state was not examined.
- The nucleophile is treated as water; hydroxide and other protonation states were not examined.
- ELF values are for α-spin electrons. Charges are Mulliken charges of the QM region only.
- The QM/MM reactant is not the same atomic arrangement as any MD frame, because the metal coordination was rebuilt.
