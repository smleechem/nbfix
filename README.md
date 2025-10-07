# NBFIX Tuning for NaCl Solutions via Osmotic Pressure Method

<img width="3212" height="1503" alt="fig1" src="https://github.com/user-attachments/assets/b06fc485-a221-429b-84b5-27980acc45a3" />


> *Figure.* (a) Slab system used for osmotic-pressure calculations (half-harmonic walls along **z**).
> (b) Osmotic pressure at 298 K for 3 M and 5 M NaCl — experiment (φ–based line), **with NBFIX**, and **without NBFIX**.

---

## Why this repository?

**0) Motivation.** Most legacy ion force fields were fit to properties in the **dilute limit** (e.g., hydration free energies, ion–water RDFs). When pushed to **highly concentrated electrolytes**, they often mis-predict thermodynamic observables — notably **osmotic pressure** — because **ion–ion** interactions become dominant. To address this, it is common to introduce **NBFIX** parameters to **modulate specific ion–ion nonbonded interactions** (beyond Lorentz–Berthelot), while leaving water/ion–water terms unchanged.

**1) Method.** This project implements the **osmotic pressure matching** strategy of Luo & Roux (J. Phys. Chem. Lett. 2010) to tune NBFIX. We simulate slab systems with half-harmonic walls and compute the osmotic pressure from wall forces, then adjust NBFIX to bring simulation in line with experiment.
Reference: [Luo & Roux, *J. Phys. Chem. Lett.* **1**, 183–189 (2010)](https://pubs.acs.org/doi/10.1021/jz900079w).

---

## Repository layout

```
nbfix/
├─ inputfiles/
│  ├─ generate_nacl.sh           # (2) Build 3M / 5M NaCl slab PDBs at 298 K
│  └─ ...                        # FF snippets, φ-tables, etc.
├─ simulation/
│  ├─ run.py                     # (3) OpenMM driver: NPT→NVT, half-harmonic walls, reporters
│  └─ ...                        # helpers, logs
├─ analysis/
│  ├─ compute_osmotic_pressure.py# (4) Read sel.pdb/sel.dcd → block-avg osmotic pressure (+ SEM)
│  ├─ plot_osmotic_pressure.py   # (4) MD vs φ-based experiment (no argparse version included)
│  └─ ...                        # plotting utilities
└─ docs/
   └─ img/
      └─ osmotic_nacl_298K.png   # figure used in README
```

---

## Quick start

### Dependencies

* **Python ≥ 3.9**
* **OpenMM ≥ 8.0**
* **MDAnalysis ≥ 2.6**
* **openmm-mdanalysis-reporter** (to write atom selections to DCD)
* NumPy, pandas, Matplotlib

Install (example):

```bash
pip install openmm MDAnalysis openmm-mdanalysis-reporter numpy pandas matplotlib
```

### 2) Generate inputs (PDBs for 3 M and 5 M)

From `inputfiles/`:

```bash
bash generate_nacl.sh
```

This creates NaCl slab systems (Na⁺/Cl⁻ + water) with **Lx = 2.4 nm, Ly = 2.4 nm, Lz = 9.6 nm** (i.e., 2.4×2.4×9.6 nm³), targeting **3 M** and **5 M** at 298 K.

### 3) Run simulations (OpenMM, half-harmonic walls)

From `simulation/`:

```bash
python run.py
```

What it does:

* NPT equilibration → NVT production at 298 K.
* Adds a **half-harmonic wall** along **z**:
  [
  U_\text{wall}(z) = \tfrac{k_z}{2},\Big(\max(0, z - z_\text{max})^2 + \min(0, z - z_\text{min})^2\Big)
  ]
* Applies the wall **only** to ions (resnames `SOD`, `CLA`).
* Writes reporters:

  * `eq_npt.dcd / .log / .chk` (NPT)
  * `md_nvt.dcd / .log / .chk` (NVT, full system)
  * **`md_nvt_sel.pdb` & `md_nvt_sel.dcd`** — *ion-only* coordinates via **MDAReporter** (requires `openmm-mdanalysis-reporter`).
* During NVT, computes **cumulative osmotic pressure** every 10 ps and logs to `md_nvt_osmo.log` (mean & SEM ready for `plt.errorbar`).

### 4) Analyze osmotic pressure and plot

From `analysis/`:

1. **Compute osmotic pressure from `sel.pdb`/`sel.dcd`** using the same wall-based recipe as the run-time logger, but offline and with **block averaging**:

```bash
python compute_osmotic_pressure.py \
  --top ../simulation/md_nvt_sel.pdb \
  --traj ../simulation/md_nvt_sel.dcd \
  --zmin 24 --zmax 72 --kz 10 \
  --frame-ps 1.0 --block-ps 100.0 \
  --out nacl5m_298k.csv
```

The output CSV contains:

```
block,time_ps,block_mean_bar,cum_mean_bar,cum_sem_bar
...
```

Use `cum_mean_bar ± cum_sem_bar` for publication-quality error bars.

2. **Plot MD vs experiment** using an internal φ-based curve at 25 °C:

* **With NBFIX** CSV: `nacl3m_298k.csv`, `nacl5m_298k.csv`
* **Without NBFIX** CSV: `nacl3m_wonbfix_298k.csv`, `nacl5m_wonbfix_298k.csv`

```bash
python plot_osmotic_pressure.py
```

The script:

* Parses filenames to detect **molarity** and **WONBFIX** vs **BASE** (with NBFIX).
* Reads the **last** cumulative mean/SEM from each CSV.
* Builds an **experimental** curve from a φ(M) table at 25 °C:
  [
  \pi = \phi,i,M,R,T \quad\text{with}\quad i=2,; R=0.08314,\text{L·bar·mol}^{-1}\text{·K}^{-1},; T=298.15,\text{K}
  ]
* Produces a figure similar to the one above.

> You can replace the default φ(M) with your preferred dataset; linear interpolation is used between tabulated M.

---

## How NBFIX enters

* **NBFIX** overrides default mixing rules for specific ion–ion pairs (e.g., Na⁺–Cl⁻ σ/ε in Lennard–Jones) to **soften or stiffen** the short-range part of the effective interaction while leaving ion–water parameters intact.
* The workflow here adjusts NBFIX so that the **simulated osmotic pressure** at 298 K for **3 M and 5 M** tracks experiment. Once fitted, you can cross-check at **1 M, 2 M, 4 M** and other salts.

---

## Practical notes

* **Block averaging.** For robust SEM, choose the **block length ≥ the autocorrelation time** of your observable. For osmotic pressure in concentrated salt, a safe starting point is **100 ps** blocks and then scan 50–200 ps to find a plateau in SEM.
* **Units.** The wall constant `k_z` is provided in **kcal·mol⁻¹·Å⁻²**, `z_min/max` in **Å**, and the area conversion uses **Å² → m² (×1e−20)**. The factor **69.5×10⁻¹²** converts **kcal·mol⁻¹·Å⁻¹ → N** for the force surrogate in the Luo–Roux recipe; dividing by area and **1.01×10⁵** yields **bar**.
* **Selections.** Only ions (`resname SOD CLA`) feel the wall and are written to `md_nvt_sel.dcd`. Water and other species remain unaffected by the external wall potential.
* **Reproducibility.** The NVT code can re-compute cumulative osmotic pressure from the saved ion-only DCD (post-run) to guard against runtime logging artifacts.

---

## End-to-end workflow (TL;DR)

1. **Build** inputs (3 M & 5 M)
   `bash inputfiles/generate_nacl.sh`
2. **Simulate** with current NBFIX guess
   `python simulation/run.py`
3. **Analyze** ion-only trajectory
   `python analysis/compute_osmotic_blockavg.py ...`
4. **Compare** to experiment (φ-based line)
   `python analysis/plot_osmotic_series_phi.py`
5. **Adjust NBFIX** → repeat steps 2–4 until MD ≈ experiment.

---

## Reference

* **Osmotic-pressure method**:
  Luo, Y.; Roux, B. *J. Phys. Chem. Lett.* **2010**, 1, 183–189.
  DOI: 10.1021/jz900079w

Please cite the above if you use this method, and cite this repository if the scripts help your work.

---

## License

MIT (see `LICENSE`).

---

## Acknowledgements

Thanks to the OpenMM and MDAnalysis communities. The MD selection writer uses **openmm-mdanalysis-reporter** to dump `sel.dcd` efficiently.

---

### Maintainer

* **S.M. Lee** — issues and PRs welcome at [https://github.com/smleechem/nbfix/](https://github.com/smleechem/nbfix/).
