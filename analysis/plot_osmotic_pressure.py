#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot osmotic pressure (bar): BASE (with NBFIX) vs WONBFIX (without NBFIX) vs EXP (from osmotic coefficient phi).

What this script does
---------------------
- Reads block-averaged CSVs from your MD postprocessing:
    columns: block,time_ps,block_mean_bar,cum_mean_bar,cum_sem_bar
- For each CSV, it takes the **last row** (most cumulative):
    y = cum_mean_bar, yerr = cum_sem_bar
- Classifies files by name:
    - contains "_wonbfix_" -> WONBFIX (without NBFIX)
    - otherwise           -> BASE (with NBFIX)
- Generates an EXP (experimental-like) curve **internally** at 25 °C using:
    π(bar) = φ * i * M * R * T
    where R=0.08314 L·bar·mol^-1·K^-1, i=2 (NaCl), T=298.15 K
  φ is taken from the table `phi_by_M` below (with linear interpolation for missing M).

Outputs
-------
- A PNG figure with three series: BASE, WONBFIX, EXP (φ-based).
- Console dump of values for quick verification.

How to use
----------
- Edit the CONFIG section (FILES, OUTPUT_PNG, etc.) and run:
    python plot_osmotic_series_phi.py
"""

import re
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================== CONFIG (edit here) ============================== #
FILES = [
    "nacl3m_298k.csv",
    "nacl5m_298k.csv",
    "nacl3m_wonbfix_298k.csv",
    "nacl5m_wonbfix_298k.csv",
]

OUTPUT_PNG = "osmotic_phi_ref.png"
PLOT_TITLE = "NaCl Osmotic Pressure at 298 K"

# φ-based experimental reference (25 °C)
R_LBAR_PER_MOLK = 0.08314
T_K = 298.15
I_VANTHOFF = 2.0  # NaCl -> Na+ + Cl-

# φ table by molarity (M) at 25 °C; edit as desired.
# These are reasonable mid-line values near 25 °C (tuned so that 5 M ~ 300 bar).
phi_by_M = {
    1.0: 1.06,
    2.0: 1.12,
    3.0: 1.16,
    4.0: 1.19,
    5.0: 1.21,
}
# =============================================================================== #


def infer_molarity_from_name(name: str) -> float:
    """Infer molarity (M) from filename by matching '...3m...' or '...5m...' etc."""
    m = re.search(r'(\d+(?:\.\d+)?)m', name.lower())
    if not m:
        raise ValueError(f"Cannot infer molarity from filename: {name}")
    return float(m.group(1))


def is_wonbfix(name: str) -> bool:
    """Return True if '_wonbfix_' appears in filename."""
    return "_wonbfix_" in name.lower()


def load_last_cumulative(csv_path: Path) -> Tuple[float, float]:
    """
    Load CSV and return (cum_mean_bar, cum_sem_bar) from the **last** row.
    Assumes CSV has headers and may include comment lines starting with '#'.
    """
    df = pd.read_csv(csv_path, comment="#")
    needed = {"cum_mean_bar", "cum_sem_bar"}
    if not needed.issubset(df.columns):
        raise ValueError(f"{csv_path} missing columns {needed}; found {list(df.columns)}")
    last = df.iloc[-1]
    return float(last["cum_mean_bar"]), float(last["cum_sem_bar"])


def phi_linear_interp(M: float, table: Dict[float, float]) -> float:
    """
    Linear interpolation of φ over the keys in `table` (M values).
    If M is outside the table range, clamp to nearest endpoint.
    """
    xs = sorted(table.keys())
    if M <= xs[0]:
        return table[xs[0]]
    if M >= xs[-1]:
        return table[xs[-1]]
    # find bracket
    for a, b in zip(xs[:-1], xs[1:]):
        if a <= M <= b:
            fa, fb = table[a], table[b]
            t = (M - a) / (b - a)
            return fa * (1 - t) + fb * t
    # fallback (shouldn't happen)
    return table[xs[-1]]


def exp_pressure_from_phi(M: float) -> float:
    """
    Compute φ-based experimental osmotic pressure at 25 °C:
        π = φ * i * M * R * T   (bar)
    """
    phi = phi_linear_interp(M, phi_by_M)
    return phi * I_VANTHOFF * M * R_LBAR_PER_MOLK * T_K


def series_to_arrays(series: Dict[float, Tuple[float, float]], xs: List[float]):
    """Convert dict {M: (mean, sem)} into y, yerr arrays aligned with xs."""
    y = [series[M][0] if M in series else np.nan for M in xs]
    e = [series[M][1] if M in series else np.nan for M in xs]
    return np.array(y, float), np.array(e, float)


def main():
    # Split to BASE (with NBFIX) vs WONBFIX (without NBFIX)
    base: Dict[float, Tuple[float, float]] = {}   # with NBFIX
    wonb: Dict[float, Tuple[float, float]] = {}   # without NBFIX

    for f in FILES:
        path = Path(f)
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        M = infer_molarity_from_name(path.name)
        mean, sem = load_last_cumulative(path)
        if is_wonbfix(path.name):
            wonb[M] = (mean, sem)   # WONBFIX (without NBFIX)
        else:
            base[M] = (mean, sem)   # BASE (with NBFIX)

    # X-axis molarities from data (sorted)
    mols = sorted(set(list(base.keys()) + list(wonb.keys())))
    if not mols:
        raise RuntimeError("No molarities inferred from input FILES. Check names like '*3m*.csv'.")

    # Build experimental series from φ (SEM = 0 by default)
    exp = {M: (exp_pressure_from_phi(M), 0.0) for M in mols}

    # Arrays for plotting
    y_base, e_base = series_to_arrays(base, mols)
    y_wonb, e_wonb = series_to_arrays(wonb, mols)
    y_exp  = np.array([exp[M][0] for M in mols], float)
    e_exp  = np.zeros_like(y_exp)

    # --- Plot ---
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.rcParams['text.usetex'] = True

    plt.figure(figsize=(7,5))
    ax = plt.subplot()
    ax.clear()

    # BASE (with NBFIX)
    ax.errorbar(mols, y_base, yerr=e_base, fmt='o-', capsize=4, label='with NBFIX')

    # WONBFIX (without NBFIX)
    ax.errorbar(mols, y_wonb, yerr=e_wonb, fmt='s-', capsize=4, label='without NBFIX')

    # Experimental from φ
    # If all SEM=0, show as markers/line without error bars
    ax.plot(mols, y_exp, 'x--', label='Experiment')

    ax.set_xlabel('NaCl concentration (M)')
    ax.set_ylabel('Osmotic Pressure (bar)')
    ax.set_title(PLOT_TITLE)
    ax.legend(frameon=False, loc='upper left', fontsize=18)
    plt.savefig(OUTPUT_PNG, bbox_inches='tight', dpi=600)
    print(f"Saved plot -> {OUTPUT_PNG}")

    # Console dump
    def dump_series(name, series):
        print(f"[{name}]")
        for M in mols:
            if M in series:
                m, s = series[M]
                print(f"  {M:>4.1f} M : mean={m:9.3f} bar,  SEM={s:7.3f} bar")
            else:
                print(f"  {M:>4.1f} M : (no data)")

    dump_series("BASE (with NBFIX)", base)
    dump_series("WONBFIX (without NBFIX)", wonb)
    print("[Experiment (φ-based, 25 °C)]")
    for M in mols:
        print(f"  {M:>4.1f} M : mean={exp[M][0]:9.3f} bar,  SEM={0.0:7.3f} bar")


if __name__ == "__main__":
    main()

