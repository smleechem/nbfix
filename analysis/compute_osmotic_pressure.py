#!/usr/bin/env python3
"""
Compute osmotic pressure (bar) from md_nvt_sel.dcd using block averaging.

Defaults:
    zmin = 24 Å, zmax = 72 Å, k_z = 10 kcal/mol/Å^2
Assumptions:
    - md_nvt_sel.dcd contains only the selected ions (e.g., SOD/CLA).
    - Frame stride is 1 ps by default (can be changed with --frame-ps).
    - Box is read from the trajectory (Lx, Ly in Å each frame).

Per-frame recipe:
    Upper = k_z * max(0, z - z_max)      # kcal/mol/Å
    Lower = k_z * max(0, z_min - z)      # kcal/mol/Å
    W     = 0.5 * (sum Upper + sum Lower)
    P(bar) = (W * 69.5e-12) / (Area_m2) / 1.01e5
where Area_m2 = (Lx * Ly * 1e-20)

Block averaging:
    - Split frames into contiguous blocks of size B frames (B = block_ps / frame_ps).
    - For each block j, compute block mean pressure m_j.
    - Cumulative mean after b blocks: mean(m_1..m_b)
    - Cumulative SEM after b blocks: std(m_1..m_b, ddof=1) / sqrt(b)

Output (CSV-like):
    block  time_ps  block_mean_bar  cum_mean_bar  cum_sem_bar

Usage:
    python compute_osmotic_blockavg.py \
        --top md_nvt_sel.pdb \
        --traj md_nvt_sel.dcd \
        --zmin 24 --zmax 72 --kz 10 \
        --frame-ps 1.0 --block-ps 50.0 \
        -out md_nvt_osmo_blockavg.csv
"""

import argparse
import sys
import numpy as np

def parse_args():
    p = argparse.ArgumentParser(description="Block-averaged osmotic pressure from md_nvt_sel.dcd")
    p.add_argument("--top",   required=True, help="Subset topology (e.g., md_nvt_sel.pdb)")
    p.add_argument("--traj",  required=True, help="Subset trajectory (e.g., md_nvt_sel.dcd)")
    p.add_argument("--zmin",  type=float, default=24.0, help="z_min in Å (default: 24.0)")
    p.add_argument("--zmax",  type=float, default=72.0, help="z_max in Å (default: 72.0)")
    p.add_argument("--kz",    type=float, default=10.0, help="k_z in kcal/mol/Å^2 (default: 10.0)")
    p.add_argument("--frame-ps", type=float, default=1.0, help="time per frame in ps (default: 1.0)")
    p.add_argument("--block-ps", type=float, default=100.0, help="block length in ps (default: 10.0)")
    p.add_argument("-out",   default="md_nvt_osmo_blockavg.csv", help="output file (default: md_nvt_osmo_blockavg.csv)")
    return p.parse_args()

def per_frame_pressure_bar(zA, Lx, Ly, zmin, zmax, kz_kcal_per_mol_A2):
    """
    Compute per-frame osmotic pressure (bar) from z positions and box lengths (Å).
    """
    # Upper/Lower magnitudes per ion (kcal/mol/Å)
    upper = kz_kcal_per_mol_A2 * np.maximum(0.0, zA - zmax)
    lower = kz_kcal_per_mol_A2 * np.maximum(0.0, zmin - zA)
    W = 0.5 * (np.sum(upper) + np.sum(lower))  # kcal/mol/Å

    # Area (Å^2 -> m^2)
    area_m2 = float(Lx) * float(Ly) * 1e-20
    if area_m2 <= 0.0:
        return np.nan

    # kcal/mol/Å -> N (× 69.5e-12), then /m^2 -> Pa, then /1.01e5 -> bar
    p_bar = (W * 69.5e-12) / area_m2 / 1.01e5
    return p_bar

def main():
    args = parse_args()

    # Lazy import MDAnalysis
    try:
        import MDAnalysis as mda
    except Exception:
        print("ERROR: MDAnalysis is required. Install with: pip install MDAnalysis", file=sys.stderr)
        sys.exit(1)

    u = mda.Universe(args.top, args.traj)

    zmin  = float(args.zmin)
    zmax  = float(args.zmax)
    kz    = float(args.kz)
    dt_ps = float(args.frame_ps)
    block_ps = float(args.block_ps)

    if dt_ps <= 0.0 or block_ps <= 0.0:
        print("ERROR: --frame-ps and --block-ps must be positive.", file=sys.stderr)
        sys.exit(1)

    B = int(round(block_ps / dt_ps))  # frames per block
    if B < 1:
        B = 1

    # Collect per-frame pressures
    p_frames = []
    for ts in u.trajectory:
        # z positions in Å
        zA = u.atoms.positions[:, 2].astype(float)
        # Box lengths in Å from MDAnalysis (Lx, Ly, Lz, alpha, beta, gamma)
        Lx, Ly, Lz, alpha, beta, gamma = ts.dimensions
        p_bar = per_frame_pressure_bar(zA, Lx, Ly, zmin, zmax, kz)
        p_frames.append(p_bar)

    p_frames = np.asarray(p_frames, dtype=float)
    n_frames = p_frames.size
    if n_frames == 0:
        print("No frames found.", file=sys.stderr)
        sys.exit(1)

    # Drop NaN frames if any
    good = np.isfinite(p_frames)
    if not np.all(good):
        n_bad = np.count_nonzero(~good)
        print(f"[WARN] Dropping {n_bad} invalid frames.", file=sys.stderr)
        p_frames = p_frames[good]
        n_frames = p_frames.size
        if n_frames == 0:
            print("All frames invalid.", file=sys.stderr)
            sys.exit(1)

    # Partition into blocks of length B (truncate tail if incomplete)
    n_blocks = n_frames // B
    if n_blocks == 0:
        # If fewer than B frames, treat all as one block
        n_blocks = 1
        B = n_frames

    p_frames = p_frames[: n_blocks * B]
    blocks = p_frames.reshape(n_blocks, B)

    # Block means
    block_means = np.mean(blocks, axis=1)  # shape (n_blocks,)

    # Cumulative stats across blocks: mean and SEM
    cum_means = np.empty(n_blocks, dtype=float)
    cum_sems  = np.empty(n_blocks, dtype=float)
    for b in range(1, n_blocks + 1):
        mb = block_means[:b]
        cum_means[b-1] = np.mean(mb)
        if b > 1:
            s = np.std(mb, ddof=1)
            cum_sems[b-1] = s / np.sqrt(b)
        else:
            cum_sems[b-1] = 0.0

    # Time at the end of each block (ps)
    times_ps = (np.arange(1, n_blocks + 1) * B * dt_ps).astype(float)

    # Write CSV-like output for direct plt.errorbar use
    with open(args.out, "w") as f:
        f.write("# Osmotic pressure (block averaging)\n")
        f.write(f"# zmin_A,{zmin}\n")
        f.write(f"# zmax_A,{zmax}\n")
        f.write(f"# k_z_kcal_per_mol_per_A2,{kz}\n")
        f.write(f"# frame_ps,{dt_ps}\n")
        f.write(f"# block_ps,{block_ps}\n")
        f.write("block,time_ps,block_mean_bar,cum_mean_bar,cum_sem_bar\n")
        for idx in range(n_blocks):
            f.write(f"{idx+1},{times_ps[idx]:.6f},{block_means[idx]:.8f},{cum_means[idx]:.8f},{cum_sems[idx]:.8f}\n")

    # Final global mean ± SEM from block means
    if n_blocks > 1:
        s_blocks = np.std(block_means, ddof=1)
        sem_global = s_blocks / np.sqrt(n_blocks)
    else:
        sem_global = 0.0
    mean_global = float(np.mean(block_means))

    print("--------------------------------------------------")
    print(f"Frames processed : {n_frames}  (block size: {B} frames, {n_blocks} blocks)")
    print(f"z_min / z_max    : {zmin:.3f} Å / {zmax:.3f} Å")
    print(f"k_z              : {kz:.6f} kcal/mol/Å^2")
    print(f"Global mean      : {mean_global:.6f} bar")
    print(f"Global SEM       : {sem_global:.6f} bar  (block averaged)")
    print(f"Saved            : {args.out}")
    print("--------------------------------------------------")

    # Quick hint for plotting:
    # import pandas as pd; df = pd.read_csv(args.out, comment='#')
    # plt.errorbar(df.time_ps, df.cum_mean_bar, yerr=df.cum_sem_bar, fmt='o-')

if __name__ == "__main__":
    main()

