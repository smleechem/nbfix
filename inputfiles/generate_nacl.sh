#!/usr/bin/env bash
set -euo pipefail

# ========== User-configurable ==========
# Box in Å
LX_A=48.0
LY_A=48.0
# Regions (z in Å): [0,24] slab | (24,72] NaCl solution | (72,92] slab
Z1_MIN=0.0;  Z1_MAX=24.0     # bottom slab
Z2_MIN=24.0; Z2_MAX=72.0     # NaCl solution region (48 Å thick)
Z3_MIN=72.0; Z3_MAX=96.0     # top slab

# NaCl concentrations to generate (M)
CONCS="1 3 5"

# Water density and molar mass (adjust if needed)
WATER_DENSITY=0.997          # g/mL @ 25 °C
M_WATER=18.01528             # g/mol

# Templates
TEMPLATES_DIR="../templates"  # must contain: tip3p.pdb, na.pdb, cl.pdb

# Packmol options
TOL=2.0
FILETYPE="pdb"

# Output directory
OUTDIR="./nacl"
# ======================================

# --- Preflight ---
command -v packmol >/dev/null 2>&1 || { echo "Error: packmol not found in PATH"; exit 1; }
mkdir -p "${OUTDIR}"

# 1) Compute counts (pure-water slabs + NaCl solution region per concentration)
#    - Slab water counts from density
#    - NaCl pairs and solution water from exact molarity & density (48×48×48 Å region)
#    - Write Packmol inputs for each target concentration
/usr/bin/env python3 - "$LX_A" "$LY_A" "$Z1_MIN" "$Z1_MAX" "$Z2_MIN" "$Z2_MAX" "$Z3_MIN" "$Z3_MAX" \
                      "$WATER_DENSITY" "$M_WATER" "$TEMPLATES_DIR" "$TOL" "$FILETYPE" "$OUTDIR" $CONCS <<'PY'
import sys
from pathlib import Path
from typing import Dict, List

# ---- Constants ----
AVOGADRO = 6.022_140_76e23        # 1/mol
M_NACL   = 58.442769              # g/mol

# ---- Args from bash ----
(LX, LY, Z1_MIN, Z1_MAX, Z2_MIN, Z2_MAX, Z3_MIN, Z3_MAX,
 RHO_W, M_WATER, TDIR, TOL, FTYPE, OUTDIR, *CONCS) = sys.argv[1:]

LX = float(LX); LY = float(LY)
Z1_MIN = float(Z1_MIN); Z1_MAX = float(Z1_MAX)
Z2_MIN = float(Z2_MIN); Z2_MAX = float(Z2_MAX)
Z3_MIN = float(Z3_MIN); Z3_MAX = float(Z3_MAX)
RHO_W = float(RHO_W); M_WATER = float(M_WATER)
TOL = float(TOL); OUTDIR = Path(OUTDIR)
CONCS = [float(c) for c in CONCS]
TDIR = Path(TDIR)

# Total box (for CRYST1)
LZ_TOTAL = Z3_MAX - Z1_MIN

# Density table for NaCl(aq) @ 25 °C (g/mL). Replace with your reference if desired.
DENSITY_NACL_25C: Dict[float, float] = {
    1.0: 1.035,
    3.0: 1.106,
    5.0: 1.173,
}

def vol_L(ax: float, ay: float, az: float) -> float:
    # Volume in liters; 1 Å^3 = 1e-27 L
    return ax * ay * az * 1e-27

def water_count_from_density(rho_g_mL: float, V_L: float) -> int:
    # N = rho * 1000 * V / M * N_A  (rounded)
    moles = (rho_g_mL * 1000.0 * V_L) / M_WATER
    return int(round(moles * AVOGADRO))

def nacl_pairs_from_molarity(c_M: float, V_L: float) -> int:
    return int(round(c_M * V_L * AVOGADRO))

def nacl_solution_water_for_density(Npair: int, rho_g_mL: float, V_L: float) -> int:
    # mass_target = rho * 1000 * V
    mass_target_g = rho_g_mL * 1000.0 * V_L
    # Nwater = (mass_target*N_A - Npair*M_NACL) / M_WATER
    Nw = int(round((mass_target_g * AVOGADRO - Npair * M_NACL) / M_WATER))
    if Nw < 0:
        raise SystemExit("Negative water count for NaCl solution region; adjust inputs.")
    return Nw

# --- Volumes ---
V1 = vol_L(LX, LY, (Z1_MAX - Z1_MIN))  # bottom slab
V2 = vol_L(LX, LY, (Z2_MAX - Z2_MIN))  # NaCl solution region (48 Å)
V3 = vol_L(LX, LY, (Z3_MAX - Z3_MIN))  # top slab

# --- Pure-water counts for slabs ---
Nw_slab_bottom = water_count_from_density(RHO_W, V1)
Nw_slab_top    = water_count_from_density(RHO_W, V3)

# --- Packmol template (layered) ---
PKT = """tolerance {tolerance}

filetype {filetype}
output {out_pdb}

# Bottom water slab (z = {z1min} to {z1max})
structure {tip3p}
  number {nwater_bottom}
  inside box 0.0 0.0 {z1min} {lx} {ly} {z1max}
end structure

# NaCl solution region (z = {z2min} to {z2max})
structure {tip3p}
  number {nwater_sol}
  inside box 0.0 0.0 {z2min} {lx} {ly} {z2max}
end structure

structure {na}
  number {npair}
  inside box 0.0 0.0 {z2min} {lx} {ly} {z2max}
end structure

structure {cl}
  number {npair}
  inside box 0.0 0.0 {z2min} {lx} {ly} {z2max}
end structure

# Top water slab (z = {z3min} to {z3max})
structure {tip3p}
  number {nwater_top}
  inside box 0.0 0.0 {z3min} {lx} {ly} {z3max}
end structure
"""

# --- Generate one Packmol input per concentration ---
OUTDIR.mkdir(parents=True, exist_ok=True)
summary_rows: List[dict] = []

for c in CONCS:
    if c not in DENSITY_NACL_25C:
        raise SystemExit(f"No density available for {c} M in DENSITY_NACL_25C.")
    rho = DENSITY_NACL_25C[c]

    # NaCl solution counts for the middle region (48×48×48 Å)
    npair = nacl_pairs_from_molarity(c, V2)
    nwater_sol = nacl_solution_water_for_density(npair, rho, V2)

    # Compose Packmol input
    out_pdb = f"slab_nacl{int(c)}m.pdb"
    text = PKT.format(
        tolerance=TOL, filetype=FTYPE, out_pdb=out_pdb,
        tip3p=str(TDIR / "tip3p.pdb"),
        na=str(TDIR / "na.pdb"),
        cl=str(TDIR / "cl.pdb"),
        lx=LX, ly=LY,
        z1min=Z1_MIN, z1max=Z1_MAX,
        z2min=Z2_MIN, z2max=Z2_MAX,
        z3min=Z3_MIN, z3max=Z3_MAX,
        nwater_bottom=Nw_slab_bottom,
        nwater_sol=nwater_sol,
        nwater_top=Nw_slab_top,
        npair=npair,
    )
    (OUTDIR / f"slab_nacl{int(c)}m.inp").write_text(text)

    # Save summary
    summary_rows.append({
        "M": c,
        "NaCl_pairs": npair,
        "Na+": npair, "Cl-": npair,
        "Water_bottom": Nw_slab_bottom,
        "Water_solution": nwater_sol,
        "Water_top": Nw_slab_top,
        "Box_A": f"{LX} x {LY} x {LZ_TOTAL}",
        "Regions": f"[{Z1_MIN},{Z1_MAX}] | [{Z2_MIN},{Z2_MAX}] | [{Z3_MIN},{Z3_MAX}]",
    })

# Write CSV summary
import pandas as pd
pd.DataFrame(summary_rows).sort_values("M").to_csv(OUTDIR / "slab_nacl_counts_summary.csv", index=False)

# Write a run script
run = """#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
for M in """ + " ".join(str(int(c)) for c in CONCS) + """; do
  echo "Running Packmol for slab_nacl${M}m ..."
  packmol < "slab_nacl${M}m.inp"
done
echo "Done."
"""
(OUTDIR / "run_packmol_layers.sh").write_text(run)

# Emit final CRYST1 string for the full box
print(f"CRYST1_LINE={LX:9.3f}{LY:9.3f}{(Z3_MAX-Z1_MIN):9.3f}")
PY

echo "[Info] Generated Packmol inputs in: ${OUTDIR}"
echo "[Info] Summary CSV: ${OUTDIR}/slab_nacl_counts_summary.csv"

# 2) Run Packmol for all layered systems
chmod +x "${OUTDIR}/run_packmol_layers.sh"
( cd "${OUTDIR}" && ./run_packmol_layers.sh )

# 3) Post-process: insert CRYST1  (48, 48, 92) AFTER REMARK, BEFORE first ATOM/HETATM
#    Use exactly one CRYST1 line; replace if present, insert otherwise.
CRYST1_LINE=$(printf 'CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1' \
              "${LX_A}" "${LY_A}" "$(awk -v z1="${Z1_MIN}" -v z3="${Z3_MAX}" 'BEGIN{printf "%.3f", z3 - z1}')")

(
  cd "${OUTDIR}"
  shopt -s nullglob
  for f in slab_nacl*m.pdb; do
    awk -v line="$CRYST1_LINE" '
      BEGIN { inserted = 0 }
      /^CRYST1/ { if (!inserted) { print line; inserted = 1 } ; next }
      /^(ATOM  |HETATM)/ { if (!inserted) { print line; inserted = 1 } ; print; next }
      { print }
      END { if (!inserted) print line }
    ' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
    echo "[Info] Inserted/updated CRYST1 in: $f"
  done
)

echo "[Done] Outputs in: ${OUTDIR}"

