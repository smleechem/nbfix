from pathlib import Path
import numpy as np

try:
    # OpenMM 8.x style
    from openmm import (
        Platform, NonbondedForce, CustomNonbondedForce, CustomExternalForce,
        MonteCarloAnisotropicBarostat, LangevinIntegrator
    )
    from openmm.app import (
        Element, PDBFile, ForceField, Simulation,
        DCDReporter, StateDataReporter, CheckpointReporter, PME, HAngles, Topology
    )
except Exception:
    # Fallback for older OpenMM installs
    from simtk.openmm import (
        Platform, NonbondedForce, CustomNonbondedForce, CustomExternalForce,
        MonteCarloAnisotropicBarostat, LangevinIntegrator
    )
    from simtk.openmm.app import (
        Element, PDBFile, ForceField, Simulation,
        DCDReporter, StateDataReporter, CheckpointReporter, PME, HAngles, Topology
    )

# Try importing MDAReporter (pip install openmm-mdanalysis-reporter)
try:
    from mdareporter import MDAReporter  # MDAnalysis-based reporter
    HAS_MDA = True
except Exception:
    HAS_MDA = False

from simtk.unit import (
    atmosphere, kelvin, nanometer, angstroms, kilocalories_per_mole,
    microseconds, nanoseconds, picoseconds, femtoseconds, kilojoule_per_mole
)

# =========================== User-configurable I/O =========================== #
STR_INP = Path("../../inputfiles/nacl")
STR_FF  = Path("../../forcefields")
STR_OUT = Path("./")

PDB_IN   = STR_INP / "slab_nacl3m.pdb"
FF_XML   = STR_FF  / "water.xml"

# NPT phase outputs
BEFOREMIN_PDB = STR_OUT / "beforemin.pdb"
MIN_PDB       = STR_OUT / "min.pdb"
EQ_PDB        = STR_OUT / "eq_npt.pdb"
DCD_FILE      = STR_OUT / "eq_npt.dcd"
LOG_FILE      = STR_OUT / "eq_npt.log"
CHK_FILE      = STR_OUT / "eq_npt.chk"
ENER_LOG_TXT  = STR_OUT / "eq_npt_ener.log"

# NPT equilibration schedule: blocks of 1000 steps totaling ~1 ns
EQ_TIME         = 1 * nanoseconds

# Production length (NVT total time)
PROD_TIME               = 20 * nanoseconds

# MD block runtime for NVT (e.g., 1 ps). Steps per block are derived from timestep.
MD_BLOCK_TIME           = 1 * picoseconds

# Cumulative osmotic pressure reporting cadence (e.g., every 100 ps)
REPORT_CADENCE_TIME     = 1 * nanoseconds

# Optional: recompute osmotic pressure after the run by reading the subset DCD
RECOMPUTE_OSMOTIC_FROM_DCD = False  # set True to post-process md_nvt_sel.dcd

# Checkpointing cadence
CHKPOINT_STEPS          = 10000

# NVT phase outputs
NVT_DCD_FILE            = STR_OUT / "md_nvt.dcd"          # full-system DCD
NVT_SEL_DCD_FILE        = STR_OUT / "md_nvt_sel.dcd"      # subset DCD via MDAReporter (1 ps)
NVT_SEL_PDB             = STR_OUT / "md_nvt_sel.pdb"      # subset PDB (saved before subset DCD)
NVT_LOG_FILE            = STR_OUT / "md_nvt.log"
NVT_CHK_FILE            = STR_OUT / "md_nvt.chk"
NVT_ENER_LOG            = STR_OUT / "md_nvt_ener.log"
NVT_OSMO_LOG            = STR_OUT / "md_nvt_osmo.log"
NVT_PDB                 = STR_OUT / "md_nvt.pdb"

# ============================= Physical parameters =========================== #
pressure    = 1 * atmosphere
temperature = 298 * kelvin
cutoff      = 1.2 * nanometer
freq        = 1 / picoseconds           # Langevin friction (ps^-1)
timestep    = 1 * femtoseconds

# Wall potential globals
kz   = 10 * kilocalories_per_mole / angstroms**2
zmax = 72 * angstroms
zmin = 24 * angstroms

# Target atoms: residue names SOD and CLA
allowed_resnames = {"SOD", "CLA"}

def get_force(system, cls):
    """Return the first force of a given class from the system, or raise if not found."""
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if isinstance(f, cls):
            return f
    raise RuntimeError(f"No force of type {cls.__name__} found in the system.")


def build_subset_topology(src_top: Topology, atom_indices):
    """
    Build a new Topology containing only atoms in `atom_indices` (order preserved),
    copying chains/residues as needed. Returns (new_topology, ordered_old_indices).
    """
    atom_indices = set(int(i) for i in atom_indices)
    new_top = Topology()
    chain_map = {}
    residue_map = {}
    atom_map = {}

    for chain in src_top.chains():
        new_chain = None
        for res in chain.residues():
            keep_atoms = [a for a in res.atoms() if a.index in atom_indices]
            if not keep_atoms:
                continue
            if new_chain is None:
                new_chain = new_top.addChain(id=chain.id)
                chain_map[chain] = new_chain
            new_res = new_top.addResidue(res.name, chain_map[chain], id=res.id)
            residue_map[res] = new_res
            for a in keep_atoms:
                new_atom = new_top.addAtom(a.name, a.element, new_res, id=str(a.id) if hasattr(a, "id") else None)
                atom_map[a.index] = new_atom

    for bond in src_top.bonds():
        i = bond[0].index
        j = bond[1].index
        if i in atom_map and j in atom_map:
            new_top.addBond(atom_map[i], atom_map[j])

    ordered = [a.index for a in src_top.atoms() if a.index in atom_map]
    return new_top, ordered


def main():
    # ------------------------------- Load input ------------------------------ #
    pdb = PDBFile(str(PDB_IN))
    ff  = ForceField(str(FF_XML))

    # ------------------------------ Build system ---------------------------- #
    system = ff.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=cutoff,
        constraints=HAngles,
        ignoreExternalBonds=True,
        rigidWater=True,
    )

    # Langevin integrator for NPT equilibration
    integ_eq = LangevinIntegrator(temperature, freq, timestep)

    # -------------------------- Configure nonbonded -------------------------- #
    # Electrostatics (PME)
    electrostatics_force = get_force(system, NonbondedForce)
    electrostatics_force.setNonbondedMethod(NonbondedForce.PME)
    electrostatics_force.setCutoffDistance(cutoff)

    # Sterics (assume first CustomNonbondedForce holds LJ)
    sterics_force = get_force(system, CustomNonbondedForce)
    sterics_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    sterics_force.setCutoffDistance(cutoff)

    # Attach an anisotropic barostat (NPT)
    barostat = MonteCarloAnisotropicBarostat([0, 0, pressure], temperature, False, False, True, 100)
    barostat_index = system.addForce(barostat)

    # --------------------------- External Z-wall ----------------------------- #
    z_wall = CustomExternalForce(
        "(k_z/2.0) * ( (max(0, z - z_max))^2 + (min(0, z - z_min))^2 )"
    )
    z_wall.addGlobalParameter("k_z",   kz)
    z_wall.addGlobalParameter("z_max", zmax)
    z_wall.addGlobalParameter("z_min", zmin)

    # Apply the external wall only to atoms in residues SOD or CLA; collect selection indices
    sel_indices = []
    for atom in pdb.topology.atoms():
        rname = atom.residue.name.upper()
        if rname in allowed_resnames:
            z_wall.addParticle(atom.index)
            sel_indices.append(atom.index)

    system.addForce(z_wall)

    # Assign each force to a unique group so we can report component energies
    for i in range(system.getNumForces()):
        system.getForce(i).setForceGroup(i)

    # ----------------------------- Make simulation -------------------------- #
    platform   = Platform.getPlatformByName("CUDA")
    properties = {"CudaPrecision": "mixed", "DeviceIndex": "0"}

    sim = Simulation(pdb.topology, system, integ_eq, platform, properties)
    sim.context.setPositions(pdb.positions)

    # Save pre-minimization coordinates
    state = sim.context.getState(getEnergy=True, getForces=True, getPositions=True)
    with open(BEFOREMIN_PDB, "w") as f:
        PDBFile.writeFile(sim.topology, state.getPositions(), f)

    # ------------------------------- Minimize ------------------------------- #
    print("Wrote initial positions. Minimizing...")
    sim.minimizeEnergy(maxIterations=10_000_000)
    print("Minimization finished!")

    # Energy summary post-minimization
    state = sim.context.getState(getEnergy=True, getVelocities=True, getPositions=True)
    print("KE:", state.getKineticEnergy())
    print("PE:", state.getPotentialEnergy())
    for i in range(system.getNumForces()):
        e_i = sim.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()
        print(f"Group {i}: {e_i}")

    # Save minimized coordinates
    with open(MIN_PDB, "w") as f:
        PDBFile.writeFile(sim.topology, state.getPositions(), f)

    # ---------------------------- NPT equilibration -------------------------- #
    print("Equilibrating (NPT)...")
    sim.context.setVelocitiesToTemperature(temperature)

    STEPS_PER_BLOCK = 1000
    EQ_BLOCKS       = int((EQ_TIME / (STEPS_PER_BLOCK * timestep)) + 0.5)  # round to nearest int

    # NPT reporters
    dcd_reporter  = DCDReporter(str(DCD_FILE), STEPS_PER_BLOCK)
    data_reporter = StateDataReporter(
        str(LOG_FILE), STEPS_PER_BLOCK,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, density=True, speed=True
    )
    chk_reporter  = CheckpointReporter(str(CHK_FILE), 50_000)

    sim.reporters.extend([dcd_reporter, data_reporter, chk_reporter])

    # Prime reporters with the current state
    dcd_reporter.report(sim, state)
    data_reporter.report(sim, state)
    print("Simulating NPT...")

    # Per-force-group energy log for NPT
    with open(ENER_LOG_TXT, "w") as flog:
        flog.write("# Energy log file (NPT)\n# x1 : block index\n")
        for j in range(system.getNumForces()):
            flog.write(f"# x{j+2} : group {j} {type(system.getForce(j))} (kJ/mol)\n")

        for i in range(1, EQ_BLOCKS + 1):
            sim.step(STEPS_PER_BLOCK)
            state = sim.context.getState(getEnergy=True, getPositions=True, getVelocities=True)
            print(f"NPT block {i}/{EQ_BLOCKS} | KE {state.getKineticEnergy()} | PE {state.getPotentialEnergy()}")

            flog.write(f"{i}")
            for j in range(system.getNumForces()):
                e_comp = sim.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()
                flog.write(f"  {e_comp.value_in_unit(kilojoule_per_mole):.6f}")
            flog.write("\n")
            flog.flush()

    print("NPT equilibration done!")

    # ----------------------------- Final NPT outputs ------------------------- #
    state_npt = sim.context.getState(getEnergy=True, getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    sim.topology.setPeriodicBoxVectors(state_npt.getPeriodicBoxVectors())
    with open(EQ_PDB, "w") as f:
        PDBFile.writeFile(sim.topology, state_npt.getPositions(), f)
    print(f"Wrote NPT output coordinate file: {EQ_PDB}")

    # =========================== NVT production (100 ns) ===================== #
    # Take final NPT positions/velocities/box to start NVT seamlessly
    pos_npt = state_npt.getPositions()
    vel_npt = state_npt.getVelocities()
    box_npt = state_npt.getPeriodicBoxVectors()

    # Remove the barostat to switch to NVT
    system.removeForce(barostat_index)

    # NVT integrator (same T/step)
    integ_md = LangevinIntegrator(temperature, freq, timestep)

    # New Simulation for NVT with the same system/platform
    simmd = Simulation(pdb.topology, system, integ_md, platform, properties)

    # Reset the MD clock (optional)
    simmd.context.setTime(0.0)

    # Initialize periodic box, positions, and velocities from NPT
    simmd.context.setPeriodicBoxVectors(box_npt[0], box_npt[1], box_npt[2])
    simmd.context.setPositions(pos_npt)
    simmd.context.setVelocities(vel_npt)

    # Steps per block and report cadence (in blocks)
    MD_BLOCK_STEPS = max(1, int((MD_BLOCK_TIME / timestep) + 0.5))
    REPORT_CADENCE_BLOCKS = max(1, int((REPORT_CADENCE_TIME / MD_BLOCK_TIME) + 0.5))

    ## Scalar parameters for osmotic calculation (in requested units)
    k_wall_val = kz.value_in_unit(kilocalories_per_mole / angstroms**2)  # kcal/mol/Å^2
    z_max_val  = zmax.value_in_unit(angstroms)                            # Å
    z_min_val  = zmin.value_in_unit(angstroms)                            # Å

    # Create osmotic log BEFORE NVT starts
    with open(NVT_OSMO_LOG, "w") as fosmo:
        fosmo.write("# Osmotic pressure log (cumulative; every {} blocks = {})\n"
                    .format(REPORT_CADENCE_BLOCKS, REPORT_CADENCE_TIME))
        fosmo.write("# Columns: block_index  time_ps  pressure_bar\n")

    # Cumulative counters (since t=0)
    frames_accum      = 0                         # number of reported blocks accumulated
    upper_sum_total   = 0.0                       # sum over frames of sum_i k*(z_i - z_max) for z_i>z_max
    lower_sum_total   = 0.0                       # sum over frames of sum_i k*(z_min - z_i) for z_i<z_min
    # (units of these sums: kcal/mol/Å)

    # ------------------------------- NVT reporters --------------------------- #
    # Full-system DCD every block
    simmd.reporters.append(DCDReporter(str(NVT_DCD_FILE), MD_BLOCK_STEPS, enforcePeriodicBox=True))

    # StateData every block
    simmd.reporters.append(StateDataReporter(
        str(NVT_LOG_FILE), MD_BLOCK_STEPS,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, density=True, speed=True
    ))

    # Checkpoints
    simmd.reporters.append(CheckpointReporter(str(NVT_CHK_FILE), CHKPOINT_STEPS))

    # Prime reporters with the current state (frame 0)
    state0 = simmd.context.getState(getEnergy=True, getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    simmd.reporters[0].report(simmd, state0)  # DCD
    simmd.reporters[1].report(simmd, state0)  # StateData

    # Selected-atom PDB (SOD/CLA) at the start of NVT
    sub_top, ordered_sel = build_subset_topology(simmd.topology, sel_indices)

    # Convert positions to a numeric (N,3) ndarray in a length unit, then slice
    pos0_nm = np.array(state0.getPositions(asNumpy=True).value_in_unit(nanometer))  # (N,3) float
    sub_pos_nm = pos0_nm[ordered_sel, :]                                            # (n_sel,3) float

    # Reattach units so PDBFile.writeFile receives a Quantity array
    sub_pos = sub_pos_nm * nanometer

    with open(NVT_SEL_PDB, "w") as f:
        PDBFile.writeFile(sub_top, sub_pos, f)

    print(f"Wrote selected-atom PDB: {NVT_SEL_PDB}")

    # MDAnalysis-based subset DCD every 1 ps (selection = SOD/CLA)
    if HAS_MDA:
        selection_str = " or ".join(f"resname {r}" for r in sorted(allowed_resnames))
        ONE_PS_STEPS = max(1, int((1 * picoseconds) / timestep + 0.5))
        simmd.reporters.append(
            MDAReporter(str(NVT_SEL_DCD_FILE), ONE_PS_STEPS, selection=selection_str)
        )
    else:
        print("[WARN] MDAReporter not available. Install with: pip install openmm-mdanalysis-reporter")

    print("Starting NVT production...")

    # ---------------- Run NVT with per-block progress + cumulative osmotic calc --------------- #

    prod_steps = int((PROD_TIME / timestep) + 0.5)
    md_blocks  = prod_steps // MD_BLOCK_STEPS  # no remainder section

    # Cumulative containers since t=0
    accum_posz_all = []  # list of 1D arrays (Å) appended each block for selected atoms
    osmo_values    = []  # cumulative osmotic values (bar) at report times

    with open(NVT_ENER_LOG, "w") as fener, open(NVT_OSMO_LOG, "a") as fosmo:
        # Energy log header
        fener.write("# Energy log file (NVT)\n")
        fener.write("# x1 : block index (each block = MD_BLOCK_TIME)\n")
        for j in range(system.getNumForces()):
            f = system.getForce(j)
            fener.write(f"# x{j+2} : group {j} {type(f)} (kJ/mol)\n")

        for i in range(1, md_blocks + 1):
            simmd.step(MD_BLOCK_STEPS)

            # State for energies + positions + box
            state_i = simmd.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
            print(f"NVT block {i}/{md_blocks} | KE {state_i.getKineticEnergy()} | PE {state_i.getPotentialEnergy()}")

            # ---- component energies ----
            fener.write(f"{i}")
            for j in range(system.getNumForces()):
                e_comp = simmd.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy()
                fener.write(f"  {e_comp.value_in_unit(kilojoule_per_mole):.6f}")
            fener.write("\n")
            fener.flush()

            # ---- accumulate selected-atom z positions (Å) cumulatively ----
            posA = state_i.getPositions(asNumpy=True).value_in_unit(angstroms)   # (N,3) Å
            zsel = posA[ordered_sel, 2]                                          # (N_sel,) Å

            # Upper-wall contribution: k * max(0, z - z_max)
            # Lower-wall contribution: k * max(0, z_min - z)
            # (both return non-negative "magnitudes" per ion per frame)
            upper_sum_total += float(np.sum(k_wall_val * np.maximum(0.0, zsel - z_max_val)))
            lower_sum_total += float(np.sum(k_wall_val * np.maximum(0.0, z_min_val - zsel)))
            frames_accum    += 1

            # ---- every REPORT_CADENCE_BLOCKS, compute cumulative osmotic pressure ----
            if (i % REPORT_CADENCE_BLOCKS) == 0:
                # Per-frame averaged magnitudes (kcal/mol/Å)
                Fupper = upper_sum_total / frames_accum
                Flower = lower_sum_total / frames_accum

                # Convert kcal/mol/Å → Newtons using your factor
                Fwall_N = (Fupper + Flower) / 2.0 * 69.5e-12

                # Instantaneous box Lx, Ly (Å) → Area (m^2)
                a, b, c = state_i.getPeriodicBoxVectors()
                aA = np.array(a.value_in_unit(angstroms))
                bA = np.array(b.value_in_unit(angstroms))
                Lx = float(np.linalg.norm(aA))
                Ly = float(np.linalg.norm(bA))
                Area_m2 = Lx * Ly * 1e-20

                # Osmotic pressure (bar)
                Osmotic_bar = Fwall_N / Area_m2 / 1.01e5

                time_ps = i * MD_BLOCK_TIME.value_in_unit(picoseconds)
                with open(NVT_OSMO_LOG, "a") as fosmo:
                    fosmo.write(f"{i}  {time_ps:.3f}  {Osmotic_bar:.6f}\n")

    # Save final NVT structure
    state_final = simmd.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=True)
    simmd.topology.setPeriodicBoxVectors(state_final.getPeriodicBoxVectors())
    with open(NVT_PDB, "w") as f:
        PDBFile.writeFile(simmd.topology, state_final.getPositions(), f)
    print(f"Wrote NVT output coordinate file: {NVT_PDB}")

    # ------------------ Optional post-run recomputation from DCD ------------------ #
    if RECOMPUTE_OSMOTIC_FROM_DCD and HAS_MDA:
        try:
            import MDAnalysis as mda
            print("Recomputing cumulative osmotic pressure from md_nvt_sel.dcd ...")

            # Universe built from the subset topology and DCD (1 ps per frame from MDAReporter)
            u = mda.Universe(str(NVT_SEL_PDB), str(NVT_SEL_DCD_FILE))

            # Scalars (same as above)
            k_wall_val = kz.value_in_unit(kilocalories_per_mole / angstroms**2)  # kcal/mol/Å^2
            z_max_val  = zmax.value_in_unit(angstroms)                            # Å
            z_min_val  = zmin.value_in_unit(angstroms)                            # Å

            cum_posz = []
            with open(NVT_OSMO_LOG, "a") as fosmo:
                fosmo.write("# Post-run recomputation from DCD (cumulative)\n")
                for fi, ts in enumerate(u.trajectory, start=1):
                    zA = u.atoms.positions[:, 2].astype(float)  # Å
                    cum_posz.append(zA)

                    if (fi % REPORT_CADENCE_BLOCKS) == 0:
                        pos_z = np.concatenate(cum_posz).astype(float)
                        n_frames = len(pos_z)

                        # Box lengths from MDAnalysis (Å, Å, Å, deg, deg, deg)
                        Lx, Ly, Lz, alpha, beta, gamma = ts.dimensions
                        Area = float(Lx) * float(Ly) * 1e-20  # m^2

                        i_z_max = pos_z >= z_max_val
                        i_z_min = pos_z <= z_min_val
                        UpperWall = k_wall_val * np.abs(pos_z[i_z_max] - z_max_val)
                        LowerWall = k_wall_val * np.abs(pos_z[i_z_min] - z_min_val)
                        Fupper = float(np.sum(UpperWall) / n_frames)
                        Flower = float(np.sum(LowerWall) / n_frames)
                        Fwall  = (Fupper + Flower) / 2.0 * 69.5e-12  # N
                        Osmotic = Fwall / Area / 1.01e5

                        time_ps = fi * MD_BLOCK_TIME.value_in_unit(picoseconds)  # 1 frame == 1 ps in our MDAReporter
                        fosmo.write(f"{fi}  {time_ps:.3f}  {Osmotic:.6f}\n")
                        fosmo.flush()

            print("Post-run DCD recomputation completed.")
        except Exception as e:
            print(f"[WARN] DCD post-processing failed: {e}")


if __name__ == "__main__":
    main()

