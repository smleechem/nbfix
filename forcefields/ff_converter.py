import math

def charmm_to_openmm_sigma_epsilon(Rmin_over_2_ang, Emin_kcal_per_mol):
    """
    Convert CHARMM-style LJ parameters (Rmin/2 in Å, Emin in kcal/mol)
    to OpenMM-style sigma (nm) and epsilon (kJ/mol).

    Parameters
    ----------
    Rmin_over_2_ang : float
        CHARMM Rmin/2 in Angstrom
    Emin_kcal_per_mol : float
        Energy at the minimum (typically negative) in kcal/mol

    Returns
    -------
    sigma_nm : float
        OpenMM sigma in nm
    epsilon_kJ_per_mol : float
        OpenMM epsilon (well depth, positive) in kJ/mol
    """
    Rmin_ang = 2.0 * Rmin_over_2_ang                  # Å
    sigma_ang = Rmin_ang / (2.0 ** (1.0 / 6.0))       # Å
    sigma_nm = sigma_ang * 0.1                        # nm
    epsilon_kJ_per_mol = abs(Emin_kcal_per_mol) * 4.184
    return sigma_nm, epsilon_kJ_per_mol


# ---- Example with your values ----
Rmin_over_2_ang = 1.41075       # Å
Emin_kcal_per_mol = -0.0469     # kcal/mol

sigma_nm, epsilon_kJ = charmm_to_openmm_sigma_epsilon(Rmin_over_2_ang, Emin_kcal_per_mol)
print(f"sigma (nm)   = {sigma_nm:.9f}")
print(f"epsilon (kJ/mol) = {epsilon_kJ:.6f}")

