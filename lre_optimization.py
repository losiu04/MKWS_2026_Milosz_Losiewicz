"""
lre_optimization.py — N2O/Ethanol Liquid Rocket Engine (LRE) Optimization
=======================================================================

Thermochemical equilibrium analysis and performance optimization for a
bipropellant liquid rocket engine using Nitrous Oxide (N2O) as oxidizer
and Ethanol (C2H5OH) as fuel.

Computational approach:
  1. Load a reduced Cantera mechanism (n2o_ethanol.yaml) containing only
     the relevant species: N2O, C2H5OH, CO2, H2O, N2, O2, CO, H2, OH, H,
     O, NO.
  2. Sweep over a range of Oxidizer-to-Fuel (O/F) mass ratios.
  3. For each O/F ratio, compute the HP (constant enthalpy, constant
     pressure) equilibrium to obtain:
       - Adiabatic flame temperature, T_c
       - Mean molar mass of product gases, M_mean
       - Specific heat ratio, gamma = cp/cv
       - Equilibrium product composition
  4. Compute rocket performance parameters:
       - Characteristic velocity:     c* = sqrt(R * T_c / (gamma * M_mean)) / Gamma
       - Thrust coefficient:          Cf = f(gamma, expansion ratio, pressure ratio)
       - Specific impulse:            Isp = c* * Cf / g0
  5. Determine the optimal O/F ratio that maximizes Isp.
  6. At the optimal O/F, compute the nozzle expansion ratio (epsilon) for
     ideal expansion (p_exit = p_ambient).

Key equations (Vandenkerckhove):
  Gamma = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
  c* = sqrt(R * T_c / (gamma * M_mean)) / Gamma

References:
  - Sutton, G.P. & Biblarz, O., "Rocket Propulsion Elements", 9th ed.
  - Cantera: https://cantera.org

Author: Auto-generated for MKWS project
Course: Metody komputerowe w spalaniu (Computer Methods in Combustion)
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

# ============================================================================
# Configuration
# ============================================================================

# Chamber and ambient conditions
P_CHAMBER = 45.0e5          # Chamber pressure [Pa] (45 bar)
T_INLET   = 300.0            # Inlet temperature [K]
P_AMBIENT = 1.01325e5        # Ambient pressure [Pa] (1 atm)

# O/F sweep parameters
OF_MIN  = 1.0                # Minimum O/F ratio [-]
OF_MAX  = 15.0               # Maximum O/F ratio [-]
OF_STEP = 0.05               # O/F step size [-]

# Mechanism file (relative to this script's location)
SCRIPT_DIR = Path(__file__).resolve().parent
MECH_FILE = str(SCRIPT_DIR / 'mechanisms' / 'n2o_ethanol.yaml')
PHASE_NAME = 'n2o_ethanol_gas'

# Output directory for figures
FIGURES_DIR = SCRIPT_DIR / 'figures'

# Standard gravity
G0 = 9.80665                 # [m/s^2]

# Universal gas constant
R_UNIV = 8314.462618         # [J/(kmol*K)]

# ============================================================================
# Helper: Vandenkerckhove function
# ============================================================================

def compute_gamma_star(gamma):
    """
    Vandenkerckhove's Gamma(gamma) function.

    Gamma = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

    This function appears in the definition of characteristic velocity c*
    and relates the specific heat ratio to the sonic flow conditions at
    the throat of a converging-diverging nozzle.

    Parameters
    ----------
    gamma : float or ndarray
        Specific heat ratio cp/cv [-].

    Returns
    -------
    float or ndarray
        Vandenkerckhove Gamma [-].
    """
    return (np.sqrt(gamma)
            * (2.0 / (gamma + 1.0))
            ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))))


# ============================================================================
# Helper: Characteristic velocity
# ============================================================================

def compute_c_star(T_c, M_mean, gamma):
    """
    Characteristic velocity c* [m/s].

    c* = sqrt(R * T_c / (gamma * M_mean)) / Gamma(gamma)

    where:
      R      = universal gas constant [J/(kmol*K)]
      T_c    = combustion temperature [K]
      M_mean = mean molar mass of products [kg/kmol]
      gamma  = specific heat ratio [-]

    Parameters
    ----------
    T_c    : float
        Combustion temperature [K].
    M_mean : float
        Mean molar mass of product gases [kg/kmol].
    gamma  : float
        Specific heat ratio [-].

    Returns
    -------
    float
        Characteristic velocity [m/s].
    """
    Gamma = compute_gamma_star(gamma)
    a = np.sqrt(R_UNIV * T_c / (gamma * M_mean))
    return a / Gamma


# ============================================================================
# Helper: Thrust coefficient
# ============================================================================

def compute_thrust_coefficient(gamma, p_c, p_e, p_a, M_e=None, epsilon=None):
    """
    Thrust coefficient Cf [-].

    For a diverging nozzle with known exit Mach number and expansion ratio:
      Cf = Gamma * sqrt(2*gamma/(gamma-1) * (1 - (p_e/p_c)^((gamma-1)/gamma)))
           + (p_e - p_a)/p_c * epsilon

    Parameters
    ----------
    gamma   : float
        Specific heat ratio [-].
    p_c     : float
        Chamber pressure [Pa].
    p_e     : float
        Nozzle exit pressure [Pa] (= p_a for ideal expansion).
    p_a     : float
        Ambient pressure [Pa].
    M_e     : float, optional
        Exit Mach number [-].
    epsilon : float, optional
        Nozzle expansion ratio A_e/A_t [-].

    Returns
    -------
    float
        Thrust coefficient [-].
    """
    if M_e is None and epsilon is None:
        raise ValueError("Either M_e or epsilon must be provided.")

    # Pressure ratio
    pr = p_e / p_c

    # Ideal thrust coefficient (vacuum contribution excluded)
    Gamma = compute_gamma_star(gamma)
    Cf_ideal = Gamma * np.sqrt(
        2.0 * gamma / (gamma - 1.0)
        * (1.0 - pr ** ((gamma - 1.0) / gamma))
    )

    # Pressure correction term (p_e - p_a)/p_c * epsilon
    if epsilon is not None:
        Cf = Cf_ideal + (p_e - p_a) / p_c * epsilon
    else:
        Cf = Cf_ideal

    return Cf


# ============================================================================
# Helper: Exit Mach number for ideal expansion (p_e = p_a)
# ============================================================================

def compute_exhaust_mach(gamma, p_c, p_a):
    """
    Exit Mach number M_e for isentropic expansion to p_e = p_a.

    From isentropic relation:
      p_c / p_e = (1 + (gamma-1)/2 * M_e^2) ^ (gamma/(gamma-1))

    Rearranged:
      M_e = sqrt(2/(gamma-1) * ((p_c/p_e)^((gamma-1)/gamma) - 1))

    Parameters
    ----------
    gamma : float
        Specific heat ratio [-].
    p_c   : float
        Chamber pressure [Pa].
    p_a   : float
        Ambient pressure [Pa] (= p_e for ideal expansion).

    Returns
    -------
    float
        Exit Mach number [-].
    """
    pr = p_c / p_a
    M_e = np.sqrt(
        2.0 / (gamma - 1.0) * (pr ** ((gamma - 1.0) / gamma) - 1.0)
    )
    return M_e


# ============================================================================
# Helper: Nozzle expansion ratio
# ============================================================================

def compute_expansion_ratio(gamma, M_e):
    """
    Nozzle expansion ratio epsilon = A_e / A_t [-].

    From isentropic area-Mach relation:
      epsilon = (1/M_e) * ((1 + (gamma-1)/2 * M_e^2)
                / ((gamma+1)/2)) ^ ((gamma+1)/(2*(gamma-1)))

    Parameters
    ----------
    gamma : float
        Specific heat ratio [-].
    M_e   : float
        Exit Mach number [-].

    Returns
    -------
    float
        Nozzle expansion ratio A_e / A_t [-].
    """
    term = (1.0 + (gamma - 1.0) / 2.0 * M_e ** 2) / ((gamma + 1.0) / 2.0)
    epsilon = (1.0 / M_e) * term ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return epsilon


# ============================================================================
# Helper: Specific impulse
# ============================================================================

def compute_specific_impulse(c_star, Cf):
    """
    Specific impulse Isp [s].

    Isp = c* * Cf / g0

    Parameters
    ----------
    c_star : float
        Characteristic velocity [m/s].
    Cf     : float
        Thrust coefficient [-].

    Returns
    -------
    float
        Specific impulse [seconds].
    """
    return c_star * Cf / G0


# ============================================================================
# Main sweep: O/F ratio scan
# ============================================================================

def run_of_sweep(gas):
    """
    Sweep over O/F ratios and compute combustion/performance parameters.

    For each O/F value:
      1. Compute the molar composition:
           n_oxid = O/F * M_ethanol / M_n2o
           n_fuel = 1.0
      2. Set gas state: T_INLET, P_CHAMBER, {C2H5OH, N2O}
      3. Equilibrate at constant H, P (HP)
      4. Extract: T_c, mean_M, gamma, product mole fractions
      5. Compute: c*, Cf (for ideal expansion at p_e=p_a), Isp

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas object with the N2O/Ethanol mechanism loaded.

    Returns
    -------
    dict
        Dictionary containing arrays of O/F and corresponding results.
    """
    # Pre-compute molar mass ratio for O/F mass-to-mole conversion
    # O/F (mass) = (m_oxid / m_fuel)
    # n_oxid / n_fuel = O/F * M_fuel / M_oxid
    M_ethanol = gas.molecular_weights[gas.species_index('C2H5OH')]
    M_n2o     = gas.molecular_weights[gas.species_index('N2O')]
    mass_ratio = M_ethanol / M_n2o  # converts O/F mass ratio to mole ratio

    print(f'M(C2H5OH) = {M_ethanol:.4f} kg/kmol')
    print(f'M(N2O)    = {M_n2o:.4f} kg/kmol')
    print(f'Mass ratio factor (M_fuel/M_oxid) = {mass_ratio:.4f}')

    # Generate O/F array
    of_values = np.arange(OF_MIN, OF_MAX + OF_STEP / 2, OF_STEP)
    n_points = len(of_values)
    print(f'Sweeping O/F from {OF_MIN:.1f} to {OF_MAX:.1f} ({n_points} points)')

    # Pre-allocate result arrays
    results = {
        'OF':           of_values,
        'Tc':           np.full(n_points, np.nan),
        'M_mean':       np.full(n_points, np.nan),
        'gamma':        np.full(n_points, np.nan),
        'c_star':       np.full(n_points, np.nan),
        'Cf':           np.full(n_points, np.nan),
        'Isp':          np.full(n_points, np.nan),
        'Me':           np.full(n_points, np.nan),
        'epsilon':      np.full(n_points, np.nan),
        'species':      {},  # mole fractions for each species at each O/F
    }

    # Get the list of all species for tracking composition
    all_species = gas.species_names
    for sp_name in all_species:
        results['species'][sp_name] = np.full(n_points, np.nan)

    # ------------------------------------------------------------------
    # Warm-up: one equilibration to ensure the gas object is initialized
    # ------------------------------------------------------------------
    print('  Warm-up equilibration ...')
    gas.TPX = T_INLET, P_CHAMBER, {'C2H5OH': 1.0, 'N2O': OF_MIN * mass_ratio}
    gas.equilibrate('HP')
    print(f'  Warm-up complete: T = {gas.T:.1f} K')

    # ------------------------------------------------------------------
    # Main sweep loop
    # ------------------------------------------------------------------
    for i, of in enumerate(of_values):
        # Convert mass-based O/F to mole-based reactant ratio
        n_oxid = of * mass_ratio  # kmol N2O per kmol C2H5OH
        n_fuel = 1.0              # kmol C2H5OH

        # Set gas state
        gas.TPX = T_INLET, P_CHAMBER, {'C2H5OH': n_fuel, 'N2O': n_oxid}

        # Attempt HP equilibrium; skip if it fails (e.g., too rich/lean)
        try:
            gas.equilibrate('HP')
        except Exception as exc:
            warnings.warn(f'O/F = {of:.2f}: equilibrate failed: {exc}')
            continue

        # Extract equilibrium properties
        T_c          = gas.T
        mean_M       = gas.mean_molecular_weight  # [kg/kmol]
        gamma        = gas.cp / gas.cv
        p_exit       = P_AMBIENT  # ideal expansion: p_e = p_a

        # Compute performance parameters
        c_star_val = compute_c_star(T_c, mean_M, gamma)

        # Compute exhaust Mach and expansion ratio for ideal expansion
        M_e_val = compute_exhaust_mach(gamma, P_CHAMBER, p_exit)
        eps_val = compute_expansion_ratio(gamma, M_e_val)

        # Thrust coefficient with pressure correction
        Cf_val = compute_thrust_coefficient(
            gamma, P_CHAMBER, p_exit, P_AMBIENT,
            M_e=M_e_val, epsilon=eps_val
        )

        Isp_val = compute_specific_impulse(c_star_val, Cf_val)

        # Store results
        results['Tc'][i]      = T_c
        results['M_mean'][i]  = mean_M
        results['gamma'][i]   = gamma
        results['c_star'][i]  = c_star_val
        results['Cf'][i]      = Cf_val
        results['Isp'][i]     = Isp_val
        results['Me'][i]      = M_e_val
        results['epsilon'][i] = eps_val

        # Store species mole fractions
        # Use species_index to avoid reference issues with gas[sp_name].X
        for sp_name in all_species:
            k = gas.species_index(sp_name)
            results['species'][sp_name][i] = gas.X[k]

        # Progress indicator
        if (i + 1) % 200 == 0 or i == 0 or i == n_points - 1:
            print(f'  O/F = {of:6.2f}: Tc = {T_c:7.1f} K, '
                  f'M = {mean_M:.3f}, gamma = {gamma:.4f}, '
                  f'Isp = {Isp_val:.1f} s')

    return results


# ============================================================================
# Nozzle optimization at optimal O/F
# ============================================================================

def optimize_nozzle(gas, of_opt, p_c=P_CHAMBER, p_a=P_AMBIENT):
    """
    Perform detailed nozzle optimization at the optimal O/F ratio.

    At the optimal O/F, compute:
      - Exit Mach number (isentropic)
      - Expansion ratio epsilon = A_e / A_t
      - Thrust coefficient Cf

    Parameters
    ----------
    gas    : ct.Solution
        Cantera gas object.
    of_opt : float
        Optimal O/F ratio.
    p_c    : float
        Chamber pressure [Pa].
    p_a    : float
        Ambient pressure [Pa].

    Returns
    -------
    dict
        Nozzle optimization results.
    """
    M_ethanol = gas.molecular_weights[gas.species_index('C2H5OH')]
    M_n2o     = gas.molecular_weights[gas.species_index('N2O')]
    mass_ratio = M_ethanol / M_n2o

    # Set state at optimal O/F and equilibrate
    gas.TPX = T_INLET, p_c, {'C2H5OH': 1.0, 'N2O': of_opt * mass_ratio}
    gas.equilibrate('HP')

    T_c    = gas.T
    M_mean = gas.mean_molecular_weight
    gamma  = gas.cp / gas.cv

    # Isentropic expansion
    M_e = compute_exhaust_mach(gamma, p_c, p_a)
    eps = compute_expansion_ratio(gamma, M_e)
    c_star_val = compute_c_star(T_c, M_mean, gamma)
    Cf_val = compute_thrust_coefficient(
        gamma, p_c, p_a, p_a, M_e=M_e, epsilon=eps
    )
    Isp_val = compute_specific_impulse(c_star_val, Cf_val)

    nozzle_results = {
        'OF_opt':   of_opt,
        'Tc':       T_c,
        'M_mean':   M_mean,
        'gamma':    gamma,
        'c_star':   c_star_val,
        'M_e':      M_e,
        'epsilon':  eps,
        'Cf':       Cf_val,
        'Isp':      Isp_val,
    }

    return nozzle_results


# ============================================================================
# Plotting: Main performance curves
# ============================================================================

def plot_results(results):
    """
    Generate a 2x2 figure showing Isp, c*, Tc, and mean_M vs O/F.

    All plots share the O/F x-axis for easy comparison. The optimal O/F
    (where Isp is maximized) is marked with a vertical dashed line.

    Parameters
    ----------
    results : dict
        Dictionary of result arrays from run_of_sweep().
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    of = results['OF']

    # Find the O/F that maximizes Isp
    idx_opt = np.nanargmax(results['Isp'])
    of_opt = of[idx_opt]

    # Use a clean plotting style
    plt.style.use('seaborn-v0_8-whitegrid')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # --- Panel 1: Specific Impulse Isp ---
    ax = axes[0, 0]
    ax.plot(of, results['Isp'], 'b-', linewidth=2, label='Isp')
    ax.axvline(of_opt, color='r', linestyle='--', linewidth=1.2,
               label=f'Optimum O/F = {of_opt:.2f}')
    ax.plot(of_opt, results['Isp'][idx_opt], 'ro', markersize=8)
    ax.set_xlabel('O/F Ratio [-]', fontsize=12)
    ax.set_ylabel('Isp [s]', fontsize=12)
    ax.set_title('Specific Impulse vs O/F Ratio', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.set_xlim(OF_MIN, OF_MAX)
    ax.tick_params(axis='both', labelsize=10)

    # --- Panel 2: Characteristic Velocity c* ---
    ax = axes[0, 1]
    ax.plot(of, results['c_star'] / 1e3, 'g-', linewidth=2, label='c*')
    ax.axvline(of_opt, color='r', linestyle='--', linewidth=1.2)
    ax.plot(of_opt, results['c_star'][idx_opt] / 1e3, 'ro', markersize=8)
    ax.set_xlabel('O/F Ratio [-]', fontsize=12)
    ax.set_ylabel('c* [km/s]', fontsize=12)
    ax.set_title('Characteristic Velocity vs O/F Ratio', fontsize=13, fontweight='bold')
    ax.set_xlim(OF_MIN, OF_MAX)
    ax.tick_params(axis='both', labelsize=10)

    # --- Panel 3: Combustion Temperature Tc ---
    ax = axes[1, 0]
    ax.plot(of, results['Tc'] / 1e3, 'm-', linewidth=2, label='Tc')
    ax.axvline(of_opt, color='r', linestyle='--', linewidth=1.2)
    ax.plot(of_opt, results['Tc'][idx_opt] / 1e3, 'ro', markersize=8)
    ax.set_xlabel('O/F Ratio [-]', fontsize=12)
    ax.set_ylabel('Tc [kK]', fontsize=12)
    ax.set_title('Adiabatic Flame Temperature vs O/F Ratio',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(OF_MIN, OF_MAX)
    ax.tick_params(axis='both', labelsize=10)

    # --- Panel 4: Mean Molar Mass ---
    ax = axes[1, 1]
    ax.plot(of, results['M_mean'], 'c-', linewidth=2, label='M_mean')
    ax.axvline(of_opt, color='r', linestyle='--', linewidth=1.2)
    ax.plot(of_opt, results['M_mean'][idx_opt], 'ro', markersize=8)
    ax.set_xlabel('O/F Ratio [-]', fontsize=12)
    ax.set_ylabel('Mean Molar Mass [kg/kmol]', fontsize=12)
    ax.set_title('Mean Molecular Weight of Products vs O/F Ratio',
                 fontsize=13, fontweight='bold')
    ax.set_xlim(OF_MIN, OF_MAX)
    ax.tick_params(axis='both', labelsize=10)

    plt.tight_layout()
    fig_path = FIGURES_DIR / 'lre_performance.png'
    plt.savefig(str(fig_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Performance plot saved: {fig_path}')


# ============================================================================
# Plotting: Product composition
# ============================================================================

def plot_composition(results):
    """
    Generate a figure showing equilibrium product mole fractions vs O/F.

    The major species (CO2, H2O, CO, H2, N2, O2, OH, H, O, NO) are shown
    in a single plot.

    Parameters
    ----------
    results : dict
        Dictionary of result arrays from run_of_sweep().
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    of = results['OF']
    species = results['species']

    # Major species to plot
    major_species = ['CO2', 'H2O', 'CO', 'H2', 'N2', 'O2', 'OH', 'H', 'O', 'NO']
    colors = plt.cm.tab10(np.linspace(0, 1, len(major_species)))

    plt.style.use('seaborn-v0_8-whitegrid')

    fig, ax = plt.subplots(figsize=(14, 8))

    for i, sp_name in enumerate(major_species):
        ax.plot(of, species[sp_name], color=colors[i], linewidth=2,
                label=sp_name)

    ax.set_xlabel('O/F Ratio [-]', fontsize=12)
    ax.set_ylabel('Mole Fraction X [-]', fontsize=12)
    ax.set_title('Equilibrium Product Composition vs O/F Ratio\n'
                 f'(N2O/Ethanol, p_chamber = {P_CHAMBER/1e5:.0f} bar)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, ncol=2, loc='upper right')
    ax.set_xlim(OF_MIN, OF_MAX)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', labelsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = FIGURES_DIR / 'lre_composition.png'
    plt.savefig(str(fig_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Composition plot saved: {fig_path}')


# ============================================================================
# Plotting: Isp vs Expansion Ratio
# ============================================================================

def plot_isp_vs_epsilon(results):
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    of = results['OF']
    idx_opt = np.nanargmax(results['Isp'])
    of_opt = of[idx_opt]

    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.plot(results['epsilon'], results['Isp'], 'b-', linewidth=2, label='Isp')

    eps_opt = results['epsilon'][idx_opt]
    isp_opt = results['Isp'][idx_opt]
    ax.plot(eps_opt, isp_opt, 'ro', markersize=8)
    ax.axvline(eps_opt, color='r', linestyle='--', linewidth=1.2)
    ax.axhline(isp_opt, color='r', linestyle='--', linewidth=1.2)

    ax.set_xlabel('Expansion Ratio $\\varepsilon = A_e/A_t$ [-]', fontsize=12)
    ax.set_ylabel('Specific Impulse Isp [s]', fontsize=12)
    ax.set_title('Specific Impulse vs Nozzle Expansion Ratio', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.tick_params(axis='both', labelsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = FIGURES_DIR / 'lre_isp_vs_epsilon.png'
    plt.savefig(str(fig_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f'  Isp vs epsilon plot saved: {fig_path}')


# ============================================================================
# Print summary table
# ============================================================================

def print_summary(results, nozzle_opt):
    """
    Print a formatted summary table of optimization results.

    Parameters
    ----------
    results    : dict
        Results from run_of_sweep().
    nozzle_opt : dict
        Results from optimize_nozzle().
    """
    of = results['OF']
    idx_opt = np.nanargmax(results['Isp'])
    of_opt = of[idx_opt]

    print('\n' + '='*70)
    print('  N2O/Ethanol LRE OPTIMIZATION SUMMARY')
    print('='*70)
    print(f'  Chamber pressure:         {P_CHAMBER/1e5:.1f} bar')
    print(f'  Ambient pressure:         {P_AMBIENT/1e5:.4f} bar')
    print(f'  Inlet temperature:         {T_INLET:.1f} K')
    print(f'  O/F sweep range:          {OF_MIN:.1f} - {OF_MAX:.1f}')
    print(f'  O/F sweep step:           {OF_STEP:.2f}')
    print('-'*70)
    print(f'  Optimal O/F ratio:        {of_opt:.2f}')
    print(f'  Maximum Isp:              {results["Isp"][idx_opt]:.2f} s')
    print(f'  c* at optimum:            {results["c_star"][idx_opt]:.1f} m/s')
    print(f'  Tc at optimum:            {results["Tc"][idx_opt]:.1f} K')
    print(f'  Mean M at optimum:        {results["M_mean"][idx_opt]:.4f} kg/kmol')
    print(f'  Gamma at optimum:         {results["gamma"][idx_opt]:.4f}')
    print('-'*70)
    print('  NOZZLE OPTIMIZATION (ideal expansion, p_e = p_a)')
    print(f'  Exit Mach number Me:      {nozzle_opt["M_e"]:.4f}')
    print(f'  Expansion ratio epsilon:  {nozzle_opt["epsilon"]:.4f}')
    print(f'  Thrust coefficient Cf:    {nozzle_opt["Cf"]:.4f}')
    print(f'  Isp at optimum (nozzle):  {nozzle_opt["Isp"]:.2f} s')
    print('='*70)

    # Additional metrics: O/F at max Tc for comparison
    idx_max_tc = np.nanargmax(results['Tc'])
    print(f'\n  Note: Maximum Tc = {results["Tc"][idx_max_tc]:.1f} K '
          f'occurs at O/F = {of[idx_max_tc]:.2f}')
    print(f'  (Demonstrates that peak Isp and peak Tc do NOT coincide.)')


# ============================================================================
# Main
# ============================================================================

def main():
    """
    Main execution routine.

    Steps:
      1. Verify that the mechanism file exists.
      2. Load the Cantera Solution.
      3. Run the O/F sweep.
      4. Find the optimal Isp.
      5. Run nozzle optimization.
      6. Generate plots.
      7. Print summary.
    """
    print('='*70)
    print('  N2O/Ethanol LRE — Thermochemical Optimization')
    print('='*70)

    # ------------------------------------------------------------------
    # Step 1: Check mechanism file
    # ------------------------------------------------------------------
    mech_path = Path(MECH_FILE)
    if not mech_path.exists():
        raise FileNotFoundError(
            f'Mechanism file not found: {MECH_FILE}\n'
            f'Run generate_yaml.py first to create it.'
        )
    print(f'\nMechanism file: {mech_path.resolve()}')

    # ------------------------------------------------------------------
    # Step 2: Load Cantera Solution
    # ------------------------------------------------------------------
    print('Loading Cantera gas phase ...')
    gas = ct.Solution(str(mech_path), PHASE_NAME)
    print(f'  Loaded: {gas.n_species} species')
    print(f'  Species: {gas.species_names}')

    # ------------------------------------------------------------------
    # Step 3: Run O/F sweep
    # ------------------------------------------------------------------
    print('\n' + '-'*70)
    print('O/F SWEEP')
    print('-'*70)
    results = run_of_sweep(gas)

    # ------------------------------------------------------------------
    # Step 4: Find optimum
    # ------------------------------------------------------------------
    valid = ~np.isnan(results['Isp'])
    if not np.any(valid):
        raise RuntimeError('No valid equilibrium solutions found!')

    idx_opt = np.nanargmax(results['Isp'])
    of_opt = results['OF'][idx_opt]
    print(f'\nOptimal O/F = {of_opt:.2f} (Isp = {results["Isp"][idx_opt]:.1f} s)')

    # ------------------------------------------------------------------
    # Step 5: Nozzle optimization at optimal O/F
    # ------------------------------------------------------------------
    print('\n' + '-'*70)
    print('NOZZLE OPTIMIZATION')
    print('-'*70)
    nozzle_results = optimize_nozzle(gas, of_opt)

    # ------------------------------------------------------------------
    # Step 6: Generate plots
    # ------------------------------------------------------------------
    print('\n' + '-'*70)
    print('GENERATING PLOTS')
    print('-'*70)
    plot_results(results)
    plot_composition(results)
    plot_isp_vs_epsilon(results)

    # ------------------------------------------------------------------
    # Step 7: Print summary
    # ------------------------------------------------------------------
    print_summary(results, nozzle_results)

    print('\nDone. All results computed and figures saved.')


# ============================================================================
# Entry point
# ============================================================================

if __name__ == '__main__':
    main()
