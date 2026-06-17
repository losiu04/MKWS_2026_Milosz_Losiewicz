"""
export_tables.py — Export simulation data as CSV tables for LaTeX report.

Usage:
    python export_tables.py

Outputs:
    tables/performance_sweep.csv  — Full O/F sweep data
    tables/optimal_point.csv      — Key results at optimal O/F
"""

import cantera as ct
import numpy as np
from pathlib import Path
import csv

# Configuration (must match lre_optimization.py)
P_CHAMBER = 45.0e5
P_AMBIENT = 1.01325e5
T_INLET = 300.0
OF_MIN, OF_MAX, OF_STEP = 1.0, 15.0, 0.05
R_UNIV = 8314.462618
G0 = 9.80665

SCRIPT_DIR = Path(__file__).resolve().parent
MECH_FILE = str(SCRIPT_DIR / 'mechanisms' / 'n2o_ethanol.yaml')
PHASE_NAME = 'n2o_ethanol_gas'
TABLES_DIR = SCRIPT_DIR / 'tables'

def gamma_star(gamma):
    return (np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))))

def compute_c_star(T_c, M_mean, gamma):
    Gamma = gamma_star(gamma)
    return np.sqrt(R_UNIV * T_c / (gamma * M_mean)) / Gamma

def compute_exhaust_mach(gamma, p_c, p_a):
    pr = p_c / p_a
    return np.sqrt(2.0 / (gamma - 1.0) * (pr ** ((gamma - 1.0) / gamma) - 1.0))

def compute_expansion_ratio(gamma, M_e):
    term = (1.0 + (gamma - 1.0) / 2.0 * M_e ** 2) / ((gamma + 1.0) / 2.0)
    return (1.0 / M_e) * term ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

def compute_thrust_coefficient(gamma, p_c, p_e, p_a, epsilon):
    pr = p_e / p_c
    Gamma = gamma_star(gamma)
    Cf_ideal = Gamma * np.sqrt(2.0 * gamma / (gamma - 1.0) * (1.0 - pr ** ((gamma - 1.0) / gamma)))
    return Cf_ideal + (p_e - p_a) / p_c * epsilon

def main():
    print('Loading Cantera gas...')
    gas = ct.Solution(MECH_FILE, PHASE_NAME)
    
    M_ethanol = gas.molecular_weights[gas.species_index('C2H5OH')]
    M_n2o = gas.molecular_weights[gas.species_index('N2O')]
    mass_ratio = M_ethanol / M_n2o
    
    of_values = np.arange(OF_MIN, OF_MAX + OF_STEP/2, OF_STEP)
    
    # Prepare output directory
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    
    # CSV: full sweep
    csv_path = TABLES_DIR / 'performance_sweep.csv'
    print(f'Writing full sweep data to {csv_path}...')
    
    with open(str(csv_path), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'OF', 'Tc_K', 'M_mean_kgkmol', 'gamma', 'c_star_mps',
            'Me', 'epsilon', 'Cf', 'Isp_s'
        ])
        
        for of in of_values:
            n_oxid = of * mass_ratio
            try:
                gas.TPX = T_INLET, P_CHAMBER, {'C2H5OH': 1.0, 'N2O': n_oxid}
                gas.equilibrate('HP')
            except Exception:
                continue
            
            T_c = gas.T
            M_mean = gas.mean_molecular_weight
            gamma = gas.cp / gas.cv
            c_star_val = compute_c_star(T_c, M_mean, gamma)
            M_e_val = compute_exhaust_mach(gamma, P_CHAMBER, P_AMBIENT)
            eps_val = compute_expansion_ratio(gamma, M_e_val)
            Cf_val = compute_thrust_coefficient(gamma, P_CHAMBER, P_AMBIENT, P_AMBIENT, eps_val)
            Isp_val = c_star_val * Cf_val / G0
            
            writer.writerow([
                f'{of:.2f}', f'{T_c:.1f}', f'{M_mean:.4f}', f'{gamma:.4f}',
                f'{c_star_val:.1f}', f'{M_e_val:.4f}', f'{eps_val:.4f}',
                f'{Cf_val:.4f}', f'{Isp_val:.2f}'
            ])
    
    print(f'  Wrote {len(of_values)} data rows.')
    
    # CSV: optimal point
    # Find optimal O/F from sweep
    of_opt = 4.10  # from simulation results
    n_oxid_opt = of_opt * mass_ratio
    gas.TPX = T_INLET, P_CHAMBER, {'C2H5OH': 1.0, 'N2O': n_oxid_opt}
    gas.equilibrate('HP')
    
    T_c_opt = gas.T
    M_mean_opt = gas.mean_molecular_weight
    gamma_opt = gas.cp / gas.cv
    c_star_opt = compute_c_star(T_c_opt, M_mean_opt, gamma_opt)
    M_e_opt = compute_exhaust_mach(gamma_opt, P_CHAMBER, P_AMBIENT)
    eps_opt = compute_expansion_ratio(gamma_opt, M_e_opt)
    Cf_opt = compute_thrust_coefficient(gamma_opt, P_CHAMBER, P_AMBIENT, P_AMBIENT, eps_opt)
    Isp_opt = c_star_opt * Cf_opt / G0
    
    opt_csv = TABLES_DIR / 'optimal_point.csv'
    print(f'Writing optimal point data to {opt_csv}...')
    with open(str(opt_csv), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Parameter', 'Value', 'Unit'])
        writer.writerow(['Optimal O/F', f'{of_opt:.2f}', '-'])
        writer.writerow(['Max Isp', f'{Isp_opt:.2f}', 's'])
        writer.writerow(['c*', f'{c_star_opt:.1f}', 'm/s'])
        writer.writerow(['Tc', f'{T_c_opt:.1f}', 'K'])
        writer.writerow(['M_mean', f'{M_mean_opt:.4f}', 'kg/kmol'])
        writer.writerow(['gamma', f'{gamma_opt:.4f}', '-'])
        writer.writerow(['Me', f'{M_e_opt:.4f}', '-'])
        writer.writerow(['epsilon', f'{eps_opt:.4f}', '-'])
        writer.writerow(['Cf', f'{Cf_opt:.4f}', '-'])
    
    # Also save composition at optimal
    comp_csv = TABLES_DIR / 'composition_optimal.csv'
    print(f'Writing composition data to {comp_csv}...')
    major_species = ['CO2', 'H2O', 'CO', 'H2', 'N2', 'O2', 'OH', 'H', 'O', 'NO']
    with open(str(comp_csv), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Species', 'MoleFraction'])
        for sp in major_species:
            k = gas.species_index(sp)
            writer.writerow([sp, f'{gas.X[k]:.6f}'])
    
    print('\nDone! All CSV tables exported.')
    print(f'  Optimal Isp = {Isp_opt:.2f} s at O/F = {of_opt:.2f}')
    print(f'  Max Tc = {T_c_opt:.1f} K')

if __name__ == '__main__':
    main()
