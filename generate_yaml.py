"""
generate_yaml.py — Generate the n2o_ethanol.yaml mechanism file for Cantera.

Strategy:
  Read the GRI-Mech 3.0 thermo data directly from the gri30.yaml file,
  extract the species we need, add C2H5OH (Ethanol) with Burcat NASA7 data,
  and write a clean YAML output file.
"""

import cantera as ct
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR = Path('mechanisms')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
MECH_OUTPUT = OUTPUT_DIR / 'n2o_ethanol.yaml'

# Species we need from gri30.yaml
GRI30_SPECIES = ['N2O', 'CO2', 'H2O', 'N2', 'O2', 'CO', 'H2', 'OH', 'H', 'O', 'NO']

# NASA7 coefficients for C2H5OH (Ethanol) from Burcat database
# Low  T (200-1000 K):  [a1, a2, a3, a4, a5, a6, a7]
# High T (1000-6000 K): [a1, a2, a3, a4, a5, a6, a7]
C2H5OH_LOW  = [4.858, -0.00319, 1.25e-4, -1.38e-7, 5.90e-11, -28500, 2.70]
C2H5OH_HIGH = [7.257, 0.00758, -2.92e-6, 5.31e-10, -3.52e-14, -31300, -6.52]

# ---------------------------------------------------------------------------
# Step 1: Load gri30.yaml and extract species thermo data
# ---------------------------------------------------------------------------
gas = ct.Solution('gri30.yaml')

# Collect the YAML species entries for our target species
species_yaml_entries = {}
for name in GRI30_SPECIES:
    sp = gas.species(name)
    thermo_data = sp.thermo.input_data
    entry = {
        'name': name,
        'composition': dict(sp.composition),
        'thermo': thermo_data
    }
    species_yaml_entries[name] = entry
    print(f'  Collected: {name} (composition: {entry["composition"]})')

print(f'\nTotal species collected: {len(species_yaml_entries)}')

# ---------------------------------------------------------------------------
# Step 2: Create C2H5OH species entry
# ---------------------------------------------------------------------------
ethanol_entry = {
    'name': 'C2H5OH',
    'composition': {'C': 2, 'H': 6, 'O': 1},
    'thermo': {
        'model': 'NASA7',
        'temperature-ranges': [200.0, 1000.0, 6000.0],
        'data': [C2H5OH_LOW, C2H5OH_HIGH]
    }
}
species_yaml_entries['C2H5OH'] = ethanol_entry
print(f'  Added: C2H5OH with Burcat NASA7 data')

# ---------------------------------------------------------------------------
# Step 3: Build the full YAML document as a string
# ---------------------------------------------------------------------------
def compose_species_yaml(entry):
    """Compose a YAML species entry string."""
    name = entry['name']
    comp = entry['composition']
    thermo = entry['thermo']
    
    lines = []
    lines.append(f"  - name: {name}")
    
    # Composition
    comp_str = ', '.join(f'{k}: {v}' for k, v in comp.items())
    lines.append(f"    composition: {{{comp_str}}}")
    
    # Thermo
    lines.append(f"    thermo:")
    lines.append(f"      model: {thermo['model']}")
    tr = thermo['temperature-ranges']
    lines.append(f"      temperature-ranges: [{tr[0]}, {tr[1]}, {tr[2]}]")
    lines.append(f"      data:")
    for coeffs in thermo['data']:
        coeff_str = ', '.join(f'{c:.14g}' for c in coeffs)
        lines.append(f"      - [{coeff_str}]")
    
    return '\n'.join(lines)

# Build the full YAML content
yaml_lines = []
yaml_lines.append('description: |-')
yaml_lines.append('  Reduced N2O/Ethanol mechanism for LRE (Liquid Rocket Engine)')
yaml_lines.append('  thermochemical equilibrium analysis.')
yaml_lines.append('  Species extracted from GRI-Mech 3.0, with C2H5OH (Ethanol)')
yaml_lines.append('  added from Burcat NASA7 database.')
yaml_lines.append('generator: generate_yaml.py')
yaml_lines.append('cantera-version: 3.2.0')
yaml_lines.append('')
yaml_lines.append('phases:')
yaml_lines.append('- name: n2o_ethanol_gas')
yaml_lines.append('  thermo: ideal-gas')
yaml_lines.append('  elements: [O, H, C, N]')

# Format species list as a YAML flow sequence (single line with continuations)
all_names = sorted(species_yaml_entries.keys())
species_line = '  species: [' + ', '.join(all_names) + ']'
yaml_lines.append(species_line)

yaml_lines.append('  skip-undeclared-elements: true')
yaml_lines.append('  state: {T: 300.0, P: 1 atm}')
yaml_lines.append('')

# Species section
yaml_lines.append('species:')
for name in all_names:
    entry = species_yaml_entries[name]
    yaml_lines.append(compose_species_yaml(entry))
    yaml_lines.append('')

yaml_content = '\n'.join(yaml_lines)

# ---------------------------------------------------------------------------
# Step 4: Write the YAML file
# ---------------------------------------------------------------------------
print(f'\nWriting mechanism to: {MECH_OUTPUT}')
with open(str(MECH_OUTPUT), 'w') as f:
    f.write(yaml_content)
print(f'  File written ({MECH_OUTPUT.stat().st_size} bytes)')

# Preview first 30 lines
print('\n--- Preview (first 30 lines) ---')
for i, line in enumerate(yaml_content.split('\n')[:30]):
    print(f'{i+1:3d}: {line}')

# ---------------------------------------------------------------------------
# Step 5: Verify the generated mechanism
# ---------------------------------------------------------------------------
print('\n' + '='*60)
print('Verifying the generated mechanism ...')
print('='*60)

try:
    gas_test = ct.Solution(str(MECH_OUTPUT), 'n2o_ethanol_gas')
    print(f'  Loaded: {gas_test.n_species} species')
    print(f'  Species: {gas_test.species_names}')

    # Verify C2H5OH
    if 'C2H5OH' in gas_test.species_names:
        print('  C2H5OH: PRESENT')
    else:
        raise RuntimeError('C2H5OH not found!')

    # HP equilibrium test
    print('\n  HP equilibrium test:')
    gas_test.TPX = 300.0, 20e5, {'C2H5OH': 1.0, 'N2O': 6.0}
    gas_test.equilibrate('HP')
    print(f'  T_combustion = {gas_test.T:.1f} K')
    print(f'  Mean molar mass = {gas_test.mean_molecular_weight:.4f} kg/kmol')
    print(f'  gamma = {gas_test.cp / gas_test.cv:.4f}')
    print('  Verification PASSED')

except Exception as e:
    print(f'  ERROR: {e}')
    import traceback
    traceback.print_exc()
    print('  Verification FAILED')
    raise

print('\n' + '='*60)
print('Mechanism file is ready!')
print('='*60)
