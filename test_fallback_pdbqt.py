import sys
sys.path.append('.')
from scripts.prep_structures import _simple_pdb_to_pdbqt

input_pdb = 'test_protein.pdb'
output_pdbqt = 'test_protein.pdbqt'

_simple_pdb_to_pdbqt(input_pdb, output_pdbqt)

with open(output_pdbqt) as f:
    lines = f.readlines()

plus_bug = any('+0.000' in line for line in lines)
charge_format_ok = all(line.rstrip().endswith((' C',' N',' O',' S',' P',' H')) for line in lines if line.startswith(('ATOM','HETATM')))

print(f"+0.000 bug present: {plus_bug}")
print(f"All atom lines end with valid AutoDock types: {charge_format_ok}")

# Print a few atom lines for manual inspection
print('Sample atom lines:')
for line in lines:
    if line.startswith(('ATOM','HETATM')):
        print(line.strip()) 