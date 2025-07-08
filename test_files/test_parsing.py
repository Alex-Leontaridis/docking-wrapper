from scripts.docking_results_parser import DockingResultsParser

print('--- Combined Parsing Test (generate_summary) ---')
parser = DockingResultsParser('.', '.')
df = parser.generate_summary('test_ligand')
print(df)

print('--- Vina Parsing Test ---')
parser = DockingResultsParser('.', '.')
vina_results = parser.parse_vina_output('test_vina_output.txt')
print(f'Vina results: {len(vina_results)} poses found')
if vina_results:
    print(f'First pose: {vina_results[0]}')

print('\n--- GNINA Parsing Test ---')
gnina_results = parser.parse_gnina_output('test_gnina_output.txt')
print(f'GNINA results: {len(gnina_results)} poses found')
if gnina_results:
    print(f'First pose: {gnina_results[0]}')

print('\n--- DiffDock Parsing Test ---')
diffdock_results = parser.parse_diffdock_output('test_diffdock_output')
print(f'DiffDock results: {len(diffdock_results)} poses found')
if diffdock_results:
    print(f'First pose: {diffdock_results[0]}') 