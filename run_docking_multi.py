#!/usr/bin/env python3
import argparse
import os
import sys
import logging
from pathlib import Path
import subprocess
import json
import time
import numpy as np
import shutil

# --- Constants for output directories ---
VINA_OUT = 'vina_output'
GNINA_OUT = 'gnina_output'
DIFFDOCK_OUT = 'diffdock_output'
LOGS_DIR = 'logs'
LOG_FILE = os.path.join(LOGS_DIR, 'docking_run.log')

# --- Setup logging ---
def setup_logging():
    os.makedirs(LOGS_DIR, exist_ok=True)
    logging.basicConfig(
        filename=LOG_FILE,
        filemode='a',
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

# --- Ensure output directories exist ---
def ensure_output_dirs(base_dir):
    for d in [VINA_OUT, GNINA_OUT, DIFFDOCK_OUT, LOGS_DIR]:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)

# --- Validate input files ---
def validate_file(path, desc):
    if not path or not os.path.isfile(path):
        logging.error(f"{desc} file '{path}' does not exist or is not a file.")
        sys.exit(1)

def extract_box_from_protein(protein_path):
    """
    Extract docking box parameters from a protein PDBQT file using multiple strategies:
    1. Bound ligand detection (original method)
    2. Largest cavity detection using protein geometry
    3. Geometric center of protein as fallback
    4. Conservative default box size
    
    Returns (center_x, center_y, center_z, size_x, size_y, size_z)
    """
    heavy_atoms_set = set([
        'C', 'N', 'O', 'P', 'S', 'F', 'CL', 'BR', 'I', 'B', 'SE', 'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'V', 'W', 'MO', 'CD', 'HG', 'SR', 'K', 'NA', 'CS', 'BA', 'AL', 'CR', 'TI', 'PB', 'SB', 'AS', 'SN', 'AG', 'AU', 'GA', 'IN', 'TL', 'PT', 'RB', 'LI', 'Y', 'Zr', 'Nb', 'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Ta', 'Bi', 'U', 'Th', 'Pa', 'Ac', 'Ra', 'Fr', 'Po', 'At', 'Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He'
    ])
    
    ligand_atoms = []
    protein_atoms = []
    ligand_resname = None
    
    with open(protein_path, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip().upper() if len(line) >= 78 else atom_name[0].upper()
                
                # Exclude water, ions, and hydrogens
                if resname in ['HOH', 'WAT']:
                    continue
                if element in ['H', 'D', 'T']:
                    continue
                if element not in heavy_atoms_set:
                    continue
                    
                # Use the first organic ligand found
                if ligand_resname is None:
                    ligand_resname = (resname, chain, resseq)
                if (resname, chain, resseq) == ligand_resname:
                    ligand_atoms.append((x, y, z))
                elif len(ligand_atoms) > 0:
                    break  # Only use the first ligand
                    
            elif line.startswith('ATOM'):
                # Collect protein atoms for cavity detection
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                protein_atoms.append((x, y, z))
    
    # Strategy 1: Use bound ligand if found
    if len(ligand_atoms) > 0 and len(ligand_atoms) <= 50:
        coords = np.array(ligand_atoms)
        center = coords.mean(axis=0)
        min_xyz = coords.min(axis=0)
        max_xyz = coords.max(axis=0)
        size = (max_xyz - min_xyz) + 8.0  # 8 Å padding
        logging.info(f'Strategy 1: Using bound ligand ({len(ligand_atoms)} atoms)')
        return tuple(center) + tuple(size)
    
    # Strategy 2: Cavity detection using protein geometry
    if len(protein_atoms) > 0:
        protein_coords = np.array(protein_atoms)
        
        # Find potential cavities by looking for regions with low atom density
        # Use a grid-based approach to find cavities
        cavities = find_protein_cavities(protein_coords)
        
        if cavities:
            # Use the largest cavity
            best_cavity = max(cavities, key=lambda c: c['volume'])
            center = best_cavity['center']
            size = best_cavity['size']
            logging.info(f'Strategy 2: Using largest detected cavity (volume: {best_cavity["volume"]:.1f} Ų)')
            return tuple(center) + tuple(size)
    
    # Strategy 3: Geometric center of protein as fallback
    if len(protein_atoms) > 0:
        protein_coords = np.array(protein_atoms)
        center = protein_coords.mean(axis=0)
        
        # Calculate conservative box size based on protein dimensions
        min_xyz = protein_coords.min(axis=0)
        max_xyz = protein_coords.max(axis=0)
        protein_span = max_xyz - min_xyz
        
        # Use 40% of protein span or minimum 20Å, maximum 30Å per dimension
        size = np.clip(protein_span * 0.4, 20.0, 30.0)
        
        logging.info(f'Strategy 3: Using protein geometric center with conservative box size')
        return tuple(center) + tuple(size)
    
    # Strategy 4: Last resort - default values
    logging.warning('Strategy 4: No protein atoms found, using default center and box size')
    return (0.0, 0.0, 0.0, 25.0, 25.0, 25.0)

def find_protein_cavities(protein_coords, grid_spacing=2.0, probe_radius=1.4):
    """
    Find potential binding cavities in protein using a grid-based approach.
    Returns list of cavities sorted by volume.
    """
    cavities = []
    
    # Define bounding box around protein
    min_xyz = protein_coords.min(axis=0) - 5.0  # 5Å padding
    max_xyz = protein_coords.max(axis=0) + 5.0
    
    # Create 3D grid
    x_range = np.arange(min_xyz[0], max_xyz[0], grid_spacing)
    y_range = np.arange(min_xyz[1], max_xyz[1], grid_spacing)
    z_range = np.arange(min_xyz[2], max_xyz[2], grid_spacing)
    
    # Find grid points that are in cavities (not too close to protein atoms)
    cavity_points = []
    
    for x in x_range:
        for y in y_range:
            for z in z_range:
                grid_point = np.array([x, y, z])
                
                # Calculate distance to nearest protein atom
                distances = np.linalg.norm(protein_coords - grid_point, axis=1)
                min_distance = distances.min()
                
                # Point is in cavity if it's not too close to protein atoms
                # but also not too far (we want interior cavities)
                if probe_radius + 1.0 < min_distance < 8.0:
                    cavity_points.append(grid_point)
    
    if not cavity_points:
        return cavities
    
    cavity_coords = np.array(cavity_points)
    
    # Cluster cavity points to find distinct cavities
    from sklearn.cluster import DBSCAN
    
    try:
        clustering = DBSCAN(eps=grid_spacing * 2, min_samples=3).fit(cavity_coords)
        labels = clustering.labels_
        
        # Process each cluster (cavity)
        for cluster_id in set(labels):
            if cluster_id == -1:  # Skip noise points
                continue
                
            cluster_points = cavity_coords[labels == cluster_id]
            if len(cluster_points) < 5:  # Skip very small cavities
                continue
            
            # Calculate cavity properties
            center = cluster_points.mean(axis=0)
            min_xyz = cluster_points.min(axis=0)
            max_xyz = cluster_points.max(axis=0)
            
            # Box size with padding
            size = (max_xyz - min_xyz) + 10.0  # 10Å padding for cavities
            
            # Ensure minimum box size
            size = np.maximum(size, 15.0)
            
            # Cavity volume estimation
            volume = np.prod(size)
            
            cavities.append({
                'center': center,
                'size': size,
                'volume': volume,
                'points': len(cluster_points)
            })
    
    except ImportError:
        # Fallback if sklearn not available: use simple geometric clustering
        logging.warning('sklearn not available, using simplified cavity detection')
        
        # Simple approach: find the most central cavity region
        if len(cavity_coords) >= 5:
            center = cavity_coords.mean(axis=0)
            # Find points within reasonable distance of center
            distances = np.linalg.norm(cavity_coords - center, axis=1)
            central_points = cavity_coords[distances <= 8.0]
            
            if len(central_points) >= 5:
                center = central_points.mean(axis=0)
                min_xyz = central_points.min(axis=0)
                max_xyz = central_points.max(axis=0)
                size = np.maximum((max_xyz - min_xyz) + 10.0, 15.0)
                
                cavities.append({
                    'center': center,
                    'size': size,
                    'volume': np.prod(size),
                    'points': len(central_points)
                })
    
    # Sort by volume (largest first)
    cavities.sort(key=lambda c: c['volume'], reverse=True)
    
    return cavities

def run_vina(protein, ligand, output_dir, box_params):
    """Run AutoDock Vina CLI."""
    vina_bin = 'vina'  # Assumes vina is in PATH
    vina_out = os.path.join(output_dir, VINA_OUT, 'vina_out.pdbqt')
    status = {'success': False, 'error': None, 'time': None}
    start = time.time()
    try:
        if not all(box_params):
            raise ValueError('All box parameters (center_x, center_y, center_z, size_x, size_y, size_z) are required for Vina.')
        cmd = [
            vina_bin,
            '--receptor', protein,
            '--ligand', ligand,
            '--center_x', str(box_params[0]),
            '--center_y', str(box_params[1]),
            '--center_z', str(box_params[2]),
            '--size_x', str(box_params[3]),
            '--size_y', str(box_params[4]),
            '--size_z', str(box_params[5]),
            '--out', vina_out
        ]
        logging.info(f'Running Vina: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)
        if not os.path.isfile(vina_out):
            raise RuntimeError('Vina did not produce output file.')
        status['success'] = True
        logging.info('Vina docking completed successfully.')
    except Exception as e:
        status['error'] = str(e)
        logging.error(f'[Vina] Docking failed: {e}', exc_info=True)
    status['time'] = round(time.time() - start, 2)
    return status

def run_gnina(protein, ligand, output_dir, use_gpu=False):
    """Run GNINA CLI for docking and scoring."""
    # Try different possible locations for gnina binary
    gnina_paths = [
        'gnina',  # In PATH
        './gnina/build/gnina',  # Local build
        os.path.expanduser('~/gnina/build/gnina'),  # User build
    ]
    
    gnina_bin = None
    for path in gnina_paths:
        if os.path.isfile(path) or (path == 'gnina' and shutil.which('gnina')):
            gnina_bin = path
            break
    
    if not gnina_bin:
        return {
            'success': False, 
            'error': 'GNINA binary not found. Please install GNINA or build it from source.',
            'time': 0.0
        }
    
    gnina_out = os.path.join(output_dir, GNINA_OUT, 'gnina_out.pdbqt')
    gnina_log = os.path.join(output_dir, GNINA_OUT, 'gnina.log')
    gnina_score = os.path.join(output_dir, GNINA_OUT, 'gnina_scores.txt')
    status = {'success': False, 'error': None, 'time': None}
    start = time.time()
    try:
        cmd = [
            gnina_bin,
            '--receptor', protein,
            '--ligand', ligand,
            '--out', gnina_out,
            '--log', gnina_log,
            '--score_only',
            '--scorefile', gnina_score
        ]
        if use_gpu:
            cmd.append('--gpu')
        logging.info(f'Running GNINA: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)
        if not os.path.isfile(gnina_out):
            raise RuntimeError('GNINA did not produce output file.')
        status['success'] = True
        logging.info('GNINA docking completed successfully.')
    except Exception as e:
        status['error'] = str(e)
        logging.error(f'[GNINA] Docking failed: {e}', exc_info=True)
    status['time'] = round(time.time() - start, 2)
    return status

def run_diffdock(protein, ligand, output_dir):
    """Run DiffDock in blind docking mode."""
    # Check for DiffDock inference script
    diffdock_script = os.path.join('DiffDock', 'inference.py')
    if not os.path.isfile(diffdock_script):
        return {
            'success': False,
            'error': f'DiffDock script not found at {diffdock_script}. Please clone DiffDock repository.',
            'time': 0.0
        }
    
    diffdock_out = os.path.join(output_dir, DIFFDOCK_OUT, 'diffdock_out.sdf')
    status = {'success': False, 'error': None, 'time': None}
    start = time.time()
    try:
        # Create a simple config for DiffDock
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write('protein_path,ligand_description,protein_sequence\n')
            f.write(f'{os.path.abspath(protein)},ligand,\n')
            csv_file = f.name
        
        cmd = [
            'python3', diffdock_script,
            '--protein_ligand_csv', csv_file,
            '--out_dir', os.path.join(output_dir, DIFFDOCK_OUT),
            '--inference_steps', '20',
            '--samples_per_complex', '10',
            '--batch_size', '10'
        ]
        logging.info(f'Running DiffDock: {" ".join(cmd)}')
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Clean up temp file
        os.unlink(csv_file)
        
        # Check for output (DiffDock creates complex directory structure)
        output_found = False
        diffdock_output_dir = os.path.join(output_dir, DIFFDOCK_OUT)
        for root, dirs, files in os.walk(diffdock_output_dir):
            if any(f.endswith('.sdf') for f in files):
                output_found = True
                break
        
        if not output_found:
            raise RuntimeError('DiffDock did not produce output files.')
        
        status['success'] = True
        logging.info('DiffDock docking completed successfully.')
    except Exception as e:
        status['error'] = str(e)
        logging.error(f'[DiffDock] Docking failed: {e}', exc_info=True)
        # Clean up temp file if it exists
        try:
            if 'csv_file' in locals():
                os.unlink(csv_file)
        except:
            pass
    status['time'] = round(time.time() - start, 2)
    return status

def run_batch_docking(protein_files, ligand_files, output_base_dir, use_gnina=False, use_diffdock=False, box_params=None):
    """
    Run docking in batch mode for multiple protein-ligand combinations.
    
    Args:
        protein_files: List of protein file paths
        ligand_files: List of ligand file paths  
        output_base_dir: Base directory for outputs
        use_gnina: Enable GNINA docking
        use_diffdock: Enable DiffDock docking
        box_params: Box parameters [cx, cy, cz, sx, sy, sz] or None for auto-detection
    
    Returns:
        Dictionary with batch results
    """
    batch_results = {}
    total_combinations = len(protein_files) * len(ligand_files)
    current_combo = 0
    
    logging.info(f"Starting batch docking: {len(protein_files)} proteins × {len(ligand_files)} ligands = {total_combinations} combinations")
    
    for protein_file in protein_files:
        protein_name = os.path.splitext(os.path.basename(protein_file))[0]
        
        for ligand_file in ligand_files:
            current_combo += 1
            ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
            combo_name = f"{protein_name}_{ligand_name}"
            
            logging.info(f"Processing combination {current_combo}/{total_combinations}: {combo_name}")
            
            # Create unique output directory for this combination
            combo_output_dir = os.path.join(output_base_dir, combo_name)
            ensure_output_dirs(combo_output_dir)
            
            # Determine box parameters for this protein
            if box_params is None:
                try:
                    cx, cy, cz, sx, sy, sz = extract_box_from_protein(protein_file)
                    current_box_params = [cx, cy, cz, sx, sy, sz]
                except Exception as e:
                    logging.error(f"Failed to extract box parameters for {protein_file}: {e}")
                    batch_results[combo_name] = {'error': f'Box extraction failed: {e}'}
                    continue
            else:
                current_box_params = box_params
            
            # Run docking for this combination
            backend_status = {}
            failed_runs = {}
            
            # Run Vina
            vina_status = run_vina(protein_file, ligand_file, combo_output_dir, current_box_params)
            backend_status['vina'] = vina_status
            if not vina_status['success']:
                failed_runs['vina'] = vina_status['error']
            
            # Optionally run GNINA
            if use_gnina:
                gnina_status = run_gnina(protein_file, ligand_file, combo_output_dir)
                backend_status['gnina'] = gnina_status
                if not gnina_status['success']:
                    failed_runs['gnina'] = gnina_status['error']
            
            # Optionally run DiffDock
            if use_diffdock:
                diffdock_status = run_diffdock(protein_file, ligand_file, combo_output_dir)
                backend_status['diffdock'] = diffdock_status
                if not diffdock_status['success']:
                    failed_runs['diffdock'] = diffdock_status['error']
            
            # Save results for this combination
            batch_results[combo_name] = {
                'protein': protein_file,
                'ligand': ligand_file,
                'output_dir': combo_output_dir,
                'backend_status': backend_status,
                'failed_runs': failed_runs
            }
            
            # Save failed runs if any
            if failed_runs:
                failed_json = os.path.join(combo_output_dir, LOGS_DIR, 'failed_runs.json')
                with open(failed_json, 'w') as f:
                    json.dump(failed_runs, f, indent=2)
    
    logging.info(f"Batch docking completed: {current_combo} combinations processed")
    return batch_results


# --- Main CLI ---
def main():
    parser = argparse.ArgumentParser(description='Docking Wrapper for Vina, GNINA, and DiffDock')
    parser.add_argument('--protein', required=True, help='Path to preprocessed .pdbqt protein file')
    parser.add_argument('--ligand', required=True, help='Path to preprocessed .pdbqt ligand file')
    parser.add_argument('--batch_proteins', nargs='+', help='Multiple protein files for batch processing')
    parser.add_argument('--batch_ligands', nargs='+', help='Multiple ligand files for batch processing')
    parser.add_argument('--use_gnina', action='store_true', help='Enable GNINA docking')
    parser.add_argument('--use_diffdock', action='store_true', help='Enable DiffDock docking')
    parser.add_argument('--output_dir', default='.', help='Base output directory (default: current directory)')
    parser.add_argument('--center_x', type=float, help='Center X coordinate for Vina')
    parser.add_argument('--center_y', type=float, help='Center Y coordinate for Vina')
    parser.add_argument('--center_z', type=float, help='Center Z coordinate for Vina')
    parser.add_argument('--size_x', type=float, help='Box size X for Vina')
    parser.add_argument('--size_y', type=float, help='Box size Y for Vina')
    parser.add_argument('--size_z', type=float, help='Box size Z for Vina')
    args = parser.parse_args()

    # Setup logging and output dirs
    ensure_output_dirs(args.output_dir)
    setup_logging()
    logging.info('Starting docking wrapper')

    # Determine if we're in batch mode
    if args.batch_proteins or args.batch_ligands:
        # Batch mode processing
        protein_files = args.batch_proteins if args.batch_proteins else [args.protein]
        ligand_files = args.batch_ligands if args.batch_ligands else [args.ligand]
        
        # Validate all batch files
        for protein_file in protein_files:
            validate_file(protein_file, 'Protein')
        for ligand_file in ligand_files:
            validate_file(ligand_file, 'Ligand')
        
        # Determine box parameters
        box_params = None
        if all(p is not None for p in [args.center_x, args.center_y, args.center_z, args.size_x, args.size_y, args.size_z]):
            box_params = [args.center_x, args.center_y, args.center_z, args.size_x, args.size_y, args.size_z]
            logging.info('Using provided box parameters for all combinations')
        else:
            logging.info('Box parameters will be auto-detected for each protein')
        
        # Run batch docking
        batch_results = run_batch_docking(
            protein_files, ligand_files, args.output_dir,
            use_gnina=args.use_gnina, use_diffdock=args.use_diffdock,
            box_params=box_params
        )
        
        # Print batch summary
        print(f"\n=== Batch Docking Summary ({len(batch_results)} combinations) ===")
        successful_combos = 0
        for combo_name, result in batch_results.items():
            if 'error' in result:
                print(f"{combo_name}: FAILED - {result['error']}")
            else:
                backend_successes = sum(1 for status in result['backend_status'].values() if status['success'])
                total_backends = len(result['backend_status'])
                if backend_successes > 0:
                    successful_combos += 1
                print(f"{combo_name}: {backend_successes}/{total_backends} backends successful")
        print(f"Overall: {successful_combos}/{len(batch_results)} combinations had successful runs")
        print("======================\n")
        
    else:
        # Single docking mode (original behavior)
        validate_file(args.protein, 'Protein')
        validate_file(args.ligand, 'Ligand')

        # Extract or validate Vina box params
        box_params = [args.center_x, args.center_y, args.center_z, args.size_x, args.size_y, args.size_z]
        if not all(p is not None for p in box_params):
            try:
                logging.info('Box parameters not provided. Attempting to extract from protein file...')
                cx, cy, cz, sx, sy, sz = extract_box_from_protein(args.protein)
                args.center_x, args.center_y, args.center_z = cx, cy, cz
                args.size_x, args.size_y, args.size_z = sx, sy, sz
                box_params = [cx, cy, cz, sx, sy, sz]
                logging.info(f'Extracted box center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), size: ({sx:.2f}, {sy:.2f}, {sz:.2f})')
            except Exception as e:
                logging.error(f'Failed to extract box parameters: {e}')
                print(f'ERROR: {e}')
                sys.exit(1)

        backend_status = {}
        failed_runs = {}

        # Only run Vina if all box params are provided
        if all(p is not None for p in [args.center_x, args.center_y, args.center_z, args.size_x, args.size_y, args.size_z]):
            vina_status = run_vina(args.protein, args.ligand, args.output_dir, [args.center_x, args.center_y, args.center_z, args.size_x, args.size_y, args.size_z])
            backend_status['vina'] = vina_status
            if not vina_status['success']:
                failed_runs['vina'] = vina_status['error']
        else:
            msg = 'Vina skipped: All box parameters (--center_x, --center_y, --center_z, --size_x, --size_y, --size_z) must be provided.'
            logging.warning(msg)
            backend_status['vina'] = {'success': False, 'error': msg, 'time': 0.0}

        # Optionally run GNINA
        if args.use_gnina:
            gnina_status = run_gnina(args.protein, args.ligand, args.output_dir)
            backend_status['gnina'] = gnina_status
            if not gnina_status['success']:
                failed_runs['gnina'] = gnina_status['error']

        # Optionally run DiffDock
        if args.use_diffdock:
            diffdock_status = run_diffdock(args.protein, args.ligand, args.output_dir)
            backend_status['diffdock'] = diffdock_status
            if not diffdock_status['success']:
                failed_runs['diffdock'] = diffdock_status['error']

        # Save failed runs if any
        if failed_runs:
            failed_json = os.path.join(args.output_dir, LOGS_DIR, 'failed_runs.json')
            with open(failed_json, 'w') as f:
                json.dump(failed_runs, f, indent=2)
            logging.info(f'Failed runs logged in {failed_json}')
        else:
            logging.info('All selected backends completed successfully.')

        # Print final summary to console
        print("\n=== Docking Summary ===")
        for backend, stat in backend_status.items():
            if stat['success']:
                print(f"{backend.upper()}: SUCCESS (time: {stat['time']}s)")
            else:
                print(f"{backend.upper()}: SKIPPED/FAILED (time: {stat['time']}s)")
                print(f"  Reason: {stat['error']}")
        print("======================\n")

if __name__ == '__main__':
    main() 