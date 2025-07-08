#!/usr/bin/env python3
import argparse
import os
import sys
import logging
import subprocess
import json
import time
import numpy as np
import shutil
import platform
from pathlib import Path

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Setup logging first
from utils.logging import setup_logging as setup_docking_logging
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir

# Import config after logging is set up
from config import config

# --- Constants for output directories ---
VINA_OUT = 'vina_output'
GNINA_OUT = 'gnina_output'
DIFFDOCK_OUT = 'diffdock_output'
LOGS_DIR = 'logs'
LOG_FILE = os.path.join(LOGS_DIR, 'docking_run.log')

# --- Setup logging ---
from utils.logging import setup_logging as setup_docking_logging

# Create logs directory if it doesn't exist
ensure_dir(LOGS_DIR)

# Configure logging
logger = setup_docking_logging(__name__)

# --- Ensure output directories exist ---
def ensure_output_dirs(base_dir):
    for d in [VINA_OUT, GNINA_OUT, DIFFDOCK_OUT, LOGS_DIR]:
        os.makedirs(os.path.join(base_dir, d), exist_ok=True)

# --- Validate input files ---
def validate_file(path, desc):
    if not path or not os.path.isfile(path):
        logger.error(f"{desc} file '{path}' does not exist or is not a file.")
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
        logger.info(f'Strategy 1: Using bound ligand ({len(ligand_atoms)} atoms)')
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
            logger.info(f'Strategy 2: Using largest detected cavity (volume: {best_cavity["volume"]:.1f} A^3)')
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
        
        logger.info(f'Strategy 3: Using protein geometric center with conservative box size')
        return tuple(center) + tuple(size)
    
    # Strategy 4: Last resort - default values
    logger.warning('Strategy 4: No protein atoms found, using default center and box size')
    return (0.0, 0.0, 0.0, 25.0, 25.0, 25.0)

# Cluster cavity points to find distinct cavities
try:
    from sklearn.cluster import DBSCAN
    SKLEARN_AVAILABLE = True
except ImportError:
    logger.error("scikit-learn is required for cavity clustering but is not installed. Install it with: pip install scikit-learn")
    SKLEARN_AVAILABLE = False

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
    
    # Replace the direct use of DBSCAN with a check for SKLEARN_AVAILABLE
    if SKLEARN_AVAILABLE:
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
        except Exception as e:
            logger.warning(f'sklearn DBSCAN failed, using simplified cavity detection: {e}')
            # fallback below
    else:
        # Fallback if sklearn not available: use simple geometric clustering
        logger.warning('sklearn not available, using simplified cavity detection')
        
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

def find_binary(binary_name, env_var=None, config_path=None):
    """
    Search for a binary in the following order:
    1. Current working directory (highest priority)
    2. Config path (if provided)
    3. Environment variable (if provided)
    4. PATH
    5. For Windows: Check WSL if binary is Linux-only (like GNINA)
    Returns the path to the binary or None if not found.
    """
    import shutil
    import os
    import platform
    
    cwd = os.getcwd()
    # 1. Current working directory (highest priority for local binaries)
    local_path = os.path.join(cwd, binary_name)
    if os.path.isfile(local_path):
        logger.info(f"Found {binary_name} in current directory: {local_path}")
        return local_path
    
    # Check for .bat files on Windows (for dummy scripts)
    if platform.system().lower() == 'windows':
        base_name = os.path.splitext(binary_name)[0]  # Remove .exe if present
        bat_path = os.path.join(cwd, f'{base_name}.bat')
        if os.path.isfile(bat_path):
            logger.info(f"Found {base_name}.bat in current directory: {bat_path}")
            return bat_path
    
    # Special case: gnina.exe or gnina in current directory
    if binary_name == 'gnina':
        for ext in ('', '.exe', '.bat'):
            gnina_local = os.path.join(cwd, f'gnina{ext}')
            if os.path.isfile(gnina_local):
                logger.info(f"Found gnina in current directory: {gnina_local}")
                return gnina_local
    # 2. Config path
    if config_path and os.path.isfile(config_path):
        logger.info(f"Found {binary_name} via config: {config_path}")
        return config_path
    # 3. Environment variable
    if env_var and os.environ.get(env_var):
        env_path = os.environ[env_var]
        if os.path.isfile(env_path):
            logger.info(f"Found {binary_name} via environment variable {env_var}: {env_path}")
            return env_path
    # 4. PATH
    which_path = shutil.which(binary_name)
    if which_path:
        logger.info(f"Found {binary_name} in PATH: {which_path}")
        return which_path
    # 5. For Windows: Check WSL for Linux-only binaries
    if platform.system().lower() == 'windows':
        try:
            result = subprocess.run(['wsl', 'which', binary_name], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                wsl_path = result.stdout.strip()
                if wsl_path:
                    logger.info(f"Found {binary_name} in WSL: {wsl_path}")
                    return f"wsl:{wsl_path}"  # Mark as WSL binary
        except:
            pass
    logger.warning(f"Binary {binary_name} not found in any location")
    return None

# Platform-specific engine availability
system = platform.system().lower()

# GNINA: Linux only, but allow on Windows/Mac if WSL is available
if system in ['windows', 'darwin']:  # Windows or Mac
    # Check if WSL is available
    try:
        result = subprocess.run(['wsl', '--list', '--quiet'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and result.stdout.strip():
            logger.info(f"WSL detected on {system.capitalize()}. GNINA will be available via WSL.")
            # Keep GNINA enabled if WSL is available
        else:
            # Only disable GNINA if config has the attribute
            if hasattr(config, 'engines') and 'gnina' in config.engines:
                config.engines["gnina"]["enabled"] = False
            logger.info(f"GNINA is not available on {system.capitalize()}. GNINA features will be disabled.")
    except:
        # Only disable GNINA if config has the attribute
        if hasattr(config, 'engines') and 'gnina' in config.engines:
            config.engines["gnina"]["enabled"] = False
        logger.info(f"GNINA is not available on {system.capitalize()}. GNINA features will be disabled.")

# UMol: Linux only  
if system in ['windows', 'darwin']:  # Windows or Mac
    # Only disable UMol if config has the attribute
    if hasattr(config, 'ml_models') and 'umol' in config.ml_models:
        config.ml_models["umol"]["enabled"] = False
    logger.info(f"UMol is not available on {system.capitalize()}. UMol features will be disabled.")
elif system == 'linux':
    logger.info("UMol is available on Linux. UMol features will be enabled if configured.")

# Log platform detection
logger.info(f"Detected platform: {system.capitalize()}")
logger.info(f"Available engines: Vina (all platforms), GNINA (Linux/WSL), DiffDock (all platforms)")
logger.info(f"Available ML models: EquiBind (all platforms), NeuralPLexer (all platforms), UMol (Linux only)")

def run_vina(protein, ligand, output_dir, box_params):
    """Run AutoDock Vina CLI."""
    vina_bin = find_binary('vina', env_var='VINA_PATH', config_path=getattr(config, 'vina_path', None))
    
    if not vina_bin:
        return {
            'success': False,
            'error': 'Vina binary not found. Please ensure vina.exe is in the current directory or set VINA_PATH.',
            'time': 0.0
        }
    
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
        logger.info(f'Running Vina: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)
        if not os.path.isfile(vina_out):
            raise RuntimeError('Vina did not produce output file.')
        status['success'] = True
        logger.info('Vina docking completed successfully.')
    except Exception as e:
        status['error'] = str(e)
        logger.error(f'[Vina] Docking failed: {e}', exc_info=True)
    status['time'] = round(time.time() - start, 2)
    return status

def load_backend_config():
    """Load backend configuration from installation."""
    config_file = Path.cwd() / "backend_config.json"
    if config_file.exists():
        try:
            with open(config_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Failed to load backend config: {e}")
    return None

def find_gnina_binary():
    """Find GNINA binary using environment variables and platform detection."""
    # Check environment variable first
    gnina_path = os.environ.get('GNINA_PATH')
    if gnina_path and os.path.isfile(gnina_path):
        logger.info(f"Found GNINA binary from environment: {gnina_path}")
        return gnina_path
    
    # Check if gnina is in PATH
    if shutil.which('gnina'):
        logger.info("Found GNINA binary in PATH")
        return 'gnina'
    
    # Platform-specific common locations
    system = platform.system().lower()
    if system == 'windows':
        possible_paths = [
            'gnina.exe',  # In PATH
            './gnina.exe',  # Current directory
            './gnina.bat',  # Current directory
            './gnina',  # Current directory
        ]
    elif system == 'darwin':  # macOS
        possible_paths = [
            '/usr/local/bin/gnina',
            '/opt/homebrew/bin/gnina',
            os.path.expanduser('~/gnina'),
            './gnina',
        ]
    else:  # Linux
        possible_paths = [
            '/usr/local/bin/gnina',
            '/usr/bin/gnina',
            '/opt/conda/envs/docking/bin/gnina',
            os.path.expanduser('~/gnina'),
            './gnina',
        ]
    
    for path in possible_paths:
        if shutil.which(path):
            logger.info(f"Found GNINA binary at: {path}")
            return path
    
    logger.warning("GNINA binary not found. Please set GNINA_PATH environment variable or install GNINA.")
    return None

def _is_dummy_diffdock_script(script_path):
    """
    Check if a DiffDock script is the dummy placeholder with comprehensive detection.
    
    Args:
        script_path: Path to the DiffDock script
        
    Returns:
        bool: True if it's a dummy script, False if it's a real script
    """
    try:
        # Check file size first (dummy scripts are usually small)
        file_size = os.path.getsize(script_path)
        if file_size < 1000:  # Less than 1KB is suspicious
            logger.warning(f"DiffDock script is very small ({file_size} bytes), likely a dummy")
            return True
        
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # Check for dummy script indicators
        dummy_indicators = [
            'DIFFDOCK NOT INSTALLED',
            'This is a dummy script',
            'DiffDock not installed',
            'placeholder script',
            'DUMMY_SCRIPT',
            'raise NotImplementedError',
            'print("DiffDock not available")',
            'sys.exit(1)  # DiffDock not installed',
            'return False  # DiffDock not available'
        ]
        
        # Check if any dummy indicators are present
        for indicator in dummy_indicators:
            if indicator.lower() in content.lower():
                logger.warning(f"Found dummy script indicator: '{indicator}' in {script_path}")
                return True
        
        # Check for real DiffDock indicators
        real_indicators = [
            'import torch',
            'import torch_geometric',
            'from diffdock',
            'class DiffDock',
            'def inference',
            'def predict',
            'torch.load',
            'model.eval()',
            'protein_ligand_csv',
            'out_dir'
        ]
        
        # Count real indicators
        real_count = sum(1 for indicator in real_indicators if indicator.lower() in content.lower())
        
        # If we have very few real indicators, it's likely a dummy
        if real_count < 3:
            logger.warning(f"DiffDock script has only {real_count} real indicators, likely a dummy")
            return True
        
        # Check for proper imports and structure
        if 'import' not in content or 'def ' not in content:
            logger.warning("DiffDock script lacks proper Python structure")
            return True
        
        # If we get here, it looks like a real script
        return False
        
    except Exception as e:
        logger.warning(f"Could not analyze DiffDock script {script_path}: {e}")
        # If we can't read the file, assume it's not a dummy to be safe
        return False

def find_diffdock_script():
    """
    Find DiffDock inference script using environment variables and platform detection.
    Returns:
        str or None: Path to the DiffDock script, or None if not found or dummy detected
    """
    logger = setup_docking_logging('DiffDock')
    # Check environment variable first
    diffdock_path = os.environ.get('DIFFDOCK_PATH')
    if diffdock_path:
        script_path = os.path.join(diffdock_path, 'inference.py')
        logger.info(f"Checking DIFFDOCK_PATH: {script_path}")
        if os.path.isfile(script_path):
            if _is_dummy_diffdock_script(script_path):
                logger.warning(f"Found DiffDock dummy script from environment: {script_path}")
                logger.warning("This is a placeholder script. Please install the real DiffDock.")
                logger.info("Install DiffDock with: git clone https://github.com/gcorso/DiffDock.git")
                return None
            else:
                logger.info(f"Found DiffDock script from environment: {script_path}")
                return script_path
        else:
            logger.warning(f"DIFFDOCK_PATH set to {diffdock_path} but inference.py not found")
    # Check current working directory for DiffDock directory
    cwd = os.getcwd()
    local_diffdock = os.path.join(cwd, 'DiffDock')
    logger.info(f"Checking local DiffDock directory: {local_diffdock}")
    if os.path.isdir(local_diffdock):
        script_path = os.path.join(local_diffdock, 'inference.py')
        logger.info(f"Checking for inference.py: {script_path}")
        if os.path.isfile(script_path):
            if _is_dummy_diffdock_script(script_path):
                logger.warning(f"Found DiffDock dummy script at: {script_path}")
                logger.warning("This is a placeholder script. Please install the real DiffDock.")
                logger.info("Install DiffDock with: git clone https://github.com/gcorso/DiffDock.git")
                return None
            else:
                logger.info(f"Found DiffDock script at: {script_path}")
                return script_path
    # Also check lowercase 'diffdock' directory
    local_diffdock_lower = os.path.join(cwd, 'diffdock')
    logger.info(f"Checking local diffdock directory: {local_diffdock_lower}")
    if os.path.isdir(local_diffdock_lower):
        script_path = os.path.join(local_diffdock_lower, 'inference.py')
        logger.info(f"Checking for inference.py: {script_path}")
        if os.path.isfile(script_path):
            if _is_dummy_diffdock_script(script_path):
                logger.warning(f"Found DiffDock dummy script at: {script_path}")
                logger.warning("This is a placeholder script. Please install the real DiffDock.")
                logger.info("Install DiffDock with: git clone https://github.com/gcorso/DiffDock.git")
                return None
            else:
                logger.info(f"Found DiffDock script at: {script_path}")
                return script_path
    # Platform-specific common locations (without hardcoded paths)
    system = platform.system().lower()
    if system == 'windows':
        possible_paths = [
            os.path.join(os.path.expanduser('~'), 'DiffDock', 'inference.py'),  # User home
            os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'DiffDock', 'inference.py'),  # Program Files
        ]
    else:  # Linux/macOS
        possible_paths = [
            os.path.join(cwd, 'DiffDock', 'inference.py'),  # Current directory
            os.path.join(cwd, 'diffdock', 'inference.py'),  # Current directory
            os.path.join(os.path.expanduser('~'), 'DiffDock', 'inference.py'),  # User home
            '/opt/DiffDock/inference.py',  # Docker installation
            '/usr/local/DiffDock/inference.py',  # System installation
            '/usr/share/DiffDock/inference.py',  # Alternative system location
        ]
    for path in possible_paths:
        logger.info(f"Checking possible path: {path}")
        if os.path.isfile(path):
            if _is_dummy_diffdock_script(path):
                logger.warning(f"Found DiffDock dummy script at: {path}")
                logger.warning("This is a placeholder script. Please install the real DiffDock.")
                logger.info("Install DiffDock with: git clone https://github.com/gcorso/DiffDock.git")
                return None
            else:
                logger.info(f"Found DiffDock script at: {path}")
                return path
    logger.warning("DiffDock inference script not found. Please set DIFFDOCK_PATH environment variable or install DiffDock.")
    logger.info("Install DiffDock with: git clone https://github.com/gcorso/DiffDock.git")
    return None

def run_gnina(protein, ligand, output_dir, use_gpu=False):
    """Run GNINA CLI for docking and scoring."""
    gnina_bin = find_binary('gnina', env_var='GNINA_PATH', config_path=getattr(config, 'gnina_path', None))
    
    if not gnina_bin:
        return {
            'success': False,
            'error': 'GNINA binary not found. Please ensure gnina is in the current directory or set GNINA_PATH.',
            'time': 0.0
        }
    
    # Check if GNINA is a WSL binary
    if gnina_bin.startswith('wsl:'):
        # Extract the WSL path
        wsl_path = gnina_bin[4:]  # Remove 'wsl:' prefix
        gnina_bin = ['wsl', wsl_path]
    else:
        # Check if this is a batch file on Windows
        if platform.system().lower() == 'windows' and gnina_bin.endswith('.bat'):
            # It's a Windows batch file, run directly
            gnina_bin = [gnina_bin]
        else:
            # Check if this is a Linux binary on Windows
            if platform.system().lower() == 'windows':
                try:
                    # Try to execute the binary to see if it's a valid Windows executable
                    result = subprocess.run([gnina_bin, '--version'], 
                                          capture_output=True, text=True, timeout=10)
                    if result.returncode != 0:
                        # If it fails, it might be a Linux binary - try WSL
                        logger.info(f"GNINA appears to be a Linux binary. Attempting to run via WSL.")
                        gnina_bin = ['wsl', gnina_bin]
                except (subprocess.TimeoutExpired, OSError, FileNotFoundError):
                    # If execution fails, try WSL
                    logger.info(f"GNINA execution failed. Attempting to run via WSL.")
                    gnina_bin = ['wsl', gnina_bin]
            else:
                gnina_bin = [gnina_bin]
    
    gnina_out_dir = os.path.join(output_dir, GNINA_OUT)
    ensure_dir(gnina_out_dir)
    
    status = {'success': False, 'error': None, 'time': None}
    start = time.time()
    
    try:
        # Build GNINA command
        cmd = gnina_bin + [
            '--receptor', protein,
            '--ligand', ligand,
            '--out', os.path.join(gnina_out_dir, 'gnina_out.sdf'),
            '--num_modes', '9',
            '--exhaustiveness', '8'
        ]
        
        if use_gpu:
            cmd.append('--gpu')
        
        logger.info(f'Running GNINA: {" ".join(cmd)}')
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Check if output was produced
        output_files = list(Path(gnina_out_dir).glob('*.sdf'))
        if not output_files:
            raise RuntimeError('GNINA did not produce any SDF output files.')
        
        status['success'] = True
        logger.info(f'GNINA completed successfully. Generated {len(output_files)} poses.')
        
    except subprocess.CalledProcessError as e:
        status['error'] = f"GNINA failed with exit code {e.returncode}: {e.stderr}"
        logger.error(f'[GNINA] Docking failed: {status["error"]}')
    except Exception as e:
        status['error'] = str(e)
        logger.error(f'[GNINA] Docking failed: {e}', exc_info=True)
    
    status['time'] = round(time.time() - start, 2)
    return status

def run_diffdock(protein, ligand, output_dir):
    """Run DiffDock for pose prediction."""
    diffdock_script = find_binary('inference.py', env_var='DIFFDOCK_PATH', config_path=getattr(config, 'diffdock_path', None))
    if not os.path.isfile(diffdock_script):
        return {
            'success': False,
            'error': 'DiffDock script not found at expected location. Please set DIFFDOCK_PATH or install DiffDock.',
            'time': 0.0
        }
    diffdock_out_dir = os.path.join(output_dir, DIFFDOCK_OUT)
    ensure_dir(diffdock_out_dir)
    status = {'success': False, 'error': None, 'time': None}
    start = time.time()
    
    try:
        # Create temporary CSV file for DiffDock input
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            # DiffDock expects: protein_path,ligand_description,complex_name
            f.write('protein_path,ligand_description,complex_name\n')
            
            # Convert ligand to SMILES if it's a file
            if ligand.endswith('.pdbqt'):
                # Extract SMILES from PDBQT (simplified)
                ligand_name = os.path.splitext(os.path.basename(ligand))[0]
                # For now, use a placeholder SMILES - in production you'd convert properly
                smiles = "CCO"  # Placeholder
            else:
                smiles = ligand  # Assume it's already SMILES
                ligand_name = "ligand"
            
            f.write(f'{protein},{smiles},{ligand_name}\n')
            csv_file = f.name
        
        # Run DiffDock
        cmd = [
            'python3', diffdock_script,
            '--protein_ligand_csv', csv_file,
            '--out_dir', diffdock_out_dir,
            '--inference_steps', '20',
            '--samples_per_complex', '10',
            '--batch_size', '10'
        ]
        
        logger.info(f'Running DiffDock: {" ".join(cmd)}')
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Clean up temporary file
        os.unlink(csv_file)
        
        # Check if output was produced
        output_files = list(Path(diffdock_out_dir).glob('**/*.sdf'))
        if not output_files:
            raise RuntimeError('DiffDock did not produce any SDF output files.')
        
        status['success'] = True
        logger.info(f'DiffDock completed successfully. Generated {len(output_files)} poses.')
        
    except subprocess.CalledProcessError as e:
        status['error'] = f"DiffDock failed with exit code {e.returncode}: {e.stderr}"
        logger.error(f'[DiffDock] Docking failed: {status["error"]}')
        if 'csv_file' in locals():
            try:
                os.unlink(csv_file)
            except:
                pass
    except Exception as e:
        status['error'] = str(e)
        logger.error(f'[DiffDock] Docking failed: {e}', exc_info=True)
        if 'csv_file' in locals():
            try:
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
    
    logger.info(f"Starting batch docking: {len(protein_files)} proteins × {len(ligand_files)} ligands = {total_combinations} combinations")
    
    for protein_file in protein_files:
        protein_name = os.path.splitext(os.path.basename(protein_file))[0]
        
        for ligand_file in ligand_files:
            current_combo += 1
            ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
            combo_name = f"{protein_name}_{ligand_name}"
            
            logger.info(f"Processing combination {current_combo}/{total_combinations}: {combo_name}")
            
            # Create unique output directory for this combination
            combo_output_dir = os.path.join(output_base_dir, combo_name)
            ensure_output_dirs(combo_output_dir)
            
            # Determine box parameters for this protein
            if box_params is None:
                try:
                    cx, cy, cz, sx, sy, sz = extract_box_from_protein(protein_file)
                    current_box_params = [cx, cy, cz, sx, sy, sz]
                except Exception as e:
                    logger.error(f"Failed to extract box parameters for {protein_file}: {e}")
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
                batch_results[combo_name] = {
                    'protein': protein_file,
                    'ligand': ligand_file,
                    'output_dir': combo_output_dir,
                    'backend_status': backend_status,
                    'failed_runs': failed_runs
                }
                continue
            
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
    
    logger.info(f"Batch docking completed: {current_combo} combinations processed")
    return batch_results


# --- Main CLI ---
def main():
    # Check for required external binaries
    required_binaries = ["vina"]
    if '--use_gnina' in sys.argv or '--use_diffdock' in sys.argv:
        if '--use_gnina' in sys.argv:
            required_binaries.append("gnina")
        if '--use_diffdock' in sys.argv:
            required_binaries.append("diffdock")
    missing = [b for b in required_binaries if shutil.which(b) is None]
    if missing:
        print(f"ERROR: Missing required external binaries: {', '.join(missing)}")
        print("Please install them and ensure they are in your PATH.")
        sys.exit(1)

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
    # Logging is already set up at module level
    logger.info('Starting docking wrapper')

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
            logger.info('Using provided box parameters for all combinations')
        else:
            logger.info('Box parameters will be auto-detected for each protein')
        
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
                logger.info('Box parameters not provided. Attempting to extract from protein file...')
                cx, cy, cz, sx, sy, sz = extract_box_from_protein(args.protein)
                args.center_x, args.center_y, args.center_z = cx, cy, cz
                args.size_x, args.size_y, args.size_z = sx, sy, sz
                box_params = [cx, cy, cz, sx, sy, sz]
                logger.info(f'Extracted box center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), size: ({sx:.2f}, {sy:.2f}, {sz:.2f})')
            except Exception as e:
                logger.error(f'Failed to extract box parameters: {e}')
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
                batch_results = {args.protein: {'error': vina_status['error']}}
                print(f"Vina failed: {vina_status['error']}")
                sys.exit(1)
        else:
            msg = 'Vina skipped: All box parameters (--center_x, --center_y, --center_z, --size_x, --size_y, --size_z) must be provided.'
            logger.warning(msg)
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
            logger.info(f'Failed runs logged in {failed_json}')
        else:
            logger.info('All selected backends completed successfully.')

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