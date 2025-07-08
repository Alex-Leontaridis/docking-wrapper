#!/usr/bin/env python3
"""
UMol Docking Runner

Runs UMol docking for protein-ligand complexes.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import subprocess
import time
from pathlib import Path
from typing import Optional, Dict, Any

# Add the project root to the path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from utils.logging import setup_logging

def run_umol(protein: str, ligand: str, output: str, 
             num_poses: int = 10, exhaustiveness: int = 8,
             energy_range: int = 3, **kwargs) -> Dict[str, Any]:
    """
    Run UMol docking.
    
    Args:
        protein: Path to protein PDB file
        ligand: Path to ligand file (SDF, MOL, etc.)
        output: Path to output file
        num_poses: Number of poses to generate
        exhaustiveness: Exhaustiveness parameter
        energy_range: Energy range parameter
        **kwargs: Additional UMol parameters
    
    Returns:
        Dictionary with results
    """
    # Setup centralized logging
    logger = setup_logging('UMol')
    
    # Validate inputs
    if not os.path.exists(protein):
        logger.error(f"Protein file not found: {protein}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'Protein file not found: {protein}'}
    
    if not os.path.exists(ligand):
        logger.error(f"Ligand file not found: {ligand}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'Ligand file not found: {ligand}'}
    
    # Build UMol command
    cmd = [
        'umol',
        '--receptor', protein,
        '--ligand', ligand,
        '--out', output,
        '--num_poses', str(num_poses),
        '--exhaustiveness', str(exhaustiveness),
        '--energy_range', str(energy_range)
    ]
    
    # Add additional parameters
    for key, value in kwargs.items():
        if value is not None:
            cmd.extend([f'--{key}', str(value)])
    
    # Run UMol
    start_time = time.time()
    try:
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        
        if result.returncode == 0:
            elapsed = time.time() - start_time
            logger.info(f"UMol completed in {elapsed:.2f}s")
            logger.info(f"Output written to {output}")
            return {
                'success': True,
                'output_file': output,
                'elapsed_time': elapsed,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
        else:
            logger.error(f"UMol failed: {result.stderr}", error_code='UMOL_FAILED')
            return {
                'success': False,
                'error': result.stderr,
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
    except subprocess.TimeoutExpired:
        logger.error("UMol timed out after 1 hour", error_code='TIMEOUT_ERROR')
        return {'success': False, 'error': 'UMol timed out after 1 hour'}
    except FileNotFoundError:
        logger.error("UMol not found. Install UMol with: python install_real_tools.py", error_code='ENGINE_NOT_FOUND')
        return {'success': False, 'error': 'UMol not found. Install UMol with: python install_real_tools.py'}
    except Exception as e:
        logger.error(f"UMol failed: {e}", error_code='SYSTEM_ERROR')
        return {'success': False, 'error': str(e)}

def main():
    """Main function for command line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run UMol docking')
    parser.add_argument('--protein', required=True, help='Protein PDB file')
    parser.add_argument('--ligand', required=True, help='Ligand file')
    parser.add_argument('--output', required=True, help='Output file')
    parser.add_argument('--num_poses', type=int, default=10, help='Number of poses')
    parser.add_argument('--exhaustiveness', type=int, default=8, help='Exhaustiveness')
    parser.add_argument('--energy_range', type=int, default=3, help='Energy range')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging('UMol')
    
    # Run UMol
    result = run_umol(
        protein=args.protein,
        ligand=args.ligand,
        output=args.output,
        num_poses=args.num_poses,
        exhaustiveness=args.exhaustiveness,
        energy_range=args.energy_range
    )
    
    if result['success']:
        logger.info("UMol completed successfully")
        sys.exit(0)
    else:
        logger.error(f"UMol failed: {result['error']}")
        sys.exit(1)

if __name__ == '__main__':
    main() 