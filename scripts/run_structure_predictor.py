#!/usr/bin/env python3
"""
Structure Prediction Runner

Runs various structure prediction tools (ColabFold, OpenFold, ESMFold).
"""

import os
import sys
import subprocess
import time
from pathlib import Path
from typing import Optional, Dict, Any

# Add the project root to the path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from utils.logging import setup_logging

def run_colabfold(fasta: str, output_dir: str, **kwargs) -> Dict[str, Any]:
    """Run ColabFold structure prediction."""
    logger = setup_logging('ColabFold')
    
    # Validate input
    if not os.path.exists(fasta):
        logger.error(f"FASTA file not found: {fasta}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'FASTA file not found: {fasta}'}
    
    # Build command
    cmd = [
        'colabfold_batch',
        '--type', 'auto',
        '--num-recycle', '3',
        '--num-models', '5',
        '--amber',
        '--templates',
        '--use-gpu-relax',
        fasta,
        output_dir
    ]
    
    # Add additional parameters
    for key, value in kwargs.items():
        if value is not None:
            cmd.extend([f'--{key}', str(value)])
    
    # Run ColabFold
    start_time = time.time()
    try:
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        if result.returncode == 0:
            elapsed = time.time() - start_time
            logger.info(f"ColabFold completed in {elapsed:.2f}s")
            return {
                'success': True,
                'output_dir': output_dir,
                'elapsed_time': elapsed,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
        else:
            logger.error(f"ColabFold failed: {result.stderr}", error_code='COLABFOLD_FAILED')
            return {
                'success': False,
                'error': result.stderr,
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
    except subprocess.TimeoutExpired:
        logger.error("ColabFold timed out after 2 hours", error_code='TIMEOUT_ERROR')
        return {'success': False, 'error': 'ColabFold timed out after 2 hours'}
    except FileNotFoundError:
        logger.error("ColabFold not found. Install with: pip install colabfold", error_code='ENGINE_NOT_FOUND')
        return {'success': False, 'error': 'ColabFold not found. Install with: pip install colabfold'}
    except Exception as e:
        logger.error(f"ColabFold failed: {e}", error_code='SYSTEM_ERROR')
        return {'success': False, 'error': str(e)}

def run_openfold(fasta: str, output_dir: str, **kwargs) -> Dict[str, Any]:
    """Run OpenFold structure prediction."""
    logger = setup_logging('OpenFold')
    
    # Validate input
    if not os.path.exists(fasta):
        logger.error(f"FASTA file not found: {fasta}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'FASTA file not found: {fasta}'}
    
    # Build command
    cmd = [
        'openfold',
        '--fasta_paths', fasta,
        '--output_dir', output_dir,
        '--model_device', 'cuda:0',
        '--config_preset', 'model_1_ptm'
    ]
    
    # Add additional parameters
    for key, value in kwargs.items():
        if value is not None:
            cmd.extend([f'--{key}', str(value)])
    
    # Run OpenFold
    start_time = time.time()
    try:
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        if result.returncode == 0:
            elapsed = time.time() - start_time
            logger.info(f"OpenFold completed in {elapsed:.2f}s")
            return {
                'success': True,
                'output_dir': output_dir,
                'elapsed_time': elapsed,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
        else:
            logger.error(f"OpenFold failed: {result.stderr}", error_code='OPENFOLD_FAILED')
            return {
                'success': False,
                'error': result.stderr,
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
    except subprocess.TimeoutExpired:
        logger.error("OpenFold timed out after 2 hours", error_code='TIMEOUT_ERROR')
        return {'success': False, 'error': 'OpenFold timed out after 2 hours'}
    except FileNotFoundError:
        logger.error("OpenFold not found. Install with: pip install openfold", error_code='ENGINE_NOT_FOUND')
        return {'success': False, 'error': 'OpenFold not found. Install with: pip install openfold'}
    except Exception as e:
        logger.error(f"OpenFold failed: {e}", error_code='SYSTEM_ERROR')
        return {'success': False, 'error': str(e)}

def run_esmfold(fasta: str, output_dir: str, **kwargs) -> Dict[str, Any]:
    """Run ESMFold structure prediction."""
    logger = setup_logging('ESMFold')
    
    # Validate input
    if not os.path.exists(fasta):
        logger.error(f"FASTA file not found: {fasta}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'FASTA file not found: {fasta}'}
    
    # Build command
    cmd = [
        'esm-fold',
        '--fasta', fasta,
        '--output-dir', output_dir,
        '--num-recycles', '4'
    ]
    
    # Add additional parameters
    for key, value in kwargs.items():
        if value is not None:
            cmd.extend([f'--{key}', str(value)])
    
    # Run ESMFold
    start_time = time.time()
    try:
        logger.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        if result.returncode == 0:
            elapsed = time.time() - start_time
            logger.info(f"ESMFold completed in {elapsed:.2f}s")
            return {
                'success': True,
                'output_dir': output_dir,
                'elapsed_time': elapsed,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
        else:
            logger.error(f"ESMFold failed: {result.stderr}", error_code='ESMFOLD_FAILED')
            return {
                'success': False,
                'error': result.stderr,
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
    except subprocess.TimeoutExpired:
        logger.error("ESMFold timed out after 2 hours", error_code='TIMEOUT_ERROR')
        return {'success': False, 'error': 'ESMFold timed out after 2 hours'}
    except FileNotFoundError:
        logger.error("ESMFold not found. Install with: pip install esm", error_code='ENGINE_NOT_FOUND')
        return {'success': False, 'error': 'ESMFold not found. Install with: pip install esm'}
    except Exception as e:
        logger.error(f"ESMFold failed: {e}", error_code='SYSTEM_ERROR')
        return {'success': False, 'error': str(e)}

def predict_structure(fasta: str, output_dir: str, model: str = 'auto', **kwargs) -> Dict[str, Any]:
    """
    Predict protein structure using the specified model.
    
    Args:
        fasta: Path to FASTA file
        output_dir: Output directory
        model: Model to use ('colabfold', 'openfold', 'esmfold', 'auto')
        **kwargs: Additional parameters
    
    Returns:
        Dictionary with results
    """
    logger = setup_logging('StructurePredictor')
    
    # Validate input
    if not os.path.exists(fasta):
        logger.error(f"FASTA file not found: {fasta}", error_code='FILE_NOT_FOUND')
        return {'success': False, 'error': f'FASTA file not found: {fasta}'}
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check for existing structure
    pdb_path = os.path.join(output_dir, "structure.pdb")
    log_path = os.path.join(output_dir, "prediction.log")
    
    if os.path.exists(pdb_path):
        logger.info(f"Structure already exists: {pdb_path}")
        logger.info(f"Log: {log_path}")
        return {
            'success': True,
            'output_file': pdb_path,
            'cached': True,
            'message': 'Structure already exists'
        }
    
    # Try different models if auto mode
    models_to_try = []
    if model == 'auto':
        models_to_try = ['colabfold', 'openfold', 'esmfold']
    else:
        models_to_try = [model]
    
    for model_name in models_to_try:
        logger.info(f"Trying {model_name}...")
        
        if model_name == 'colabfold':
            result = run_colabfold(fasta, output_dir, **kwargs)
        elif model_name == 'openfold':
            result = run_openfold(fasta, output_dir, **kwargs)
        elif model_name == 'esmfold':
            result = run_esmfold(fasta, output_dir, **kwargs)
        else:
            logger.error(f"Unknown model: {model_name}", error_code='INVALID_VALUE')
            continue
        
        if result['success']:
            logger.info(f"{model_name} prediction complete. Output: {pdb_path}")
            return result
        else:
            logger.warning(f"{model_name} failed: {result['error']}")
    
    # All models failed
    logger.error(f"All structure predictors failed. See log: {log_path}", error_code='ALL_MODELS_FAILED')
    logger.info("Try installing tools with: python install_real_tools.py")
    return {'success': False, 'error': 'All structure predictors failed'}

def main():
    """Main function for command line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Run structure prediction')
    parser.add_argument('--fasta', required=True, help='FASTA file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--model', default='auto', choices=['colabfold', 'openfold', 'esmfold', 'auto'], 
                       help='Model to use')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging('StructurePredictor')
    
    # Run prediction
    result = predict_structure(
        fasta=args.fasta,
        output_dir=args.output,
        model=args.model
    )
    
    if result['success']:
        logger.info("Structure prediction completed successfully")
        sys.exit(0)
    else:
        logger.error(f"Structure prediction failed: {result['error']}")
        sys.exit(1)

if __name__ == '__main__':
    main() 