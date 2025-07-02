#!/usr/bin/env python3
"""
Results Parsing & Pose Evaluation Script
Task 3: Parse and score docking results from Vina, GNINA, and DiffDock
Outputs: summary.csv and failed_runs.json
"""

import os
import sys
import re
import json
import logging
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np

# Optional imports for RMSD calculation
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign, rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False

# --- Constants ---
VINA_OUT_DIR = 'vina_output'
GNINA_OUT_DIR = 'gnina_output'
DIFFDOCK_OUT_DIR = 'diffdock_output'
LOGS_DIR = 'logs'
OUTPUT_CSV = 'summary.csv'
FAILED_JSON = 'failed_parsing.json'

class DockingResultsParser:
    """Main class for parsing docking results from multiple software packages."""
    
    def __init__(self, base_dir: str = '.', output_dir: str = '.'):
        self.base_dir = Path(base_dir)
        self.output_dir = Path(output_dir)
        self.results = []
        self.failed_runs = {}
        self.setup_logging()
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_file = self.output_dir / 'parsing.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def parse_vina_results(self, ligand_name: str = "ligand") -> List[Dict]:
        """Parse Vina PDBQT output files."""
        vina_results = []
        vina_dir = self.base_dir / VINA_OUT_DIR
        
        if not vina_dir.exists():
            self.logger.warning(f"Vina output directory not found: {vina_dir}")
            return vina_results
            
        # Look for vina output files (avoid duplicates)
        vina_files = list(vina_dir.glob("*.pdbqt"))
        # Add subdirectory files only if not already found
        subdir_files = list(vina_dir.glob("**/vina_out.pdbqt"))
        for sf in subdir_files:
            if sf not in vina_files:
                vina_files.append(sf)
        
        if not vina_files:
            self.logger.warning(f"No Vina PDBQT files found in {vina_dir}")
            self.failed_runs['vina'] = "No output files found"
            return vina_results
            
        for vina_file in vina_files:
            try:
                poses = self._parse_vina_pdbqt(vina_file)
                for i, pose in enumerate(poses):
                    result = {
                        'ligand_name': ligand_name,
                        'method': 'Vina',
                        'affinity_kcal_mol': pose['score'],
                        'rmsd_lb': pose.get('rmsd_lb', None),
                        'rmsd_ub': pose.get('rmsd_ub', None),
                        'RMSD': pose.get('rmsd_lb', None),  # Use lower bound RMSD as default
                        'pose_rank': i + 1,
                        'pose_path': str(vina_file),
                        'model_number': pose['model']
                    }
                    vina_results.append(result)
                    
                self.logger.info(f"Parsed {len(poses)} Vina poses from {vina_file}")
                
            except Exception as e:
                error_msg = f"Failed to parse Vina file {vina_file}: {str(e)}"
                self.logger.error(error_msg)
                self.failed_runs[f'vina_{vina_file.name}'] = error_msg
                
        return vina_results
        
    def _parse_vina_pdbqt(self, file_path: Path) -> List[Dict]:
        """Parse individual Vina PDBQT file."""
        poses = []
        current_model = None
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    current_model = int(line.split()[1])
                elif line.startswith('REMARK VINA RESULT:'):
                    # Parse: REMARK VINA RESULT:    -4.145      0.000      0.000
                    parts = line.split()
                    if len(parts) >= 6:
                        score = float(parts[3])
                        rmsd_lb = float(parts[4]) if parts[4] != '0.000' else None
                        rmsd_ub = float(parts[5]) if parts[5] != '0.000' else None
                        
                        poses.append({
                            'model': current_model,
                            'score': score,
                            'rmsd_lb': rmsd_lb,
                            'rmsd_ub': rmsd_ub
                        })
                        
        return poses
        
    def parse_gnina_results(self, ligand_name: str = "ligand") -> List[Dict]:
        """Parse GNINA output files."""
        gnina_results = []
        gnina_dir = self.base_dir / GNINA_OUT_DIR
        
        if not gnina_dir.exists():
            self.logger.warning(f"GNINA output directory not found: {gnina_dir}")
            return gnina_results
            
        # Look for GNINA score file
        score_file = gnina_dir / 'gnina_scores.txt'
        pdbqt_file = gnina_dir / 'gnina_out.pdbqt'
        log_file = gnina_dir / 'gnina.log'
        
        try:
            if score_file.exists():
                # Parse score file (preferred method)
                poses = self._parse_gnina_scores(score_file)
            elif pdbqt_file.exists():
                # Parse PDBQT file as fallback
                poses = self._parse_gnina_pdbqt(pdbqt_file)
            elif log_file.exists():
                # Parse log file as last resort
                poses = self._parse_gnina_log(log_file)
            else:
                raise FileNotFoundError("No GNINA output files found")
                
            for i, pose in enumerate(poses):
                result = {
                    'ligand_name': ligand_name,
                    'method': 'GNINA',
                    'affinity_kcal_mol': pose.get('cnn_score', pose.get('score')),
                    'RMSD': pose.get('rmsd'),
                    'pose_rank': i + 1,
                    'pose_path': str(pdbqt_file) if pdbqt_file.exists() else str(score_file),
                    'cnn_score': pose.get('cnn_score'),
                    'cnn_affinity': pose.get('cnn_affinity'),
                    'vina_score': pose.get('vina_score')
                }
                gnina_results.append(result)
                
            self.logger.info(f"Parsed {len(poses)} GNINA poses")
            
        except Exception as e:
            error_msg = f"Failed to parse GNINA results: {str(e)}"
            self.logger.error(error_msg)
            self.failed_runs['gnina'] = error_msg
            
        return gnina_results
        
    def _parse_gnina_scores(self, score_file: Path) -> List[Dict]:
        """Parse GNINA scores.txt file."""
        poses = []
        
        with open(score_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    # GNINA score format varies, try to parse common patterns
                    parts = line.split()
                    if len(parts) >= 2:
                        pose = {}
                        try:
                            # Common format: cnn_score cnn_affinity [other scores...]
                            pose['cnn_score'] = float(parts[0])
                            if len(parts) > 1:
                                pose['cnn_affinity'] = float(parts[1])
                            if len(parts) > 2:
                                pose['vina_score'] = float(parts[2])
                            poses.append(pose)
                        except ValueError:
                            continue
                            
        return poses
        
    def _parse_gnina_pdbqt(self, pdbqt_file: Path) -> List[Dict]:
        """Parse GNINA PDBQT output file."""
        poses = []
        current_model = None
        
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    current_model = int(line.split()[1])
                elif line.startswith('REMARK'):
                    # Look for scoring information in remarks
                    if 'CNN' in line or 'score' in line.lower():
                        # Try to extract scores from remark lines
                        numbers = re.findall(r'-?\d+\.?\d*', line)
                        if numbers:
                            pose = {'model': current_model}
                            pose['score'] = float(numbers[0])
                            poses.append(pose)
                            
        return poses
        
    def _parse_gnina_log(self, log_file: Path) -> List[Dict]:
        """Parse GNINA log file for scoring information."""
        poses = []
        
        with open(log_file, 'r') as f:
            content = f.read()
            # Look for scoring patterns in log
            score_patterns = [
                r'CNN\s+score:\s*(-?\d+\.?\d*)',
                r'Affinity:\s*(-?\d+\.?\d*)',
                r'Score:\s*(-?\d+\.?\d*)'
            ]
            
            for pattern in score_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    for i, score in enumerate(matches):
                        if i < len(poses):
                            poses[i]['score'] = float(score)
                        else:
                            poses.append({'score': float(score)})
                    break
                    
        return poses
        
    def parse_diffdock_results(self, ligand_name: str = "ligand") -> List[Dict]:
        """Parse DiffDock output files."""
        diffdock_results = []
        diffdock_dir = self.base_dir / DIFFDOCK_OUT_DIR
        
        if not diffdock_dir.exists():
            self.logger.warning(f"DiffDock output directory not found: {diffdock_dir}")
            return diffdock_results
            
        try:
            # DiffDock creates complex directory structure
            # Look for SDF files and confidence scores
            sdf_files = list(diffdock_dir.rglob("*.sdf"))
            
            if not sdf_files:
                raise FileNotFoundError("No DiffDock SDF files found")
                
            for sdf_file in sdf_files:
                poses = self._parse_diffdock_sdf(sdf_file)
                for i, pose in enumerate(poses):
                    result = {
                        'ligand_name': ligand_name,
                        'method': 'DiffDock',
                        'affinity_kcal_mol': pose.get('confidence_score'),
                        'RMSD': pose.get('rmsd'),
                        'pose_rank': i + 1,
                        'pose_path': str(sdf_file),
                        'confidence_score': pose.get('confidence_score'),
                        'model_number': pose.get('model', i + 1)
                    }
                    diffdock_results.append(result)
                    
            # Also look for confidence files
            conf_files = list(diffdock_dir.rglob("*confidence*"))
            if conf_files:
                self._add_diffdock_confidences(diffdock_results, conf_files)
                
            self.logger.info(f"Parsed {len(diffdock_results)} DiffDock poses")
            
        except Exception as e:
            error_msg = f"Failed to parse DiffDock results: {str(e)}"
            self.logger.error(error_msg)
            self.failed_runs['diffdock'] = error_msg
            
        return diffdock_results
        
    def _parse_diffdock_sdf(self, sdf_file: Path) -> List[Dict]:
        """Parse DiffDock SDF file."""
        poses = []
        
        if RDKIT_AVAILABLE:
            # Use RDKit to parse SDF
            try:
                supplier = Chem.SDMolSupplier(str(sdf_file))
                for i, mol in enumerate(supplier):
                    if mol is None:
                        continue
                    pose = {'model': i + 1}
                    
                    # Extract properties from SDF
                    for prop in mol.GetPropNames():
                        if 'confidence' in prop.lower():
                            pose['confidence_score'] = float(mol.GetProp(prop))
                        elif 'score' in prop.lower():
                            pose['score'] = float(mol.GetProp(prop))
                            
                    poses.append(pose)
                    
            except Exception as e:
                self.logger.warning(f"RDKit failed to parse {sdf_file}: {e}")
                
        if not poses:
            # Fallback: manual SDF parsing
            poses = self._parse_sdf_manually(sdf_file)
            
        return poses
        
    def _parse_sdf_manually(self, sdf_file: Path) -> List[Dict]:
        """Manually parse SDF file for confidence scores."""
        poses = []
        current_pose = {}
        model_count = 0
        
        with open(sdf_file, 'r') as f:
            for line in f:
                if line.strip() == '$$$$':
                    # End of molecule
                    if current_pose:
                        current_pose['model'] = model_count + 1
                        poses.append(current_pose)
                        model_count += 1
                        current_pose = {}
                elif 'confidence' in line.lower():
                    # Extract confidence score
                    numbers = re.findall(r'-?\d+\.?\d*', line)
                    if numbers:
                        current_pose['confidence_score'] = float(numbers[0])
                        
        return poses
        
    def _add_diffdock_confidences(self, results: List[Dict], conf_files: List[Path]):
        """Add confidence scores from separate confidence files."""
        for conf_file in conf_files:
            try:
                with open(conf_file, 'r') as f:
                    confidences = [float(line.strip()) for line in f if line.strip()]
                    
                # Match confidences to results
                for i, confidence in enumerate(confidences):
                    if i < len(results):
                        if results[i]['affinity_kcal_mol'] is None:
                            results[i]['affinity_kcal_mol'] = confidence
                        results[i]['confidence_score'] = confidence
                        
            except Exception as e:
                self.logger.warning(f"Failed to parse confidence file {conf_file}: {e}")
                
    def calculate_rmsd(self, reference_file: str, pose_file: str, 
                      method: str = 'auto') -> Optional[float]:
        """Calculate RMSD between reference and pose structures."""
        if not RDKIT_AVAILABLE and not MDANALYSIS_AVAILABLE:
            self.logger.warning("No RMSD calculation libraries available")
            return None
            
        try:
            if method == 'rdkit' and RDKIT_AVAILABLE:
                return self._calculate_rmsd_rdkit(reference_file, pose_file)
            elif method == 'mdanalysis' and MDANALYSIS_AVAILABLE:
                return self._calculate_rmsd_mda(reference_file, pose_file)
            elif method == 'auto':
                # Try RDKit first, fallback to MDAnalysis
                if RDKIT_AVAILABLE:
                    return self._calculate_rmsd_rdkit(reference_file, pose_file)
                elif MDANALYSIS_AVAILABLE:
                    return self._calculate_rmsd_mda(reference_file, pose_file)
                    
        except Exception as e:
            self.logger.warning(f"RMSD calculation failed: {e}")
            
        return None
        
    def _calculate_rmsd_rdkit(self, ref_file: str, pose_file: str) -> Optional[float]:
        """Calculate RMSD using RDKit."""
        try:
            if ref_file.endswith('.sdf'):
                ref_mol = Chem.MolFromMolFile(ref_file)
            else:
                # Convert PDBQT to MOL format (simplified)
                ref_mol = None
                
            if pose_file.endswith('.sdf'):
                pose_mol = Chem.MolFromMolFile(pose_file)
            else:
                pose_mol = None
                
            if ref_mol and pose_mol:
                # Align molecules and calculate RMSD
                rmsd = rdMolAlign.AlignMol(pose_mol, ref_mol)
                return rmsd
                
        except Exception as e:
            self.logger.warning(f"RDKit RMSD calculation failed: {e}")
            
        return None
        
    def _calculate_rmsd_mda(self, ref_file: str, pose_file: str) -> Optional[float]:
        """Calculate RMSD using MDAnalysis."""
        try:
            ref_u = mda.Universe(ref_file)
            pose_u = mda.Universe(pose_file)
            
            # Select non-hydrogen atoms
            ref_atoms = ref_u.select_atoms("not name H*")
            pose_atoms = pose_u.select_atoms("not name H*")
            
            if len(ref_atoms) == len(pose_atoms):
                rmsd_val = rms.rmsd(ref_atoms.positions, pose_atoms.positions)
                return rmsd_val
                
        except Exception as e:
            self.logger.warning(f"MDAnalysis RMSD calculation failed: {e}")
            
        return None
        
    def generate_summary(self, ligand_name: str = "ligand", 
                        reference_structure: Optional[str] = None) -> pd.DataFrame:
        """Generate comprehensive summary of all docking results."""
        self.logger.info("Starting results parsing...")
        
        # Parse results from all methods
        vina_results = self.parse_vina_results(ligand_name)
        gnina_results = self.parse_gnina_results(ligand_name)
        diffdock_results = self.parse_diffdock_results(ligand_name)
        
        # Combine all results
        all_results = vina_results + gnina_results + diffdock_results
        
        if not all_results:
            self.logger.warning("No docking results found!")
            return pd.DataFrame()
            
        # Calculate additional RMSD if reference provided
        if reference_structure and os.path.exists(reference_structure):
            self.logger.info(f"Calculating RMSD against reference: {reference_structure}")
            for result in all_results:
                if result.get('pose_path') and os.path.exists(result['pose_path']):
                    rmsd = self.calculate_rmsd(reference_structure, result['pose_path'])
                    if rmsd is not None:
                        result['RMSD_calculated'] = rmsd
                        
        # Create DataFrame
        df = pd.DataFrame(all_results)
        
        # Ensure required columns exist
        required_columns = ['ligand_name', 'method', 'affinity_kcal_mol', 'RMSD', 'pose_path']
        for col in required_columns:
            if col not in df.columns:
                df[col] = None
                
        # Sort by affinity (best scores first)
        if 'affinity_kcal_mol' in df.columns and len(df) > 0:
            # Handle NaN values by putting them at the end
            df = df.sort_values('affinity_kcal_mol', ascending=True).reset_index(drop=True)
            
        self.logger.info(f"Generated summary with {len(df)} total poses")
        return df
        
    def save_results(self, df: pd.DataFrame, output_csv: str = OUTPUT_CSV, 
                    failed_json: str = FAILED_JSON):
        """Save results to CSV and failed runs to JSON."""
        # Save main results
        output_path = self.output_dir / output_csv
        df.to_csv(output_path, index=False)
        self.logger.info(f"Results saved to: {output_path}")
        
        # Save failed runs
        if self.failed_runs:
            failed_path = self.output_dir / failed_json
            with open(failed_path, 'w') as f:
                json.dump(self.failed_runs, f, indent=2)
            self.logger.info(f"Failed runs logged to: {failed_path}")
            
        # Print summary
        self.print_summary(df)
        
    def print_summary(self, df: pd.DataFrame):
        """Print summary statistics."""
        if df.empty:
            print("\nNo results to summarize.")
            return
            
        print(f"\n{'='*60}")
        print("DOCKING RESULTS SUMMARY")
        print(f"{'='*60}")
        
        print(f"Total poses: {len(df)}")
        
        if 'method' in df.columns:
            method_counts = df['method'].value_counts()
            print("\nPoses by method:")
            for method, count in method_counts.items():
                print(f"  {method}: {count}")
                
        if 'affinity_kcal_mol' in df.columns:
            valid_scores = df['affinity_kcal_mol'].dropna()
            if not valid_scores.empty:
                print(f"\nAffinity statistics (kcal/mol):")
                print(f"  Best score: {valid_scores.min():.3f}")
                print(f"  Worst score: {valid_scores.max():.3f}")
                print(f"  Mean score: {valid_scores.mean():.3f}")
                
        if 'RMSD' in df.columns:
            valid_rmsd = df['RMSD'].dropna()
            if not valid_rmsd.empty:
                print(f"\nRMSD statistics (Å):")
                print(f"  Best RMSD: {valid_rmsd.min():.3f}")
                print(f"  Worst RMSD: {valid_rmsd.max():.3f}")
                print(f"  Mean RMSD: {valid_rmsd.mean():.3f}")
                
        if self.failed_runs:
            print(f"\nFailed runs: {len(self.failed_runs)}")
            for method, error in self.failed_runs.items():
                print(f"  {method}: {error[:100]}...")
                
        print(f"{'='*60}")

def process_batch_directories(base_dirs, output_dir, reference=None, verbose=False):
    """
    Process multiple docking result directories in batch mode for scalability.
    
    Args:
        base_dirs: List of base directories containing docking outputs
        output_dir: Directory to save consolidated results
        reference: Reference structure for RMSD calculation
        verbose: Enable verbose logging
    
    Returns:
        Consolidated DataFrame with all results
    """
    all_results = []
    batch_failed_runs = {}
    
    for i, base_dir in enumerate(base_dirs):
        ligand_name = os.path.basename(os.path.abspath(base_dir))
        logging.info(f"Processing directory {i+1}/{len(base_dirs)}: {base_dir} (ligand: {ligand_name})")
        
        try:
            # Initialize parser for this directory
            parser_obj = DockingResultsParser(base_dir, output_dir)
            parser_obj.setup_logging()
            
            # Generate summary for this directory
            df = parser_obj.generate_summary(ligand_name, reference)
            
            if not df.empty:
                # Add directory source information
                df['source_directory'] = base_dir
                all_results.append(df)
                
            # Collect failed runs
            if parser_obj.failed_runs:
                batch_failed_runs[base_dir] = parser_obj.failed_runs
                
        except Exception as e:
            error_msg = f"Failed to process directory {base_dir}: {str(e)}"
            logging.error(error_msg)
            batch_failed_runs[base_dir] = {'directory_processing': error_msg}
    
    # Consolidate all results
    if all_results:
        consolidated_df = pd.concat(all_results, ignore_index=True)
        
        # Sort by affinity across all results
        if 'affinity_kcal_mol' in consolidated_df.columns:
            consolidated_df = consolidated_df.sort_values('affinity_kcal_mol', ascending=True).reset_index(drop=True)
            
        logging.info(f"Batch processing complete: {len(consolidated_df)} total poses from {len(all_results)} directories")
    else:
        consolidated_df = pd.DataFrame()
        logging.warning("No valid results found in any directory")
    
    # Save batch failed runs
    if batch_failed_runs:
        batch_failed_path = os.path.join(output_dir, 'batch_failed_parsing.json')
        with open(batch_failed_path, 'w') as f:
            json.dump(batch_failed_runs, f, indent=2)
        logging.info(f"Batch failed runs logged to: {batch_failed_path}")
    
    return consolidated_df, batch_failed_runs


def main():
    """Main CLI interface."""
    parser = argparse.ArgumentParser(
        description='Parse and score docking results from Vina, GNINA, and DiffDock',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single directory
  python3 parse_and_score_results.py
  python3 parse_and_score_results.py --ligand_name aspirin --base_dir ./results
  python3 parse_and_score_results.py --reference ligand_native.sdf --output summary_with_rmsd.csv
  
  # Batch processing multiple directories
  python3 parse_and_score_results.py --batch_dirs ./results_1 ./results_2 ./results_3 --output batch_summary.csv
        """
    )
    
    parser.add_argument('--base_dir', default='.',
                       help='Base directory containing docking output folders (default: current directory)')
    parser.add_argument('--batch_dirs', nargs='+',
                       help='Multiple base directories for batch processing')
    parser.add_argument('--output_dir', default='.',
                       help='Output directory for results (default: current directory)')
    parser.add_argument('--ligand_name', default='ligand',
                       help='Name of the ligand being analyzed (default: ligand)')
    parser.add_argument('--reference', 
                       help='Reference structure file for RMSD calculation (SDF, PDB, etc.)')
    parser.add_argument('--output_csv', default=OUTPUT_CSV,
                       help=f'Output CSV filename (default: {OUTPUT_CSV})')
    parser.add_argument('--failed_json', default=FAILED_JSON,
                       help=f'Failed runs JSON filename (default: {FAILED_JSON})')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Determine processing mode
    if args.batch_dirs:
        # Batch processing mode
        logging.info(f"Batch mode: processing {len(args.batch_dirs)} directories")
        
        # Validate all directories exist
        for batch_dir in args.batch_dirs:
            if not os.path.isdir(batch_dir):
                logging.error(f"Directory not found: {batch_dir}")
                sys.exit(1)
        
        # Process all directories
        consolidated_df, batch_failed_runs = process_batch_directories(
            args.batch_dirs, args.output_dir, args.reference, args.verbose
        )
        
        # Save consolidated results
        if not consolidated_df.empty:
            output_path = os.path.join(args.output_dir, args.output_csv)
            consolidated_df.to_csv(output_path, index=False)
            logging.info(f"Batch results saved to: {output_path}")
            
            # Print batch summary
            print(f"\n{'='*60}")
            print("BATCH PROCESSING SUMMARY")
            print(f"{'='*60}")
            print(f"Directories processed: {len(args.batch_dirs)}")
            print(f"Total poses: {len(consolidated_df)}")
            
            if 'method' in consolidated_df.columns:
                method_counts = consolidated_df['method'].value_counts()
                print("\nPoses by method:")
                for method, count in method_counts.items():
                    print(f"  {method}: {count}")
            
            if 'source_directory' in consolidated_df.columns:
                dir_counts = consolidated_df['source_directory'].value_counts()
                print("\nPoses by directory:")
                for directory, count in dir_counts.items():
                    print(f"  {os.path.basename(directory)}: {count}")
            
            print(f"{'='*60}")
        
        # Determine exit code
        if batch_failed_runs and consolidated_df.empty:
            sys.exit(1)  # All processing failed
        elif batch_failed_runs:
            sys.exit(2)  # Partial success
        else:
            sys.exit(0)  # Full success
    
    else:
        # Single directory processing mode (original behavior)
        parser_obj = DockingResultsParser(args.base_dir, args.output_dir)
        
        # Generate summary
        df = parser_obj.generate_summary(args.ligand_name, args.reference)
        
        # Save results
        parser_obj.save_results(df, args.output_csv, args.failed_json)
        
        # Return appropriate exit code
        if parser_obj.failed_runs and df.empty:
            sys.exit(1)  # All parsing failed
        elif parser_obj.failed_runs:
            sys.exit(2)  # Partial success
        else:
            sys.exit(0)  # Full success

if __name__ == '__main__':
    main()
