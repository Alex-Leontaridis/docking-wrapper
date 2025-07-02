#!/usr/bin/env python3
"""
Batch Molecular Docking Pipeline
Task 4: Orchestrates structure preparation, docking, and results parsing across multiple ligands
Supports AutoDock Vina, GNINA, and DiffDock with comprehensive logging and error handling
"""

import os
import sys
import json
import time
import logging
import argparse
import traceback
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import subprocess
import shutil
from multiprocessing import Pool, cpu_count
from functools import partial

# Add scripts directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import existing modules
from prep_structures import (
    prepare_protein, prepare_ligand_single, validate_file, 
    SUPPORTED_PROTEIN_EXT, SUPPORTED_LIGAND_EXT
)
from run_docking_multi import (
    run_vina, run_gnina, run_diffdock, extract_box_from_protein,
    ensure_output_dirs as ensure_single_output_dirs
)
from parse_and_score_results import DockingResultsParser

# --- Constants ---
DEFAULT_CONFIG = {
    "engines": {
        "vina": {"enabled": True, "exhaustiveness": 8, "num_modes": 9},
        "gnina": {"enabled": False, "use_gpu": False, "cnn_scoring": "rescore"},
        "diffdock": {"enabled": False, "inference_steps": 20, "samples_per_complex": 10}
    },
    "box": {
        "auto_detect": True,
        "default_size": [25.0, 25.0, 25.0],
        "padding": 8.0
    },
    "timeouts": {
        "preparation": 300,
        "vina": 1800,
        "gnina": 3600,
        "diffdock": 7200
    },
    "parallel": {
        "max_workers": min(4, cpu_count()),
        "chunk_size": 1
    },
    "logging": {
        "level": "INFO",
        "max_file_size_mb": 100,
        "backup_count": 5
    }
}

class BatchDockingPipeline:
    """Main batch docking pipeline orchestrator."""
    
    def __init__(self, config_path: Optional[str] = None, output_dir: str = "outputs"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Load configuration
        self.config = self._load_config(config_path)
        
        # Setup logging
        self.setup_logging()
        
        # Initialize tracking
        self.ligand_results = {}
        self.failed_ligands = {}
        self.total_ligands = 0
        self.successful_ligands = 0
        self.start_time = time.time()
        
        # Ensure base output directories
        self._ensure_output_dirs()
    
    def _load_config(self, config_path: Optional[str]) -> Dict:
        """Load and validate configuration."""
        config = DEFAULT_CONFIG.copy()
        
        if config_path and os.path.exists(config_path):
            try:
                with open(config_path, 'r') as f:
                    user_config = json.load(f)
                # Deep merge configuration
                config = self._deep_merge_config(config, user_config)
            except Exception as e:
                print(f"Warning: Failed to load config from {config_path}: {e}")
                print("Using default configuration")
        
        return config
    
    def _deep_merge_config(self, base: Dict, update: Dict) -> Dict:
        """Deep merge configuration dictionaries."""
        for key, value in update.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                base[key] = self._deep_merge_config(base[key], value)
            else:
                base[key] = value
        return base
    
    def setup_logging(self):
        """Setup comprehensive logging."""
        log_dir = self.output_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        
        # Main log file
        log_file = log_dir / "batch_log.txt"
        
        # Configure logging
        log_level = getattr(logging, self.config["logging"]["level"])
        
        # Create formatters
        detailed_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)-20s | %(message)s'
        )
        simple_formatter = logging.Formatter('%(levelname)s: %(message)s')
        
        # Configure root logger
        self.logger = logging.getLogger('BatchPipeline')
        self.logger.setLevel(log_level)
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # File handler with rotation
        from logging.handlers import RotatingFileHandler
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=self.config["logging"]["max_file_size_mb"] * 1024 * 1024,
            backupCount=self.config["logging"]["backup_count"]
        )
        file_handler.setFormatter(detailed_formatter)
        file_handler.setLevel(log_level)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(simple_formatter)
        console_handler.setLevel(log_level)
        
        # Add handlers
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def _ensure_output_dirs(self):
        """Ensure all necessary output directories exist."""
        base_dirs = [
            "logs",
            "prepared_structures",
            "docking_results", 
            "parsed_results",
            "summary_reports"
        ]
        
        for dir_name in base_dirs:
            (self.output_dir / dir_name).mkdir(exist_ok=True)
    
    def discover_ligands(self, ligand_input: str) -> List[Tuple[str, Path]]:
        """
        Discover ligands from input path (file or directory).
        Returns list of (ligand_name, ligand_path) tuples.
        """
        ligands = []
        ligand_path = Path(ligand_input)
        
        if ligand_path.is_file():
            # Single ligand file
            ligand_name = ligand_path.stem
            ligands.append((ligand_name, ligand_path))
            
        elif ligand_path.is_dir():
            # Directory of ligands
            extensions = SUPPORTED_LIGAND_EXT
            for ext in extensions:
                for ligand_file in ligand_path.glob(f"*{ext}"):
                    ligand_name = ligand_file.stem
                    ligands.append((ligand_name, ligand_file))
        else:
            raise ValueError(f"Ligand input path does not exist: {ligand_input}")
        
        if not ligands:
            raise ValueError(f"No supported ligand files found in: {ligand_input}")
        
        self.logger.info(f"Discovered {len(ligands)} ligand(s) to process")
        return ligands
    
    def prepare_protein_once(self, protein_file: str) -> str:
        """Prepare protein structure once for all ligands."""
        self.logger.info(f"Preparing protein structure: {protein_file}")
        
        # Validate protein file
        validate_file(protein_file, SUPPORTED_PROTEIN_EXT, "Protein")
        
        # Prepare protein in the prepared_structures directory
        prepared_dir = self.output_dir / "prepared_structures"
        protein_name = Path(protein_file).stem
        prepared_protein = prepared_dir / f"{protein_name}_prepared.pdbqt"
        
        # Use existing preparation function but with custom output
        original_wd = os.getcwd()
        try:
            os.chdir(prepared_dir)
            prepared_file = prepare_protein(protein_file)
            # Move to final location with proper naming
            if prepared_file != str(prepared_protein):
                # Ensure target directory exists
                prepared_protein.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(prepared_file, prepared_protein)
        finally:
            os.chdir(original_wd)
        
        self.logger.info(f"Protein prepared: {prepared_protein}")
        return str(prepared_protein)
    
    def process_single_ligand(self, ligand_info: Tuple[str, Path], 
                            prepared_protein: str) -> Dict[str, Any]:
        """Process a single ligand through the complete pipeline."""
        ligand_name, ligand_path = ligand_info
        
        result = {
            'ligand_name': ligand_name,
            'ligand_path': str(ligand_path),
            'success': False,
            'stages': {},
            'timings': {},
            'errors': {}
        }
        
        stage_start = time.time()
        
        try:
            self.logger.info(f"Processing ligand: {ligand_name}")
            
            # Create ligand-specific output directory
            ligand_output_dir = self.output_dir / "docking_results" / ligand_name
            ligand_output_dir.mkdir(exist_ok=True)
            
            # Stage 1: Ligand preparation
            self.logger.info(f"[{ligand_name}] Stage 1: Preparing ligand structure")
            stage_start = time.time()
            
            prepared_ligand_dir = self.output_dir / "prepared_structures"
            prepared_ligand = prepared_ligand_dir / f"{ligand_name}_prepared.pdbqt"
            
            try:
                prepare_ligand_single(str(ligand_path), str(prepared_ligand))
                result['stages']['preparation'] = True
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.info(f"[{ligand_name}] Ligand prepared successfully")
            except Exception as e:
                result['stages']['preparation'] = False
                result['errors']['preparation'] = str(e)
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] Ligand preparation failed: {e}")
                return result
            
            # Ensure docking output directories
            ensure_single_output_dirs(str(ligand_output_dir))
            
            # Stage 2: Docking with enabled engines
            enabled_engines = []
            if self.config["engines"]["vina"]["enabled"]:
                enabled_engines.append("vina")
            if self.config["engines"]["gnina"]["enabled"]:
                enabled_engines.append("gnina")
            if self.config["engines"]["diffdock"]["enabled"]:
                enabled_engines.append("diffdock")
            
            self.logger.info(f"[{ligand_name}] Stage 2: Running docking with engines: {enabled_engines}")
            
            # Get docking box parameters
            if self.config["box"]["auto_detect"]:
                try:
                    box_params = extract_box_from_protein(prepared_protein)
                    self.logger.info(f"[{ligand_name}] Auto-detected box parameters: {box_params}")
                except Exception as e:
                    self.logger.warning(f"[{ligand_name}] Box auto-detection failed: {e}, using defaults")
                    default_size = self.config["box"]["default_size"]
                    box_params = (0.0, 0.0, 0.0, *default_size)
            else:
                default_size = self.config["box"]["default_size"]
                box_params = (0.0, 0.0, 0.0, *default_size)
            
            # Run each enabled docking engine
            for engine in enabled_engines:
                engine_start = time.time()
                self.logger.info(f"[{ligand_name}] Running {engine.upper()} docking")
                
                try:
                    if engine == "vina":
                        status = run_vina(
                            prepared_protein, 
                            str(prepared_ligand), 
                            str(ligand_output_dir), 
                            box_params
                        )
                    elif engine == "gnina":
                        status = run_gnina(
                            prepared_protein, 
                            str(prepared_ligand), 
                            str(ligand_output_dir),
                            self.config["engines"]["gnina"]["use_gpu"]
                        )
                    elif engine == "diffdock":
                        status = run_diffdock(
                            prepared_protein, 
                            str(prepared_ligand), 
                            str(ligand_output_dir)
                        )
                    
                    result['stages'][engine] = status['success']
                    result['timings'][engine] = time.time() - engine_start
                    
                    if status['success']:
                        self.logger.info(f"[{ligand_name}] {engine.upper()} completed successfully")
                    else:
                        result['errors'][engine] = status.get('error', 'Unknown error')
                        self.logger.error(f"[{ligand_name}] {engine.upper()} failed: {status.get('error')}")
                        
                except Exception as e:
                    result['stages'][engine] = False
                    result['errors'][engine] = str(e)
                    result['timings'][engine] = time.time() - engine_start
                    self.logger.error(f"[{ligand_name}] {engine.upper()} failed with exception: {e}")
            
            # Stage 3: Results parsing
            self.logger.info(f"[{ligand_name}] Stage 3: Parsing docking results")
            parsing_start = time.time()
            
            try:
                parser = DockingResultsParser(
                    base_dir=str(ligand_output_dir),
                    output_dir=str(self.output_dir / "parsed_results" / ligand_name)
                )
                
                summary_df = parser.generate_summary(ligand_name=ligand_name)
                
                if not summary_df.empty:
                    # Save ligand-specific results
                    ligand_results_dir = self.output_dir / "parsed_results" / ligand_name
                    ligand_results_dir.mkdir(exist_ok=True)
                    
                    summary_file = ligand_results_dir / f"{ligand_name}_summary.csv"
                    summary_df.to_csv(summary_file, index=False)
                    
                    # Save failed runs
                    if parser.failed_runs:
                        failed_file = ligand_results_dir / f"{ligand_name}_failed.json"
                        with open(failed_file, 'w') as f:
                            json.dump(parser.failed_runs, f, indent=2)
                    
                    result['stages']['parsing'] = True
                    result['summary_file'] = str(summary_file)
                    result['poses_found'] = len(summary_df)
                    self.logger.info(f"[{ligand_name}] Parsing completed, found {len(summary_df)} poses")
                else:
                    result['stages']['parsing'] = False
                    result['errors']['parsing'] = "No poses found in results"
                    self.logger.warning(f"[{ligand_name}] Parsing completed but no poses found")
                
                result['timings']['parsing'] = time.time() - parsing_start
                
            except Exception as e:
                result['stages']['parsing'] = False
                result['errors']['parsing'] = str(e)
                result['timings']['parsing'] = time.time() - parsing_start
                self.logger.error(f"[{ligand_name}] Results parsing failed: {e}")
            
            # Determine overall success
            any_docking_success = any(result['stages'].get(engine, False) for engine in enabled_engines)
            parsing_success = result['stages'].get('parsing', False)
            
            result['success'] = any_docking_success and parsing_success
            result['total_time'] = time.time() - stage_start
            
            if result['success']:
                self.logger.info(f"[{ligand_name}] Processing completed successfully in {result['total_time']:.1f}s")
            else:
                self.logger.warning(f"[{ligand_name}] Processing completed with issues in {result['total_time']:.1f}s")
            
        except Exception as e:
            result['success'] = False
            result['errors']['general'] = str(e)
            result['total_time'] = time.time() - stage_start
            self.logger.error(f"[{ligand_name}] Processing failed: {e}")
            self.logger.debug(f"[{ligand_name}] Traceback: {traceback.format_exc()}")
        
        return result
    
    def run_batch_pipeline(self, protein_file: str, ligand_input: str, 
                          parallel: bool = True) -> Dict[str, Any]:
        """Run the complete batch pipeline."""
        pipeline_start = time.time()
        
        self.logger.info("="*80)
        self.logger.info("STARTING BATCH MOLECULAR DOCKING PIPELINE")
        self.logger.info("="*80)
        
        try:
            # Stage 0: Discover ligands
            ligands = self.discover_ligands(ligand_input)
            self.total_ligands = len(ligands)
            
            # Stage 1: Prepare protein once
            prepared_protein = self.prepare_protein_once(protein_file)
            
            # Stage 2: Process ligands in parallel or serial
            if parallel and self.config["parallel"]["max_workers"] > 1 and len(ligands) > 1:
                self.logger.info(f"Processing {len(ligands)} ligands in parallel with {self.config['parallel']['max_workers']} workers")
                
                with Pool(processes=self.config["parallel"]["max_workers"]) as pool:
                    process_func = partial(self.process_single_ligand, prepared_protein=prepared_protein)
                    results = pool.map(process_func, ligands)
            else:
                self.logger.info(f"Processing {len(ligands)} ligands serially")
                results = []
                for ligand_info in ligands:
                    result = self.process_single_ligand(ligand_info, prepared_protein)
                    results.append(result)
            
            # Stage 3: Aggregate results
            self.logger.info("Aggregating results and generating final summary")
            
            for result in results:
                ligand_name = result['ligand_name']
                if result['success']:
                    self.ligand_results[ligand_name] = result
                    self.successful_ligands += 1
                else:
                    self.failed_ligands[ligand_name] = result
            
            # Generate final summary
            final_summary = self._generate_final_summary()
            
            # Save batch log
            self._save_batch_log(results, pipeline_start)
            
            pipeline_time = time.time() - pipeline_start
            success_rate = (self.successful_ligands / self.total_ligands) * 100
            
            self.logger.info("="*80)
            self.logger.info("BATCH PIPELINE COMPLETED")
            self.logger.info(f"Total ligands: {self.total_ligands}")
            self.logger.info(f"Successful: {self.successful_ligands}")
            self.logger.info(f"Failed: {len(self.failed_ligands)}")
            self.logger.info(f"Success rate: {success_rate:.1f}%")
            self.logger.info(f"Total time: {pipeline_time:.1f}s")
            self.logger.info("="*80)
            
            return {
                'success': True,
                'total_ligands': self.total_ligands,
                'successful_ligands': self.successful_ligands,
                'failed_ligands': len(self.failed_ligands),
                'success_rate': success_rate,
                'total_time': pipeline_time,
                'final_summary_file': final_summary
            }
            
        except Exception as e:
            self.logger.error(f"Batch pipeline failed: {e}")
            self.logger.debug(f"Traceback: {traceback.format_exc()}")
            return {
                'success': False,
                'error': str(e),
                'total_time': time.time() - pipeline_start
            }
    
    def _generate_final_summary(self) -> str:
        """Generate final aggregated summary CSV."""
        all_results = []
        
        # Collect all individual results
        for ligand_name, result in self.ligand_results.items():
            summary_file = result.get('summary_file')
            if summary_file and os.path.exists(summary_file):
                try:
                    df = pd.read_csv(summary_file)
                    all_results.append(df)
                except Exception as e:
                    self.logger.warning(f"Failed to read summary for {ligand_name}: {e}")
        
        if all_results:
            # Combine all results
            final_df = pd.concat(all_results, ignore_index=True)
            
            # Add batch-level statistics
            if not final_df.empty:
                # Sort by affinity (best first)
                if 'affinity_kcal_mol' in final_df.columns:
                    final_df = final_df.sort_values('affinity_kcal_mol', ascending=True)
                
                # Add ranking within each ligand
                final_df['global_rank'] = range(1, len(final_df) + 1)
            
            # Save final summary
            final_summary_file = self.output_dir / "final_summary.csv"
            final_df.to_csv(final_summary_file, index=False)
            
            self.logger.info(f"Final summary saved: {final_summary_file}")
            self.logger.info(f"Total poses in final summary: {len(final_df)}")
            
            return str(final_summary_file)
        else:
            self.logger.warning("No successful results to aggregate")
            return ""
    
    def _save_batch_log(self, results: List[Dict], pipeline_start: float):
        """Save detailed batch processing log."""
        batch_log = {
            'pipeline_start': pipeline_start,
            'pipeline_end': time.time(),
            'total_time': time.time() - pipeline_start,
            'config': self.config,
            'ligand_results': results,
            'summary': {
                'total_ligands': self.total_ligands,
                'successful_ligands': self.successful_ligands,
                'failed_ligands': len(self.failed_ligands),
                'success_rate': (self.successful_ligands / self.total_ligands) * 100 if self.total_ligands > 0 else 0
            }
        }
        
        batch_log_file = self.output_dir / "batch_log.json"
        with open(batch_log_file, 'w') as f:
            json.dump(batch_log, f, indent=2, default=str)
        
        self.logger.info(f"Detailed batch log saved: {batch_log_file}")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Batch Molecular Docking Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single protein, directory of ligands
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligand_dir/
  
  # With custom config
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --config pipeline_config.json
  
  # Enable all engines
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-gnina --enable-diffdock
  
  # Serial processing
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --serial
        """
    )
    
    # Required arguments
    parser.add_argument('--protein', required=True,
                       help='Path to protein file (.pdb or .pdbqt)')
    parser.add_argument('--ligands', required=True,
                       help='Path to ligand file or directory of ligand files')
    
    # Optional arguments
    parser.add_argument('--config', 
                       help='Path to pipeline configuration JSON file')
    parser.add_argument('--output-dir', default='outputs',
                       help='Output directory (default: outputs)')
    parser.add_argument('--enable-gnina', action='store_true',
                       help='Enable GNINA docking')
    parser.add_argument('--enable-diffdock', action='store_true',
                       help='Enable DiffDock docking')
    parser.add_argument('--serial', action='store_true',
                       help='Process ligands serially instead of in parallel')
    parser.add_argument('--max-workers', type=int,
                       help='Maximum number of parallel workers')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Create pipeline
    pipeline = BatchDockingPipeline(
        config_path=args.config,
        output_dir=args.output_dir
    )
    
    # Override config with CLI arguments
    if args.enable_gnina:
        pipeline.config["engines"]["gnina"]["enabled"] = True
    if args.enable_diffdock:
        pipeline.config["engines"]["diffdock"]["enabled"] = True
    if args.max_workers:
        pipeline.config["parallel"]["max_workers"] = args.max_workers
    if args.verbose:
        pipeline.config["logging"]["level"] = "DEBUG"
        pipeline.setup_logging()  # Reconfigure with new level
    
    # Run pipeline
    result = pipeline.run_batch_pipeline(
        protein_file=args.protein,
        ligand_input=args.ligands,
        parallel=not args.serial
    )
    
    # Exit with appropriate code
    if result['success']:
        success_rate = result.get('success_rate', 0)
        if success_rate == 100:
            sys.exit(0)  # Perfect success
        elif success_rate >= 50:
            sys.exit(0)  # Acceptable success
        else:
            sys.exit(2)  # Partial success
    else:
        sys.exit(1)  # Total failure


if __name__ == "__main__":
    main() 