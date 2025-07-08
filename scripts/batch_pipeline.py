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
import threading
from multiprocessing import Lock as ProcessLock
from multiprocessing.managers import SharedMemoryManager
import atexit
import signal

# Pretty output helpers
try:
    from colorama import init as colorama_init, Fore, Style
    colorama_init()
    COLOR_OK = Fore.GREEN + Style.BRIGHT
    COLOR_FAIL = Fore.RED + Style.BRIGHT
    COLOR_WARN = Fore.YELLOW + Style.BRIGHT
    COLOR_INFO = Fore.CYAN + Style.BRIGHT
    COLOR_RESET = Style.RESET_ALL
except ImportError:
    COLOR_OK = COLOR_FAIL = COLOR_WARN = COLOR_INFO = COLOR_RESET = ''

def pretty_status(success):
    return f"{COLOR_OK}OK{COLOR_RESET}" if success else f"{COLOR_FAIL}FAIL{COLOR_RESET}"

def pretty_stage(stage):
    return f"{COLOR_INFO}{stage}{COLOR_RESET}"

try:
    from tabulate import tabulate
    HAVE_TABULATE = True
except ImportError:
    HAVE_TABULATE = False

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import existing modules
from prep_structures import (
    prepare_protein, prepare_ligand_single, validate_file, 
    SUPPORTED_PROTEIN_EXT, SUPPORTED_LIGAND_EXT
)
from run_docking_multi import (
    run_vina, run_gnina, run_diffdock, extract_box_from_protein,
    ensure_output_dirs as ensure_single_output_dirs
)
from docking_results_parser import DockingResultsParser
from config import Config
from utils.logging import setup_logging, get_global_logger

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
        """
        Initialize the batch docking pipeline.
        
        Args:
            config_path: Path to configuration JSON file
            output_dir: Output directory for results
        """
        self.output_dir = Path(output_dir)
        self.config = self._load_config(config_path)
        self.setup_logging()
        self._ensure_output_dirs()
        
        # Thread safety and resource management
        self._lock = threading.Lock()
        self._process_lock = ProcessLock()
        self._shared_memory_manager = None
        self._active_processes = set()
        self._cleanup_registered = False
        
        # Register cleanup handlers
        self._register_cleanup_handlers()
        
        # Statistics tracking (thread-safe)
        self.total_ligands = 0
        self.successful_ligands = 0
        self.failed_ligands = set()
        self.ligand_results = {}
        
        # Check external dependencies
        self._check_external_binaries()
    
    def _check_external_binaries(self):
        """Check for required external binaries and exit if missing."""
        required_binaries = ["vina", "gnina", "diffdock"]
        missing = []
        for binary in required_binaries:
            if shutil.which(binary) is None:
                missing.append(binary)
        # Check for MGLTools (optional, but warn if missing)
        self.config_obj = Config()
        if not self.config_obj.validate_mgltools():
            print(f"{COLOR_WARN}Warning: MGLTools not found. Protein preparation may be limited.{COLOR_RESET}")
            print(f"To install MGLTools, visit: http://mgltools.scripps.edu/downloads/downloads/tools/downloads")
            print(f"Or set MGLTOOLS_PATH environment variable to point to your installation.")
        if missing:
            print(f"{COLOR_FAIL}Missing required external binaries: {', '.join(missing)}{COLOR_RESET}")
            print(f"Please install them and ensure they are in your PATH.")
            sys.exit(1)
    
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
        
        # Convert protein file to absolute path
        abs_protein_file = os.path.abspath(protein_file)
        
        # Prepare protein in the prepared_structures directory
        prepared_dir = self.output_dir / "prepared_structures"
        protein_name = Path(protein_file).stem
        prepared_protein = prepared_dir / f"{protein_name}_prepared.pdbqt"
        
        # Use custom preparation with absolute paths to avoid directory changes
        prepared_file = self._prepare_protein_custom(abs_protein_file, str(prepared_protein))
        
        self.logger.info(f"Protein prepared: {prepared_protein}")
        return str(prepared_protein)
    
    def _prepare_protein_custom(self, protein_file: str, output_file: str) -> str:
        """Custom protein preparation that handles absolute paths correctly."""
        import subprocess
        
        ext = os.path.splitext(protein_file)[1].lower()
        
        if ext == '.pdb':
            try:
                self.logger.info("Preparing protein: cleaning, adding hydrogens, assigning Gasteiger charges, converting to PDBQT...")
                
                # Ensure output directory exists
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                
                mgltools_pythonsh = self.config_obj.get_mgltools_pythonsh()
                prepare_script = self.config_obj.get_mgltools_prepare_script()
                
                # Check if MGLTools is available
                if mgltools_pythonsh and prepare_script and os.path.exists(mgltools_pythonsh) and os.path.exists(prepare_script):
                    cmd = [
                        mgltools_pythonsh, prepare_script,
                        '-r', protein_file,
                        '-o', output_file,
                        '-A', 'hydrogens',  # Add hydrogens
                        '-U', 'waters',     # Remove waters
                    ]
                    
                    self.logger.info(f"Running MGLTools: {' '.join(cmd)}")
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    
                    if result.returncode != 0:
                        self.logger.error(f"MGLTools failed: {result.stderr}\n{result.stdout}")
                        raise Exception(f"MGLTools preparation failed: {result.stderr}")
                else:
                    # Fallback for Docker/Linux environments
                    self.logger.info("MGLTools not found. Using simple PDB to PDBQT conversion...")
                    self._simple_pdb_to_pdbqt(protein_file, output_file)
                    # Skip cleaning step for simple conversion - it's already in proper format
                    self.logger.info(f"Basic PDB to PDBQT conversion completed: {output_file}")
                    return output_file
                
                # Clean up PDBQT formatting issues (only for MGLTools output)
                try:
                    self.logger.info("Cleaning PDBQT formatting...")
                    self._clean_pdbqt_formatting(output_file)
                except Exception as e:
                    self.logger.warning(f"PDBQT cleaning failed, but continuing: {e}")
                    # Don't raise - cleaning is optional
                
                self.logger.info(f"Protein cleaned, hydrogens added, Gasteiger charges assigned, and saved as: {output_file}")
                return output_file
                
            except Exception as e:
                self.logger.error(f"Error during protein preparation: {e}")
                raise
                
        elif ext == '.pdbqt':
            # If already PDBQT, just copy to target location
            if protein_file != output_file:
                os.makedirs(os.path.dirname(output_file), exist_ok=True)
                shutil.copy2(protein_file, output_file)
            self.logger.info(f"PDBQT file copied to: {output_file}")
            return output_file
        else:
            raise ValueError(f"Unsupported protein file format: {ext}")
    
    def _clean_pdbqt_formatting(self, pdbqt_file: str):
        """Simplified PDBQT cleaning for batch pipeline."""
        import tempfile
        
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pdbqt')
        
        try:
            with open(pdbqt_file, 'r') as infile:
                for line in infile:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Basic cleanup - skip problematic alternate conformations
                        alt_loc = line[16:17].strip()
                        if alt_loc and alt_loc not in ['', ' ', 'A']:
                            continue
                    temp_file.write(line)
            
            temp_file.close()
            shutil.move(temp_file.name, pdbqt_file)
            self.logger.info(f"Cleaned PDBQT formatting in {pdbqt_file}")
            
        except Exception as e:
            self.logger.error(f"Error cleaning PDBQT file: {e}")
            if os.path.exists(temp_file.name):
                os.unlink(temp_file.name)
            # Don't raise - cleaning is optional
    
    def _simple_pdb_to_pdbqt(self, pdb_file: str, pdbqt_file: str):
        """Simple PDB to PDBQT conversion for Docker environments."""
        try:
            with open(pdb_file, 'r') as f:
                pdb_lines = f.readlines()
            
            pdbqt_lines = []
            for line in pdb_lines:
                if line.startswith(('ATOM', 'HETATM')):
                    # Convert PDB line to basic PDBQT format
                    if len(line) >= 54:  # Minimum length for coordinates
                        # Take only the standard PDB part (up to column 78) but remove element column
                        # PDB format: columns 77-78 contain element symbol which we need to replace
                        pdb_part = line[:76].rstrip()  # Stop before element column
                        
                        # Determine AutoDock atom type based on atom name
                        atom_name = line[12:16].strip()
                        if atom_name.startswith('C'):
                            autodock_type = "C"
                        elif atom_name.startswith('N'):
                            autodock_type = "N"
                        elif atom_name.startswith('O'):
                            autodock_type = "O"
                        elif atom_name.startswith('S'):
                            autodock_type = "S"
                        elif atom_name.startswith('P'):
                            autodock_type = "P"
                        elif atom_name.startswith('H'):
                            autodock_type = "H"
                        else:
                            autodock_type = "C"  # Default to carbon
                        
                        # Format: PDB_part + spaces + charge + space + atom_type
                        # Ensure proper spacing to column 79
                        padding_needed = 76 - len(pdb_part)
                        if padding_needed > 0:
                            pdb_part += ' ' * padding_needed
                        
                        # Use proper AutoDock charge format (6.3f) instead of hardcoded +0.000
                        charge = 0.000  # Default charge, could be calculated from atom type
                        pdbqt_line = f"{pdb_part}  {charge:>6.3f} {autodock_type}"
                        pdbqt_lines.append(pdbqt_line + '\n')
                else:
                    # Skip header and other PDB-specific lines that Vina doesn't need
                    # Only keep essential structural information
                    if line.startswith(('REMARK', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                        pdbqt_lines.append(line)
            
            with open(pdbqt_file, 'w') as f:
                f.writelines(pdbqt_lines)
            
            self.logger.info(f"Basic PDB to PDBQT conversion completed: {pdbqt_file}")
            
        except Exception as e:
            self.logger.error(f"Simple PDB to PDBQT conversion failed: {e}")
            raise
    
    def process_single_ligand(self, ligand_info: Tuple[str, Path], 
                            prepared_protein: str) -> Dict[str, Any]:
        """Process a single ligand through the complete pipeline with strict error handling."""
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
            ligand_output_dir = self.output_dir / "docking_results" / ligand_name
            ligand_output_dir.mkdir(exist_ok=True)
            # Stage 1: Ligand preparation
            self.logger.info(f"[{ligand_name}] Stage 1: Preparing ligand structure")
            stage_start = time.time()
            if not os.path.exists(ligand_path):
                result['stages']['preparation'] = False
                result['errors']['preparation'] = f"Ligand file not found: {ligand_path}"
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] Ligand file not found: {ligand_path}")
                return result  # Early exit
            prepared_ligand_dir = self.output_dir / "prepared_structures"
            prepared_ligand = prepared_ligand_dir / f"{ligand_name}_prepared.pdbqt"
            try:
                prepare_ligand_single(str(ligand_path), str(prepared_ligand))
                if not os.path.exists(prepared_ligand):
                    raise Exception("Ligand preparation did not create output file")
                result['stages']['preparation'] = True
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.info(f"[{ligand_name}] Ligand prepared successfully")
            except Exception as e:
                result['stages']['preparation'] = False
                result['errors']['preparation'] = str(e)
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] Ligand preparation failed: {e}")
                return result  # Early exit
            ensure_single_output_dirs(str(ligand_output_dir))
            # Stage 2: Docking with enabled engines
            enabled_engines = []
            if self.config["engines"]["vina"]["enabled"]:
                enabled_engines.append("vina")
            if self.config["engines"]["gnina"]["enabled"]:
                enabled_engines.append("gnina")
            if self.config["engines"]["diffdock"]["enabled"]:
                enabled_engines.append("diffdock")
            if not enabled_engines:
                result['stages']['docking'] = False
                result['errors']['docking'] = "No docking engines enabled in configuration"
                result['timings']['docking'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] No docking engines enabled")
                return result  # Early exit
            if not os.path.exists(prepared_protein):
                result['stages']['docking'] = False
                result['errors']['docking'] = f"Prepared protein file not found: {prepared_protein}"
                result['timings']['docking'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] Prepared protein file not found: {prepared_protein}")
                return result  # Early exit
            self.logger.info(f"[{ligand_name}] Stage 2: Running docking with engines: {enabled_engines}")
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
            # Run each enabled docking engine, fail early if any engine fails
            engine_failures = []
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
                    
                    result['timings'][engine] = time.time() - engine_start
                    
                    if status['success']:
                        output_files_exist = self._validate_engine_output(engine, str(ligand_output_dir))
                        if not output_files_exist:
                            error_msg = f"{engine.upper()} reported success but no output files found"
                            result['stages'][engine] = False
                            result['errors'][engine] = error_msg
                            engine_failures.append(f"{engine}: {error_msg}")
                            self.logger.error(f"[{ligand_name}] {error_msg}")
                        else:
                            result['stages'][engine] = True
                            self.logger.info(f"[{ligand_name}] {engine.upper()} completed successfully")
                    else:
                        error_msg = status.get('error', 'Unknown error')
                        result['stages'][engine] = False
                        result['errors'][engine] = error_msg
                        engine_failures.append(f"{engine}: {error_msg}")
                        self.logger.error(f"[{ligand_name}] {engine.upper()} failed: {error_msg}")
                        
                except Exception as e:
                    error_msg = str(e)
                    result['stages'][engine] = False
                    result['errors'][engine] = error_msg
                    result['timings'][engine] = time.time() - engine_start
                    engine_failures.append(f"{engine}: {error_msg}")
                    self.logger.error(f"[{ligand_name}] {engine.upper()} failed with exception: {e}")
            
            # Check if any engines failed and decide whether to continue
            if engine_failures:
                if self.config.get("strict_mode", True):  # Default to strict mode
                    result['stages']['docking'] = False
                    result['errors']['docking'] = f"One or more docking engines failed: {'; '.join(engine_failures)}"
                    result['timings']['docking'] = time.time() - stage_start
                    self.logger.error(f"[{ligand_name}] Docking failed due to engine failures: {'; '.join(engine_failures)}")
                    return result  # Early exit in strict mode
                else:
                    # Non-strict mode: continue with successful engines only
                    successful_engines = [eng for eng in enabled_engines if result['stages'].get(eng, False)]
                    if not successful_engines:
                        result['stages']['docking'] = False
                        result['errors']['docking'] = f"All docking engines failed: {'; '.join(engine_failures)}"
                        result['timings']['docking'] = time.time() - stage_start
                        self.logger.error(f"[{ligand_name}] All docking engines failed: {'; '.join(engine_failures)}")
                        return result  # Early exit if all engines failed
                    else:
                        result['stages']['docking'] = True
                        self.logger.warning(f"[{ligand_name}] Some engines failed but continuing with successful ones: {successful_engines}")
            else:
                result['stages']['docking'] = True
            # Stage 3: Results parsing
            self.logger.info(f"[{ligand_name}] Stage 3: Parsing docking results")
            parsing_start = time.time()
            try:
                parser_output_dir = self.output_dir / "parsed_results" / ligand_name
                parser_output_dir.mkdir(parents=True, exist_ok=True)
                parser = DockingResultsParser(
                    base_dir=str(ligand_output_dir),
                    output_dir=str(parser_output_dir)
                )
                summary_df = parser.generate_summary(ligand_name=ligand_name)
                if not summary_df.empty:
                    ligand_results_dir = self.output_dir / "parsed_results" / ligand_name
                    ligand_results_dir.mkdir(exist_ok=True)
                    summary_file = ligand_results_dir / f"{ligand_name}_summary.csv"
                    summary_df.to_csv(summary_file, index=False)
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
                    return result  # Early exit
                result['timings']['parsing'] = time.time() - parsing_start
            except Exception as e:
                result['stages']['parsing'] = False
                result['errors']['parsing'] = str(e)
                result['timings']['parsing'] = time.time() - parsing_start
                self.logger.error(f"[{ligand_name}] Results parsing failed: {e}")
                return result  # Early exit
            # If we reach here, all steps succeeded
            result['success'] = True
            result['total_time'] = time.time() - stage_start
            self.logger.info(f"[{ligand_name}] Processing completed successfully in {result['total_time']:.1f}s")
        except Exception as e:
            result['success'] = False
            result['errors']['general'] = str(e)
            result['total_time'] = time.time() - stage_start
            self.logger.error(f"[{ligand_name}] Processing failed: {e}")
            self.logger.debug(f"[{ligand_name}] Traceback: {traceback.format_exc()}")
        return result
    
    def _validate_engine_output(self, engine: str, output_dir: str) -> bool:
        """Validate that a docking engine produced expected output files."""
        output_path = Path(output_dir) / f"{engine}_output"
        
        if not output_path.exists():
            return False
        
        # Check for expected output files based on engine
        if engine == "vina":
            vina_out = output_path / "vina_out.pdbqt"
            return vina_out.exists() and vina_out.stat().st_size > 0
        elif engine == "gnina":
            # GNINA typically produces .sdf files
            sdf_files = list(output_path.glob("*.sdf"))
            return len(sdf_files) > 0
        elif engine == "diffdock":
            # DiffDock produces .sdf files and confidence.txt
            sdf_files = list(output_path.glob("*.sdf"))
            confidence_file = output_path / "diffdock_confidence.txt"
            return len(sdf_files) > 0 and confidence_file.exists()
        else:
            # Unknown engine, assume success if directory exists and has files
            return len(list(output_path.glob("*"))) > 0
    
    def print_final_summary(self, results: List[Dict]):
        print(f"\n{COLOR_INFO}{'='*60}{COLOR_RESET}")
        print(f"{COLOR_INFO}BATCH PIPELINE SUMMARY{COLOR_RESET}")
        print(f"{COLOR_INFO}{'='*60}{COLOR_RESET}")
        headers = ["Ligand", "Status", "Prep", "Vina", "GNINA", "DiffDock", "Parse", "Time (s)", "Error"]
        table = []
        for r in results:
            stages = r.get('stages', {})
            row = [
                r.get('ligand_name', ''),
                pretty_status(r.get('success', False)),
                pretty_status(stages.get('preparation', False)),
                pretty_status(stages.get('vina', False)) if 'vina' in stages else '-',
                pretty_status(stages.get('gnina', False)) if 'gnina' in stages else '-',
                pretty_status(stages.get('diffdock', False)) if 'diffdock' in stages else '-',
                pretty_status(stages.get('parsing', False)),
                f"{r.get('total_time', 0):.1f}",
                (next(iter(r.get('errors', {}).values()), '')[:40] if r.get('errors') else '')
            ]
            table.append(row)
        if HAVE_TABULATE:
            try:
                print(tabulate(table, headers=headers, tablefmt="fancy_grid"))
            except UnicodeEncodeError:
                # Fallback to simple table format if Unicode issues occur
                print(tabulate(table, headers=headers, tablefmt="simple"))
        else:
            # fallback: aligned text
            col_widths = [max(len(str(x)) for x in col) for col in zip(*([headers]+table))]
            fmt = '  '.join(f'{{:<{w}}}' for w in col_widths)
            print(fmt.format(*headers))
            for row in table:
                print(fmt.format(*row))
        print(f"{COLOR_INFO}{'='*60}{COLOR_RESET}\n")
    
    def run_batch_pipeline(self, protein_file: str, ligand_input: str, parallel: bool = True) -> Dict[str, Any]:
        pipeline_start = time.time()
        print(f"\n{COLOR_INFO}{'='*80}{COLOR_RESET}")
        print(f"{COLOR_INFO}STARTING BATCH MOLECULAR DOCKING PIPELINE{COLOR_RESET}")
        print(f"{COLOR_INFO}{'='*80}{COLOR_RESET}")
        
        try:
            ligands = self.discover_ligands(ligand_input)
            self.total_ligands = len(ligands)
            prepared_protein = self.prepare_protein_once(protein_file)
            
            # Initialize shared memory manager for parallel processing
            if parallel and self.config["parallel"]["max_workers"] > 1 and len(ligands) > 1:
                self._shared_memory_manager = SharedMemoryManager()
                self._shared_memory_manager.start()
                
                print(f"{COLOR_INFO}Processing {len(ligands)} ligands in parallel with {self.config['parallel']['max_workers']} workers{COLOR_RESET}")
                
                # Use thread-safe processing with proper resource management
                with Pool(processes=self.config["parallel"]["max_workers"]) as pool:
                    # Create a thread-safe wrapper for the processing function
                    process_func = partial(self._process_single_ligand_safe, prepared_protein=prepared_protein)
                    
                    # Process ligands with proper error handling
                    results = []
                    for ligand_info in ligands:
                        try:
                            result = pool.apply_async(process_func, (ligand_info,))
                            results.append(result)
                        except Exception as e:
                            self.logger.error(f"Failed to submit ligand {ligand_info[0]} for processing: {e}")
                            results.append(None)
                    
                    # Collect results with timeout
                    processed_results = []
                    for i, result in enumerate(results):
                        try:
                            if result is not None:
                                processed_result = result.get(timeout=3600)  # 1 hour timeout
                                processed_results.append(processed_result)
                                
                                # Update statistics thread-safely
                                ligand_name = ligands[i][0]
                                success = processed_result.get('success', False)
                                self._update_statistics(ligand_name, success, processed_result)
                            else:
                                # Handle failed submission
                                ligand_name = ligands[i][0]
                                failed_result = {
                                    'ligand_name': ligand_name,
                                    'success': False,
                                    'errors': {'submission': 'Failed to submit for processing'},
                                    'total_time': 0
                                }
                                processed_results.append(failed_result)
                                self._update_statistics(ligand_name, False, failed_result)
                                
                        except Exception as e:
                            self.logger.error(f"Failed to get result for ligand {ligands[i][0]}: {e}")
                            failed_result = {
                                'ligand_name': ligands[i][0],
                                'success': False,
                                'errors': {'processing': str(e)},
                                'total_time': 0
                            }
                            processed_results.append(failed_result)
                            self._update_statistics(ligands[i][0], False, failed_result)
                
                results = processed_results
            else:
                print(f"{COLOR_INFO}Processing {len(ligands)} ligands serially{COLOR_RESET}")
                results = []
                for ligand_info in ligands:
                    try:
                        result = self.process_single_ligand(ligand_info, prepared_protein)
                        results.append(result)
                        
                        # Update statistics
                        ligand_name = ligand_info[0]
                        success = result.get('success', False)
                        self._update_statistics(ligand_name, success, result)
                        
                    except Exception as e:
                        self.logger.error(f"Failed to process ligand {ligand_info[0]}: {e}")
                        failed_result = {
                            'ligand_name': ligand_info[0],
                            'success': False,
                            'errors': {'processing': str(e)},
                            'total_time': 0
                        }
                        results.append(failed_result)
                        self._update_statistics(ligand_info[0], False, failed_result)
            
            self.print_final_summary(results)
            return {
                'success': True,
                'total_ligands': self.total_ligands,
                'successful_ligands': self.successful_ligands,
                'failed_ligands': len(self.failed_ligands),
                'success_rate': (self.successful_ligands / self.total_ligands) * 100 if self.total_ligands > 0 else 0,
                'total_time': time.time() - pipeline_start,
                'final_summary_file': self._generate_final_summary()
            }
        except Exception as e:
            print(f"{COLOR_FAIL}Batch pipeline failed: {e}{COLOR_RESET}")
            print(f"{COLOR_FAIL}{traceback.format_exc()}{COLOR_RESET}")
            return {
                'success': False,
                'error': str(e),
                'total_time': time.time() - pipeline_start
            }
        finally:
            # Ensure cleanup happens
            self._cleanup_resources()
    
    def _process_single_ligand_safe(self, ligand_info: Tuple[str, Path], prepared_protein: str) -> Dict[str, Any]:
        """
        Thread-safe wrapper for processing a single ligand.
        
        Args:
            ligand_info: Tuple of (ligand_name, ligand_path)
            prepared_protein: Path to prepared protein file
            
        Returns:
            Processing result dictionary
        """
        try:
            # Use process lock for shared resource access
            with self._process_lock:
                return self.process_single_ligand(ligand_info, prepared_protein)
        except Exception as e:
            self.logger.error(f"Error in thread-safe ligand processing: {e}")
            return {
                'ligand_name': ligand_info[0],
                'success': False,
                'errors': {'thread_safe': str(e)},
                'total_time': 0
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
    
    def _register_cleanup_handlers(self):
        """Register cleanup handlers for graceful shutdown."""
        if not self._cleanup_registered:
            atexit.register(self._cleanup_resources)
            signal.signal(signal.SIGINT, self._signal_handler)
            signal.signal(signal.SIGTERM, self._signal_handler)
            self._cleanup_registered = True
    
    def _signal_handler(self, signum, frame):
        """Handle termination signals gracefully."""
        self.logger.info(f"Received signal {signum}, cleaning up resources...")
        self._cleanup_resources()
        exit(1)
    
    def _cleanup_resources(self):
        """Clean up shared resources and active processes."""
        try:
            with self._lock:
                # Clean up shared memory manager
                if self._shared_memory_manager:
                    self._shared_memory_manager.shutdown()
                    self._shared_memory_manager = None
                
                # Terminate active processes
                for process in self._active_processes.copy():
                    if process.is_alive():
                        process.terminate()
                        process.join(timeout=5)
                        if process.is_alive():
                            process.kill()
                
                self._active_processes.clear()
                
            self.logger.info("Resource cleanup completed")
        except Exception as e:
            self.logger.error(f"Error during cleanup: {e}")
    
    def _update_statistics(self, ligand_name: str, success: bool, result: dict):
        """Thread-safe update of pipeline statistics."""
        with self._lock:
            if success:
                self.successful_ligands += 1
                self.ligand_results[ligand_name] = result
            else:
                self.failed_ligands.add(ligand_name)
    
    def _register_process(self, process):
        """Register a process for cleanup tracking."""
        with self._lock:
            self._active_processes.add(process)
    
    def _unregister_process(self, process):
        """Unregister a process from cleanup tracking."""
        with self._lock:
            self._active_processes.discard(process)


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