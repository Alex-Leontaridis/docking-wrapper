#!/usr/bin/env python3
"""
Batch Molecular Docking Pipeline
Task 4: Orchestrates structure preparation, docking, and results parsing across multiple ligands
Supports AutoDock Vina, GNINA, DiffDock, and enhanced ML models with comprehensive logging and error handling
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
from multiprocessing import Pool, cpu_count, Manager, Lock
from functools import partial
import atexit
import signal
import threading
from multiprocessing.managers import SharedMemoryManager
from multiprocessing import Lock as ProcessLock
import platform

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import path manager and logging
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context, setup_global_error_handler, DockingLogger

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

# Import existing modules
from prep_structures import (
    prepare_protein, prepare_ligand_single, validate_file, 
    SUPPORTED_PROTEIN_EXT, SUPPORTED_LIGAND_EXT
)
from run_docking_multi import (
    run_vina, run_gnina, run_diffdock, extract_box_from_protein,
    ensure_output_dirs as ensure_single_output_dirs, find_binary
)
from docking_results_parser import DockingResultsParser
from config import DockingConfig

# Import enhanced ML/Analysis modules
try:
    from run_equibind import run_equibind_inference
    from run_neuralplexer import run_neuralplexer_inference
    from run_umol import run_umol
    from run_structure_predictor import predict_structure
    from run_boltz2 import run_boltz2_prediction
    from extract_interactions import run_plip
    from run_druggability import run_fpocket_analysis
#    from model_consensus import run_consensus_analysis
    from compute_confidence import compute_confidence_score
    ML_MODULES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Some ML modules not available: {e}")
    ML_MODULES_AVAILABLE = False

def process_single_ligand_standalone(ligand_info: Tuple[str, Path], prepared_protein: str, config: Dict, output_dir: str) -> Dict[str, Any]:
    """
    Standalone function for processing a single ligand in multiprocessing context.
    This function doesn't depend on instance variables and can be pickled.
    
    Args:
        ligand_info: Tuple of (ligand_name, ligand_path)
        prepared_protein: Path to prepared protein file
        config: Configuration dictionary
        output_dir: Output directory path
        
    Returns:
        Processing result dictionary
    """
    ligand_name, ligand_path = ligand_info
    result = {
        'ligand_name': ligand_name,
        'ligand_path': str(ligand_path),
        'success': False,
        'stages': {},
        'timings': {},
        'errors': {},
        'ml_results': {},
        'analysis_results': {}
    }
    
    try:
        # Import necessary functions here to avoid pickling issues
        from scripts.run_docking_multi import run_vina, run_gnina, run_diffdock, extract_box_from_protein
        from scripts.prep_structures import prepare_ligand_single, validate_file, SUPPORTED_LIGAND_EXT
        from utils.path_manager import get_path_manager
        
        # Setup logging for this process
        import logging
        logger = logging.getLogger(f"process_{ligand_name}")
        
        # Create output directory for this ligand
        ligand_output_dir = Path(output_dir) / "docking_results" / ligand_name
        ligand_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Stage 1: Ligand preparation
        logger.info(f"Processing ligand: {ligand_name}")
        stage_start = time.time()
        
        if not os.path.exists(ligand_path):
            result['stages']['preparation'] = False
            result['errors']['preparation'] = f"Ligand file not found: {ligand_path}"
            result['timings']['preparation'] = time.time() - stage_start
            return result
        
        # Prepare ligand
        prepared_dir = Path(output_dir) / "prepared_structures"
        prepared_ligand = prepared_dir / f"{ligand_name}_prepared.pdbqt"
        try:
            prepare_ligand_single(str(ligand_path), str(prepared_ligand))
            if not os.path.exists(prepared_ligand):
                raise Exception("Ligand preparation did not create output file")
            result['stages']['preparation'] = True
            result['timings']['preparation'] = time.time() - stage_start
        except Exception as e:
            result['stages']['preparation'] = False
            result['errors']['preparation'] = str(e)
            result['timings']['preparation'] = time.time() - stage_start
            return result
        
        # Stage 2: Traditional Docking
        enabled_engines = []
        if config["engines"]["vina"]["enabled"]:
            enabled_engines.append("vina")
        if config["engines"]["gnina"]["enabled"]:
            enabled_engines.append("gnina")
        if config["engines"]["diffdock"]["enabled"]:
            enabled_engines.append("diffdock")
        
        if enabled_engines:
            logger.info(f"Running traditional docking with engines: {enabled_engines}")
            
            # Get box parameters
            if config["box"]["auto_detect"]:
                try:
                    box_params = extract_box_from_protein(prepared_protein)
                except Exception as e:
                    logger.warning(f"Box auto-detection failed: {e}, using defaults")
                    default_size = config["box"]["default_size"]
                    box_params = (0.0, 0.0, 0.0, *default_size)
            else:
                default_size = config["box"]["default_size"]
                box_params = (0.0, 0.0, 0.0, *default_size)
            
            # Run each enabled docking engine
            for engine in enabled_engines:
                engine_start = time.time()
                logger.info(f"Running {engine.upper()} docking")
                try:
                    if engine == "vina":
                        status = run_vina(prepared_protein, str(prepared_ligand), str(ligand_output_dir), box_params)
                    elif engine == "gnina":
                        status = run_gnina(prepared_protein, str(prepared_ligand), str(ligand_output_dir), config["engines"]["gnina"]["use_gpu"])
                    elif engine == "diffdock":
                        status = run_diffdock(prepared_protein, str(prepared_ligand), str(ligand_output_dir))
                    
                    result['timings'][engine] = time.time() - engine_start
                    
                    if status['success']:
                        result['stages'][engine] = True
                        logger.info(f"{engine.upper()} completed successfully")
                    else:
                        result['stages'][engine] = False
                        result['errors'][engine] = status.get('error', 'Unknown error')
                        logger.error(f"{engine.upper()} failed: {result['errors'][engine]}")
                        
                except Exception as e:
                    result['stages'][engine] = False
                    result['errors'][engine] = str(e)
                    result['timings'][engine] = time.time() - engine_start
                    logger.error(f"{engine.upper()} failed with exception: {e}")
        
        # Mark as successful if at least one stage succeeded
        successful_stages = [stage for stage, success in result['stages'].items() if success]
        if successful_stages:
            result['success'] = True
            logger.info(f"Processing completed successfully. Successful stages: {successful_stages}")
        else:
            result['success'] = False
            logger.error("All stages failed")
        
        result['total_time'] = time.time() - stage_start
        return result
        
    except Exception as e:
        logger.error(f"Unexpected error in ligand processing: {e}")
        result['success'] = False
        result['errors']['unexpected'] = str(e)
        result['total_time'] = time.time() - stage_start
        return result

# --- Constants ---
DEFAULT_CONFIG = {
    "engines": {
        "vina": {"enabled": True, "exhaustiveness": 8, "num_modes": 9},
        "gnina": {"enabled": False, "use_gpu": False, "cnn_scoring": "rescore"},
        "diffdock": {"enabled": False, "inference_steps": 20, "samples_per_complex": 10}
    },
    "ml_models": {
        "equibind": {"enabled": False, "timeout": 600, "use_gpu": False},
        "neuralplexer": {"enabled": False, "timeout": 600, "use_gpu": False},
        "umol": {"enabled": False, "timeout": 600, "use_gpu": False},
        "structure_predictor": {"enabled": False, "timeout": 1800, "method": "colabfold"}
    },
    "analysis": {
        "boltz2": {"enabled": False, "timeout": 300},
        "interactions": {"enabled": False, "timeout": 300, "method": "plip"},
        "druggability": {"enabled": False, "timeout": 300},
        "consensus": {"enabled": False, "rmsd_threshold": 2.0},
        "confidence": {"enabled": False}
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
        "diffdock": 7200,
        "parsing": 300,
        "overall_per_ligand": 10800
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
    
    def __init__(self, config_path: Optional[str] = None, output_dir: str = None):
        """
        Initialize the batch docking pipeline.
        
        Args:
            config_path: Path to configuration JSON file
            output_dir: Output directory for results (optional, uses path manager if not provided)
        """
        # Setup global error handler
        setup_global_error_handler()
        
        # Initialize path manager
        self.path_manager = get_path_manager()
        
        # Setup logging first
        self.logger = setup_logging('BatchPipeline')
        log_startup('BatchDockingPipeline', '1.0.0')
        
        # Create individual tool loggers
        self.vina_logger = setup_logging('Vina', auto_log_file=False)
        self.gnina_logger = setup_logging('Gnina', auto_log_file=False)
        self.diffdock_logger = setup_logging('DiffDock', auto_log_file=False)
        
        # Set output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            output_dir_path = self.path_manager.get_absolute_path("outputs")
            if output_dir_path:
                self.output_dir = output_dir_path
            else:
                self.output_dir = Path("outputs")
        
        self.config = self._load_config(config_path)
        
        # Ensure output and log directories exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        log_dir = self.output_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        
        # Platform-specific availability check
        system = platform.system().lower()
        
        # GNINA: Linux only, but allow on Windows/Mac if WSL is available
        if system in ['windows', 'darwin']:  # Windows or Mac
            # Check if WSL is available
            try:
                result = subprocess.run(['wsl', '--list', '--quiet'], 
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0 and result.stdout.strip():
                    self.logger.info(f"WSL detected on {system.capitalize()}. GNINA will be available via WSL.")
                    # Keep GNINA enabled if WSL is available
                else:
                    self.config["engines"]["gnina"]["enabled"] = False
                    self.logger.warning(f"GNINA is not available on {system.capitalize()} (no WSL detected). GNINA features will be disabled.")
            except:
                self.config["engines"]["gnina"]["enabled"] = False
                self.logger.warning(f"GNINA is not available on {system.capitalize()} (WSL not available). GNINA features will be disabled.")
        
        # UMol: Linux only
        if system in ['windows', 'darwin']:  # Windows or Mac
            self.config["ml_models"]["umol"]["enabled"] = False
            self.logger.warning(f"UMol is not available on {system.capitalize()}. UMol features will be disabled.")
        elif system == 'linux':
            self.logger.info("UMol is available on Linux. UMol features will be enabled if configured.")
        
        # Log platform detection
        self.logger.info(f"Detected platform: {system.capitalize()}")
        self.logger.info(f"Available engines: Vina (all platforms), GNINA (Linux/WSL), DiffDock (all platforms)")
        self.logger.info(f"Available ML models: EquiBind (all platforms), NeuralPLexer (all platforms), UMol (Linux only)")
        
        self._ensure_output_dirs()
        
        # Process safety and resource management (using multiprocessing-safe locks)
        from multiprocessing import Lock, Manager
        self._lock = Lock()  # Multiprocessing-safe lock
        self._process_lock = Lock()  # Multiprocessing-safe lock
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
        import platform
        
        # Platform-specific binary names
        system = platform.system().lower()
        if system == 'windows':
            required_binaries = ["vina.exe", "gnina"]
        else:  # Linux/macOS
            required_binaries = ["vina", "gnina"]
        
        missing = []
        
        for binary in required_binaries:
            found_path = find_binary(binary)
            if found_path is None:
                missing.append(binary)
            else:
                self.logger.info(f"Found {binary}: {found_path}")
        
        # Check for DiffDock script (only if enabled in config)
        if self.config["engines"]["diffdock"]["enabled"]:
            diffdock_path = self.path_manager.get_path("diffdock")
            if not diffdock_path:
                missing.append('DiffDock/inference.py')
                self.logger.warning("DiffDock is enabled in config but not found. Disabling DiffDock.")
                self.config["engines"]["diffdock"]["enabled"] = False
            else:
                # Check if it's a file (inference.py) or directory containing inference.py
                if os.path.isfile(diffdock_path):
                    self.logger.info(f"Found DiffDock script: {diffdock_path}")
                elif os.path.isdir(diffdock_path):
                    inference_script = os.path.join(diffdock_path, 'inference.py')
                    if os.path.isfile(inference_script):
                        self.logger.info(f"Found DiffDock directory: {diffdock_path}")
                    else:
                        missing.append('DiffDock/inference.py')
                        self.logger.warning("DiffDock directory found but inference.py missing. Disabling DiffDock.")
                        self.config["engines"]["diffdock"]["enabled"] = False
                else:
                    missing.append('DiffDock/inference.py')
                    self.logger.warning("DiffDock path exists but is neither file nor directory. Disabling DiffDock.")
                    self.config["engines"]["diffdock"]["enabled"] = False
        else:
            self.logger.info("DiffDock is disabled in configuration - skipping DiffDock checks")
        
        # Check for MGLTools (optional, but warn if missing)
        self.config_obj = DockingConfig()
        if not self.config_obj.validate_mgltools():
            self.logger.warning("MGLTools not found. Protein preparation may be limited.")
            self.logger.info("To install MGLTools, visit: http://mgltools.scripps.edu/downloads/downloads/tools/downloads")
            self.logger.info("Or set MGLTOOLS_PATH environment variable to point to your installation.")
        
        if missing:
            self.logger.error(f"Missing required external binaries: {', '.join(missing)}")
            self.logger.error("Please install them and ensure they are in your PATH or current directory.")
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
                self.logger.warning(f"Failed to load config from {config_path}: {e}")
                self.logger.info("Using default configuration")
        
        return config
    
    def _deep_merge_config(self, base: Dict, update: Dict) -> Dict:
        """Deep merge configuration dictionaries."""
        for key, value in update.items():
            if key in base and isinstance(base[key], dict) and isinstance(value, dict):
                base[key] = self._deep_merge_config(base[key], value)
            else:
                base[key] = value
        return base
    
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
    
    def _setup_ligand_logging(self, ligand_name: str) -> Dict[str, DockingLogger]:
        """Setup individual tool loggers for a specific ligand."""
        ligand_log_dir = self.output_dir / "docking_results" / ligand_name / "logs"
        ligand_log_dir.mkdir(parents=True, exist_ok=True)
        
        # Create individual tool loggers with specific log files
        tool_loggers = {}
        
        # Vina logger
        vina_log_file = ligand_log_dir / "vina.log"
        tool_loggers['vina'] = setup_logging(
            f'Vina_{ligand_name}', 
            str(vina_log_file),
            auto_log_file=False
        )
        
        # Gnina logger
        gnina_log_file = ligand_log_dir / "gnina.log"
        tool_loggers['gnina'] = setup_logging(
            f'Gnina_{ligand_name}', 
            str(gnina_log_file),
            auto_log_file=False
        )
        
        # DiffDock logger
        diffdock_log_file = ligand_log_dir / "diffdock.log"
        tool_loggers['diffdock'] = setup_logging(
            f'DiffDock_{ligand_name}', 
            str(diffdock_log_file),
            auto_log_file=False
        )
        
        # Parsing logger
        parsing_log_file = ligand_log_dir / "parsing.log"
        tool_loggers['parsing'] = setup_logging(
            f'Parsing_{ligand_name}', 
            str(parsing_log_file),
            auto_log_file=False
        )
        
        return tool_loggers
    
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
            'errors': {},
            'ml_results': {},
            'analysis_results': {}
        }
        stage_start = time.time()
        try:
            with self.logger.timer(f"ligand_processing:{ligand_name}"):
                self.logger.info(f"Processing ligand: {ligand_name}")
                ligand_output_dir = self.output_dir / "docking_results" / ligand_name
                ligand_output_dir.mkdir(exist_ok=True)
                
                # Setup individual tool loggers for this ligand
                tool_loggers = self._setup_ligand_logging(ligand_name)
            
            # Stage 1: Ligand preparation
            self.logger.subsection(f"Stage 1: Preparing ligand structure")
            stage_start = time.time()
            if not os.path.exists(ligand_path):
                result['stages']['preparation'] = False
                result['errors']['preparation'] = f"Ligand file not found: {ligand_path}"
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.error(f"Ligand file not found: {ligand_path}", error_code='FILE_NOT_FOUND')
                return result  # Early exit
                
            prepared_ligand_dir = self.output_dir / "prepared_structures"
            prepared_ligand = prepared_ligand_dir / f"{ligand_name}_prepared.pdbqt"
            try:
                with self.logger.timer(f"ligand_preparation:{ligand_name}"):
                    prepare_ligand_single(str(ligand_path), str(prepared_ligand))
                    if not os.path.exists(prepared_ligand):
                        raise Exception("Ligand preparation did not create output file")
                result['stages']['preparation'] = True
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.info(f"Ligand prepared successfully")
            except Exception as e:
                result['stages']['preparation'] = False
                result['errors']['preparation'] = str(e)
                result['timings']['preparation'] = time.time() - stage_start
                self.logger.error(f"Ligand preparation failed: {e}", error_code='LIGAND_PREP_FAILED')
                return result  # Early exit
                
            ensure_single_output_dirs(str(ligand_output_dir))
            
            # Stage 2: Traditional Docking with enabled engines
            enabled_engines = []
            if self.config["engines"]["vina"]["enabled"]:
                enabled_engines.append("vina")
            if self.config["engines"]["gnina"]["enabled"]:
                enabled_engines.append("gnina")
            if self.config["engines"]["diffdock"]["enabled"]:
                enabled_engines.append("diffdock")
                
            if not enabled_engines and not any(self.config["ml_models"].values()):
                result['stages']['docking'] = False
                result['errors']['docking'] = "No docking engines or ML models enabled in configuration"
                result['timings']['docking'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] No docking engines or ML models enabled")
                return result  # Early exit
                
            if not os.path.exists(prepared_protein):
                result['stages']['docking'] = False
                result['errors']['docking'] = f"Prepared protein file not found: {prepared_protein}"
                result['timings']['docking'] = time.time() - stage_start
                self.logger.error(f"[{ligand_name}] Prepared protein file not found: {prepared_protein}")
                return result  # Early exit
                
            # Run traditional docking engines if enabled
            if enabled_engines:
                self.logger.info(f"[{ligand_name}] Stage 2: Running traditional docking with engines: {enabled_engines}")
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
                engine_failures = []
                for engine in enabled_engines:
                    engine_start = time.time()
                    self.logger.info(f"[{ligand_name}] Running {engine.upper()} docking")
                    try:
                        if engine == "vina":
                            # Log to individual vina log file
                            tool_loggers['vina'].info(f"Starting VINA docking for {ligand_name}")
                            tool_loggers['vina'].info(f"Protein: {prepared_protein}")
                            tool_loggers['vina'].info(f"Ligand: {prepared_ligand}")
                            tool_loggers['vina'].info(f"Box parameters: {box_params}")
                            
                            status = run_vina(
                                prepared_protein, 
                                str(prepared_ligand), 
                                str(ligand_output_dir), 
                                box_params
                            )
                            
                            if status['success']:
                                tool_loggers['vina'].info(f"VINA docking completed successfully in {status.get('time', 0):.2f}s")
                            else:
                                tool_loggers['vina'].error(f"VINA docking failed: {status.get('error', 'Unknown error')}")
                                
                        elif engine == "gnina":
                            # Log to individual gnina log file
                            tool_loggers['gnina'].info(f"Starting GNINA docking for {ligand_name}")
                            tool_loggers['gnina'].info(f"Protein: {prepared_protein}")
                            tool_loggers['gnina'].info(f"Ligand: {prepared_ligand}")
                            tool_loggers['gnina'].info(f"GPU enabled: {self.config['engines']['gnina']['use_gpu']}")
                            
                            status = run_gnina(
                                prepared_protein, 
                                str(prepared_ligand), 
                                str(ligand_output_dir),
                                self.config["engines"]["gnina"]["use_gpu"]
                            )
                            
                            if status['success']:
                                tool_loggers['gnina'].info(f"GNINA docking completed successfully in {status.get('time', 0):.2f}s")
                            else:
                                tool_loggers['gnina'].error(f"GNINA docking failed: {status.get('error', 'Unknown error')}")
                                
                        elif engine == "diffdock":
                            # Log to individual diffdock log file
                            tool_loggers['diffdock'].info(f"Starting DiffDock docking for {ligand_name}")
                            tool_loggers['diffdock'].info(f"Protein: {prepared_protein}")
                            tool_loggers['diffdock'].info(f"Ligand: {prepared_ligand}")
                            
                            status = run_diffdock(
                                prepared_protein, 
                                str(prepared_ligand), 
                                str(ligand_output_dir)
                            )
                            
                            if status['success']:
                                tool_loggers['diffdock'].info(f"DiffDock docking completed successfully in {status.get('time', 0):.2f}s")
                            else:
                                tool_loggers['diffdock'].error(f"DiffDock docking failed: {status.get('error', 'Unknown error')}")
                        
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
            else:
                result['stages']['docking'] = True  # Skip traditional docking if no engines enabled
            
            # Stage 3: ML Model Docking (if enabled and available)
            ml_enabled = any(model_config.get("enabled", False) for model_config in self.config["ml_models"].values())
            self.logger.info(f"[{ligand_name}] ML_MODULES_AVAILABLE: {ML_MODULES_AVAILABLE}, ml_enabled: {ml_enabled}")
            if ML_MODULES_AVAILABLE and ml_enabled:
                self.logger.info(f"[{ligand_name}] Stage 3: Running ML model docking")
                ml_start = time.time()
                
                # Create ML output directory
                ml_output_dir = ligand_output_dir / "ml_models"
                ml_output_dir.mkdir(exist_ok=True)
                
                # Run enabled ML models
                ml_poses = []
                ml_models_found = False
                for model_name, model_config in self.config["ml_models"].items():
                    if not model_config.get("enabled", False):
                        continue

                    # Check if model directory exists
                    model_dir = self.path_manager.get_path("models", model_name)
                    if not model_dir or not os.path.exists(model_dir):
                        self.logger.warning(f"[{ligand_name}] {model_name} model directory not found: {model_dir}. Skipping {model_name}.")
                        result['stages'][f"ml_{model_name}"] = False
                        result['errors'][f"ml_{model_name}"] = f"{model_name} model directory not found: {model_dir}"
                        continue
                    ml_models_found = True

                    model_start = time.time()
                    self.logger.info(f"[{ligand_name}] Running {model_name.upper()} ML docking")
                    
                    try:
                        model_output_dir = ml_output_dir / model_name
                        model_output_dir.mkdir(exist_ok=True)
                        
                        if model_name == "equibind":
                            success, output_file = run_equibind_inference(
                                prepared_protein, str(ligand_path), str(model_output_dir), model_config
                            )
                        elif model_name == "neuralplexer":
                            success, output_file = run_neuralplexer_inference(
                                prepared_protein, str(ligand_path), str(model_output_dir)
                            )
                        elif model_name == "umol":
                            umol_result = run_umol(
                                prepared_protein, str(ligand_path), str(model_output_dir)
                            )
                            success = umol_result.get('success', False)
                            output_file = umol_result.get('output_file') if success else None
                        elif model_name == "structure_predictor":
                            # Assume the input is a FASTA file for structure prediction
                            fasta_path = prepared_protein  # Or adjust as needed
                            struct_result = predict_structure(
                                fasta=fasta_path,
                                output_dir=str(model_output_dir),
                                model='auto'
                            )
                            success = struct_result.get('success', False)
                            output_file = struct_result.get('output_file') if success else None
                        else:
                            continue
                            
                        result['timings'][f"ml_{model_name}"] = time.time() - model_start
                        
                        if success and output_file:
                            result['stages'][f"ml_{model_name}"] = True
                            result['ml_results'][model_name] = output_file
                            ml_poses.append(output_file)
                            self.logger.info(f"[{ligand_name}] {model_name.upper()} completed successfully")
                        else:
                            result['stages'][f"ml_{model_name}"] = False
                            result['errors'][f"ml_{model_name}"] = f"{model_name} failed to produce output"
                            self.logger.warning(f"[{ligand_name}] {model_name.upper()} failed")
                            
                    except Exception as e:
                        result['stages'][f"ml_{model_name}"] = False
                        result['errors'][f"ml_{model_name}"] = str(e)
                        result['timings'][f"ml_{model_name}"] = time.time() - model_start
                        self.logger.error(f"[{ligand_name}] {model_name.upper()} failed with exception: {e}")
                if not ml_models_found:
                    self.logger.warning(f"[{ligand_name}] No ML model directories found. Skipping ML docking.")
                result['timings']['ml_docking'] = time.time() - ml_start
                result['ml_poses'] = ml_poses
            
            # Stage 4: Analysis (if enabled and available)
            if ML_MODULES_AVAILABLE and any(self.config["analysis"].values()):
                self.logger.info(f"[{ligand_name}] Stage 4: Running analysis tools")
                analysis_start = time.time()
                
                # Create analysis output directory
                analysis_output_dir = ligand_output_dir / "analysis"
                analysis_output_dir.mkdir(exist_ok=True)
                
                # Run enabled analysis tools
                for analysis_name, analysis_config in self.config["analysis"].items():
                    if not analysis_config.get("enabled", False):
                        continue
                        
                    analysis_tool_start = time.time()
                    self.logger.info(f"[{ligand_name}] Running {analysis_name.upper()} analysis")
                    
                    try:
                        tool_output_dir = analysis_output_dir / analysis_name
                        tool_output_dir.mkdir(exist_ok=True)
                        
                        if analysis_name == "boltz2":
                            success, output_file = run_boltz2_prediction(
                                prepared_protein, str(ligand_path), str(tool_output_dir)
                            )
                        elif analysis_name == "interactions":
                            # Use best pose for interaction analysis
                            best_pose = self._get_best_pose(ligand_output_dir, ml_poses)
                            if best_pose:
                                output_json = str(tool_output_dir / f"{ligand_name}_interactions.json")
                                ok, msg = run_plip(best_pose, output_json)
                                success = ok
                                output_file = output_json if ok else None
                            else:
                                success, output_file = False, None
                        elif analysis_name == "druggability":
                            # Use protein for druggability analysis
                            try:
                                result_fpocket = run_fpocket_analysis(prepared_protein, str(tool_output_dir))
                                success = result_fpocket.get('success', False)
                                output_file = result_fpocket.get('output_file') if success else None
                            except Exception as e:
                                success = False
                                output_file = None
                        elif analysis_name == "consensus":
                            # Use CLI of model_consensus.py for consensus analysis
                            try:
                                consensus_output_dir = ligand_output_dir / "consensus"
                                consensus_output_dir.mkdir(exist_ok=True)
                                consensus_file = consensus_output_dir / f"{ligand_name}_consensus.json"
                                pose_files = result.get('ml_poses', [])
                                if len(pose_files) < 2:
                                    success = False
                                    output_file = None
                                else:
                                    cmd = [
                                        sys.executable, "scripts/model_consensus.py",
                                        "--poses", *pose_files,
                                        "--output", str(consensus_file),
                                        "--ligand_id", ligand_name,
                                        "--rmsd_threshold", str(self.config["analysis"]["consensus"].get("rmsd_threshold", 2.0))
                                    ]
                                    proc = subprocess.run(cmd, capture_output=True, text=True)
                                    success = proc.returncode == 0 and os.path.exists(consensus_file)
                                    output_file = str(consensus_file) if success else None
                                    if not success:
                                        self.logger.error(f"Consensus analysis failed: {proc.stderr}")
                            except Exception as e:
                                success = False
                                output_file = None
                        elif analysis_name == "confidence":
                            # Use CLI of compute_confidence.py for confidence scoring
                            try:
                                confidence_output_dir = ligand_output_dir / "confidence"
                                confidence_output_dir.mkdir(exist_ok=True)
                                confidence_file = confidence_output_dir / f"{ligand_name}_confidence.json"
                                consensus_json = result.get('consensus_file')
                                druggability_json = result['analysis_results'].get('druggability') if result.get('analysis_results') else None
                                affinity_json = result['analysis_results'].get('boltz2') if result.get('analysis_results') else None
                                interaction_json = result['analysis_results'].get('interactions') if result.get('analysis_results') else None
                                if not all([consensus_json, druggability_json, affinity_json, interaction_json]):
                                    success = False
                                    output_file = None
                                else:
                                    cmd = [
                                        sys.executable, "scripts/compute_confidence.py",
                                        "--consensus_json", consensus_json,
                                        "--druggability_json", druggability_json,
                                        "--affinity_json", affinity_json,
                                        "--interaction_json", interaction_json,
                                        "--output", str(confidence_file),
                                        "--ligand_id", ligand_name
                                    ]
                                    proc = subprocess.run(cmd, capture_output=True, text=True)
                                    success = proc.returncode == 0 and os.path.exists(confidence_file)
                                    output_file = str(confidence_file) if success else None
                                    if not success:
                                        self.logger.error(f"Confidence scoring failed: {proc.stderr}")
                            except Exception as e:
                                success = False
                                output_file = None
                        else:
                            continue
                            
                        result['timings'][f"analysis_{analysis_name}"] = time.time() - analysis_tool_start
                        
                        if success and output_file:
                            result['stages'][f"analysis_{analysis_name}"] = True
                            result['analysis_results'][analysis_name] = output_file
                            self.logger.info(f"[{ligand_name}] {analysis_name.upper()} completed successfully")
                        else:
                            result['stages'][f"analysis_{analysis_name}"] = False
                            result['errors'][f"analysis_{analysis_name}"] = f"{analysis_name} failed to produce output"
                            self.logger.warning(f"[{ligand_name}] {analysis_name.upper()} failed")
                            
                    except Exception as e:
                        result['stages'][f"analysis_{analysis_name}"] = False
                        result['errors'][f"analysis_{analysis_name}"] = str(e)
                        result['timings'][f"analysis_{analysis_name}"] = time.time() - analysis_tool_start
                        self.logger.error(f"[{ligand_name}] {analysis_name.upper()} failed with exception: {e}")
                
                result['timings']['analysis'] = time.time() - analysis_start
            
            # Stage 5: Consensus Analysis (if enabled and poses available)
            if (ML_MODULES_AVAILABLE and 
                self.config["analysis"]["consensus"]["enabled"] and 
                len(result.get('ml_poses', [])) > 1):
                
                self.logger.info(f"[{ligand_name}] Stage 5: Running consensus analysis")
                consensus_start = time.time()
                
                try:
                    consensus_output_dir = ligand_output_dir / "consensus"
                    consensus_output_dir.mkdir(exist_ok=True)
                    
                    consensus_file = consensus_output_dir / f"{ligand_name}_consensus.json"
                    # run_consensus_analysis is not available on this platform
                    self.logger.warning(f"[{ligand_name}] Consensus analysis is not available (run_consensus_analysis not imported)")
                    success, output_file = False, None
                    result['timings']['consensus'] = time.time() - consensus_start
                    result['stages']['consensus'] = False
                    result['errors']['consensus'] = 'Consensus analysis not available on this platform.'
                except Exception as e:
                    result['timings']['consensus'] = time.time() - consensus_start
                    result['stages']['consensus'] = False
                    result['errors']['consensus'] = str(e)
            
            # Stage 6: Confidence Scoring (if enabled)
            if (ML_MODULES_AVAILABLE and 
                self.config["analysis"]["confidence"]["enabled"] and
                result.get('consensus_file') and
                result.get('analysis_results')):
                
                self.logger.info(f"[{ligand_name}] Stage 6: Computing confidence score")
                confidence_start = time.time()
                
                try:
                    confidence_output_dir = ligand_output_dir / "confidence"
                    confidence_output_dir.mkdir(exist_ok=True)
                    
                    confidence_file = confidence_output_dir / f"{ligand_name}_confidence.json"
                    success, output_file = compute_confidence_score(
                        result['consensus_file'],
                        result['analysis_results'].get('druggability'),
                        result['analysis_results'].get('boltz2'),
                        result['analysis_results'].get('interactions'),
                        str(confidence_file),
                        ligand_name
                    )
                    
                    result['timings']['confidence'] = time.time() - confidence_start
                    
                    if success and output_file:
                        result['stages']['confidence'] = True
                        result['confidence_file'] = output_file
                        self.logger.info(f"[{ligand_name}] Confidence scoring completed successfully")
                    else:
                        result['stages']['confidence'] = False
                        result['errors']['confidence'] = "Confidence scoring failed"
                        self.logger.warning(f"[{ligand_name}] Confidence scoring failed")
                        
                except Exception as e:
                    result['stages']['confidence'] = False
                    result['errors']['confidence'] = str(e)
                    result['timings']['confidence'] = time.time() - confidence_start
                    self.logger.error(f"[{ligand_name}] Confidence scoring failed with exception: {e}")
            
            # Stage 7: Results parsing (traditional engines)
            if enabled_engines:
                self.logger.info(f"[{ligand_name}] Stage 7: Parsing traditional docking results")
                parsing_start = time.time()
                try:
                    # Log to individual parsing log file
                    tool_loggers['parsing'].info(f"Starting results parsing for {ligand_name}")
                    tool_loggers['parsing'].info(f"Output directory: {ligand_output_dir}")
                    
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
                        
                        # Log parsing results to individual log
                        tool_loggers['parsing'].info(f"Parsing completed successfully")
                        tool_loggers['parsing'].info(f"Found {len(summary_df)} poses")
                        tool_loggers['parsing'].info(f"Summary saved to: {summary_file}")
                        if parser.failed_runs:
                            tool_loggers['parsing'].warning(f"Failed runs: {len(parser.failed_runs)}")
                        
                        self.logger.info(f"[{ligand_name}] Parsing completed, found {len(summary_df)} poses")
                    else:
                        result['stages']['parsing'] = False
                        result['errors']['parsing'] = "No poses found in results"
                        tool_loggers['parsing'].warning("No poses found in results")
                        self.logger.warning(f"[{ligand_name}] Parsing completed but no poses found")
                        return result  # Early exit
                    result['timings']['parsing'] = time.time() - parsing_start
                except Exception as e:
                    result['stages']['parsing'] = False
                    result['errors']['parsing'] = str(e)
                    result['timings']['parsing'] = time.time() - parsing_start
                    tool_loggers['parsing'].error(f"Results parsing failed: {e}")
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
    
    def _get_best_pose(self, ligand_output_dir: Path, ml_poses: List[str]) -> Optional[str]:
        """Get the best pose for analysis, prioritizing ML poses over traditional docking."""
        # First try to get the best ML pose
        if ml_poses:
            # For now, return the first ML pose (could be enhanced with scoring)
            return ml_poses[0]
        
        # Fallback to traditional docking poses
        vina_out = ligand_output_dir / "vina_output" / "vina_out.pdbqt"
        if vina_out.exists():
            return str(vina_out)
        
        gnina_out = ligand_output_dir / "gnina_output"
        if gnina_out.exists():
            gnina_files = list(gnina_out.glob("*.pdbqt")) + list(gnina_out.glob("*.sdf"))
            if gnina_files:
                return str(gnina_files[0])
        
        diffdock_out = ligand_output_dir / "diffdock_output"
        if diffdock_out.exists():
            diffdock_files = list(diffdock_out.glob("*.sdf"))
            if diffdock_files:
                return str(diffdock_files[0])
        
        return None
    
    def print_final_summary(self, results: List[Dict]):
        print(f"\n{COLOR_INFO}{'='*80}{COLOR_RESET}")
        print(f"{COLOR_INFO}ENHANCED BATCH PIPELINE SUMMARY{COLOR_RESET}")
        print(f"{COLOR_INFO}{'='*80}{COLOR_RESET}")
        
        # Determine which columns to show based on what was actually run
        all_stages = set()
        for r in results:
            all_stages.update(r.get('stages', {}).keys())
        
        # Build headers dynamically
        headers = ["Ligand", "Status", "Prep"]
        
        # Traditional docking engines
        if any(stage in all_stages for stage in ['vina', 'gnina', 'diffdock']):
            headers.extend(["Vina", "GNINA", "DiffDock"])
        
        # ML models
        ml_models = ['ml_equibind', 'ml_neuralplexer', 'ml_umol', 'ml_structure_predictor']
        if any(stage in all_stages for stage in ml_models):
            headers.extend(["EquiBind", "NeuralPLexer", "UMol", "StructPred"])
        
        # Analysis tools
        analysis_tools = ['analysis_boltz2', 'analysis_interactions', 'analysis_druggability']
        if any(stage in all_stages for stage in analysis_tools):
            headers.extend(["Boltz2", "Interactions", "Druggability"])
        
        # Consensus and confidence
        if 'consensus' in all_stages:
            headers.append("Consensus")
        if 'confidence' in all_stages:
            headers.append("Confidence")
        
        # Traditional parsing
        if 'parsing' in all_stages:
            headers.append("Parse")
        
        headers.extend(["Time (s)", "Error"])
        
        table = []
        for r in results:
            stages = r.get('stages', {})
            row = [
                r.get('ligand_name', ''),
                pretty_status(r.get('success', False)),
                pretty_status(stages.get('preparation', False))
            ]
            
            # Traditional docking engines
            if any(stage in all_stages for stage in ['vina', 'gnina', 'diffdock']):
                row.extend([
                    pretty_status(stages.get('vina', False)) if 'vina' in stages else '-',
                    pretty_status(stages.get('gnina', False)) if 'gnina' in stages else '-',
                    pretty_status(stages.get('diffdock', False)) if 'diffdock' in stages else '-'
                ])
            
            # ML models
            if any(stage in all_stages for stage in ml_models):
                row.extend([
                    pretty_status(stages.get('ml_equibind', False)) if 'ml_equibind' in stages else '-',
                    pretty_status(stages.get('ml_neuralplexer', False)) if 'ml_neuralplexer' in stages else '-',
                    pretty_status(stages.get('ml_umol', False)) if 'ml_umol' in stages else '-',
                    pretty_status(stages.get('ml_structure_predictor', False)) if 'ml_structure_predictor' in stages else '-'
                ])
            
            # Analysis tools
            if any(stage in all_stages for stage in analysis_tools):
                row.extend([
                    pretty_status(stages.get('analysis_boltz2', False)) if 'analysis_boltz2' in stages else '-',
                    pretty_status(stages.get('analysis_interactions', False)) if 'analysis_interactions' in stages else '-',
                    pretty_status(stages.get('analysis_druggability', False)) if 'analysis_druggability' in stages else '-'
                ])
            
            # Consensus and confidence
            if 'consensus' in all_stages:
                row.append(pretty_status(stages.get('consensus', False)) if 'consensus' in stages else '-')
            if 'confidence' in all_stages:
                row.append(pretty_status(stages.get('confidence', False)) if 'confidence' in stages else '-')
            
            # Traditional parsing
            if 'parsing' in all_stages:
                row.append(pretty_status(stages.get('parsing', False)) if 'parsing' in stages else '-')
            
            row.extend([
                f"{r.get('total_time', 0):.1f}",
                (next(iter(r.get('errors', {}).values()), '')[:40] if r.get('errors') else '')
            ])
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
        
        # Print ML/Analysis summary if any were run
        ml_results = [r for r in results if r.get('ml_results') or r.get('analysis_results')]
        if ml_results:
            print(f"\n{COLOR_INFO}ML/Analysis Results Summary:{COLOR_RESET}")
            for r in ml_results:
                ligand_name = r.get('ligand_name', '')
                ml_count = len(r.get('ml_results', {}))
                analysis_count = len(r.get('analysis_results', {}))
                if ml_count > 0 or analysis_count > 0:
                    print(f"  {ligand_name}: {ml_count} ML poses, {analysis_count} analysis results")
        
        print(f"{COLOR_INFO}{'='*80}{COLOR_RESET}\n")
    
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
                
                # Use process-safe processing with proper resource management
                with Pool(processes=self.config["parallel"]["max_workers"]) as pool:
                    # Create a process-safe wrapper for the processing function
                    # Pass necessary data as arguments instead of using instance variables
                    process_func = partial(process_single_ligand_standalone, 
                                         prepared_protein=prepared_protein,
                                         config=self.config,
                                         output_dir=str(self.output_dir))
                    
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
                                
                                # Update statistics (no need for thread safety in main process)
                                ligand_name = ligands[i][0]
                                success = processed_result.get('success', False)
                                if success:
                                    self.successful_ligands += 1
                                    self.ligand_results[ligand_name] = processed_result
                                else:
                                    self.failed_ligands.add(ligand_name)
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
                                self.failed_ligands.add(ligand_name)
                                
                        except Exception as e:
                            self.logger.error(f"Failed to get result for ligand {ligands[i][0]}: {e}")
                            failed_result = {
                                'ligand_name': ligands[i][0],
                                'success': False,
                                'errors': {'processing': str(e)},
                                'total_time': 0
                            }
                            processed_results.append(failed_result)
                            self.failed_ligands.add(ligands[i][0])
                
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
                        if success:
                            self.successful_ligands += 1
                            self.ligand_results[ligand_name] = result
                        else:
                            self.failed_ligands.add(ligand_name)
                        
                    except Exception as e:
                        self.logger.error(f"Failed to process ligand {ligand_info[0]}: {e}")
                        failed_result = {
                            'ligand_name': ligand_info[0],
                            'success': False,
                            'errors': {'processing': str(e)},
                            'total_time': 0
                        }
                        results.append(failed_result)
                        self.failed_ligands.add(ligand_info[0])
            
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
        description="Enhanced Batch Molecular Docking Pipeline with ML/Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Traditional docking only
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligand_dir/
  
  # With custom config
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --config pipeline_config.json
  
  # Enable traditional engines
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-gnina --enable-diffdock
  
  # Enable ML models
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-equibind --enable-neuralplexer
  
  # Enable analysis tools
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-analysis --enable-consensus
  
  # Full enhanced pipeline
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-all-ml --enable-all-analysis
  
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
    
    # Traditional docking engines
    parser.add_argument('--enable-gnina', action='store_true',
                       help='Enable GNINA docking')
    parser.add_argument('--enable-diffdock', action='store_true',
                       help='Enable DiffDock docking')
    
    # ML models
    parser.add_argument('--enable-equibind', action='store_true',
                       help='Enable EquiBind ML pose prediction')
    parser.add_argument('--enable-neuralplexer', action='store_true',
                       help='Enable NeuralPLexer ML pose prediction')
    parser.add_argument('--enable-umol', action='store_true',
                       help='Enable UMol structure-free pose prediction')
    parser.add_argument('--enable-structure-predictor', action='store_true',
                       help='Enable protein structure prediction')
    parser.add_argument('--enable-all-ml', action='store_true',
                       help='Enable all ML models')
    
    # Analysis tools
    parser.add_argument('--enable-boltz2', action='store_true',
                       help='Enable Boltz2 binding affinity prediction')
    parser.add_argument('--enable-interactions', action='store_true',
                       help='Enable protein-ligand interaction analysis')
    parser.add_argument('--enable-druggability', action='store_true',
                       help='Enable binding site druggability analysis')
    parser.add_argument('--enable-all-analysis', action='store_true',
                       help='Enable all analysis tools')
    
    # Consensus and confidence
    parser.add_argument('--enable-consensus', action='store_true',
                       help='Enable pose consensus analysis')
    parser.add_argument('--enable-confidence', action='store_true',
                       help='Enable confidence scoring')
    
    # Processing options
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
    # Traditional engines
    if args.enable_gnina:
        pipeline.config["engines"]["gnina"]["enabled"] = True
    if args.enable_diffdock:
        pipeline.config["engines"]["diffdock"]["enabled"] = True
    
    # ML models
    if args.enable_equibind or args.enable_all_ml:
        pipeline.config["ml_models"]["equibind"]["enabled"] = True
    if args.enable_neuralplexer or args.enable_all_ml:
        pipeline.config["ml_models"]["neuralplexer"]["enabled"] = True
    if args.enable_umol or args.enable_all_ml:
        pipeline.config["ml_models"]["umol"]["enabled"] = True
    if args.enable_structure_predictor or args.enable_all_ml:
        pipeline.config["ml_models"]["structure_predictor"]["enabled"] = True
    
    # Analysis tools
    if args.enable_boltz2 or args.enable_all_analysis:
        pipeline.config["analysis"]["boltz2"]["enabled"] = True
    if args.enable_interactions or args.enable_all_analysis:
        pipeline.config["analysis"]["interactions"]["enabled"] = True
    if args.enable_druggability or args.enable_all_analysis:
        pipeline.config["analysis"]["druggability"]["enabled"] = True
    
    # Consensus and confidence
    if args.enable_consensus or args.enable_all_analysis:
        pipeline.config["analysis"]["consensus"]["enabled"] = True
    if args.enable_confidence or args.enable_all_analysis:
        pipeline.config["analysis"]["confidence"]["enabled"] = True
    
    # Processing options
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