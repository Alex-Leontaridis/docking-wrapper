#!/usr/bin/env python3
"""
Traditional Docking Engines Implementation

This module implements concrete docking engines that inherit from the DockingEngine interface.
"""

import os
import sys
import subprocess
import time
import shutil
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
import logging
import json
import platform
import structlog
from utils.metrics import DOCKING_DURATION, INFERENCE_DURATION, PIPELINE_ERRORS
import sentry_sdk

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .interfaces import DockingEngine, DockingResult, ComponentStatus
from config import config
from utils.path_manager import get_path_manager, ensure_dir
from utils.logging import setup_logging
from utils.metrics import start_metrics_server
from utils.sentry import init_sentry
from config import SENTRY_DSN, METRICS_PORT

# Initialize observability
setup_logging()
if SENTRY_DSN:
    init_sentry(SENTRY_DSN)
start_metrics_server(METRICS_PORT)


class AutoDockVina(DockingEngine):
    """AutoDock Vina docking engine implementation."""
    
    def __init__(self, name: str = "vina", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.vina_path = None
        self.mgltools_path = None
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to Vina and MGLTools."""
        # Get paths from config
        self.vina_path = config.vina_path
        self.mgltools_path = config.mgltools_path
        
        # Fallback to path manager
        if not self.vina_path:
            path_manager = get_path_manager()
            self.vina_path = path_manager.get_binary_path('vina')
    
    def check_availability(self) -> ComponentStatus:
        """Check if Vina is available."""
        if not self.vina_path or not os.path.exists(self.vina_path):
            self.logger.warning("Vina executable not found")
            return ComponentStatus.UNAVAILABLE
        
        # Test if vina is executable
        try:
            result = subprocess.run([self.vina_path, "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing Vina: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get Vina capabilities."""
        return {
            "name": "AutoDock Vina",
            "version": "1.2.0",
            "supports_gpu": False,
            "supports_multiple_poses": True,
            "max_poses": 9,
            "supports_flexible_docking": False,
            "supports_ensemble_docking": False
        }
    
    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get supported file formats."""
        return {
            "receptor": [".pdbqt"],
            "ligand": [".pdbqt", ".mol2"]
        }
    
    def prepare_inputs(self, receptor: Path, ligand: Path, **kwargs) -> Dict[str, Path]:
        """Prepare input files for Vina docking."""
        # Vina requires PDBQT format
        prepared_receptor = receptor
        prepared_ligand = ligand
        
        # Convert receptor to PDBQT if needed
        if receptor.suffix.lower() != '.pdbqt':
            prepared_receptor = self._convert_to_pdbqt(receptor, is_receptor=True)
        
        # Convert ligand to PDBQT if needed
        if ligand.suffix.lower() != '.pdbqt':
            prepared_ligand = self._convert_to_pdbqt(ligand, is_receptor=False)
        
        return {
            "receptor": prepared_receptor,
            "ligand": prepared_ligand
        }
    
    def _convert_to_pdbqt(self, input_file: Path, is_receptor: bool = True) -> Path:
        """Convert file to PDBQT format using MGLTools."""
        if not self.mgltools_path:
            raise RuntimeError("MGLTools not available for PDBQT conversion")
        
        output_file = input_file.with_suffix('.pdbqt')
        
        if is_receptor:
            script = os.path.join(self.mgltools_path, 'MGLToolsPckgs', 'AutoDockTools', 
                                'Utilities24', 'prepare_receptor4.py')
            cmd = [
                os.path.join(self.mgltools_path, 'bin', 'pythonsh'),
                script,
                '-r', str(input_file),
                '-o', str(output_file),
                '-A', 'bonds'
            ]
        else:
            script = os.path.join(self.mgltools_path, 'MGLToolsPckgs', 'AutoDockTools', 
                                'Utilities24', 'prepare_ligand4.py')
            cmd = [
                os.path.join(self.mgltools_path, 'bin', 'pythonsh'),
                script,
                '-l', str(input_file),
                '-o', str(output_file),
                '-A', 'bonds'
            ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode != 0:
                raise RuntimeError(f"PDBQT conversion failed: {result.stderr}")
            return output_file
        except Exception as e:
            raise RuntimeError(f"Error converting to PDBQT: {e}")
    
    def run(self, receptor: Path, ligand: Path, **kwargs) -> DockingResult:
        logger = structlog.get_logger(module=__name__, engine=self.name)
        start_time = time.time()
        with DOCKING_DURATION.labels(engine=self.name).time():
            try:
                self.set_status(ComponentStatus.RUNNING)
                logger.info("Starting docking", receptor=str(receptor), ligand=str(ligand))
                # Validate inputs
                if not receptor.exists():
                    PIPELINE_ERRORS.labels(stage="vina_input").inc()
                    logger.error("Receptor file not found", receptor=str(receptor))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Receptor file not found: {receptor}"
                    )
                if not ligand.exists():
                    PIPELINE_ERRORS.labels(stage="vina_input").inc()
                    logger.error("Ligand file not found", ligand=str(ligand))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Ligand file not found: {ligand}"
                    )
                # Prepare inputs
                prepared_inputs = self.prepare_inputs(receptor, ligand, **kwargs)
                # Get docking parameters
                center = kwargs.get('center', (0, 0, 0))
                size = kwargs.get('size', (25, 25, 25))
                exhaustiveness = kwargs.get('exhaustiveness', 8)
                num_poses = kwargs.get('num_poses', 9)
                energy_range = kwargs.get('energy_range', 3)
                # Create output directory
                output_dir = kwargs.get('output_dir', Path.cwd() / 'vina_output')
                ensure_dir(output_dir)
                # Prepare Vina command
                output_file = output_dir / f"{ligand.stem}_vina_out.pdbqt"
                log_file = output_dir / f"{ligand.stem}_vina.log"
                cmd = [
                    self.vina_path,
                    '--receptor', str(prepared_inputs['receptor']),
                    '--ligand', str(prepared_inputs['ligand']),
                    '--out', str(output_file),
                    '--log', str(log_file),
                    '--center_x', str(center[0]),
                    '--center_y', str(center[1]),
                    '--center_z', str(center[2]),
                    '--size_x', str(size[0]),
                    '--size_y', str(size[1]),
                    '--size_z', str(size[2]),
                    '--exhaustiveness', str(exhaustiveness),
                    '--num_modes', str(num_poses),
                    '--energy_range', str(energy_range)
                ]
                logger.info("Running Vina", command=cmd)
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                if result.returncode != 0:
                    PIPELINE_ERRORS.labels(stage="vina_run").inc()
                    logger.error("Vina failed", stderr=result.stderr)
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Vina failed: {result.stderr}",
                        execution_time=time.time() - start_time
                    )
                # Parse results
                poses, scores = self._parse_vina_results(output_file, log_file)
                self.set_status(ComponentStatus.COMPLETED)
                logger.info("Docking completed", poses=len(poses), scores=scores)
                return DockingResult(
                    success=True,
                    poses=poses,
                    scores=scores,
                    metadata={
                        "engine": "AutoDock Vina",
                        "exhaustiveness": exhaustiveness,
                        "num_poses": num_poses,
                        "center": center,
                        "size": size
                    },
                    execution_time=time.time() - start_time
                )
            except Exception as e:
                self.set_status(ComponentStatus.FAILED)
                PIPELINE_ERRORS.labels(stage="vina_exception").inc()
                logger.exception("Exception in Vina run", error=str(e))
                sentry_sdk.capture_exception(e)
                return DockingResult(
                    success=False,
                    poses=[],
                    scores=[],
                    error_message=str(e),
                    execution_time=time.time() - start_time
                )
    
    def _parse_vina_results(self, output_file: Path, log_file: Path) -> Tuple[List[Path], List[float]]:
        """Parse Vina output files to extract poses and scores."""
        poses = []
        scores = []
        
        if not output_file.exists():
            return poses, scores
        
        # Read scores from log file
        if log_file.exists():
            with open(log_file, 'r') as f:
                for line in f:
                    if line.startswith('-----'):
                        continue
                    if line.strip() and not line.startswith('REMARK'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                score = float(parts[1])
                                scores.append(score)
                            except ValueError:
                                continue
        
        # Split output file into individual poses
        if output_file.exists():
            pose_dir = output_file.parent / 'poses'
            ensure_dir(pose_dir)
            
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Split by MODEL/ENDMDL
            pose_blocks = content.split('MODEL')
            
            for i, block in enumerate(pose_blocks[1:], 1):  # Skip first empty block
                if block.strip():
                    pose_file = pose_dir / f"pose_{i}.pdbqt"
                    with open(pose_file, 'w') as f:
                        f.write(f"MODEL{block}")
                    poses.append(pose_file)
        
        return poses, scores


class GNINA(DockingEngine):
    """GNINA docking engine implementation."""
    
    def __init__(self, name: str = "gnina", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.gnina_path = None
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to GNINA."""
        # Get path from config
        self.gnina_path = config.gnina_path
        
        # Fallback to path manager
        if not self.gnina_path:
            path_manager = get_path_manager()
            self.gnina_path = path_manager.get_binary_path('gnina')
    
    def check_availability(self) -> ComponentStatus:
        """Check if GNINA is available."""
        if not self.gnina_path or not os.path.exists(self.gnina_path):
            self.logger.warning("GNINA executable not found")
            return ComponentStatus.UNAVAILABLE
        
        # Test if gnina is executable
        try:
            result = subprocess.run([self.gnina_path, "--help"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing GNINA: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get GNINA capabilities."""
        return {
            "name": "GNINA",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_poses": True,
            "max_poses": 10,
            "supports_flexible_docking": True,
            "supports_ensemble_docking": True,
            "supports_cnn_scoring": True
        }
    
    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get supported file formats."""
        return {
            "receptor": [".pdb", ".pdbqt"],
            "ligand": [".sdf", ".mol2", ".pdbqt"]
        }
    
    def run(self, receptor: Path, ligand: Path, **kwargs) -> DockingResult:
        logger = structlog.get_logger(module=__name__, engine=self.name)
        start_time = time.time()
        with DOCKING_DURATION.labels(engine=self.name).time():
            try:
                self.set_status(ComponentStatus.RUNNING)
                logger.info("Starting docking", receptor=str(receptor), ligand=str(ligand))
                # Validate inputs
                if not receptor.exists():
                    PIPELINE_ERRORS.labels(stage="gnina_input").inc()
                    logger.error("Receptor file not found", receptor=str(receptor))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Receptor file not found: {receptor}"
                    )
                if not ligand.exists():
                    PIPELINE_ERRORS.labels(stage="gnina_input").inc()
                    logger.error("Ligand file not found", ligand=str(ligand))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Ligand file not found: {ligand}"
                    )
                # Get docking parameters
                center = kwargs.get('center', (0, 0, 0))
                size = kwargs.get('size', (25, 25, 25))
                exhaustiveness = kwargs.get('exhaustiveness', 8)
                num_poses = kwargs.get('num_poses', 10)
                use_gpu = kwargs.get('use_gpu', False)
                use_cnn = kwargs.get('use_cnn', True)
                # Create output directory
                output_dir = kwargs.get('output_dir', Path.cwd() / 'gnina_output')
                ensure_dir(output_dir)
                # Prepare GNINA command
                output_file = output_dir / f"{ligand.stem}_gnina_out.sdf"
                log_file = output_dir / f"{ligand.stem}_gnina.log"
                cmd = [
                    self.gnina_path,
                    '--receptor', str(receptor),
                    '--ligand', str(ligand),
                    '--out', str(output_file),
                    '--log', str(log_file),
                    '--center_x', str(center[0]),
                    '--center_y', str(center[1]),
                    '--center_z', str(center[2]),
                    '--size_x', str(size[0]),
                    '--size_y', str(size[1]),
                    '--size_z', str(size[2]),
                    '--exhaustiveness', str(exhaustiveness),
                    '--num_modes', str(num_poses)
                ]
                if use_gpu:
                    cmd.append('--gpu')
                if use_cnn:
                    cmd.append('--cnn')
                logger.info("Running GNINA", command=cmd)
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                if result.returncode != 0:
                    PIPELINE_ERRORS.labels(stage="gnina_run").inc()
                    logger.error("GNINA failed", stderr=result.stderr)
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"GNINA failed: {result.stderr}",
                        execution_time=time.time() - start_time
                    )
                # Parse results
                poses, scores = self._parse_gnina_results(output_file, log_file)
                self.set_status(ComponentStatus.COMPLETED)
                logger.info("Docking completed", poses=len(poses), scores=scores)
                return DockingResult(
                    success=True,
                    poses=poses,
                    scores=scores,
                    metadata={
                        "engine": "GNINA",
                        "exhaustiveness": exhaustiveness,
                        "num_poses": num_poses,
                        "use_gpu": use_gpu,
                        "use_cnn": use_cnn,
                        "center": center,
                        "size": size
                    },
                    execution_time=time.time() - start_time
                )
            except Exception as e:
                self.set_status(ComponentStatus.FAILED)
                PIPELINE_ERRORS.labels(stage="gnina_exception").inc()
                logger.exception("Exception in GNINA run", error=str(e))
                sentry_sdk.capture_exception(e)
                return DockingResult(
                    success=False,
                    poses=[],
                    scores=[],
                    error_message=str(e),
                    execution_time=time.time() - start_time
                )
    
    def _parse_gnina_results(self, output_file: Path, log_file: Path) -> Tuple[List[Path], List[float]]:
        """Parse GNINA output files to extract poses and scores."""
        poses = []
        scores = []
        
        if not output_file.exists():
            return poses, scores
        
        # Read scores from log file
        if log_file.exists():
            with open(log_file, 'r') as f:
                for line in f:
                    if line.startswith('-----'):
                        continue
                    if line.strip() and not line.startswith('REMARK'):
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                score = float(parts[1])
                                scores.append(score)
                            except ValueError:
                                continue
        
        # Split output file into individual poses
        if output_file.exists():
            pose_dir = output_file.parent / 'poses'
            ensure_dir(pose_dir)
            
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Split by $$$$ (SDF separator)
            pose_blocks = content.split('$$$$')
            
            for i, block in enumerate(pose_blocks[:-1], 1):  # Skip last empty block
                if block.strip():
                    pose_file = pose_dir / f"pose_{i}.sdf"
                    with open(pose_file, 'w') as f:
                        f.write(block.strip() + '\n$$$$\n')
                    poses.append(pose_file)
        
        return poses, scores


class DiffDock(DockingEngine):
    """DiffDock docking engine implementation."""
    
    def __init__(self, name: str = "diffdock", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.diffdock_path = None
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to DiffDock."""
        # Get path from config
        self.diffdock_path = config.diffdock_path
        
        # Fallback to path manager
        if not self.diffdock_path:
            path_manager = get_path_manager()
            self.diffdock_path = path_manager.get_binary_path('diffdock')
    
    def check_availability(self) -> ComponentStatus:
        """Check if DiffDock is available."""
        if not self.diffdock_path:
            self.logger.warning("DiffDock not found")
            return ComponentStatus.UNAVAILABLE
        
        # Check if it's a script or directory
        if os.path.isfile(self.diffdock_path):
            # Check if it's executable
            if not os.access(self.diffdock_path, os.X_OK):
                self.logger.warning("DiffDock script is not executable")
                return ComponentStatus.UNAVAILABLE
        elif os.path.isdir(self.diffdock_path):
            # Check if directory contains inference.py
            inference_script = os.path.join(self.diffdock_path, 'inference.py')
            if not os.path.exists(inference_script):
                self.logger.warning("DiffDock directory does not contain inference.py")
                return ComponentStatus.UNAVAILABLE
        
        # Test if diffdock is working
        try:
            if os.path.isfile(self.diffdock_path):
                result = subprocess.run([sys.executable, self.diffdock_path, '--help'], 
                                      capture_output=True, text=True, timeout=30)
            else:
                inference_script = os.path.join(self.diffdock_path, 'inference.py')
                result = subprocess.run([sys.executable, inference_script, '--help'], 
                                      capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing DiffDock: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get DiffDock capabilities."""
        return {
            "name": "DiffDock",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_poses": True,
            "max_poses": 10,
            "supports_flexible_docking": False,
            "supports_ensemble_docking": False,
            "supports_cnn_scoring": False,
            "uses_diffusion": True
        }
    
    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get supported file formats."""
        return {
            "receptor": [".pdb"],
            "ligand": [".sdf", ".mol2"]
        }
    
    def run(self, receptor: Path, ligand: Path, **kwargs) -> DockingResult:
        logger = structlog.get_logger(module=__name__, engine=self.name)
        start_time = time.time()
        with DOCKING_DURATION.labels(engine=self.name).time():
            try:
                self.set_status(ComponentStatus.RUNNING)
                logger.info("Starting docking", receptor=str(receptor), ligand=str(ligand))
                # Validate inputs
                if not receptor.exists():
                    PIPELINE_ERRORS.labels(stage="diffdock_input").inc()
                    logger.error("Receptor file not found", receptor=str(receptor))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Receptor file not found: {receptor}"
                    )
                if not ligand.exists():
                    PIPELINE_ERRORS.labels(stage="diffdock_input").inc()
                    logger.error("Ligand file not found", ligand=str(ligand))
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"Ligand file not found: {ligand}"
                    )
                # Get docking parameters
                num_poses = kwargs.get('num_poses', 10)
                confidence = kwargs.get('confidence', 0.5)
                use_gpu = kwargs.get('use_gpu', True)
                # Create output directory
                output_dir = kwargs.get('output_dir', Path.cwd() / 'diffdock_output')
                ensure_dir(output_dir)
                # Prepare DiffDock command
                output_file = output_dir / f"{ligand.stem}_diffdock_out.sdf"
                confidence_file = output_dir / f"{ligand.stem}_confidence.txt"
                
                if os.path.isfile(self.diffdock_path):
                    cmd = [
                        sys.executable, self.diffdock_path,
                        '--protein_path', str(receptor),
                        '--ligand', str(ligand),
                        '--out_dir', str(output_dir),
                        '--inference_steps', '20',
                        '--samples_per_complex', str(num_poses),
                        '--confidence_cutoff', str(confidence)
                    ]
                else:
                    inference_script = os.path.join(self.diffdock_path, 'inference.py')
                    cmd = [
                        sys.executable, inference_script,
                        '--protein_path', str(receptor),
                        '--ligand', str(ligand),
                        '--out_dir', str(output_dir),
                        '--inference_steps', '20',
                        '--samples_per_complex', str(num_poses),
                        '--confidence_cutoff', str(confidence)
                    ]
                
                if use_gpu:
                    cmd.append('--gpu')
                
                logger.info("Running DiffDock", command=cmd)
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)  # 2 hours timeout
                
                if result.returncode != 0:
                    PIPELINE_ERRORS.labels(stage="diffdock_run").inc()
                    logger.error("DiffDock failed", stderr=result.stderr)
                    return DockingResult(
                        success=False,
                        poses=[],
                        scores=[],
                        error_message=f"DiffDock failed: {result.stderr}",
                        execution_time=time.time() - start_time
                    )
                # Parse results
                poses, scores = self._parse_diffdock_results(output_dir, confidence_file)
                self.set_status(ComponentStatus.COMPLETED)
                logger.info("Docking completed", poses=len(poses), scores=scores)
                return DockingResult(
                    success=True,
                    poses=poses,
                    scores=scores,
                    metadata={
                        "engine": "DiffDock",
                        "num_poses": num_poses,
                        "confidence_cutoff": confidence,
                        "use_gpu": use_gpu
                    },
                    execution_time=time.time() - start_time
                )
            except Exception as e:
                self.set_status(ComponentStatus.FAILED)
                PIPELINE_ERRORS.labels(stage="diffdock_exception").inc()
                logger.exception("Exception in DiffDock run", error=str(e))
                sentry_sdk.capture_exception(e)
                return DockingResult(
                    success=False,
                    poses=[],
                    scores=[],
                    error_message=str(e),
                    execution_time=time.time() - start_time
                )
    
    def _parse_diffdock_results(self, output_dir: Path, confidence_file: Path) -> Tuple[List[Path], List[float]]:
        """Parse DiffDock output files to extract poses and scores."""
        poses = []
        scores = []
        
        # Read confidence scores
        confidence_scores = []
        if confidence_file.exists():
            with open(confidence_file, 'r') as f:
                for line in f:
                    try:
                        score = float(line.strip())
                        confidence_scores.append(score)
                    except ValueError:
                        continue
        
        # Find pose files
        pose_files = list(output_dir.glob('pose_*.sdf'))
        pose_files.sort(key=lambda x: int(x.stem.split('_')[1]))
        
        for i, pose_file in enumerate(pose_files):
            poses.append(pose_file)
            if i < len(confidence_scores):
                scores.append(confidence_scores[i])
            else:
                scores.append(0.0)  # Default score if confidence not available
        
        return poses, scores 