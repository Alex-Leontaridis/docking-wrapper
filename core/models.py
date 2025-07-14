#!/usr/bin/env python3
"""
Machine Learning Models Implementation

This module implements concrete ML models that inherit from the MLModel interface.
"""

import os
import sys
import subprocess
import time
import json
import tempfile
from pathlib import Path
from typing import Dict, Any, Optional, List
import logging
import requests
import platform

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .interfaces import MLModel, PredictionResult, ComponentStatus
from config import config
from utils.path_manager import get_path_manager, ensure_dir
from utils.logging import setup_logging


class Boltz2Model(MLModel):
    """Boltz2 API-based model implementation."""
    
    def __init__(self, name: str = "boltz2", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.api_url = config.get('api_url', 'https://api.boltz2.com') if config else 'https://api.boltz2.com'
        self.api_key = config.get('api_key', None) if config else None
        self.model_loaded = False
    
    def check_availability(self) -> ComponentStatus:
        """Check if Boltz2 API is available."""
        if not self.api_key:
            self.logger.warning("Boltz2 API key not configured")
            return ComponentStatus.UNAVAILABLE
        
        # Test API connectivity
        try:
            headers = {'Authorization': f'Bearer {self.api_key}'}
            response = requests.get(f"{self.api_url}/health", headers=headers, timeout=10)
            if response.status_code == 200:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing Boltz2 API: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get Boltz2 capabilities."""
        return {
            "name": "Boltz2",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_predictions": True,
            "max_predictions": 100,
            "supports_batch_processing": True,
            "api_based": True
        }
    
    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load the Boltz2 model (API-based, no local loading needed)."""
        if self.check_availability() == ComponentStatus.AVAILABLE:
            self.model_loaded = True
            self.set_status(ComponentStatus.AVAILABLE)
            return True
        return False
    
    def predict(self, complex_pdb: Path, **kwargs) -> PredictionResult:
        """Run Boltz2 prediction on protein-ligand complex."""
        start_time = time.time()
        
        try:
            self.set_status(ComponentStatus.RUNNING)
            
            # Validate inputs
            if not complex_pdb.exists():
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"Complex PDB file not found: {complex_pdb}"
                )
            
            if not self.model_loaded:
                if not self.load_model():
                    return PredictionResult(
                        success=False,
                        predictions={},
                        error_message="Boltz2 model not loaded"
                    )
            
            # Prepare API request
            headers = {
                'Authorization': f'Bearer {self.api_key}',
                'Content-Type': 'application/json'
            }
            
            # Read PDB file
            with open(complex_pdb, 'r') as f:
                pdb_content = f.read()
            
            # Prepare request data
            request_data = {
                'pdb_content': pdb_content,
                'prediction_type': kwargs.get('prediction_type', 'binding_affinity'),
                'confidence_threshold': kwargs.get('confidence_threshold', 0.5)
            }
            
            # Make API request
            response = requests.post(
                f"{self.api_url}/predict",
                headers=headers,
                json=request_data,
                timeout=300
            )
            
            if response.status_code != 200:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"Boltz2 API error: {response.status_code} - {response.text}",
                    execution_time=time.time() - start_time
                )
            
            # Parse response
            result_data = response.json()
            
            self.set_status(ComponentStatus.COMPLETED)
            
            return PredictionResult(
                success=True,
                predictions=result_data.get('predictions', {}),
                confidence=result_data.get('confidence'),
                metadata={
                    "model": "Boltz2",
                    "api_version": result_data.get('version'),
                    "prediction_type": request_data['prediction_type']
                },
                execution_time=time.time() - start_time
            )
            
        except Exception as e:
            self.set_status(ComponentStatus.FAILED)
            return PredictionResult(
                success=False,
                predictions={},
                error_message=str(e),
                execution_time=time.time() - start_time
            )


class EquiBindModel(MLModel):
    """EquiBind PyTorch model implementation."""
    
    def __init__(self, name: str = "equibind", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.equibind_path = None
        self.model_path = None
        self.conda_env = config.get('conda_env', 'equibind') if config else 'equibind'
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to EquiBind."""
        # Check for EquiBind in current directory
        current_dir = Path.cwd()
        possible_paths = [
            current_dir / 'EquiBind',
            current_dir / 'equibind',
            Path.home() / 'EquiBind',
            Path('/opt/EquiBind')
        ]
        
        for path in possible_paths:
            if path.exists() and (path / 'inference.py').exists():
                self.equibind_path = path
                break
        
        # Look for model checkpoint
        if self.equibind_path:
            model_paths = [
                self.equibind_path / 'runs' / 'flexible_self_docking' / 'best_checkpoint.pt',
                self.equibind_path / 'runs' / 'rigid_redocking' / 'best_checkpoint.pt',
                self.equibind_path / 'models' / 'best_checkpoint.pt'
            ]
            
            for path in model_paths:
                if path.exists():
                    self.model_path = path
                    break
    
    def check_availability(self) -> ComponentStatus:
        """Check if EquiBind is available."""
        if not self.equibind_path:
            self.logger.warning("EquiBind not found")
            return ComponentStatus.UNAVAILABLE
        
        # Check if conda environment exists
        try:
            result = subprocess.run(['conda', 'env', 'list'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0 and self.conda_env in result.stdout:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error checking conda environment: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get EquiBind capabilities."""
        return {
            "name": "EquiBind",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_predictions": True,
            "max_predictions": 10,
            "supports_batch_processing": True,
            "api_based": False,
            "requires_conda": True
        }
    
    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load the EquiBind model."""
        if model_path:
            self.model_path = model_path
        
        if not self.model_path or not self.model_path.exists():
            self.logger.warning("EquiBind model checkpoint not found")
            return False
        
        # Model is loaded when running inference
        self.model_loaded = True
        self.set_status(ComponentStatus.AVAILABLE)
        return True
    
    def predict(self, complex_pdb: Path, **kwargs) -> PredictionResult:
        """Run EquiBind prediction on protein-ligand complex."""
        start_time = time.time()
        
        try:
            self.set_status(ComponentStatus.RUNNING)
            
            # Validate inputs
            if not complex_pdb.exists():
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"Complex PDB file not found: {complex_pdb}"
                )
            
            if not self.equibind_path:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message="EquiBind not found"
                )
            
            # Create output directory
            output_dir = kwargs.get('output_dir', Path.cwd() / 'equibind_output')
            ensure_dir(output_dir)
            
            # Prepare EquiBind command
            inference_script = self.equibind_path / 'inference.py'
            output_file = output_dir / f"{complex_pdb.stem}_equibind_out.pdb"
            
            cmd = [
                'conda', 'run', '-n', self.conda_env,
                'python', str(inference_script),
                '--protein_path', str(complex_pdb),
                '--out_dir', str(output_dir),
                '--model_path', str(self.model_path) if self.model_path else '',
                '--num_predictions', str(kwargs.get('num_predictions', 1))
            ]
            
            # Run EquiBind
            self.logger.info(f"Running EquiBind with command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30 minutes
            
            if result.returncode != 0:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"EquiBind failed: {result.stderr}",
                    execution_time=time.time() - start_time
                )
            
            # Parse results
            predictions = self._parse_equibind_results(output_dir, output_file)
            
            self.set_status(ComponentStatus.COMPLETED)
            
            return PredictionResult(
                success=True,
                predictions=predictions,
                metadata={
                    "model": "EquiBind",
                    "model_path": str(self.model_path) if self.model_path else None,
                    "conda_env": self.conda_env
                },
                execution_time=time.time() - start_time
            )
            
        except Exception as e:
            self.set_status(ComponentStatus.FAILED)
            return PredictionResult(
                success=False,
                predictions={},
                error_message=str(e),
                execution_time=time.time() - start_time
            )
    
    def _parse_equibind_results(self, output_dir: Path, output_file: Path) -> Dict[str, Any]:
        """Parse EquiBind output files."""
        predictions = {}
        
        if output_file.exists():
            predictions['output_pdb'] = str(output_file)
        
        # Look for additional output files
        log_files = list(output_dir.glob('*.log'))
        if log_files:
            predictions['log_files'] = [str(f) for f in log_files]
        
        return predictions


class UMolModel(MLModel):
    """UMol binary model implementation."""
    
    def __init__(self, name: str = "umol", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.umol_path = None
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to UMol."""
        # Check for UMol in current directory
        current_dir = Path.cwd()
        possible_paths = [
            current_dir / 'Umol',
            current_dir / 'umol',
            Path.home() / 'Umol',
            Path('/opt/Umol')
        ]
        
        for path in possible_paths:
            if path.exists() and (path / 'predict.sh').exists():
                self.umol_path = path
                break
    
    def check_availability(self) -> ComponentStatus:
        """Check if UMol is available."""
        if not self.umol_path:
            self.logger.warning("UMol not found")
            return ComponentStatus.UNAVAILABLE
        
        # Check if predict script is executable
        predict_script = self.umol_path / 'predict.sh'
        if not predict_script.exists() or not os.access(predict_script, os.X_OK):
            self.logger.warning("UMol predict script not found or not executable")
            return ComponentStatus.UNAVAILABLE
        
        # Test if UMol is working
        try:
            result = subprocess.run([str(predict_script), '--help'], 
                                  capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing UMol: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get UMol capabilities."""
        return {
            "name": "UMol",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_predictions": True,
            "max_predictions": 5,
            "supports_batch_processing": False,
            "api_based": False,
            "binary_based": True
        }
    
    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load the UMol model (binary-based, no explicit loading needed)."""
        if self.check_availability() == ComponentStatus.AVAILABLE:
            self.model_loaded = True
            return True
        return False
    
    def predict(self, complex_pdb: Path, **kwargs) -> PredictionResult:
        """Run UMol prediction on protein-ligand complex."""
        start_time = time.time()
        
        try:
            self.set_status(ComponentStatus.RUNNING)
            
            # Validate inputs
            if not complex_pdb.exists():
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"Complex PDB file not found: {complex_pdb}"
                )
            
            if not self.umol_path:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message="UMol not found"
                )
            
            # Create output directory
            output_dir = kwargs.get('output_dir', Path.cwd() / 'umol_output')
            ensure_dir(output_dir)
            
            # Prepare UMol command
            predict_script = self.umol_path / 'predict.sh'
            output_file = output_dir / f"{complex_pdb.stem}_umol_out.pdb"
            
            cmd = [
                str(predict_script),
                '--input', str(complex_pdb),
                '--output', str(output_file),
                '--num_predictions', str(kwargs.get('num_predictions', 1))
            ]
            
            # Run UMol
            self.logger.info(f"Running UMol with command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30 minutes
            
            if result.returncode != 0:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"UMol failed: {result.stderr}",
                    execution_time=time.time() - start_time
                )
            
            # Parse results
            predictions = self._parse_umol_results(output_dir, output_file)
            
            self.set_status(ComponentStatus.COMPLETED)
            
            return PredictionResult(
                success=True,
                predictions=predictions,
                metadata={
                    "model": "UMol",
                    "umol_path": str(self.umol_path)
                },
                execution_time=time.time() - start_time
            )
            
        except Exception as e:
            self.set_status(ComponentStatus.FAILED)
            return PredictionResult(
                success=False,
                predictions={},
                error_message=str(e),
                execution_time=time.time() - start_time
            )
    
    def _parse_umol_results(self, output_dir: Path, output_file: Path) -> Dict[str, Any]:
        """Parse UMol output files."""
        predictions = {}
        
        if output_file.exists():
            predictions['output_pdb'] = str(output_file)
        
        # Look for additional output files
        log_files = list(output_dir.glob('*.log'))
        if log_files:
            predictions['log_files'] = [str(f) for f in log_files]
        
        return predictions


class NeuralPLexerModel(MLModel):
    """NeuralPLexer CLI model implementation."""
    
    def __init__(self, name: str = "neuralplexer", config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.neuralplexer_path = None
        self._setup_paths()
    
    def _setup_paths(self):
        """Setup paths to NeuralPLexer."""
        # Check for NeuralPLexer in current directory
        current_dir = Path.cwd()
        possible_paths = [
            current_dir / 'NeuralPLexer',
            current_dir / 'neuralplexer',
            Path.home() / 'NeuralPLexer',
            Path('/opt/NeuralPLexer')
        ]
        
        for path in possible_paths:
            if path.exists() and (path / 'predict.py').exists():
                self.neuralplexer_path = path
                break
        
        # Also check if it's installed via pip
        if not self.neuralplexer_path:
            try:
                result = subprocess.run(['neuralplexer', '--help'], 
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    self.neuralplexer_path = "installed"
            except Exception:
                pass
    
    def check_availability(self) -> ComponentStatus:
        """Check if NeuralPLexer is available."""
        if not self.neuralplexer_path:
            self.logger.warning("NeuralPLexer not found")
            return ComponentStatus.UNAVAILABLE
        
        # Test if NeuralPLexer is working
        try:
            if self.neuralplexer_path == "installed":
                result = subprocess.run(['neuralplexer', '--help'], 
                                      capture_output=True, text=True, timeout=30)
            else:
                predict_script = self.neuralplexer_path / 'predict.py'
                result = subprocess.run([sys.executable, str(predict_script), '--help'], 
                                      capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                self.set_status(ComponentStatus.AVAILABLE)
                return ComponentStatus.AVAILABLE
        except Exception as e:
            self.logger.error(f"Error testing NeuralPLexer: {e}")
        
        return ComponentStatus.UNAVAILABLE
    
    def get_capabilities(self) -> Dict[str, Any]:
        """Get NeuralPLexer capabilities."""
        return {
            "name": "NeuralPLexer",
            "version": "1.0",
            "supports_gpu": True,
            "supports_multiple_predictions": True,
            "max_predictions": 10,
            "supports_batch_processing": True,
            "api_based": False,
            "cli_based": True
        }
    
    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load the NeuralPLexer model (CLI-based, no explicit loading needed)."""
        if self.check_availability() == ComponentStatus.AVAILABLE:
            self.model_loaded = True
            return True
        return False
    
    def predict(self, complex_pdb: Path, **kwargs) -> PredictionResult:
        """Run NeuralPLexer prediction on protein-ligand complex."""
        start_time = time.time()
        
        try:
            self.set_status(ComponentStatus.RUNNING)
            
            # Validate inputs
            if not complex_pdb.exists():
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"Complex PDB file not found: {complex_pdb}"
                )
            
            if not self.neuralplexer_path:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message="NeuralPLexer not found"
                )
            
            # Create output directory
            output_dir = kwargs.get('output_dir', Path.cwd() / 'neuralplexer_output')
            ensure_dir(output_dir)
            
            # Prepare NeuralPLexer command
            output_file = output_dir / f"{complex_pdb.stem}_neuralplexer_out.pdb"
            
            if self.neuralplexer_path == "installed":
                cmd = [
                    'neuralplexer',
                    '--input', str(complex_pdb),
                    '--output', str(output_file),
                    '--num_predictions', str(kwargs.get('num_predictions', 1))
                ]
            else:
                predict_script = self.neuralplexer_path / 'predict.py'
                cmd = [
                    sys.executable, str(predict_script),
                    '--input', str(complex_pdb),
                    '--output', str(output_file),
                    '--num_predictions', str(kwargs.get('num_predictions', 1))
                ]
            
            # Run NeuralPLexer
            self.logger.info(f"Running NeuralPLexer with command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)  # 1 hour
            
            if result.returncode != 0:
                return PredictionResult(
                    success=False,
                    predictions={},
                    error_message=f"NeuralPLexer failed: {result.stderr}",
                    execution_time=time.time() - start_time
                )
            
            # Parse results
            predictions = self._parse_neuralplexer_results(output_dir, output_file)
            
            self.set_status(ComponentStatus.COMPLETED)
            
            return PredictionResult(
                success=True,
                predictions=predictions,
                metadata={
                    "model": "NeuralPLexer",
                    "neuralplexer_path": str(self.neuralplexer_path)
                },
                execution_time=time.time() - start_time
            )
            
        except Exception as e:
            self.set_status(ComponentStatus.FAILED)
            return PredictionResult(
                success=False,
                predictions={},
                error_message=str(e),
                execution_time=time.time() - start_time
            )
    
    def _parse_neuralplexer_results(self, output_dir: Path, output_file: Path) -> Dict[str, Any]:
        """Parse NeuralPLexer output files."""
        predictions = {}
        
        if output_file.exists():
            predictions['output_pdb'] = str(output_file)
        
        # Look for additional output files
        log_files = list(output_dir.glob('*.log'))
        if log_files:
            predictions['log_files'] = [str(f) for f in log_files]
        
        return predictions 