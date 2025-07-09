#!/usr/bin/env python3
"""
Centralized Configuration System

Manages paths to external tools and provides platform-specific detection.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import platform
import logging
from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context
from pathlib import Path
from typing import Optional, Dict, Any

class Config:
    """Centralized configuration for external tools and paths."""
    
    def __init__(self):
        """Initialize configuration with automatic tool detection."""
        self.logger = logging.getLogger(__name__)
        
        # Initialize tool paths
        self.mgltools_path = self.get_mgltools_path()
        self.vina_path = self._find_vina()
        self.gnina_path = self._find_gnina()
        self.diffdock_path = self._find_diffdock()
        
        # Output directories
        self.output_dir = os.environ.get('OUTPUT_DIR', 'outputs')
        
        # Log detected tools
        self._log_tool_status()
    
    def get_mgltools_path(self):
        """Get MGLTools installation path."""
        # Check environment variable first
        mgltools_path = os.environ.get('MGLTOOLS_PATH')
        if mgltools_path and os.path.exists(mgltools_path):
            return mgltools_path
        
        # Platform-specific common locations
        system = platform.system().lower()
        if system == 'windows':
            possible_paths = [
                os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'MGLTools'),
                os.path.join(os.environ.get('PROGRAMFILES(X86)', 'C:\\Program Files (x86)'), 'MGLTools'),
                os.path.join(os.path.expanduser('~'), 'MGLTools'),
                os.path.join(os.path.expanduser('~'), 'mgltools'),
            ]
        else:  # Linux/macOS
            possible_paths = [
                '/usr/local/mgltools',
                '/opt/mgltools',
                os.path.expanduser('~/mgltools'),
                os.path.expanduser('~/MGLTools'),
            ]
        
        for path in possible_paths:
            if os.path.exists(path):
                return path
        
        return None
    
    def _find_vina(self) -> Optional[str]:
        """Find Vina executable."""
        # Check environment variable first
        vina_path = os.environ.get('VINA_PATH')
        if vina_path and os.path.exists(vina_path):
            return vina_path
        
        # Check if vina is in PATH (try both variants)
        import shutil
        system = platform.system().lower()
        
        # Try platform-appropriate name first
        if system == 'windows':
            vina_names = ['vina.exe', 'vina']
        else:  # Linux/macOS
            vina_names = ['vina', 'vina.exe']
        
        for vina_name in vina_names:
            if shutil.which(vina_name):
                return vina_name
        
        # Platform-specific common locations
        if system == 'windows':
            possible_paths = [
                'vina.exe',  # In PATH
                './vina.exe',  # Current directory
                './vina.bat',  # Current directory
                './vina',  # Current directory
                'bin/vina.bat',  # Local bin directory
                'bin/vina',  # Local bin directory
            ]
        elif system == 'darwin':  # macOS
            possible_paths = [
                '/usr/local/bin/vina',
                '/opt/homebrew/bin/vina',
                os.path.expanduser('~/vina'),
                './vina',
            ]
        else:  # Linux
            possible_paths = [
                '/usr/local/bin/vina',
                '/usr/bin/vina',
                '/opt/conda/envs/docking/bin/vina',
                os.path.expanduser('~/vina'),
                './vina',
            ]
        
        for path in possible_paths:
            if shutil.which(path):
                return path
        
        return None
    
    def _find_gnina(self) -> Optional[str]:
        """Find GNINA executable."""
        # Check environment variable first
        gnina_path = os.environ.get('GNINA_PATH')
        if gnina_path and os.path.exists(gnina_path):
            return gnina_path
        
        # Check if gnina is in PATH
        import shutil
        if shutil.which('gnina'):
            return 'gnina'
        
        # Platform-specific common locations
        system = platform.system().lower()
        if system == 'windows':
            possible_paths = [
                'gnina.exe',  # In PATH
                './gnina.exe',  # Current directory
                './gnina.bat',  # Current directory
                './gnina',  # Current directory
                'bin/gnina.bat',  # Local bin directory
                'bin/gnina',  # Local bin directory
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
                return path
        
        return None
    
    def _find_diffdock(self) -> Optional[str]:
        """Find DiffDock installation."""
        # Check environment variable first
        diffdock_path = os.environ.get('DIFFDOCK_PATH')
        if diffdock_path and os.path.exists(diffdock_path):
            # Check if it's a directory containing inference.py
            inference_script = os.path.join(diffdock_path, 'inference.py')
            if os.path.isfile(inference_script):
                # Check if it's a dummy script
                if self._is_dummy_diffdock_script(inference_script):
                    return None  # Don't return dummy script
                return inference_script
            # If it's a directory but no inference.py, return the directory
            elif os.path.isdir(diffdock_path):
                return diffdock_path
        
        # Check current working directory for DiffDock
        cwd = os.getcwd()
        local_paths = [
            os.path.join(cwd, 'DiffDock', 'inference.py'),
            os.path.join(cwd, 'diffdock', 'inference.py'),
            os.path.join(cwd, 'DiffDock'),
            os.path.join(cwd, 'diffdock'),
        ]
        
        for path in local_paths:
            if os.path.isfile(path):
                # Check if it's a dummy script
                if self._is_dummy_diffdock_script(path):
                    continue  # Skip dummy scripts
                return path
            elif os.path.isdir(path):
                # Check if directory contains inference.py
                inference_script = os.path.join(path, 'inference.py')
                if os.path.isfile(inference_script):
                    # Check if it's a dummy script
                    if self._is_dummy_diffdock_script(inference_script):
                        continue  # Skip dummy scripts
                    return inference_script
        
        # Platform-specific common locations
        system = platform.system().lower()
        if system == 'windows':
            possible_paths = [
                os.path.join(os.path.expanduser('~'), 'DiffDock', 'inference.py'),
                os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'DiffDock', 'inference.py'),
                os.path.join(os.path.expanduser('~'), 'DiffDock'),
                os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'DiffDock'),
            ]
        else:  # Linux/macOS
            possible_paths = [
                os.path.join(cwd, 'DiffDock', 'inference.py'),
                os.path.join(cwd, 'diffdock', 'inference.py'),
                os.path.join(os.path.expanduser('~'), 'DiffDock', 'inference.py'),
                '/opt/DiffDock/inference.py',
                '/usr/local/DiffDock/inference.py',
                '/usr/share/DiffDock/inference.py',
                os.path.join(cwd, 'DiffDock'),
                os.path.join(cwd, 'diffdock'),
                os.path.join(os.path.expanduser('~'), 'DiffDock'),
                '/opt/DiffDock',
                '/usr/local/DiffDock',
                '/usr/share/DiffDock',
            ]
        
        for path in possible_paths:
            if os.path.isfile(path):
                # Check if it's a dummy script
                if self._is_dummy_diffdock_script(path):
                    continue  # Skip dummy scripts
                return path
            elif os.path.isdir(path):
                # Check if directory contains inference.py
                inference_script = os.path.join(path, 'inference.py')
                if os.path.isfile(inference_script):
                    # Check if it's a dummy script
                    if self._is_dummy_diffdock_script(inference_script):
                        continue  # Skip dummy scripts
                    return inference_script
        
        return None

    def _is_dummy_diffdock_script(self, script_path: str) -> bool:
        """Check if a DiffDock script is a dummy script."""
        try:
            with open(script_path, 'r', encoding='utf-8') as f:
                content = f.read()
                # Check for dummy script indicators
                dummy_indicators = [
                    "DIFFDOCK NOT INSTALLED",
                    "This is a dummy script",
                    "DiffDock Dummy Script",
                    "placeholder script"
                ]
                return any(indicator in content for indicator in dummy_indicators)
        except Exception:
            return False
    
    def _log_tool_status(self):
        """Log the status of detected tools."""
        self.logger.info("Tool detection results:")
        self.logger.info(f"  MGLTools: {'✓' if self.mgltools_path else '✗'} {self.mgltools_path or 'Not found'}")
        self.logger.info(f"  Vina: {'✓' if self.vina_path else '✗'} {self.vina_path}")
        self.logger.info(f"  GNINA: {'✓' if self.gnina_path else '✗'} {self.gnina_path}")
        self.logger.info(f"  DiffDock: {'✓' if self.diffdock_path else '✗'} {self.diffdock_path}")
    
    def get_mgltools_pythonsh(self) -> Optional[str]:
        """Get path to MGLTools pythonsh executable."""
        if self.mgltools_path:
            pythonsh_path = os.path.join(self.mgltools_path, 'bin', 'pythonsh')
            if os.path.exists(pythonsh_path):
                return pythonsh_path
        return None
    
    def get_mgltools_prepare_script(self) -> Optional[str]:
        """Get path to MGLTools prepare_receptor4.py script."""
        if self.mgltools_path:
            script_path = os.path.join(
                self.mgltools_path, 
                'MGLToolsPckgs', 
                'AutoDockTools', 
                'Utilities24', 
                'prepare_receptor4.py'
            )
            if os.path.exists(script_path):
                return script_path
        return None
    
    def validate_mgltools(self) -> bool:
        """Validate that MGLTools is properly installed and accessible."""
        pythonsh = self.get_mgltools_pythonsh()
        prepare_script = self.get_mgltools_prepare_script()
        
        if not pythonsh:
            self.logger.warning("MGLTools pythonsh not found")
            return False
        
        if not prepare_script:
            self.logger.warning("MGLTools prepare_receptor4.py script not found")
            return False
        
        return True
    
    def get_config_dict(self) -> Dict[str, Any]:
        """Get configuration as a dictionary."""
        return {
            'mgltools_path': self.mgltools_path,
            'vina_path': self.vina_path,
            'gnina_path': self.gnina_path,
            'diffdock_path': self.diffdock_path,
            'output_dir': self.output_dir,
            'mgltools_available': self.validate_mgltools()
        }

# Global configuration instance
config = Config()
DockingConfig = Config 