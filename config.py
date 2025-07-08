#!/usr/bin/env python3
"""
Centralized Configuration System

Manages paths to external tools and provides platform-specific detection.
"""

import os
import platform
import logging
from pathlib import Path
from typing import Optional, Dict, Any

class Config:
    """Centralized configuration for external tools and paths."""
    
    def __init__(self):
        """Initialize configuration with automatic tool detection."""
        self.logger = logging.getLogger(__name__)
        
        # Initialize tool paths
        self.mgltools_path = self._find_mgltools()
        self.vina_path = self._find_vina()
        self.gnina_path = self._find_gnina()
        self.diffdock_path = self._find_diffdock()
        
        # Output directories
        self.output_dir = os.environ.get('OUTPUT_DIR', 'outputs')
        
        # Log detected tools
        self._log_tool_status()
    
    def _find_mgltools(self) -> Optional[str]:
        """Find MGLTools installation."""
        # Check environment variable first
        mgltools_path = os.environ.get('MGLTOOLS_PATH')
        if mgltools_path and os.path.exists(mgltools_path):
            return mgltools_path
        
        # Platform-specific common locations
        system = platform.system().lower()
        if system == 'darwin':  # macOS
            common_paths = [
                os.path.expanduser('~/mgltools_1.5.7_MacOS-X'),
                '/opt/mgltools',
                '/usr/local/mgltools',
                '/Applications/mgltools'
            ]
        elif system == 'linux':
            common_paths = [
                '/opt/mgltools',
                '/usr/local/mgltools',
                os.path.expanduser('~/mgltools'),
                '/usr/share/mgltools'
            ]
        elif system == 'windows':
            common_paths = [
                'C:\\Program Files\\MGLTools',
                'C:\\mgltools',
                os.path.expanduser('~\\mgltools')
            ]
        else:
            common_paths = []
        
        # Check each path for the pythonsh executable
        for path in common_paths:
            pythonsh_path = os.path.join(path, 'bin', 'pythonsh')
            if os.path.exists(pythonsh_path):
                return path
        
        return None
    
    def _find_vina(self) -> Optional[str]:
        """Find Vina executable."""
        vina_path = os.environ.get('VINA_PATH')
        if vina_path and os.path.exists(vina_path):
            return vina_path
        return None
    
    def _find_gnina(self) -> Optional[str]:
        """Find GNINA executable."""
        gnina_path = os.environ.get('GNINA_PATH')
        if gnina_path and os.path.exists(gnina_path):
            return gnina_path
        return None
    
    def _find_diffdock(self) -> Optional[str]:
        """Find DiffDock installation."""
        diffdock_path = os.environ.get('DIFFDOCK_PATH')
        if diffdock_path and os.path.exists(diffdock_path):
            return diffdock_path
        return None
    
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