#!/usr/bin/env python3
"""
Centralized Path Management System

Provides platform-agnostic path handling and eliminates hardcoded paths.
Supports environment variables, configuration files, and automatic detection.
"""

import os
import sys
import platform
import json
from pathlib import Path
from typing import Optional, Dict, Any, List, Union
import shutil

class PathManager:
    """Centralized path management for the docking wrapper."""
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize the path manager.
        
        Args:
            config_file: Optional path to configuration file
        """
        self.project_root = self._get_project_root()
        self.config_file = config_file or self._get_default_config_path()
        self.paths = self._load_paths()
        self._validate_paths()
    
    def _get_project_root(self) -> Path:
        """Get the project root directory."""
        # Try to find the project root by looking for key files
        current = Path(__file__).resolve()
        while current.parent != current:
            # Check if this is the project root by looking for key files/directories
            if (current / "scripts" / "batch_pipeline.py").exists():
                return current
            if (current / "requirements.txt").exists():
                return current
            if (current / "config.py").exists():
                return current
            if (current / "DiffDock").exists() and (current / "scripts").exists():
                return current
            current = current.parent
        
        # Fallback to current working directory
        cwd = Path.cwd()
        
        # Check if current working directory looks like the project root
        if (cwd / "scripts" / "batch_pipeline.py").exists():
            return cwd
        if (cwd / "utils" / "path_manager.py").exists():
            return cwd
        if (cwd / "config.py").exists():
            return cwd
        
        # If we still can't find it, use current working directory but warn
        print(f"Warning: Could not determine project root. Using current directory: {cwd}")
        return cwd
    
    def _get_default_config_path(self) -> Path:
        """Get the default configuration file path."""
        return self.project_root / "config" / "paths.json"
    
    def _load_paths(self) -> Dict[str, Any]:
        """Load paths from configuration file or create defaults."""
        default_paths = self._get_default_paths()
        
        if self.config_file.exists():
            try:
                with open(self.config_file, 'r') as f:
                    config_paths = json.load(f)
                
                # Check if config contains empty strings and clean them up
                if self._has_empty_paths(config_paths):
                    print("Warning: Found empty paths in config file. Regenerating with current paths.")
                    self._save_paths(default_paths)
                    return default_paths
                
                # Merge with defaults
                return self._merge_paths(default_paths, config_paths)
            except Exception as e:
                print(f"Warning: Could not load path config from {self.config_file}: {e}")
                return default_paths
        else:
            # Create default config file
            self._save_paths(default_paths)
            return default_paths
    
    def _get_default_paths(self) -> Dict[str, Any]:
        """Get default paths based on platform and environment."""
        system = platform.system().lower()
        
        # Base directories
        paths = {
            "project_root": str(self.project_root),
            "scripts": str(self.project_root / "scripts"),
            "utils": str(self.project_root / "utils"),
            "config": str(self.project_root / "config"),
            "tests": str(self.project_root / "tests"),
            "docs": str(self.project_root / "docs"),
            "docker": str(self.project_root / "docker"),
            "bin": str(self.project_root / "bin"),
        }
        
        # Output directories (configurable via environment)
        paths.update({
            "outputs": os.environ.get('DOCKING_OUTPUT_DIR', str(self.project_root / "outputs")),
            "logs": os.environ.get('DOCKING_LOG_DIR', str(self.project_root / "logs")),
            "temp": os.environ.get('DOCKING_TEMP_DIR', str(self.project_root / "temp")),
            "cache": os.environ.get('DOCKING_CACHE_DIR', str(self.project_root / "cache")),
        })
        
        # External tools (with environment variable support)
        paths.update({
            "vina": self._find_tool_path("vina", "VINA_PATH", system),
            "gnina": self._find_tool_path("gnina", "GNINA_PATH", system),
            "diffdock": self._find_diffdock_path(system),
            "mgltools": self._find_mgltools_path(system),
        })
        
        # Model directories
        paths.update({
            "models": {
                "equibind": str(self.project_root / "EquiBind"),
                "neuralplexer": str(self.project_root / "NeuralPLexer"),
                "umol": str(self.project_root / "Umol"),
                "boltz2": str(self.project_root / "Boltz2"),
            }
        })
        
        # Data directories
        paths.update({
            "data": {
                "inputs": str(self.project_root / "inputs"),
                "test_data": str(self.project_root / "test_files"),
                "examples": str(self.project_root / "examples"),
            }
        })
        
        return paths
    
    def _find_tool_path(self, tool_name: str, env_var: str, system: str) -> Optional[str]:
        """Find tool path using environment variables and platform detection."""
        # Check environment variable first
        env_path = os.environ.get(env_var)
        if env_path and os.path.exists(env_path):
            return env_path
        
        # Check if tool is in PATH
        if shutil.which(tool_name):
            return tool_name
        
        # Platform-specific common locations
        if system == 'windows':
            possible_paths = [
                f"{tool_name}.exe",
                f"./{tool_name}.exe",
                f"./{tool_name}.bat",
                f"./{tool_name}",
                f"bin/{tool_name}.bat",
                f"bin/{tool_name}",
            ]
        elif system == 'darwin':  # macOS
            possible_paths = [
                f"/usr/local/bin/{tool_name}",
                f"/opt/homebrew/bin/{tool_name}",
                f"~/{tool_name}",
                f"./{tool_name}",
            ]
        else:  # Linux
            possible_paths = [
                f"/usr/local/bin/{tool_name}",
                f"/usr/bin/{tool_name}",
                f"/opt/conda/envs/docking/bin/{tool_name}",
                f"~/{tool_name}",
                f"./{tool_name}",
            ]
        
        for path in possible_paths:
            expanded_path = os.path.expanduser(path)
            if shutil.which(expanded_path):
                return expanded_path
        
        return None
    
    def _find_diffdock_path(self, system: str) -> Optional[str]:
        """Find DiffDock installation path."""
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

    def _find_mgltools_path(self, system: str) -> Optional[str]:
        """Find MGLTools installation path."""
        # Check environment variable first
        mgltools_path = os.environ.get('MGLTOOLS_PATH')
        if mgltools_path and os.path.exists(mgltools_path):
            return mgltools_path
        
        # Platform-specific common locations
        if system == 'windows':
            possible_paths = [
                os.path.join(os.environ.get('PROGRAMFILES', 'C:\\Program Files'), 'MGLTools'),
                os.path.join(os.environ.get('PROGRAMFILES(X86)', 'C:\\Program Files (x86)'), 'MGLTools'),
                os.path.expanduser('~/MGLTools'),
                os.path.expanduser('~/mgltools'),
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
    
    def _merge_paths(self, default: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
        """Merge configuration paths with defaults."""
        merged = default.copy()
        
        for key, value in config.items():
            if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
                merged[key] = self._merge_paths(merged[key], value)
            else:
                # Only use config value if it's not empty and not None
                if value and value != "":
                    merged[key] = value
                # For empty strings or None values, keep the default
        
        return merged
    
    def _has_empty_paths(self, config: Dict[str, Any]) -> bool:
        """Check if config contains empty string paths."""
        for key, value in config.items():
            if isinstance(value, dict):
                if self._has_empty_paths(value):
                    return True
            elif isinstance(value, str) and value == "":
                return True
        return False
    
    def _save_paths(self, paths: Dict[str, Any]):
        """Save paths to configuration file."""
        try:
            self.config_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.config_file, 'w') as f:
                json.dump(paths, f, indent=2)
        except Exception as e:
            print(f"Warning: Could not save path config to {self.config_file}: {e}")
    
    def _validate_paths(self):
        """Validate that critical paths exist."""
        critical_paths = [
            ("project_root", "Project root directory"),
            ("scripts", "Scripts directory"),
            ("utils", "Utils directory"),
        ]
        
        missing_paths = []
        for path_key, description in critical_paths:
            path = self.get_path(path_key)
            if not path or not os.path.exists(path):
                missing_paths.append(f"{description} not found: {path}")
        
        if missing_paths:
            print("Warning: Some critical directories are missing:")
            for missing in missing_paths:
                print(f"  - {missing}")
            print("The script will attempt to continue, but some features may not work properly.")
            print("Make sure you're running the script from the project root directory.")
    
    def get_path(self, path_key: str, sub_key: Optional[str] = None) -> Optional[str]:
        """
        Get a path by key.
        
        Args:
            path_key: Main path key
            sub_key: Optional sub-key for nested paths
            
        Returns:
            Path string or None if not found
        """
        if path_key not in self.paths:
            return None
        
        path_value = self.paths[path_key]
        
        if sub_key:
            if isinstance(path_value, dict) and sub_key in path_value:
                return path_value[sub_key]
            return None
        
        return path_value
    
    def get_absolute_path(self, path_key: str, sub_key: Optional[str] = None) -> Optional[Path]:
        """Get an absolute Path object."""
        path_str = self.get_path(path_key, sub_key)
        if path_str:
            return Path(path_str).resolve()
        return None
    
    def ensure_dir(self, path_key: str, sub_key: Optional[str] = None) -> Optional[Path]:
        """Ensure a directory exists and return its Path."""
        path = self.get_absolute_path(path_key, sub_key)
        if path:
            path.mkdir(parents=True, exist_ok=True)
            return path
        return None
    
    def join(self, *path_parts: Union[str, Path]) -> Path:
        """Join path parts using the project root as base."""
        return self.project_root.joinpath(*path_parts)
    
    def get_output_path(self, *subdirs: str) -> Path:
        """Get an output path with optional subdirectories."""
        output_dir = Path(self.get_path("outputs"))
        return output_dir.joinpath(*subdirs)
    
    def get_log_path(self, filename: str) -> Path:
        """Get a log file path."""
        log_dir = Path(self.get_path("logs"))
        return log_dir / filename
    
    def get_temp_path(self, filename: str) -> Path:
        """Get a temporary file path."""
        temp_dir = Path(self.get_path("temp"))
        return temp_dir / filename
    
    def get_cache_path(self, filename: str) -> Path:
        """Get a cache file path."""
        cache_dir = Path(self.get_path("cache"))
        return cache_dir / filename
    
    def update_path(self, path_key: str, new_path: str, sub_key: Optional[str] = None):
        """Update a path and save to configuration."""
        if sub_key:
            if path_key not in self.paths:
                self.paths[path_key] = {}
            if not isinstance(self.paths[path_key], dict):
                self.paths[path_key] = {}
            self.paths[path_key][sub_key] = new_path
        else:
            self.paths[path_key] = new_path
        
        self._save_paths(self.paths)
    
    def list_available_tools(self) -> Dict[str, bool]:
        """List available external tools."""
        tools = {}
        for tool in ["vina", "gnina", "diffdock", "mgltools"]:
            path = self.get_path(tool)
            tools[tool] = path is not None and os.path.exists(path)
        return tools
    
    def get_platform_info(self) -> Dict[str, str]:
        """Get platform information."""
        return {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python_version": sys.version,
        }


# Global path manager instance
_path_manager = None

def get_path_manager(config_file: Optional[str] = None) -> PathManager:
    """Get the global path manager instance."""
    global _path_manager
    if _path_manager is None:
        _path_manager = PathManager(config_file)
    return _path_manager

def _reset_path_manager():
    """Reset the global path manager (for testing)."""
    global _path_manager
    _path_manager = None

def get_path(path_key: str, sub_key: Optional[str] = None) -> Optional[str]:
    """Get a path using the global path manager."""
    return get_path_manager().get_path(path_key, sub_key)

def get_absolute_path(path_key: str, sub_key: Optional[str] = None) -> Optional[Path]:
    """Get an absolute path using the global path manager."""
    return get_path_manager().get_absolute_path(path_key, sub_key)

def ensure_dir(path_key: str, sub_key: Optional[str] = None) -> Optional[Path]:
    """Ensure a directory exists using the global path manager."""
    return get_path_manager().ensure_dir(path_key, sub_key) 