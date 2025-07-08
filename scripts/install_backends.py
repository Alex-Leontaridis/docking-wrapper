#!/usr/bin/env python3
"""
Backend Installation Script

This script helps install and configure all the docking backends and dependencies
for the docking-wrapper project.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import subprocess
import json
import platform
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional

class BackendInstaller:
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.config_file = self.project_root / "backend_config.json"
        self.system = platform.system().lower()
        self.is_windows = self.system == "windows"
        
    def load_config(self) -> Dict:
        """Load the current backend configuration."""
        if self.config_file.exists():
            with open(self.config_file, 'r') as f:
                return json.load(f)
        return {
            "installation_date": "",
            "backends": {
                "gnina": {"installed": False, "path": None, "version": None, "mock": True},
                "diffdock": {"installed": False, "path": None, "version": None, "mock": True},
                "vina": {"installed": False, "path": None, "version": None, "mock": True},
                "equibind": {"installed": False, "path": None, "version": None, "mock": True},
                "neuralplexer": {"installed": False, "path": None, "version": None, "mock": True},
                "umol": {"installed": False, "path": None, "version": None, "mock": True}
            }
        }
    
    def save_config(self, config: Dict):
        """Save the backend configuration."""
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)
    
    def run_command(self, command: List[str], cwd: Optional[Path] = None) -> Tuple[bool, str]:
        """Run a command and return success status and output."""
        try:
            result = subprocess.run(
                command,
                cwd=cwd or self.project_root,
                capture_output=True,
                text=True,
                timeout=300
            )
            return result.returncode == 0, result.stdout + result.stderr
        except subprocess.TimeoutExpired:
            return False, "Command timed out"
        except Exception as e:
            return False, str(e)
    
    def check_python_package(self, package_name: str) -> bool:
        """Check if a Python package is installed."""
        try:
            __import__(package_name)
            return True
        except ImportError:
            return False
    
    def install_python_package(self, package_name: str) -> bool:
        """Install a Python package."""
        print(f"Installing {package_name}...")
        success, output = self.run_command([sys.executable, "-m", "pip", "install", package_name])
        if success:
            print(f"✅ {package_name} installed successfully")
        else:
            print(f"❌ Failed to install {package_name}: {output}")
        return success
    
    def install_gnina(self) -> bool:
        """Install GNINA docking tool."""
        print("\n=== Installing GNINA ===")
        
        # Check if conda is available
        conda_available = shutil.which("conda") is not None
        if not conda_available:
            print("❌ Conda not found. GNINA requires conda for installation.")
            print("Please install conda first: https://docs.conda.io/en/latest/miniconda.html")
            return False
        
        # Create conda environment for GNINA
        env_name = "gnina_env"
        print(f"Creating conda environment: {env_name}")
        
        success, output = self.run_command([
            "conda", "create", "-n", env_name, 
            "python=3.8", "gnina", "-c", "conda-forge", "-y"
        ])
        
        if success:
            print("✅ GNINA installed successfully")
            # Update config
            config = self.load_config()
            config["backends"]["gnina"] = {
                "installed": True,
                "path": f"conda run -n {env_name} gnina",
                "version": "latest",
                "mock": False
            }
            self.save_config(config)
            return True
        else:
            print(f"❌ Failed to install GNINA: {output}")
            return False
    
    def install_diffdock(self) -> bool:
        """Install DiffDock."""
        print("\n=== Installing DiffDock ===")
        
        diffdock_dir = self.project_root / "DiffDock"
        if diffdock_dir.exists():
            print("✅ DiffDock directory already exists")
            return True
        
        # Clone DiffDock repository
        print("Cloning DiffDock repository...")
        success, output = self.run_command([
            "git", "clone", "https://github.com/gcorso/DiffDock.git"
        ])
        
        if not success:
            print(f"❌ Failed to clone DiffDock: {output}")
            return False
        
        # Install DiffDock dependencies
        print("Installing DiffDock dependencies...")
        success, output = self.run_command([
            sys.executable, "-m", "pip", "install", "-r", "requirements.txt"
        ], cwd=diffdock_dir)
        
        if success:
            print("✅ DiffDock installed successfully")
            # Update config
            config = self.load_config()
            config["backends"]["diffdock"] = {
                "installed": True,
                "path": str(diffdock_dir / "inference.py"),
                "version": "latest",
                "mock": False
            }
            self.save_config(config)
            return True
        else:
            print(f"❌ Failed to install DiffDock dependencies: {output}")
            return False
    
    def install_vina(self) -> bool:
        """Install AutoDock Vina."""
        print("\n=== Installing AutoDock Vina ===")
        
        if self.is_windows:
            print("For Windows, please download AutoDock Vina manually:")
            print("https://vina.scripps.edu/downloads/")
            print("Extract to a directory and add to PATH")
            return False
        else:
            # Try to install via conda
            success, output = self.run_command([
                "conda", "install", "-c", "conda-forge", "autodock-vina", "-y"
            ])
            
            if success:
                print("✅ AutoDock Vina installed successfully")
                # Update config
                config = self.load_config()
                config["backends"]["vina"] = {
                    "installed": True,
                    "path": "vina",
                    "version": "latest",
                    "mock": False
                }
                self.save_config(config)
                return True
            else:
                print(f"❌ Failed to install AutoDock Vina: {output}")
                return False
    
    def install_equibind(self) -> bool:
        """Install EquiBind."""
        print("\n=== Installing EquiBind ===")
        
        equibind_dir = self.project_root / "EquiBind"
        if equibind_dir.exists():
            print("✅ EquiBind directory already exists")
            return True
        
        # Clone EquiBind repository
        print("Cloning EquiBind repository...")
        success, output = self.run_command([
            "git", "clone", "https://github.com/luost26/equibind.git", "EquiBind"
        ])
        
        if not success:
            print(f"❌ Failed to clone EquiBind: {output}")
            return False
        
        # Install EquiBind dependencies
        print("Installing EquiBind dependencies...")
        success, output = self.run_command([
            "conda", "env", "create", "-f", "environment_cpuonly.yml"
        ], cwd=equibind_dir)
        
        if success:
            print("✅ EquiBind installed successfully")
            # Update config
            config = self.load_config()
            config["backends"]["equibind"] = {
                "installed": True,
                "path": str(equibind_dir / "inference.py"),
                "version": "latest",
                "mock": False
            }
            self.save_config(config)
            return True
        else:
            print(f"❌ Failed to install EquiBind: {output}")
            return False
    
    def install_neuralplexer(self) -> bool:
        """Install NeuralPLexer."""
        print("\n=== Installing NeuralPLexer ===")
        print("⚠️  NeuralPLexer installation is complex and may require:")
        print("   - CUDA-compatible GPU")
        print("   - Specific PyTorch version")
        print("   - Multiple dependencies")
        print("\nFor now, we'll create a placeholder installation.")
        
        # Create a placeholder directory
        neuralplexer_dir = self.project_root / "NeuralPLexer"
        neuralplexer_dir.mkdir(exist_ok=True)
        
        # Create a simple inference script
        inference_script = neuralplexer_dir / "inference.py"
        with open(inference_script, 'w') as f:
            f.write('''#!/usr/bin/env python3
"""
NeuralPLexer Inference Script (Placeholder)
This is a placeholder for the actual NeuralPLexer installation.
"""

import sys
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
from pathlib import Path

def main():
    print("NeuralPLexer is not fully installed.")
    print("Please follow the official installation guide:")
    print("https://github.com/DeepMind/neuralplexer")
    sys.exit(1)

if __name__ == "__main__":
    main()
''')
        
        print("✅ NeuralPLexer placeholder created")
        # Update config
        config = self.load_config()
        config["backends"]["neuralplexer"] = {
            "installed": True,
            "path": str(inference_script),
            "version": "placeholder",
            "mock": True
        }
        self.save_config(config)
        return True
    
    def install_umol(self) -> bool:
        """Install UMol."""
        print("\n=== Installing UMol ===")
        
        umol_dir = self.project_root / "Umol"
        if umol_dir.exists():
            print("✅ UMol directory already exists")
            return True
        
        # Clone UMol repository
        print("Cloning UMol repository...")
        success, output = self.run_command([
            "git", "clone", "https://github.com/dptech-corp/Uni-Mol.git", "Umol"
        ])
        
        if not success:
            print(f"❌ Failed to clone UMol: {output}")
            return False
        
        # Install UMol dependencies
        print("Installing UMol dependencies...")
        success, output = self.run_command([
            sys.executable, "-m", "pip", "install", "-r", "requirements.txt"
        ], cwd=umol_dir)
        
        if success:
            print("✅ UMol installed successfully")
            # Update config
            config = self.load_config()
            config["backends"]["umol"] = {
                "installed": True,
                "path": str(umol_dir / "predict.sh"),
                "version": "latest",
                "mock": False
            }
            self.save_config(config)
            return True
        else:
            print(f"❌ Failed to install UMol dependencies: {output}")
            return False
    
    def install_python_dependencies(self) -> bool:
        """Install required Python packages."""
        print("\n=== Installing Python Dependencies ===")
        
        required_packages = [
            "rdkit-pypi",
            "numpy",
            "pandas", 
            "scikit-learn",
            "meeko",
            "biopython",
            "colorama",
            "tabulate"
        ]
        
        all_success = True
        for package in required_packages:
            if not self.check_python_package(package.replace("-", "_")):
                if not self.install_python_package(package):
                    all_success = False
        
        return all_success
    
    def setup_environment_variables(self):
        """Set up environment variables."""
        print("\n=== Setting up Environment Variables ===")
        
        # Create environment setup script
        if self.is_windows:
            env_script = self.project_root / "setup_env.bat"
            with open(env_script, 'w') as f:
                f.write('''@echo off
REM Environment setup for docking-wrapper

REM Set project root
set DOCKING_WRAPPER_ROOT=%~dp0

REM Add project scripts to PATH
set PATH=%DOCKING_WRAPPER_ROOT%scripts;%PATH%

REM Set Python path
set PYTHONPATH=%DOCKING_WRAPPER_ROOT%;%PYTHONPATH%

echo Environment variables set for docking-wrapper
echo Run this script before using the docking tools
''')
        else:
            env_script = self.project_root / "setup_env.sh"
            with open(env_script, 'w') as f:
                f.write('''#!/bin/bash
# Environment setup for docking-wrapper

# Set project root
export DOCKING_WRAPPER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Add project scripts to PATH
export PATH="$DOCKING_WRAPPER_ROOT/scripts:$PATH"

# Set Python path
export PYTHONPATH="$DOCKING_WRAPPER_ROOT:$PYTHONPATH"

echo "Environment variables set for docking-wrapper"
echo "Run 'source setup_env.sh' before using the docking tools"
''')
            # Make executable
            os.chmod(env_script, 0o755)
        
        print(f"✅ Environment setup script created: {env_script}")
    
    def create_configuration_files(self):
        """Create default configuration files."""
        print("\n=== Creating Configuration Files ===")
        
        # Create default batch config
        default_config = {
            "input_dir": "inputs",
            "output_dir": "outputs",
            "protein_file": "protein.pdb",
            "ligand_dir": "ligands",
            "backends": {
                "vina": {"enabled": True, "config": {}},
                "gnina": {"enabled": True, "config": {}},
                "diffdock": {"enabled": True, "config": {}},
                "equibind": {"enabled": True, "config": {}},
                "neuralplexer": {"enabled": False, "config": {}},
                "umol": {"enabled": True, "config": {}}
            },
            "analysis": {
                "compute_confidence": True,
                "extract_interactions": True,
                "generate_reports": True
            }
        }
        
        config_file = self.project_root / "default_batch_config.json"
        with open(config_file, 'w') as f:
            json.dump(default_config, f, indent=2)
        
        print(f"✅ Default configuration created: {config_file}")
    
    def run_full_installation(self):
        """Run the complete installation process."""
        print("=" * 60)
        print("DOCKING-WRAPPER BACKEND INSTALLATION")
        print("=" * 60)
        
        # Install Python dependencies first
        if not self.install_python_dependencies():
            print("❌ Failed to install some Python dependencies")
            return False
        
        # Install backends
        backends = [
            ("GNINA", self.install_gnina),
            ("DiffDock", self.install_diffdock),
            ("AutoDock Vina", self.install_vina),
            ("EquiBind", self.install_equibind),
            ("NeuralPLexer", self.install_neuralplexer),
            ("UMol", self.install_umol)
        ]
        
        installed_count = 0
        for name, installer in backends:
            try:
                if installer():
                    installed_count += 1
            except Exception as e:
                print(f"❌ Error installing {name}: {e}")
        
        # Setup environment
        self.setup_environment_variables()
        self.create_configuration_files()
        
        # Final summary
        print("\n" + "=" * 60)
        print("INSTALLATION SUMMARY")
        print("=" * 60)
        print(f"✅ Installed {installed_count}/{len(backends)} backends")
        print("\nNext steps:")
        print("1. Set up environment variables:")
        if self.is_windows:
            print("   run: setup_env.bat")
        else:
            print("   run: source setup_env.sh")
        print("2. Test installation:")
        print("   python scripts/check_dependencies.py")
        print("3. Run a test:")
        print("   python scripts/batch_pipeline.py --config default_batch_config.json")
        
        return installed_count > 0

def main():
    installer = BackendInstaller()
    success = installer.run_full_installation()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main() 