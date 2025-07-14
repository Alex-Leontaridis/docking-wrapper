#!/usr/bin/env python3
"""
Dependency Checker

Validates that all required packages and external tools are available.
Provides clear installation instructions for missing dependencies.
"""

import sys
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import subprocess
import importlib
import platform
import shutil
from typing import Dict, List, Tuple

# Required Python packages
REQUIRED_PACKAGES = {
    'rdkit': 'rdkit-pypi',
    'numpy': 'numpy',
    'pandas': 'pandas',
    'scikit-learn': 'scikit-learn',
    'meeko': 'meeko',
    'biopython': 'biopython',
    'colorama': 'colorama',
    'tabulate': 'tabulate',
    'pathlib': 'pathlib',  # Built-in in Python 3.4+
}

# Optional packages (warn if missing)
OPTIONAL_PACKAGES = {
    'netCDF4': 'netCDF4',
    'openmm': 'openmm',
    'mdtraj': 'mdtraj',
}

# External tools to check
EXTERNAL_TOOLS = {
    'vina': 'AutoDock Vina',
    'gnina': 'GNINA',
    'diffdock': 'DiffDock',
}

def check_python_package(package_name: str, install_name: str = None) -> Tuple[bool, str]:
    """Check if a Python package is available."""
    if install_name is None:
        install_name = package_name
    
    try:
        importlib.import_module(package_name)
        return True, f"✓ {package_name}"
    except ImportError:
        return False, f"✗ {package_name} (install: pip install {install_name})"

def check_external_tool(tool_name: str, display_name: str) -> Tuple[bool, str]:
    """Check if an external tool is available in PATH."""
    try:
        result = subprocess.run([tool_name, '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return True, f"✓ {display_name}"
        else:
            return False, f"✗ {display_name} (not found in PATH)"
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        return False, f"✗ {display_name} (not found in PATH)"

def check_mgltools():
    """Check if MGLTools is installed and accessible."""
    print("Checking MGLTools...")
    
    # Check environment variable first
    mgltools_path = os.environ.get('MGLTOOLS_PATH')
    if mgltools_path and os.path.exists(mgltools_path):
        print(f"✓ MGLTools found via environment variable: {mgltools_path}")
        return True
    
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
            print(f"✓ MGLTools found at: {path}")
            return True
    
    print("✗ MGLTools not found")
    print("  Install from: http://mgltools.scripps.edu/downloads/downloads/tools/downloads")
    print("  Or set MGLTOOLS_PATH environment variable")
    return False

def main():
    """Main dependency checking function."""
    print("=" * 60)
    print("DEPENDENCY CHECKER")
    print("=" * 60)
    
    all_good = True
    missing_packages = []
    missing_tools = []
    
    # Check Python packages
    print("\nPython Packages:")
    print("-" * 40)
    
    for package, install_name in REQUIRED_PACKAGES.items():
        available, message = check_python_package(package, install_name)
        print(f"  {message}")
        if not available:
            all_good = False
            missing_packages.append(install_name)
    
    # Check optional packages
    print("\nOptional Python Packages:")
    print("-" * 40)
    
    for package, install_name in OPTIONAL_PACKAGES.items():
        available, message = check_python_package(package, install_name)
        print(f"  {message}")
        if not available:
            print(f"    (Optional - install with: pip install {install_name})")
    
    # Check external tools
    print("\nExternal Tools:")
    print("-" * 40)
    
    for tool, display_name in EXTERNAL_TOOLS.items():
        available, message = check_external_tool(tool, display_name)
        print(f"  {message}")
        if not available:
            missing_tools.append(display_name)
    
    # Check MGLTools
    print("\nMGLTools:")
    print("-" * 40)
    mgltools_available = check_mgltools()
    
    # Check observability dependencies
    print("\nObservability Dependencies:")
    print("-" * 40)
    try:
        import structlog
    except ImportError:
        print("WARNING: structlog is not installed. Logging will not be structured.")
    try:
        import prometheus_client
    except ImportError:
        print("WARNING: prometheus_client is not installed. Metrics will not be available.")
    try:
        import sentry_sdk
    except ImportError:
        print("WARNING: sentry_sdk is not installed. Sentry monitoring will not be available.")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    if missing_packages:
        print(f"\n❌ Missing required Python packages:")
        for package in missing_packages:
            print(f"   pip install {package}")
        all_good = False
    
    if missing_tools:
        print(f"\n⚠️  Missing external tools:")
        for tool in missing_tools:
            print(f"   - {tool}")
        print("   (These are required for docking but may be optional depending on your use case)")
    
    if mgltools_available:
        print(f"\n✅ MGLTools is available")
    else:
        print(f"\n⚠️  MGLTools not found")
        print("   This is optional but recommended for protein preparation.")
        print("   Download from: http://mgltools.scripps.edu/downloads/downloads/tools/downloads")
        print("   Or set MGLTOOLS_PATH environment variable.")
    
    if all_good:
        print(f"\n✅ All required dependencies are available!")
        print("   You can run the docking pipeline.")
    else:
        print(f"\n❌ Some dependencies are missing.")
        print("   Please install the missing packages and tools before running the pipeline.")
    
    print("\n" + "=" * 60)
    
    return all_good

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 