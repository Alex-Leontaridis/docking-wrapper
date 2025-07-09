#!/usr/bin/env python3
"""
Test script to verify the fixes for vina and DiffDock detection.
"""

import os
import platform
import sys
from pathlib import Path

# Add the project root to the path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from config import Config
from utils.path_manager import get_path_manager

def test_vina_detection():
    """Test vina detection on different platforms."""
    print("=== Testing Vina Detection ===")
    
    config = Config()
    system = platform.system().lower()
    
    print(f"Platform: {system}")
    print(f"Vina path: {config.vina_path}")
    
    if config.vina_path:
        print("✅ Vina detected successfully!")
    else:
        print("❌ Vina not detected")
        
        # Check if vina is in PATH
        import shutil
        if system == 'windows':
            vina_names = ['vina.exe', 'vina']
        else:
            vina_names = ['vina', 'vina.exe']
            
        for name in vina_names:
            path = shutil.which(name)
            if path:
                print(f"   Found {name} in PATH: {path}")
                break
        else:
            print("   No vina binary found in PATH")

def test_diffdock_detection():
    """Test DiffDock detection."""
    print("\n=== Testing DiffDock Detection ===")
    
    config = Config()
    path_manager = get_path_manager()
    
    print(f"Config DiffDock path: {config.diffdock_path}")
    print(f"Path manager DiffDock path: {path_manager.get_path('diffdock')}")
    
    # Check if DiffDock directory exists
    diffdock_dir = project_root / "DiffDock"
    if diffdock_dir.exists():
        print(f"✅ DiffDock directory found: {diffdock_dir}")
        
        # Check for inference.py
        inference_script = diffdock_dir / "inference.py"
        if inference_script.exists():
            print(f"✅ inference.py found: {inference_script}")
        else:
            print(f"❌ inference.py not found in {diffdock_dir}")
    else:
        print(f"❌ DiffDock directory not found at {diffdock_dir}")
    
    # Check environment variable
    diffdock_env = os.environ.get('DIFFDOCK_PATH')
    if diffdock_env:
        print(f"DIFFDOCK_PATH environment variable: {diffdock_env}")
    else:
        print("No DIFFDOCK_PATH environment variable set")

def test_binary_finding():
    """Test the find_binary function."""
    print("\n=== Testing Binary Finding ===")
    
    # Import the function from run_docking_multi
    sys.path.insert(0, str(project_root / "scripts"))
    from run_docking_multi import find_binary
    
    system = platform.system().lower()
    
    # Test vina finding
    if system == 'windows':
        vina_binary = 'vina.exe'
    else:
        vina_binary = 'vina'
    
    vina_path = find_binary(vina_binary)
    print(f"find_binary('{vina_binary}'): {vina_path}")
    
    # Test gnina finding
    gnina_path = find_binary('gnina')
    print(f"find_binary('gnina'): {gnina_path}")

def main():
    """Run all tests."""
    print("Testing fixes for vina and DiffDock detection...\n")
    
    test_vina_detection()
    test_diffdock_detection()
    test_binary_finding()
    
    print("\n=== Summary ===")
    print("If you see ✅ marks, the detection is working correctly.")
    print("If you see ❌ marks, there may still be issues to resolve.")

if __name__ == "__main__":
    main() 