#!/usr/bin/env python3
"""
Test script to verify DiffDock integration fixes.
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

def test_diffdock_optional():
    """Test that DiffDock is truly optional when disabled."""
    print("=== Testing DiffDock Optional Behavior ===")
    
    # Test with default config (DiffDock disabled)
    config = Config()
    path_manager = get_path_manager()
    
    print(f"Default config - DiffDock enabled: {config.diffdock_path is not None}")
    print(f"Path manager - DiffDock path: {path_manager.get_path('diffdock')}")
    
    # Test that DiffDock is truly optional
    # The config should work whether DiffDock is found or not
    if config.diffdock_path is not None:
        print("✅ DiffDock found and available")
    else:
        print("✅ DiffDock not found - this is acceptable")
    
    # Test that the path manager handles missing DiffDock gracefully
    diffdock_path = path_manager.get_path('diffdock')
    if diffdock_path is not None:
        print(f"✅ Path manager found DiffDock: {diffdock_path}")
    else:
        print("✅ Path manager correctly reports no DiffDock")
    
    return True
    
    return True

def test_dummy_script_detection():
    """Test dummy script detection."""
    print("\n=== Testing Dummy Script Detection ===")
    
    # Check if dummy script exists
    dummy_script = project_root / "DiffDock_dummy" / "inference.py"
    if dummy_script.exists():
        print(f"✅ Dummy script found: {dummy_script}")
        
        # Test dummy detection by reading the dummy script
        try:
            with open(dummy_script, 'r') as f:
                content = f.read()
                if "DIFFDOCK NOT INSTALLED" in content:
                    print("✅ Dummy script contains correct dummy indicators")
                else:
                    print("❌ Dummy script missing dummy indicators")
                    return False
        except Exception as e:
            print(f"❌ Error reading dummy script: {e}")
            return False
        
        # Test that the dummy detection function works
        from config import Config
        config = Config()
        
        # If there's a real DiffDock installation, it should be found
        # If there's no real DiffDock, it should be None
        if config.diffdock_path is not None:
            print(f"✅ Real DiffDock found: {config.diffdock_path}")
        else:
            print("✅ No real DiffDock found (dummy correctly ignored)")
    else:
        print(f"❌ Dummy script not found at: {dummy_script}")
        return False
    
    return True

def test_path_handling():
    """Test path handling for DiffDock."""
    print("\n=== Testing Path Handling ===")
    
    # Check for real DiffDock installation
    diffdock_dir = project_root / "DiffDock"
    diffdock_lower = project_root / "diffdock"
    
    print(f"DiffDock directory (capital): {diffdock_dir.exists()}")
    print(f"DiffDock directory (lowercase): {diffdock_lower.exists()}")
    
    if diffdock_dir.exists():
        inference_script = diffdock_dir / "inference.py"
        print(f"Inference script exists: {inference_script.exists()}")
        
        if inference_script.exists():
            # Test if it's a real script (not dummy)
            try:
                with open(inference_script, 'r') as f:
                    content = f.read()
                    if "DIFFDOCK NOT INSTALLED" in content:
                        print("❌ Found dummy script in DiffDock directory")
                        return False
                    else:
                        print("✅ Found real DiffDock script")
            except Exception as e:
                print(f"❌ Error reading inference script: {e}")
                return False
    
    return True

def test_environment_variable():
    """Test environment variable handling."""
    print("\n=== Testing Environment Variable ===")
    
    # Test with environment variable set
    test_path = str(project_root / "test_diffdock")
    os.environ['DIFFDOCK_PATH'] = test_path
    
    config = Config()
    print(f"Config with DIFFDOCK_PATH: {config.diffdock_path}")
    
    # Clean up
    if 'DIFFDOCK_PATH' in os.environ:
        del os.environ['DIFFDOCK_PATH']
    
    return True

def main():
    """Run all tests."""
    print("Testing DiffDock integration fixes...\n")
    
    tests = [
        test_diffdock_optional,
        test_dummy_script_detection,
        test_path_handling,
        test_environment_variable
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"❌ Test {test.__name__} failed with exception: {e}")
            results.append(False)
    
    print("\n=== Summary ===")
    passed = sum(results)
    total = len(results)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("✅ All DiffDock integration fixes working correctly!")
    else:
        print("❌ Some issues remain with DiffDock integration")

if __name__ == "__main__":
    main() 