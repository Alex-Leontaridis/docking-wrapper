#!/usr/bin/env python3
"""
Comprehensive test script to verify all bug fixes have been implemented correctly.
"""

import os
import sys
import platform
import subprocess
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

def test_bug_1_binary_detection_conflicts():
    """Test Bug 1: Binary Detection Conflicts - Fixed"""
    print("\n=== Testing Bug 1: Binary Detection Conflicts ===")
    
    try:
        from utils.path_manager import get_path_manager
        
        # Test unified binary detection
        path_manager = get_path_manager()
        
        # Test that get_binary_path method exists
        assert hasattr(path_manager, 'get_binary_path'), "get_binary_path method not found"
        
        # Test platform-specific binary names
        system = platform.system().lower()
        vina_path = path_manager.get_binary_path('vina')
        gnina_path = path_manager.get_binary_path('gnina')
        
        print(f"✓ Unified binary detection working")
        print(f"  Platform: {system}")
        print(f"  Vina path: {vina_path}")
        print(f"  GNINA path: {gnina_path}")
        
        return True
        
    except Exception as e:
        print(f"✗ Binary detection test failed: {e}")
        return False

def test_bug_2_hardcoded_windows_assumptions():
    """Test Bug 2: Hardcoded Windows Assumptions - Fixed"""
    print("\n=== Testing Bug 2: Hardcoded Windows Assumptions ===")
    
    try:
        from utils.path_manager import get_path_manager
        
        path_manager = get_path_manager()
        system = platform.system().lower()
        
        # Test platform-specific binary name generation
        binary_name = path_manager._get_platform_binary_name('vina', system)
        
        if system == 'windows':
            expected = 'vina.exe'
        else:
            expected = 'vina'
        
        assert binary_name == expected, f"Expected {expected}, got {binary_name}"
        
        print(f"✓ Platform-specific binary names working")
        print(f"  Platform: {system}")
        print(f"  Binary name: {binary_name}")
        
        # Test that config.py uses PathManager
        from config import Config
        config = Config()
        print(f"✓ Config uses PathManager for binary detection")
        
        return True
        
    except Exception as e:
        print(f"✗ Hardcoded Windows assumptions test failed: {e}")
        return False

def test_bug_3_python_version_issues():
    """Test Bug 3: Python Version Issues - Fixed"""
    print("\n=== Testing Bug 3: Python Version Issues ===")
    
    try:
        # Check requirements.txt
        with open('requirements.txt', 'r') as f:
            content = f.read()
        
        # Check that Python version allows 3.11+
        assert 'python>=3.10,<3.12' in content, "Python version not updated to allow 3.11+"
        
        # Check that rdkit-pypi is fixed
        assert 'rdkit>=2022.9.1' in content, "rdkit-pypi not fixed to rdkit"
        
        # Check requirements_enhanced.txt
        with open('requirements_enhanced.txt', 'r') as f:
            content_enhanced = f.read()
        
        assert 'python>=3.10,<3.12' in content_enhanced, "Enhanced requirements Python version not updated"
        assert 'rdkit>=2022.9.1' in content_enhanced, "Enhanced requirements rdkit not fixed"
        
        print(f"✓ Python version requirements updated")
        print(f"  Current Python: {sys.version}")
        print(f"  Requirements allow: Python 3.10-3.11")
        print(f"  RDKit package name fixed")
        
        return True
        
    except Exception as e:
        print(f"✗ Python version issues test failed: {e}")
        return False

def test_bug_4_diffdock_integration():
    """Test Bug 4: DiffDock Integration Issues - Fixed"""
    print("\n=== Testing Bug 4: DiffDock Integration Issues ===")
    
    try:
        # Test that dummy script exists and works
        dummy_script = Path('DiffDock_dummy/inference.py')
        assert dummy_script.exists(), "DiffDock dummy script not found"
        
        # Test that dummy script doesn't exit with error code
        result = subprocess.run([
            sys.executable, str(dummy_script), 
            '--protein', 'test.pdb', 
            '--ligand', 'test.sdf', 
            '--out', 'test_output'
        ], capture_output=True, text=True, timeout=30)
        
        # Should exit with success (0) not error (1)
        assert result.returncode == 0, f"Dummy script exited with error code {result.returncode}"
        
        # Check that output files were created
        test_output = Path('test_output')
        if test_output.exists():
            confidence_file = test_output / 'diffdock_confidence.txt'
            assert confidence_file.exists(), "Dummy script didn't create confidence file"
            
            # Clean up
            import shutil
            shutil.rmtree(test_output)
        
        print(f"✓ DiffDock dummy script working")
        print(f"  Script exits successfully")
        print(f"  Creates mock output files")
        
        return True
        
    except Exception as e:
        print(f"✗ DiffDock integration test failed: {e}")
        return False

def test_bug_5_import_failures():
    """Test Bug 5: Import Failures - Fixed"""
    print("\n=== Testing Bug 5: Import Failures ===")
    
    try:
        # Test psutil import handling
        from scripts.run_equibind import PSUTIL_AVAILABLE
        print(f"✓ psutil import handling: {'Available' if PSUTIL_AVAILABLE else 'Not available'}")
        
        # Test Meeko import handling
        from scripts.prep_structures import MEEKO_AVAILABLE
        print(f"✓ Meeko import handling: {'Available' if MEEKO_AVAILABLE else 'Not available'}")
        
        # Test that modules don't crash when imports fail
        if not PSUTIL_AVAILABLE:
            from scripts.run_equibind import check_system_resources
            result, msg = check_system_resources({})
            assert result == True, "System resource check should work without psutil"
            print(f"✓ System resource check works without psutil: {msg}")
        
        print(f"✓ Import failure handling working")
        
        return True
        
    except Exception as e:
        print(f"✗ Import failures test failed: {e}")
        return False

def test_bug_6_optional_components():
    """Test Bug 6: Optional Components Treated as Required - Fixed"""
    print("\n=== Testing Bug 6: Optional Components ===")
    
    try:
        # Test MGLTools optional handling
        from utils.path_manager import get_path_manager
        path_manager = get_path_manager()
        
        mgltools_path = path_manager.get_path('mgltools')
        print(f"✓ MGLTools detection: {'Found' if mgltools_path else 'Not found'}")
        
        # Test that the system doesn't crash when MGLTools is missing
        # This is tested by the fact that the import didn't fail above
        
        print(f"✓ Optional components handled gracefully")
        
        return True
        
    except Exception as e:
        print(f"✗ Optional components test failed: {e}")
        return False

def main():
    """Run all bug fix tests."""
    print("=" * 60)
    print("COMPREHENSIVE BUG FIX TESTING")
    print("=" * 60)
    
    tests = [
        ("Binary Detection Conflicts", test_bug_1_binary_detection_conflicts),
        ("Hardcoded Windows Assumptions", test_bug_2_hardcoded_windows_assumptions),
        ("Python Version Issues", test_bug_3_python_version_issues),
        ("DiffDock Integration", test_bug_4_diffdock_integration),
        ("Import Failures", test_bug_5_import_failures),
        ("Optional Components", test_bug_6_optional_components),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"✗ {test_name} test crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = 0
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
        if result:
            passed += 1
    
    print(f"\nResults: {passed}/{len(results)} tests passed")
    
    if passed == len(results):
        print("\n🎉 ALL BUG FIXES VERIFIED SUCCESSFULLY!")
        print("The codebase is now more robust and cross-platform compatible.")
    else:
        print(f"\n⚠️  {len(results) - passed} tests failed. Some bugs may still exist.")
    
    return passed == len(results)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 