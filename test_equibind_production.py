#!/usr/bin/env python3
"""
Comprehensive test suite for production-ready EquiBind script
Tests all critical production features including validation, security, error handling, and logging.
"""

import os
import sys
import tempfile
import shutil
import json
import subprocess
import time
from pathlib import Path
# import pytest  # Optional for advanced testing features

def test_config_file_exists():
    """Test that tools_config.json exists and is valid JSON"""
    assert os.path.exists("tools_config.json"), "tools_config.json not found"
    
    with open("tools_config.json", 'r') as f:
        config = json.load(f)
    
    # Check required sections
    assert "tools" in config, "Missing 'tools' section in config"
    assert "equibind" in config["tools"], "Missing 'equibind' section in config"
    assert "logging" in config, "Missing 'logging' section in config"
    assert "security" in config, "Missing 'security' section in config"
    
    print("✅ Configuration file validation passed")

def test_input_validation():
    """Test input file validation features"""
    # Create test files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Valid PDB file
        valid_pdb = os.path.join(temp_dir, "test_protein.pdb")
        with open(valid_pdb, 'w') as f:
            f.write("""TITLE     TEST PROTEIN
ATOM      1  N   ALA A   1      27.462  24.862   5.114  1.00 20.00           N
ATOM      2  CA  ALA A   1      26.336  24.862   6.114  1.00 20.00           C
HETATM    3  C1  LIG A   2      25.336  25.862   7.114  1.00 20.00           C
END""")
        
        # Valid SDF file
        valid_sdf = os.path.join(temp_dir, "test_ligand.sdf")
        with open(valid_sdf, 'w') as f:
            f.write("""aspirin
     RDKit          3D

  7  7  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.7321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  2  7  1  0  0  0  0
M  END
$$$$""")
        
        # Invalid PDB file
        invalid_pdb = os.path.join(temp_dir, "invalid.pdb")
        with open(invalid_pdb, 'w') as f:
            f.write("This is not a valid PDB file")
        
        # Test the script with valid files
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # This should fail gracefully due to missing EquiBind installation
        # but should pass validation
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", valid_pdb,
            "--ligand", valid_sdf,
            "--output", output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should fail due to missing EquiBind, but not due to validation
        assert "Invalid PDB file" not in result.stderr, "Valid PDB file was incorrectly rejected"
        assert "Invalid ligand file" not in result.stderr, "Valid SDF file was incorrectly rejected"
        
        print("✅ Input validation tests passed")

def test_security_features():
    """Test security features like path validation and file size limits"""
    import tempfile
    import platform
    # Create a large file to test size limits
    with tempfile.TemporaryDirectory() as temp_dir:
        large_file = os.path.join(temp_dir, "large.pdb")
        with open(large_file, 'w') as f:
            f.write("A" * 200 * 1024 * 1024)  # 200MB file
        
        # Test with file outside allowed paths
        outside_file = os.path.join(tempfile.gettempdir(), "test.pdb")
        with open(outside_file, 'w') as f:
            f.write("TITLE TEST\nATOM 1 N ALA 1 0.0 0.0 0.0\nEND")
        
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Test large file (should fail)
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", large_file,
            "--ligand", "inputs/ligands/aspirin.sdf",
            "--output", output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert "too large" in result.stderr or "too large" in result.stdout, "Large file size validation failed"
        
        # Test file outside allowed paths (should fail)
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", outside_file,
            "--ligand", "inputs/ligands/aspirin.sdf",
            "--output", output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert "not allowed" in result.stderr or "not allowed" in result.stdout, "Path validation failed"
        
        # Cleanup
        os.remove(outside_file)
        
        print("✅ Security feature tests passed")

def test_logging_setup():
    """Test that logging is properly configured"""
    # Check if logs directory is created
    assert os.path.exists("logs"), "Logs directory not created"
    
    # Test that the script creates log files
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", "inputs/protein.pdb",
            "--ligand", "inputs/ligands/aspirin.sdf",
            "--output", output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check if log file was created
        log_files = [f for f in os.listdir("logs") if f.endswith(".log")]
        assert len(log_files) > 0, "No log files were created"
        
        print("✅ Logging setup tests passed")

def test_error_handling():
    """Test error handling and graceful failure"""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Test with non-existent files
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", "nonexistent.pdb",
            "--ligand", "nonexistent.sdf",
            "--output", output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode != 0, "Script should fail with non-existent files"
        assert "not found" in result.stderr, "Should report file not found"
        
        print("✅ Error handling tests passed")

def test_skip_validation():
    """Test the --skip-validation flag"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create invalid files
        invalid_pdb = os.path.join(temp_dir, "invalid.pdb")
        with open(invalid_pdb, 'w') as f:
            f.write("This is not a valid PDB file")
        
        invalid_sdf = os.path.join(temp_dir, "invalid.sdf")
        with open(invalid_sdf, 'w') as f:
            f.write("This is not a valid SDF file")
        
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)
        
        # Test with skip-validation flag
        cmd = [
            sys.executable, "scripts/run_equibind.py",
            "--protein", invalid_pdb,
            "--ligand", invalid_sdf,
            "--output", output_dir,
            "--skip-validation"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        # Should fail due to missing EquiBind, but not due to validation
        assert "Invalid PDB file" not in result.stderr, "Validation should be skipped"
        assert "Invalid ligand file" not in result.stderr, "Validation should be skipped"
        
        print("✅ Skip validation tests passed")

def test_help_output():
    """Test that help output is comprehensive"""
    cmd = [sys.executable, "scripts/run_equibind.py", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    help_text = result.stdout + result.stderr
    assert result.returncode == 0, "Help command should succeed"
    assert "production-ready" in help_text, "Help should mention production-ready"
    assert "--skip-validation" in help_text, "Help should include skip-validation option"
    print("✅ Help output tests passed")

def main():
    """Run all tests"""
    print("🧪 Testing Production-Ready EquiBind Script")
    print("=" * 50)
    
    tests = [
        test_config_file_exists,
        test_input_validation,
        test_security_features,
        test_logging_setup,
        test_error_handling,
        test_skip_validation,
        test_help_output
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"❌ {test.__name__} failed: {e}")
    
    print(f"\n{'=' * 50}")
    print(f"Test Results: {passed}/{total} passed")
    print(f"Success Rate: {(passed/total)*100:.1f}%")
    
    if passed == total:
        print("🎉 All production-ready tests passed!")
        return 0
    else:
        print("⚠️  Some tests failed. Please review the output above.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 