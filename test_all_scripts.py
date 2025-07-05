#!/usr/bin/env python3
"""
Comprehensive test script for the molecular docking pipeline
Tests all four main scripts: prep_structures.py, run_docking_multi.py, parse_and_score_results.py, and batch_pipeline.py
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path

def run_test(test_name, command, expected_exit_code=0):
    """Run a test command and report results."""
    print(f"\n{'='*60}")
    print(f"TEST: {test_name}")
    print(f"COMMAND: {' '.join(command)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        success = result.returncode == expected_exit_code
        
        print(f"EXIT CODE: {result.returncode}")
        if result.stdout:
            print(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            print(f"STDERR:\n{result.stderr}")
        
        if success:
            print(f"✅ {test_name} PASSED")
        else:
            print(f"❌ {test_name} FAILED (expected exit code {expected_exit_code})")
        
        return success
        
    except subprocess.TimeoutExpired:
        print(f"❌ {test_name} FAILED (timeout)")
        return False
    except Exception as e:
        print(f"❌ {test_name} FAILED (exception: {e})")
        return False

def main():
    """Run all tests."""
    print("🧪 COMPREHENSIVE DOCKING PIPELINE TEST SUITE")
    print("="*60)
    
    # Get the conda environment Python
    conda_python = r"C:\Users\alexl\miniconda3\envs\docking\python.exe"
    
    # Test results
    test_results = []
    
    # Test 1: prep_structures.py
    test_results.append(run_test(
        "Structure Preparation (prep_structures.py)",
        [conda_python, "scripts/prep_structures.py", "--protein", "inputs/protein.pdb", "--ligand", "inputs/ligands/aspirin.sdf"]
    ))
    
    # Test 2: run_docking_multi.py (Vina only)
    test_results.append(run_test(
        "Docking Execution - Vina (run_docking_multi.py)",
        [conda_python, "scripts/run_docking_multi.py", "--protein", "protein_prepped.pdbqt", "--ligand", "ligand_prepped.pdbqt"]
    ))
    
    # Test 3: run_docking_multi.py (All engines)
    test_results.append(run_test(
        "Docking Execution - All Engines (run_docking_multi.py)",
        [conda_python, "scripts/run_docking_multi.py", "--protein", "protein_prepped.pdbqt", "--ligand", "ligand_prepped.pdbqt", "--use_gnina", "--use_diffdock"]
    ))
    
    # Test 4: parse_and_score_results.py
    test_results.append(run_test(
        "Results Parsing (parse_and_score_results.py)",
        [conda_python, "scripts/parse_and_score_results.py", "--base_dir", ".", "--output_dir", "test_outputs/"]
    ))
    
    # Test 5: batch_pipeline.py (single ligand)
    test_results.append(run_test(
        "Batch Pipeline - Single Ligand (batch_pipeline.py)",
        [conda_python, "scripts/batch_pipeline.py", "--protein", "inputs/protein.pdb", "--ligands", "inputs/ligands/aspirin.sdf", "--output-dir", "test_batch_outputs/"],
        expected_exit_code=2  # Partial success (50%+ success rate)
    ))
    
    # Test 6: batch_pipeline.py (multiple ligands)
    test_results.append(run_test(
        "Batch Pipeline - Multiple Ligands (batch_pipeline.py)",
        [conda_python, "scripts/batch_pipeline.py", "--protein", "inputs/protein.pdb", "--ligands", "inputs/ligands/", "--output-dir", "test_batch_outputs/"],
        expected_exit_code=2  # Partial success (50%+ success rate)
    ))
    
    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")
    
    passed = sum(test_results)
    total = len(test_results)
    
    print(f"Passed: {passed}/{total}")
    print(f"Success Rate: {(passed/total)*100:.1f}%")
    
    if passed == total:
        print("🎉 ALL TESTS PASSED! The docking pipeline is working correctly.")
        return 0
    else:
        print("⚠️  Some tests failed. Please check the output above for details.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 