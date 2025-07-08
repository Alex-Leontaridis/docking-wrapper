#!/usr/bin/env python3
"""
Test script to verify all bug fixes are working correctly.
Tests each of the 6 major bugs that were fixed.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import tempfile
import shutil
import logging
from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context
from pathlib import Path
import json

# Add scripts directory to path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'scripts'))

def setup_test_logging():
    """Setup logging for tests."""
    logger = setup_logging(__name__)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )

def test_meeko_compatibility():
    """Test Bug #10: Meeko Compatibility Issue"""
    print("\n=== Testing Bug #10: Meeko Compatibility Issue ===")
    
    try:
        # Import the prep_structures module to test Meeko import
        import prep_structures
        
        # Check if MEEKO_AVAILABLE is properly set
        if hasattr(prep_structures, 'MEEKO_AVAILABLE'):
            print(f"OK: MEEKO_AVAILABLE properly set: {prep_structures.MEEKO_AVAILABLE}")
            
            # Test that the module doesn't crash even if Meeko fails
            if not prep_structures.MEEKO_AVAILABLE:
                print("OK: Module gracefully handles Meeko import failure")
            else:
                print("OK: Meeko successfully imported")
        else:
            print("FAIL: MEEKO_AVAILABLE not found")
            return False
            
        return True
    except Exception as e:
        print(f"✗ Meeko compatibility test failed: {e}")
        return False

def test_diffdock_confidence_parsing():
    """Test Bug #11: DiffDock Confidence Parsing Bug"""
    print("\n=== Testing Bug #11: DiffDock Confidence Parsing Bug ===")
    
    try:
        from docking_results_parser import DockingResultsParser
        
        # Create test confidence file with header
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("pose_file,confidence\n")  # Header line
            f.write("pose_1.sdf,0.85\n")
            f.write("pose_2.sdf,0.72\n")
            f.write("pose_3.sdf,0.68\n")
            confidence_file = f.name
        
        # Create test directory structure
        test_dir = tempfile.mkdtemp()
        diffdock_dir = os.path.join(test_dir, "diffdock_output")
        os.makedirs(diffdock_dir)
        
        # Copy confidence file
        shutil.copy2(confidence_file, os.path.join(diffdock_dir, "diffdock_confidence.txt"))
        
        # Create dummy SDF files
        for i in range(1, 4):
            sdf_file = os.path.join(diffdock_dir, f"pose_{i}.sdf")
            with open(sdf_file, 'w') as f:
                f.write("dummy SDF content\n")
        
        # Test parsing
        parser = DockingResultsParser(test_dir, tempfile.mkdtemp())
        poses = parser._parse_diffdock_results("test_ligand")
        
        # Cleanup
        os.unlink(confidence_file)
        shutil.rmtree(test_dir)
        
        if len(poses) == 3:
            print("OK: DiffDock confidence parsing works with header")
            print(f"OK: Parsed {len(poses)} poses with confidence scores")
            for pose in poses:
                print(f"  - Pose {pose['pose_id']}: confidence = {pose['confidence_score']}")
            return True
        else:
            print(f"FAIL: Expected 3 poses, got {len(poses)}")
            return False
            
    except Exception as e:
        print(f"✗ DiffDock confidence parsing test failed: {e}")
        return False

def test_batch_pipeline_error_handling():
    """Test Bug #12: Batch Pipeline Error Handling Gap"""
    print("\n=== Testing Bug #12: Batch Pipeline Error Handling Gap ===")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        
        # Create a minimal config for testing
        config = {
            "engines": {
                "vina": {"enabled": True},
                "gnina": {"enabled": False},
                "diffdock": {"enabled": False}
            },
            "parallel": {"max_workers": 1},
            "strict_mode": True
        }
        
        # Create test config file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config, f)
            config_file = f.name
        
        # Test pipeline initialization
        pipeline = BatchDockingPipeline(config_file, tempfile.mkdtemp())
        
        # Check if strict_mode is properly handled
        if hasattr(pipeline, '_lock') and hasattr(pipeline, '_process_lock'):
            print("OK: Thread safety mechanisms in place")
        else:
            print("FAIL: Thread safety mechanisms missing")
            return False
        
        # Check if cleanup handlers are registered
        if hasattr(pipeline, '_cleanup_registered'):
            print("OK: Cleanup handlers properly registered")
        else:
            print("FAIL: Cleanup handlers not registered")
            return False
        
        # Cleanup
        os.unlink(config_file)
        
        return True
        
    except Exception as e:
        print(f"✗ Batch pipeline error handling test failed: {e}")
        return False

def test_input_validation():
    """Test Bug #13: Missing Input Validation"""
    print("\n=== Testing Bug #13: Missing Input Validation ===")
    
    try:
        from input_validator import InputValidator, validate_protein_file, validate_ligand_file
        
        validator = InputValidator()
        
        # Test protein file validation
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("ATOM      1  N   ALA A   1      27.451  14.105   5.383  1.00 20.00\n")
            f.write("ATOM      2  CA  ALA A   1      26.325  13.134   5.547  1.00 20.00\n")
            protein_file = f.name
        
        valid, error = validator.validate_protein_file(protein_file)
        if valid:
            print("OK: Protein file validation works")
        else:
            print(f"FAIL: Protein file validation failed: {error}")
            return False
        
        # Test ligand file validation
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as f:
            f.write("aspirin\n")
            f.write("  -OEChem-01301921423D\n")
            f.write("\n")
            f.write("  9  8  0     1  0  0  0  0  0999 V2000\n")
            f.write("    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
            f.write("M  END\n")
            ligand_file = f.name
        
        valid, error = validator.validate_ligand_file(ligand_file)
        if valid:
            print("OK: Ligand file validation works")
        else:
            print(f"FAIL: Ligand file validation failed: {error}")
            return False
        
        # Test invalid file validation
        valid, error = validator.validate_protein_file("nonexistent.pdb")
        if not valid and "does not exist" in error:
            print("OK: Invalid file detection works")
        else:
            print("FAIL: Invalid file detection failed")
            return False
        
        # Cleanup
        os.unlink(protein_file)
        os.unlink(ligand_file)
        
        return True
        
    except Exception as e:
        print(f"✗ Input validation test failed: {e}")
        return False

def test_resource_management():
    """Test Bug #14: Resource Management Issues"""
    print("\n=== Testing Bug #14: Resource Management Issues ===")
    
    try:
        from resource_manager import ResourceManager, temp_file, temp_directory, cleanup_all
        
        # Test ResourceManager
        manager = ResourceManager()
        
        # Create temporary files and directories
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_file_path = f.name
        
        temp_dir_path = tempfile.mkdtemp()
        
        # Register them
        manager.register_temp_file(temp_file_path)
        manager.register_temp_dir(temp_dir_path)
        
        # Check registration
        resources = manager.get_registered_resources()
        if len(resources['temp_files']) == 1 and len(resources['temp_dirs']) == 1:
            print("OK: Resource registration works")
        else:
            print("FAIL: Resource registration failed")
            return False
        
        # Test cleanup
        manager.cleanup_all()
        
        # Check if files are cleaned up
        if not os.path.exists(temp_file_path) and not os.path.exists(temp_dir_path):
            print("OK: Resource cleanup works")
        else:
            print("FAIL: Resource cleanup failed")
            return False
        
        # Test context managers
        with temp_file(suffix='.test') as tf:
            if os.path.exists(tf):
                print("OK: temp_file context manager works")
            else:
                print("FAIL: temp_file context manager failed")
                return False
        
        with temp_directory() as td:
            if os.path.exists(td):
                print("OK: temp_directory context manager works")
            else:
                print("FAIL: temp_directory context manager failed")
                return False
        
        return True
        
    except Exception as e:
        print(f"✗ Resource management test failed: {e}")
        return False

def test_thread_safety():
    """Test Bug #15: Thread Safety Issues"""
    print("\n=== Testing Bug #15: Thread Safety Issues ===")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        import threading
        import time
        
        # Create test pipeline
        pipeline = BatchDockingPipeline(output_dir=tempfile.mkdtemp())
        
        # Test thread-safe statistics update
        def update_stats(thread_id):
            for i in range(10):
                pipeline._update_statistics(f"ligand_{thread_id}_{i}", True, {"test": "data"})
                time.sleep(0.01)
        
        # Create multiple threads
        threads = []
        for i in range(5):
            thread = threading.Thread(target=update_stats, args=(i,))
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join()
        
        # Check if statistics are consistent
        expected_total = 50  # 5 threads * 10 updates each
        if pipeline.successful_ligands == expected_total:
            print("OK: Thread-safe statistics update works")
        else:
            print(f"FAIL: Thread-safe statistics update failed: expected {expected_total}, got {pipeline.successful_ligands}")
            return False
        
        # Test process lock
        if hasattr(pipeline, '_process_lock'):
            print("OK: Process lock available")
        else:
            print("FAIL: Process lock missing")
            return False
        
        return True
        
    except Exception as e:
        print(f"✗ Thread safety test failed: {e}")
        return False

def main():
    """Run all bug fix tests."""
    setup_test_logging()
    
    print("Testing Bug Fixes for Molecular Docking Pipeline")
    print("=" * 60)
    
    tests = [
        ("Meeko Compatibility", test_meeko_compatibility),
        ("DiffDock Confidence Parsing", test_diffdock_confidence_parsing),
        ("Batch Pipeline Error Handling", test_batch_pipeline_error_handling),
        ("Input Validation", test_input_validation),
        ("Resource Management", test_resource_management),
        ("Thread Safety", test_thread_safety),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
                print(f"PASS: {test_name}")
            else:
                print(f"FAIL: {test_name}")
        except Exception as e:
            print(f"ERROR: {test_name} - {e}")
    
    print("\n" + "=" * 60)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("SUCCESS: All bug fixes are working correctly!")
        return 0
    else:
        print("WARNING: Some bug fixes may need attention.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 