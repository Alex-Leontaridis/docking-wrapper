#!/usr/bin/env python3
"""
Test script for the batch molecular docking pipeline
Tests basic functionality without requiring full docking runs
"""

import os
import sys
import tempfile
import json
import shutil
from pathlib import Path

def create_test_data():
    """Create minimal test data for validation."""
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)
    
    # Create test protein (minimal PDB)
    test_protein = test_dir / "test_protein.pdb"
    with open(test_protein, 'w') as f:
        f.write("""HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1      20.154  15.343  15.259  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.059  15.859  1.00 20.00           C  
ATOM      3  C   ALA A   1      19.444  16.740  17.152  1.00 20.00           C  
ATOM      4  O   ALA A   1      20.556  17.231  17.282  1.00 20.00           O  
ATOM      5  CB  ALA A   1      18.461  17.000  14.873  1.00 20.00           C  
TER
END
""")
    
    # Create test ligands directory
    ligands_dir = test_dir / "ligands"
    ligands_dir.mkdir(exist_ok=True)
    
    # Create test SMILES ligands
    ligands = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O"
    }
    
    for name, smiles in ligands.items():
        ligand_file = ligands_dir / f"{name}.smi"
        with open(ligand_file, 'w') as f:
            f.write(f"{smiles}\t{name}\n")
    
    return test_protein, ligands_dir

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    
    try:
        # Add scripts directory to path
        scripts_dir = Path(__file__).parent / "scripts"
        sys.path.insert(0, str(scripts_dir))
        
        from batch_pipeline import BatchDockingPipeline, DEFAULT_CONFIG
        from prep_structures import prepare_protein, prepare_ligand_single
        from run_docking_multi import extract_box_from_protein
        from parse_and_score_results import DockingResultsParser
        print("✓ All imports successful")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_config_loading():
    """Test configuration loading."""
    print("Testing configuration loading...")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        
        # Test default config
        pipeline = BatchDockingPipeline()
        assert "engines" in pipeline.config
        assert "vina" in pipeline.config["engines"]
        
        # Test custom config
        test_config = {"engines": {"vina": {"enabled": False}}}
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_config, f)
            config_file = f.name
        
        pipeline = BatchDockingPipeline(config_path=config_file)
        assert not pipeline.config["engines"]["vina"]["enabled"]
        
        os.unlink(config_file)
        print("✓ Configuration loading successful")
        return True
    except Exception as e:
        print(f"✗ Configuration loading failed: {e}")
        return False

def test_ligand_discovery():
    """Test ligand discovery functionality."""
    print("Testing ligand discovery...")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        
        test_protein, ligands_dir = create_test_data()
        pipeline = BatchDockingPipeline()
        
        # Test directory discovery
        ligands = pipeline.discover_ligands(str(ligands_dir))
        assert len(ligands) == 3
        assert all(isinstance(lig, tuple) and len(lig) == 2 for lig in ligands)
        
        # Test single file discovery
        single_ligand = ligands_dir / "aspirin.smi"
        ligands = pipeline.discover_ligands(str(single_ligand))
        assert len(ligands) == 1
        assert ligands[0][0] == "aspirin"
        
        print("✓ Ligand discovery successful")
        return True
    except Exception as e:
        print(f"✗ Ligand discovery failed: {e}")
        return False

def test_box_detection():
    """Test automatic box detection."""
    print("Testing box detection...")
    
    try:
        from run_docking_multi import extract_box_from_protein
        
        test_protein, _ = create_test_data()
        
        # This should fall back to default parameters since no ligand in test protein
        box_params = extract_box_from_protein(str(test_protein))
        assert len(box_params) == 6  # center_x, center_y, center_z, size_x, size_y, size_z
        assert all(isinstance(param, (int, float)) for param in box_params)
        
        print("✓ Box detection successful")
        return True
    except Exception as e:
        print(f"✗ Box detection failed: {e}")
        return False

def test_dry_run():
    """Test dry run of pipeline without actual docking."""
    print("Testing pipeline dry run...")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        
        test_protein, ligands_dir = create_test_data()
        
        # Create a test config that won't actually run docking
        test_config = {
            "engines": {
                "vina": {"enabled": False},
                "gnina": {"enabled": False},
                "diffdock": {"enabled": False}
            },
            "parallel": {"max_workers": 1}
        }
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_file = Path(tmp_dir) / "test_config.json"
            with open(config_file, 'w') as f:
                json.dump(test_config, f)
            
            pipeline = BatchDockingPipeline(
                config_path=str(config_file),
                output_dir=tmp_dir
            )
            
            # Test ligand discovery
            ligands = pipeline.discover_ligands(str(ligands_dir))
            assert len(ligands) > 0
            
            # Test protein preparation (dry run - just validate input)
            from prep_structures import validate_file, SUPPORTED_PROTEIN_EXT
            validate_file(str(test_protein), SUPPORTED_PROTEIN_EXT, "Protein")
            
        print("✓ Pipeline dry run successful")
        return True
    except Exception as e:
        print(f"✗ Pipeline dry run failed: {e}")
        return False

def test_output_structure():
    """Test output directory structure creation."""
    print("Testing output structure...")
    
    try:
        from batch_pipeline import BatchDockingPipeline
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            pipeline = BatchDockingPipeline(output_dir=tmp_dir)
            
            # Check that directories were created
            expected_dirs = [
                "logs", "prepared_structures", "docking_results", 
                "parsed_results", "summary_reports"
            ]
            
            for dir_name in expected_dirs:
                assert (Path(tmp_dir) / dir_name).exists()
            
        print("✓ Output structure test successful")
        return True
    except Exception as e:
        print(f"✗ Output structure test failed: {e}")
        return False

def cleanup_test_data():
    """Clean up test data."""
    test_dir = Path("test_data")
    if test_dir.exists():
        shutil.rmtree(test_dir)

def main():
    """Run all tests."""
    print("="*60)
    print("BATCH DOCKING PIPELINE VALIDATION TESTS")
    print("="*60)
    
    tests = [
        test_imports,
        test_config_loading,
        test_ligand_discovery,
        test_box_detection,
        test_dry_run,
        test_output_structure
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"✗ {test.__name__} failed with exception: {e}")
            failed += 1
        print()
    
    # Cleanup
    cleanup_test_data()
    
    print("="*60)
    print(f"TEST RESULTS: {passed} passed, {failed} failed")
    print("="*60)
    
    if failed == 0:
        print("🎉 All tests passed! The batch pipeline is ready to use.")
        return 0
    else:
        print("❌ Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 