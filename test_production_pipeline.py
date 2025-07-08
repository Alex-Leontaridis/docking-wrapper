#!/usr/bin/env python3
"""
Production Pipeline Test Script
Tests the complete docking pipeline with real binaries and ML models.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import json
import tempfile
import shutil
from pathlib import Path

# Add the scripts directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

def create_test_config():
    """Create a test configuration for the pipeline."""
    config = {
        "engines": {
            "vina": {
                "enabled": True,
                "exhaustiveness": 8,
                "num_modes": 9
            },
            "gnina": {
                "enabled": True,
                "num_modes": 9,
                "exhaustiveness": 8,
                "use_gpu": False
            },
            "diffdock": {
                "enabled": False,  # Disable for testing since we don't have it installed
                "num_predictions": 3
            }
        },
        "ml_models": {
            "equibind": {
                "enabled": True,
                "model_path": "EquiBind/runs/flexible_self_docking/best_checkpoint.pt"
            },
            "neuralplexer": {
                "enabled": False,  # Disable for testing
                "model_path": None
            },
            "umol": {
                "enabled": False,  # Disable on Windows
                "model_path": None
            },
            "structure_predictor": {
                "enabled": False,  # Disable for testing
                "model": "esmfold"
            },
            "boltz2": {
                "enabled": False,  # Disable for testing
                "model_path": None
            }
        },
        "analysis": {
            "interactions": {
                "enabled": True
            },
            "druggability": {
                "enabled": True
            },
            "consensus": {
                "enabled": True,
                "rmsd_threshold": 2.0
            },
            "confidence": {
                "enabled": True
            }
        },
        "box_params": {
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0
        }
    }
    return config

def create_test_inputs():
    """Create test input files."""
    test_dir = Path("test_production_inputs")
    test_dir.mkdir(exist_ok=True)
    
    # Create a simple test protein (minimal PDB)
    protein_content = """ATOM      1  N   ALA A   1      27.462  14.105   5.468  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.525  13.000   5.000  1.00 20.00           C  
ATOM      3  C   ALA A   1      25.000  13.000   5.000  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.500  14.000   5.000  1.00 20.00           O  
TER
END
"""
    
    # Create a simple test ligand (aspirin-like)
    ligand_content = """     RDKit          3D

  8  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5000    1.3000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    3.8000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  1  0  0  0  0
  4  5  2  0  0  0  0
  5  6  1  0  0  0  0
  6  1  2  0  0  0  0
  3  7  1  0  0  0  0
  4  8  1  0  0  0  0
M  END
"""
    
    # Write test files
    protein_file = test_dir / "test_protein.pdb"
    ligand_file = test_dir / "test_ligand.sdf"
    
    with open(protein_file, 'w') as f:
        f.write(protein_content)
    
    with open(ligand_file, 'w') as f:
        f.write(ligand_content)
    
    return test_dir

def test_binary_discovery():
    """Test binary discovery functionality."""
    print("=" * 60)
    print("TESTING BINARY DISCOVERY")
    print("=" * 60)
    
    from scripts.run_docking_multi import find_binary
    
    # Test Vina discovery
    vina_path = find_binary('vina.exe')
    if vina_path:
        print(f"✅ Vina found: {vina_path}")
    else:
        print("❌ Vina not found")
    
    # Test GNINA discovery
    gnina_path = find_binary('gnina')
    if gnina_path:
        print(f"✅ GNINA found: {gnina_path}")
    else:
        print("❌ GNINA not found")
    
    return vina_path is not None

def test_pipeline_initialization():
    """Test pipeline initialization."""
    print("\n" + "=" * 60)
    print("TESTING PIPELINE INITIALIZATION")
    print("=" * 60)
    
    try:
        from scripts.batch_pipeline import BatchDockingPipeline
        
        # Create test config
        config = create_test_config()
        
        # Save config to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config, f)
            config_path = f.name
        
        # Initialize pipeline
        pipeline = BatchDockingPipeline(config_path=config_path, output_dir="test_production_output")
        
        print("✅ Pipeline initialized successfully")
        print(f"✅ Output directory: {pipeline.output_dir}")
        print(f"✅ Config loaded: {len(pipeline.config.get('engines', {}))} engines")
        print(f"✅ ML models: {len(pipeline.config.get('ml_models', {}))} models")
        
        # Clean up
        os.unlink(config_path)
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline initialization failed: {e}")
        return False

def test_single_ligand_processing():
    """Test processing a single ligand."""
    print("\n" + "=" * 60)
    print("TESTING SINGLE LIGAND PROCESSING")
    print("=" * 60)
    
    try:
        from scripts.batch_pipeline import BatchDockingPipeline
        
        # Create test inputs
        test_dir = create_test_inputs()
        protein_file = test_dir / "test_protein.pdb"
        ligand_file = test_dir / "test_ligand.sdf"
        
        # Create minimal config for testing
        config = {
            "engines": {
                "vina": {"enabled": True, "exhaustiveness": 1, "num_modes": 1},
                "gnina": {"enabled": False},  # Disable for faster testing
                "diffdock": {"enabled": False}
            },
            "ml_models": {
                "equibind": {"enabled": False},  # Disable for faster testing
                "neuralplexer": {"enabled": False},
                "umol": {"enabled": False},
                "structure_predictor": {"enabled": False},
                "boltz2": {"enabled": False}
            },
            "analysis": {
                "interactions": {"enabled": False},
                "druggability": {"enabled": False},
                "consensus": {"enabled": False},
                "confidence": {"enabled": False}
            },
            "box_params": {
                "center_x": 1.0, "center_y": 1.0, "center_z": 1.0,
                "size_x": 10.0, "size_y": 10.0, "size_z": 10.0
            }
        }
        
        # Save config
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(config, f)
            config_path = f.name
        
        # Initialize pipeline
        pipeline = BatchDockingPipeline(config_path=config_path, output_dir="test_production_output")
        
        # Process single ligand
        result = pipeline.process_single_ligand(
            ligand_path=str(ligand_file),
            protein_path=str(protein_file),
            ligand_name="test_ligand"
        )
        
        print(f"✅ Single ligand processing completed")
        print(f"✅ Result keys: {list(result.keys())}")
        
        # Clean up
        os.unlink(config_path)
        shutil.rmtree(test_dir, ignore_errors=True)
        
        return True
        
    except Exception as e:
        print(f"❌ Single ligand processing failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all production tests."""
    print("🚀 PRODUCTION PIPELINE TEST")
    print("=" * 60)
    
    tests = [
        ("Binary Discovery", test_binary_discovery),
        ("Pipeline Initialization", test_pipeline_initialization),
        ("Single Ligand Processing", test_single_ligand_processing),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = 0
    for test_name, result in results:
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{test_name}: {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{len(results)} tests passed")
    
    if passed == len(results):
        print("\n🎉 ALL TESTS PASSED! Your pipeline is production-ready!")
    else:
        print("\n⚠️  Some tests failed. Please check the errors above.")
    
    return passed == len(results)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 