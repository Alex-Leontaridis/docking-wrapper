#!/usr/bin/env python3
"""
Simple test script to run the batch pipeline directly.
"""

import sys
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
from pathlib import Path

# Add scripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'scripts'))

def main():
    """Test the batch pipeline with real inputs."""
    try:
        from batch_pipeline import BatchDockingPipeline
        
        print("🚀 Testing Batch Pipeline with Real Binaries")
        print("=" * 60)
        
        # Check if test files exist
        protein_file = "inputs/protein.pdb"
        ligand_file = "inputs/ligands/aspirin.sdf"
        
        if not os.path.exists(protein_file):
            print(f"❌ Protein file not found: {protein_file}")
            return False
            
        if not os.path.exists(ligand_file):
            print(f"❌ Ligand file not found: {ligand_file}")
            return False
        
        print(f"✅ Protein file: {protein_file}")
        print(f"✅ Ligand file: {ligand_file}")
        
        # Create minimal config
        config = {
            "engines": {
                "vina": {"enabled": True, "exhaustiveness": 1, "num_modes": 1},
                "gnina": {"enabled": True, "num_modes": 1, "exhaustiveness": 1, "use_gpu": False},
                "diffdock": {"enabled": False}
            },
            "ml_models": {
                "equibind": {"enabled": False},
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
        config_file = "test_simple_config.json"
        import json
        with open(config_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        print(f"✅ Config saved: {config_file}")
        
        # Initialize pipeline
        print("\n🔧 Initializing pipeline...")
        pipeline = BatchDockingPipeline(config_path=config_file, output_dir="test_simple_output")
        
        print("✅ Pipeline initialized successfully")
        
        # Process single ligand
        print("\n🔬 Processing ligand...")
        ligand_info = ("aspirin", ligand_file)
        result = pipeline.process_single_ligand(
            ligand_info=ligand_info,
            prepared_protein=protein_file
        )
        
        print("✅ Processing completed")
        print(f"Result keys: {list(result.keys())}")
        
        # Check results
        if result.get('success', False):
            print("🎉 SUCCESS: Pipeline completed successfully!")
            
            # Check for outputs
            output_dir = Path("test_simple_output")
            if output_dir.exists():
                print(f"📁 Output directory: {output_dir}")
                for item in output_dir.rglob("*"):
                    if item.is_file():
                        print(f"  📄 {item.relative_to(output_dir)}")
        else:
            print("❌ FAILED: Pipeline did not complete successfully")
            if 'error' in result:
                print(f"Error: {result['error']}")
        
        # Clean up
        if os.path.exists(config_file):
            os.remove(config_file)
        
        return result.get('success', False)
        
    except Exception as e:
        print(f"❌ Exception occurred: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 