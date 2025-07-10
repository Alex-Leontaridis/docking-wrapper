#!/usr/bin/env python3
"""
Test script to verify Linux fixes for the docking wrapper.
"""

import os
import sys
import platform
import tempfile
import shutil
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

def test_pdbqt_conversion():
    """Test the PDBQT conversion fix."""
    print("=" * 60)
    print("TESTING PDBQT CONVERSION FIX")
    print("=" * 60)
    
    # Create a simple test PDB file
    test_pdb_content = """ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00           O  
HETATM    5  C   LIG A   2      20.000  20.000  20.000  1.00 20.00           C  
HETATM    6  O   LIG A   2      21.000  21.000  21.000  1.00 20.00           O  
END"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdb = os.path.join(temp_dir, 'test.pdb')
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdb, 'w') as f:
            f.write(test_pdb_content)
        
        # Test the conversion
        try:
            from batch_pipeline import BatchDockingPipeline
            pipeline = BatchDockingPipeline()
            pipeline._simple_pdb_to_pdbqt(test_pdb, test_pdbqt)
            
            # Check the output
            with open(test_pdbqt, 'r') as f:
                lines = f.readlines()
            
            # Check for proper AutoDock types
            valid_types = {'C', 'N', 'O', 'S', 'P', 'H'}
            atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
            
            print(f"Generated {len(atom_lines)} atom lines")
            
            # Check each atom line
            for line in atom_lines:
                parts = line.strip().split()
                if len(parts) >= 2:
                    atom_type = parts[-1]  # Last part should be the AutoDock type
                    if atom_type not in valid_types:
                        print(f"❌ Invalid atom type: {atom_type} in line: {line.strip()}")
                        return False
                    else:
                        print(f"✅ Valid atom type: {atom_type}")
            
            # Check for proper charge format (should not have +0.000)
            charge_bug = any('+0.000' in line for line in lines)
            if charge_bug:
                print("❌ Found +0.000 bug in PDBQT file")
                return False
            else:
                print("✅ No +0.000 bug found")
            
            print("✅ PDBQT conversion test passed!")
            return True
            
        except Exception as e:
            print(f"❌ PDBQT conversion test failed: {e}")
            return False

def test_gnina_detection():
    """Test GNINA detection on different platforms."""
    print("\n" + "=" * 60)
    print("TESTING GNINA DETECTION")
    print("=" * 60)
    
    system = platform.system().lower()
    print(f"Platform: {system}")
    
    try:
        from run_docking_multi import find_binary
        
        # Test GNINA detection
        gnina_path = find_binary('gnina')
        print(f"GNINA path: {gnina_path}")
        
        if gnina_path:
            if system == 'windows':
                # On Windows, should find gnina.bat
                if gnina_path.endswith('.bat'):
                    print("✅ Found Windows GNINA batch file")
                    return True
                else:
                    print("❌ Expected Windows batch file but found: " + gnina_path)
                    return False
            elif system == 'linux':
                # On Linux, should find gnina.sh
                if gnina_path.endswith('.sh') or gnina_path.endswith('gnina'):
                    print("✅ Found Linux GNINA script/binary")
                    return True
                else:
                    print("❌ Expected Linux script/binary but found: " + gnina_path)
                    return False
            else:
                print("✅ GNINA found (platform not specifically tested)")
                return True
        else:
            print("❌ GNINA not found")
            return False
            
    except Exception as e:
        print(f"❌ GNINA detection test failed: {e}")
        return False

def test_meeko_fallback():
    """Test Meeko fallback functionality."""
    print("\n" + "=" * 60)
    print("TESTING MEEKO FALLBACK")
    print("=" * 60)
    
    try:
        from prep_structures import convert_ligand_to_pdbqt
        from rdkit import Chem
        
        # Create a simple test molecule
        mol = Chem.MolFromSmiles('CCO')  # Ethanol
        if mol is None:
            print("❌ Failed to create test molecule")
            return False
        
        with tempfile.TemporaryDirectory() as temp_dir:
            test_pdbqt = os.path.join(temp_dir, 'test_ligand.pdbqt')
            
            # Test conversion (should fall back to RDKit if Meeko fails)
            try:
                convert_ligand_to_pdbqt(mol, test_pdbqt)
                
                # Check if file was created
                if os.path.exists(test_pdbqt):
                    print("✅ Ligand conversion completed successfully")
                    return True
                else:
                    print("❌ PDBQT file was not created")
                    return False
                    
            except Exception as e:
                print(f"❌ Ligand conversion failed: {e}")
                return False
                
    except Exception as e:
        print(f"❌ Meeko fallback test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("Testing Linux fixes for docking wrapper...")
    print(f"Platform: {platform.system()} {platform.release()}")
    print(f"Python: {sys.version}")
    
    tests = [
        test_pdbqt_conversion,
        test_gnina_detection,
        test_meeko_fallback
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"❌ Test {test.__name__} crashed: {e}")
            results.append(False)
    
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(results)
    total = len(results)
    
    print(f"Passed: {passed}/{total}")
    
    if passed == total:
        print("🎉 All tests passed! Linux fixes are working correctly.")
        return 0
    else:
        print("❌ Some tests failed. Please check the issues above.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 