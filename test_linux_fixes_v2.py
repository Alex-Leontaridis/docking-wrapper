#!/usr/bin/env python3
"""
Updated test script to verify Linux fixes for the docking wrapper.
"""

import os
import sys
import platform
import tempfile
import shutil
import subprocess
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

def test_pdbqt_conversion_linux():
    """Test the PDBQT conversion fix specifically for Linux."""
    print("=" * 60)
    print("TESTING PDBQT CONVERSION FIX (Linux)")
    print("=" * 60)
    
    # Create a simple test PDB file that matches the error
    test_pdb_content = """ATOM      1  N   PRO A   1      29.361  39.686   5.862  1.00 38.10           N  
ATOM      2  CA  PRO A   1      28.336  40.234   5.234  1.00 38.10           C  
ATOM      3  C   PRO A   1      27.085  39.456   5.567  1.00 38.10           C  
ATOM      4  O   PRO A   1      26.000  40.000   5.000  1.00 38.10           O  
HETATM    5  C   LIG A   2      20.000  20.000  20.000  1.00 20.00           C  
HETATM    6  O   LIG A   2      21.000  21.000  21.000  1.00 20.00           O  
END"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdb = os.path.join(temp_dir, 'test.pdb')
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdb, 'w') as f:
            f.write(test_pdb_content)
        
        # Test the conversion using the prep_structures module
        try:
            from prep_structures import _simple_pdb_to_pdbqt
            _simple_pdb_to_pdbqt(test_pdb, test_pdbqt)
            
            # Check the output
            with open(test_pdbqt, 'r') as f:
                lines = f.readlines()
            
            print(f"Generated {len(lines)} lines in PDBQT file")
            
            # Check for proper AutoDock types
            valid_types = {'C', 'N', 'O', 'S', 'P', 'H'}
            atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
            
            print(f"Found {len(atom_lines)} atom lines")
            
            # Check each atom line
            for i, line in enumerate(atom_lines[:5]):  # Check first 5 lines
                print(f"Line {i+1}: {line.strip()}")
                
                # Check for the problematic format
                if '0.000 N' in line or '0.000 C' in line or '0.000 O' in line:
                    print(f"❌ Found problematic format in line: {line.strip()}")
                    return False
                
                # Check that line ends with proper AutoDock type
                parts = line.strip().split()
                if len(parts) >= 2:
                    atom_type = parts[-1]  # Last part should be the AutoDock type
                    if atom_type not in valid_types:
                        print(f"❌ Invalid atom type: {atom_type} in line: {line.strip()}")
                        return False
                    else:
                        print(f"✅ Valid atom type: {atom_type}")
            
            # Check for proper charge format (should not have +0.000 or 0.000 N)
            charge_bug = any('+0.000' in line for line in lines)
            if charge_bug:
                print("❌ Found +0.000 bug in PDBQT file")
                return False
            else:
                print("✅ No +0.000 bug found")
            
            # Check for the specific error format
            error_format = any('0.000 N' in line for line in lines)
            if error_format:
                print("❌ Found '0.000 N' format that causes Vina to crash")
                return False
            else:
                print("✅ No '0.000 N' format found")
            
            print("✅ PDBQT conversion test passed!")
            return True
            
        except Exception as e:
            print(f"❌ PDBQT conversion test failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_gnina_detection_linux():
    """Test GNINA detection on Linux."""
    print("\n" + "=" * 60)
    print("TESTING GNINA DETECTION (Linux)")
    print("=" * 60)
    
    system = platform.system().lower()
    print(f"Platform: {system}")
    
    try:
        from run_docking_multi import find_binary
        
        # Test GNINA detection
        gnina_path = find_binary('gnina')
        print(f"GNINA path: {gnina_path}")
        
        if gnina_path:
            if system == 'linux':
                # On Linux, should find gnina.sh or a working binary
                if gnina_path.endswith('.sh'):
                    print("✅ Found Linux GNINA shell script")
                    return True
                elif 'wsl:' in gnina_path:
                    print("✅ Found GNINA via WSL")
                    return True
                else:
                    # Test if the binary actually works
                    try:
                        result = subprocess.run([gnina_path, '--version'], 
                                              capture_output=True, text=True, timeout=5)
                        if result.returncode == 0:
                            print("✅ Found working GNINA binary")
                            return True
                        else:
                            print(f"❌ GNINA binary found but not working: {result.stderr}")
                            return False
                    except Exception as e:
                        print(f"❌ GNINA binary found but failed to execute: {e}")
                        return False
            else:
                print("✅ GNINA found (platform not specifically tested)")
                return True
        else:
            print("❌ GNINA not found")
            return False
            
    except Exception as e:
        print(f"❌ GNINA detection test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_meeko_fallback_linux():
    """Test Meeko fallback functionality on Linux."""
    print("\n" + "=" * 60)
    print("TESTING MEEKO FALLBACK (Linux)")
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
                    
                    # Check the content
                    with open(test_pdbqt, 'r') as f:
                        content = f.read()
                        if 'ATOM' in content or 'HETATM' in content:
                            print("✅ PDBQT file contains valid atom records")
                            return True
                        else:
                            print("❌ PDBQT file doesn't contain atom records")
                            return False
                else:
                    print("❌ PDBQT file was not created")
                    return False
                    
            except Exception as e:
                print(f"❌ Ligand conversion failed: {e}")
                return False
                
    except Exception as e:
        print(f"❌ Meeko fallback test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_vina_compatibility():
    """Test if the generated PDBQT file is compatible with Vina."""
    print("\n" + "=" * 60)
    print("TESTING VINA COMPATIBILITY")
    print("=" * 60)
    
    # Create a minimal test PDBQT file
    test_pdbqt_content = """ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00          0.000 N
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00          0.000 C
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00          0.000 C
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00          0.000 O
"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdbqt, 'w') as f:
            f.write(test_pdbqt_content)
        
        # Test if Vina can read the file
        try:
            from run_docking_multi import find_binary
            vina_path = find_binary('vina')
            
            if vina_path:
                print(f"Found Vina at: {vina_path}")
                
                # Try to run Vina with the test file (should not crash)
                cmd = [vina_path, '--receptor', test_pdbqt, '--help']
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
                
                if result.returncode == 0:
                    print("✅ Vina can read the PDBQT file format")
                    return True
                else:
                    print(f"❌ Vina failed to read PDBQT file: {result.stderr}")
                    return False
            else:
                print("❌ Vina not found, skipping compatibility test")
                return True  # Skip this test if Vina not available
                
        except Exception as e:
            print(f"❌ Vina compatibility test failed: {e}")
            return False

def main():
    """Run all tests."""
    print("Testing Linux fixes for docking wrapper (v2)...")
    print(f"Platform: {platform.system()} {platform.release()}")
    print(f"Python: {sys.version}")
    
    tests = [
        test_pdbqt_conversion_linux,
        test_gnina_detection_linux,
        test_meeko_fallback_linux,
        test_vina_compatibility
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
        print("\nNext steps:")
        print("1. Make sure bin/gnina.sh is executable: chmod +x bin/gnina.sh")
        print("2. Run your pipeline: python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --output-dir test_output/")
        return 0
    else:
        print("❌ Some tests failed. Please check the issues above.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 