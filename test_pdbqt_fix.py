#!/usr/bin/env python3
"""
Test script to verify PDBQT conversion fix.
"""

import os
import sys
import tempfile
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

def test_pdbqt_conversion():
    """Test the PDBQT conversion fix."""
    print("=" * 60)
    print("TESTING PDBQT CONVERSION FIX")
    print("=" * 60)
    
    # Create a test PDB file that matches the error
    test_pdb_content = """ATOM      1  N   PRO A   1      29.361  39.686   5.862  1.00 38.10           N  
ATOM      2  CA  PRO A   1      28.336  40.234   5.234  1.00 38.10           C  
ATOM      3  C   PRO A   1      27.085  39.456   5.567  1.00 38.10           C  
ATOM      4  O   PRO A   1      26.000  40.000   5.000  1.00 38.10           O  
END"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdb = os.path.join(temp_dir, 'test.pdb')
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdb, 'w') as f:
            f.write(test_pdb_content)
        
        # Test the conversion using batch pipeline method
        try:
            from batch_pipeline import BatchDockingPipeline
            pipeline = BatchDockingPipeline()
            pipeline._simple_pdb_to_pdbqt(test_pdb, test_pdbqt)
            
            # Check the output
            with open(test_pdbqt, 'r') as f:
                lines = f.readlines()
            
            print(f"Generated {len(lines)} lines in PDBQT file")
            
            # Check for proper AutoDock types
            valid_types = {'C', 'N', 'O', 'S', 'P', 'H'}
            atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
            
            print(f"Found {len(atom_lines)} atom lines")
            
            # Check each atom line
            for i, line in enumerate(atom_lines):
                print(f"Line {i+1}: {line.strip()}")
                
                # Check for the problematic format (specifically the space between 0.000 and atom type)
                if ' 0.000 N' in line or ' 0.000 C' in line or ' 0.000 O' in line:
                    print(f"   ❌ Found problematic format!")
                    return False
                
                # Check that line ends with proper AutoDock type
                parts = line.strip().split()
                if len(parts) >= 2:
                    atom_type = parts[-1]  # Last part should be the AutoDock type
                    if atom_type not in valid_types:
                        print(f"   ❌ Invalid atom type: {atom_type}")
                        return False
                    else:
                        print(f"   ✅ Valid atom type: {atom_type}")
            
            # Check for proper charge format (should not have +0.000 or extra spaces)
            charge_bug = any('+0.000' in line for line in lines)
            if charge_bug:
                print("❌ Found +0.000 bug in PDBQT file")
                return False
            else:
                print("✅ No +0.000 bug found")
            
            # Check for the specific error format (extra spaces before 0.000)
            error_format = any('             0.000' in line for line in lines)
            if error_format:
                print("❌ Found extra spaces before 0.000")
                return False
            else:
                print("✅ No extra spaces found")
            
            print("✅ PDBQT conversion test passed!")
            return True
            
        except Exception as e:
            print(f"❌ PDBQT conversion test failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def test_prep_structures_conversion():
    """Test the prep_structures PDBQT conversion."""
    print("\n" + "=" * 60)
    print("TESTING PREP_STRUCTURES PDBQT CONVERSION")
    print("=" * 60)
    
    # Create a test PDB file
    test_pdb_content = """ATOM      1  N   PRO A   1      29.361  39.686   5.862  1.00 38.10           N  
ATOM      2  CA  PRO A   1      28.336  40.234   5.234  1.00 38.10           C  
ATOM      3  C   PRO A   1      27.085  39.456   5.567  1.00 38.10           C  
ATOM      4  O   PRO A   1      26.000  40.000   5.000  1.00 38.10           O  
END"""
    
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdb = os.path.join(temp_dir, 'test.pdb')
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdb, 'w') as f:
            f.write(test_pdb_content)
        
        # Test the conversion using prep_structures method
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
            for i, line in enumerate(atom_lines):
                print(f"Line {i+1}: {line.strip()}")
                
                # Check for the problematic format (specifically the space between 0.000 and atom type)
                if ' 0.000 N' in line or ' 0.000 C' in line or ' 0.000 O' in line:
                    print(f"   ❌ Found problematic format!")
                    return False
                
                # Check that line ends with proper AutoDock type
                parts = line.strip().split()
                if len(parts) >= 2:
                    atom_type = parts[-1]  # Last part should be the AutoDock type
                    if atom_type not in valid_types:
                        print(f"   ❌ Invalid atom type: {atom_type}")
                        return False
                    else:
                        print(f"   ✅ Valid atom type: {atom_type}")
            
            # Check for proper charge format
            charge_bug = any('+0.000' in line for line in lines)
            if charge_bug:
                print("❌ Found +0.000 bug in PDBQT file")
                return False
            else:
                print("✅ No +0.000 bug found")
            
            # Check for the specific error format
            error_format = any('             0.000' in line for line in lines)
            if error_format:
                print("❌ Found extra spaces before 0.000")
                return False
            else:
                print("✅ No extra spaces found")
            
            print("✅ Prep_structures PDBQT conversion test passed!")
            return True
            
        except Exception as e:
            print(f"❌ Prep_structures PDBQT conversion test failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Run all tests."""
    print("Testing PDBQT conversion fixes...")
    
    # Test batch pipeline conversion
    batch_ok = test_pdbqt_conversion()
    
    # Test prep_structures conversion
    prep_ok = test_prep_structures_conversion()
    
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    if batch_ok and prep_ok:
        print("🎉 All PDBQT conversion tests passed!")
        print("\nThe PDBQT conversion fix should now work correctly.")
        print("Run your pipeline again to test the fix.")
        return 0
    else:
        print("❌ Some PDBQT conversion tests failed:")
        if not batch_ok:
            print("   - Batch pipeline PDBQT conversion failed")
        if not prep_ok:
            print("   - Prep_structures PDBQT conversion failed")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 