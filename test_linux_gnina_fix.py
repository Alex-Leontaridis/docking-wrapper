#!/usr/bin/env python3
"""
Simple test script for GNINA detection on Linux.
Run this on your Linux VirtualBox machine.
"""

import os
import sys
import platform
import subprocess
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

def test_gnina_detection():
    """Test GNINA detection on Linux."""
    print("=" * 60)
    print("TESTING GNINA DETECTION ON LINUX")
    print("=" * 60)
    
    cwd = os.getcwd()
    system = platform.system().lower()
    
    print(f"Current working directory: {cwd}")
    print(f"Platform: {system}")
    
    # Check if we're on Linux
    if system != 'linux':
        print(f"❌ This test is for Linux, but you're on {system}")
        return False
    
    # Check if bin/gnina.sh exists and is executable
    gnina_sh_path = os.path.join(cwd, 'bin', 'gnina.sh')
    print(f"\n1. Checking bin/gnina.sh:")
    print(f"   Path: {gnina_sh_path}")
    print(f"   Exists: {os.path.exists(gnina_sh_path)}")
    
    if not os.path.exists(gnina_sh_path):
        print("❌ bin/gnina.sh does not exist!")
        print("   Please make sure you copied the file to your Linux machine.")
        return False
    
    print(f"   Is file: {os.path.isfile(gnina_sh_path)}")
    print(f"   Is executable: {os.access(gnina_sh_path, os.X_OK)}")
    
    if not os.access(gnina_sh_path, os.X_OK):
        print("❌ bin/gnina.sh is not executable!")
        print("   Run: chmod +x bin/gnina.sh")
        return False
    
    # Test the script directly
    print(f"\n2. Testing bin/gnina.sh directly:")
    try:
        result = subprocess.run([gnina_sh_path, '--help'], 
                              capture_output=True, text=True, timeout=5)
        print(f"   Return code: {result.returncode}")
        print(f"   Output: {result.stdout.strip()}")
        
        if result.returncode == 0:
            print("   ✅ Script works correctly!")
        else:
            print("   ❌ Script failed!")
            return False
            
    except Exception as e:
        print(f"   ❌ Error running script: {e}")
        return False
    
    # Test the find_binary function
    print(f"\n3. Testing find_binary function:")
    try:
        from run_docking_multi import find_binary
        gnina_path = find_binary('gnina')
        print(f"   find_binary('gnina') returned: {gnina_path}")
        
        if gnina_path:
            if 'bin/gnina.sh' in gnina_path:
                print("   ✅ Found bin/gnina.sh correctly!")
                return True
            else:
                print(f"   ❌ Found wrong path: {gnina_path}")
                return False
        else:
            print("   ❌ find_binary returned None")
            return False
            
    except Exception as e:
        print(f"   ❌ Error importing find_binary: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_pdbqt_conversion():
    """Test PDBQT conversion on Linux."""
    print("\n" + "=" * 60)
    print("TESTING PDBQT CONVERSION ON LINUX")
    print("=" * 60)
    
    # Create a simple test PDB file
    test_pdb_content = """ATOM      1  N   PRO A   1      29.361  39.686   5.862  1.00 38.10           N  
ATOM      2  CA  PRO A   1      28.336  40.234   5.234  1.00 38.10           C  
ATOM      3  C   PRO A   1      27.085  39.456   5.567  1.00 38.10           C  
ATOM      4  O   PRO A   1      26.000  40.000   5.000  1.00 38.10           O  
END"""
    
    import tempfile
    with tempfile.TemporaryDirectory() as temp_dir:
        test_pdb = os.path.join(temp_dir, 'test.pdb')
        test_pdbqt = os.path.join(temp_dir, 'test.pdbqt')
        
        with open(test_pdb, 'w') as f:
            f.write(test_pdb_content)
        
        # Test the conversion
        try:
            from prep_structures import _simple_pdb_to_pdbqt
            _simple_pdb_to_pdbqt(test_pdb, test_pdbqt)
            
            # Check the output
            with open(test_pdbqt, 'r') as f:
                lines = f.readlines()
            
            print(f"Generated {len(lines)} lines in PDBQT file")
            
            # Check for the problematic format
            atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
            print(f"Found {len(atom_lines)} atom lines")
            
            # Check each atom line
            for i, line in enumerate(atom_lines[:3]):  # Check first 3 lines
                print(f"Line {i+1}: {line.strip()}")
                
                # Check for the problematic format
                if '0.000 N' in line or '0.000 C' in line or '0.000 O' in line:
                    print(f"   ❌ Found problematic format!")
                    return False
            
            print("   ✅ PDBQT conversion looks correct!")
            return True
            
        except Exception as e:
            print(f"   ❌ PDBQT conversion failed: {e}")
            import traceback
            traceback.print_exc()
            return False

def main():
    """Run all tests."""
    print("Testing Linux fixes for docking wrapper...")
    
    # Test GNINA detection
    gnina_ok = test_gnina_detection()
    
    # Test PDBQT conversion
    pdbqt_ok = test_pdbqt_conversion()
    
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    if gnina_ok and pdbqt_ok:
        print("🎉 All tests passed! Linux fixes are working correctly.")
        print("\nNext steps:")
        print("1. Run your pipeline: python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --output-dir test_output/")
        return 0
    else:
        print("❌ Some tests failed:")
        if not gnina_ok:
            print("   - GNINA detection failed")
        if not pdbqt_ok:
            print("   - PDBQT conversion failed")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 