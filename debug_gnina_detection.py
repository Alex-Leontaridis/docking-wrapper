#!/usr/bin/env python3
"""
Debug script to test GNINA detection logic.
"""

import os
import sys
import platform
import subprocess
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

def debug_gnina_detection():
    """Debug the GNINA detection process step by step."""
    print("=" * 60)
    print("DEBUGGING GNINA DETECTION")
    print("=" * 60)
    
    cwd = os.getcwd()
    system = platform.system().lower()
    
    print(f"Current working directory: {cwd}")
    print(f"Platform: {system}")
    
    # Check if gnina.sh exists
    gnina_sh_path = os.path.join(cwd, 'bin', 'gnina.sh')
    print(f"\n1. Checking for bin/gnina.sh:")
    print(f"   Path: {gnina_sh_path}")
    print(f"   Exists: {os.path.exists(gnina_sh_path)}")
    
    if os.path.exists(gnina_sh_path):
        print(f"   Is file: {os.path.isfile(gnina_sh_path)}")
        print(f"   Is executable: {os.access(gnina_sh_path, os.X_OK)}")
        
        # Check file permissions
        import stat
        st = os.stat(gnina_sh_path)
        print(f"   Permissions: {oct(st.st_mode)[-3:]}")
        
        # Try to read the file
        try:
            with open(gnina_sh_path, 'r') as f:
                first_line = f.readline().strip()
                print(f"   First line: {first_line}")
        except Exception as e:
            print(f"   Error reading file: {e}")
    
    # Check for gnina in current directory
    print(f"\n2. Checking for gnina in current directory:")
    for ext in ('', '.sh'):
        gnina_local = os.path.join(cwd, f'gnina{ext}')
        print(f"   Checking: {gnina_local}")
        print(f"   Exists: {os.path.exists(gnina_local)}")
        if os.path.exists(gnina_local):
            print(f"   Is file: {os.path.isfile(gnina_local)}")
            print(f"   Is executable: {os.access(gnina_local, os.X_OK)}")
    
    # Check PATH for gnina
    print(f"\n3. Checking PATH for gnina:")
    import shutil
    path_binary = shutil.which('gnina')
    print(f"   Found in PATH: {path_binary}")
    
    if path_binary:
        print(f"   Testing binary...")
        try:
            result = subprocess.run([path_binary, '--version'], 
                                  capture_output=True, text=True, timeout=5)
            print(f"   Return code: {result.returncode}")
            print(f"   Stdout: {result.stdout[:100]}...")
            print(f"   Stderr: {result.stderr[:100]}...")
            
            if 'error while loading shared libraries' in result.stderr:
                print(f"   ❌ Binary has missing dependencies!")
            else:
                print(f"   ✅ Binary appears to work")
                
        except Exception as e:
            print(f"   ❌ Error testing binary: {e}")
    
    # Test the find_binary function
    print(f"\n4. Testing find_binary function:")
    try:
        from run_docking_multi import find_binary
        gnina_path = find_binary('gnina')
        print(f"   find_binary('gnina') returned: {gnina_path}")
        
        if gnina_path:
            print(f"   Testing returned binary...")
            try:
                result = subprocess.run([gnina_path, '--help'], 
                                      capture_output=True, text=True, timeout=5)
                print(f"   Return code: {result.returncode}")
                print(f"   Stdout: {result.stdout[:100]}...")
                print(f"   Stderr: {result.stderr[:100]}...")
                
            except Exception as e:
                print(f"   ❌ Error testing returned binary: {e}")
        
    except Exception as e:
        print(f"   ❌ Error importing find_binary: {e}")
        import traceback
        traceback.print_exc()
    
    # Check common paths
    print(f"\n5. Checking common Linux paths:")
    common_paths = [
        '/usr/local/bin',
        '/usr/bin',
        '/opt/conda/envs/docking/bin',
        os.path.expanduser('~/bin'),
    ]
    
    for common_path in common_paths:
        binary_path = os.path.join(common_path, 'gnina')
        print(f"   {common_path}/gnina: {os.path.exists(binary_path)}")
        if os.path.exists(binary_path):
            print(f"     Is executable: {os.access(binary_path, os.X_OK)}")
            try:
                result = subprocess.run([binary_path, '--version'], 
                                      capture_output=True, text=True, timeout=5)
                if 'error while loading shared libraries' in result.stderr:
                    print(f"     ❌ Has missing dependencies")
                else:
                    print(f"     ✅ Appears to work")
            except Exception as e:
                print(f"     ❌ Error testing: {e}")

def test_gnina_script_directly():
    """Test the gnina.sh script directly."""
    print("\n" + "=" * 60)
    print("TESTING GNINA SCRIPT DIRECTLY")
    print("=" * 60)
    
    gnina_sh_path = os.path.join(os.getcwd(), 'bin', 'gnina.sh')
    
    if not os.path.exists(gnina_sh_path):
        print("❌ bin/gnina.sh does not exist!")
        return False
    
    print(f"Testing: {gnina_sh_path}")
    
    try:
        # Test with --help
        result = subprocess.run([gnina_sh_path, '--help'], 
                              capture_output=True, text=True, timeout=5)
        print(f"Return code: {result.returncode}")
        print(f"Stdout: {result.stdout}")
        print(f"Stderr: {result.stderr}")
        
        if result.returncode == 0:
            print("✅ Script works correctly!")
            return True
        else:
            print("❌ Script failed!")
            return False
            
    except Exception as e:
        print(f"❌ Error running script: {e}")
        return False

def main():
    """Run all debug tests."""
    debug_gnina_detection()
    test_gnina_script_directly()
    
    print("\n" + "=" * 60)
    print("DEBUG SUMMARY")
    print("=" * 60)
    print("Check the output above to identify the issue with GNINA detection.")

if __name__ == "__main__":
    main() 