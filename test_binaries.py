#!/usr/bin/env python3
"""
Test script to verify that GNINA and Vina binaries are found and working.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import subprocess
import platform
from pathlib import Path

def test_binary(binary_name, test_args=None):
    """Test if a binary is found and executable."""
    print(f"\nTesting {binary_name}...")
    
    # Check current directory first
    cwd = os.getcwd()
    local_path = os.path.join(cwd, binary_name)
    
    if os.path.isfile(local_path):
        print(f"✓ Found {binary_name} in current directory: {local_path}")
        binary_path = local_path
    else:
        # Check PATH
        import shutil
        which_path = shutil.which(binary_name)
        if which_path:
            print(f"✓ Found {binary_name} in PATH: {which_path}")
            binary_path = which_path
        else:
            print(f"✗ {binary_name} not found in current directory or PATH")
            return False
    
    # Test if binary is executable
    try:
        if test_args:
            cmd = [binary_path] + test_args
        else:
            cmd = [binary_path, '--version']
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print(f"✓ {binary_name} is executable and working")
            if result.stdout:
                print(f"  Output: {result.stdout.strip()}")
            return True
        else:
            print(f"✗ {binary_name} failed to execute")
            if result.stderr:
                print(f"  Error: {result.stderr.strip()}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"✗ {binary_name} timed out")
        return False
    except Exception as e:
        print(f"✗ {binary_name} error: {e}")
        return False

def main():
    """Test all binaries."""
    print("=" * 60)
    print("BINARY AVAILABILITY TEST")
    print("=" * 60)
    
    print(f"Current directory: {os.getcwd()}")
    print(f"Platform: {platform.system()} {platform.release()}")
    
    # Test Vina
    vina_ok = test_binary('vina.exe', ['--version'])
    
    # Test GNINA
    gnina_ok = test_binary('gnina', ['--version'])
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Vina: {'✓ Working' if vina_ok else '✗ Not working'}")
    print(f"GNINA: {'✓ Working' if gnina_ok else '✗ Not working'}")
    
    if vina_ok and gnina_ok:
        print("\n🎉 All binaries are working! Your pipeline is ready for production.")
    else:
        print("\n⚠️  Some binaries are not working. Please check the installation.")
    
    return vina_ok and gnina_ok

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 