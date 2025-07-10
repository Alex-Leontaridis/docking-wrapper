#!/usr/bin/env python3
"""
Test script to verify the improved binary detection
"""

import sys
import os
sys.path.append('scripts')

from run_docking_multi import find_binary

def test_detection():
    print("Testing improved binary detection...")
    print("=" * 50)
    
    # Test vina detection
    print("Testing VINA detection:")
    vina_path = find_binary('vina')
    print(f"find_binary('vina'): {vina_path}")
    
    if vina_path:
        print(f"✅ VINA found at: {vina_path}")
        # Test if it's executable
        if os.access(vina_path, os.X_OK):
            print("✅ VINA is executable")
        else:
            print("❌ VINA is not executable")
    else:
        print("❌ VINA not found")
    
    print()
    
    # Test gnina detection
    print("Testing GNINA detection:")
    gnina_path = find_binary('gnina')
    print(f"find_binary('gnina'): {gnina_path}")
    
    if gnina_path:
        print(f"✅ GNINA found at: {gnina_path}")
        # Test if it's executable
        if os.access(gnina_path, os.X_OK):
            print("✅ GNINA is executable")
        else:
            print("❌ GNINA is not executable")
    else:
        print("❌ GNINA not found")
    
    print()
    
    # Test with environment variables
    print("Testing with environment variables:")
    os.environ['VINA_PATH'] = '/usr/local/bin/vina'
    os.environ['GNINA_PATH'] = '/usr/local/bin/gnina'
    
    vina_env = find_binary('vina', env_var='VINA_PATH')
    gnina_env = find_binary('gnina', env_var='GNINA_PATH')
    
    print(f"VINA via env var: {vina_env}")
    print(f"GNINA via env var: {gnina_env}")
    
    print()
    print("=" * 50)
    print("Detection test completed!")

if __name__ == "__main__":
    test_detection() 