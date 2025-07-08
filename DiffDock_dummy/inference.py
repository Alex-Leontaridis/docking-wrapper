#!/usr/bin/env python3
"""
DiffDock Dummy Script

This is a placeholder script that runs when DiffDock is not properly installed.
It provides clear instructions for installing DiffDock.
"""

import sys
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir

def main():
    print("=" * 60)
    print("DIFFDOCK NOT INSTALLED")
    print("=" * 60)
    print()
    print("This is a dummy script. DiffDock is not properly installed.")
    print()
    print("To install DiffDock:")
    print("1. Clone the DiffDock repository:")
    print("   git clone https://github.com/gcorso/DiffDock.git")
    print()
    print("2. Follow the installation instructions in the DiffDock README:")
    print("   cd DiffDock")
    print("   pip install -e .")
    print()
    print("3. Set the DIFFDOCK_PATH environment variable:")
    print("   export DIFFDOCK_PATH=/path/to/DiffDock")
    print()
    print("4. Or update the configuration to point to your DiffDock installation.")
    print()
    print("For more information, visit:")
    print("https://github.com/gcorso/DiffDock")
    print()
    print("=" * 60)
    
    # Exit with error code
    sys.exit(1)

if __name__ == "__main__":
    main()
