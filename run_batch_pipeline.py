#!/usr/bin/env python3
"""
Wrapper script to run batch molecular docking pipeline from the root directory.
This maintains compatibility with the new directory structure.
"""

import os
import sys
from pathlib import Path

# Add scripts directory to Python path
project_root = Path(__file__).parent
scripts_dir = project_root / "scripts"
sys.path.insert(0, str(scripts_dir))

# Import and run the batch pipeline
if __name__ == "__main__":
    from batch_pipeline import main
    main() 