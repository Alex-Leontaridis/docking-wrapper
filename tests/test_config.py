#!/usr/bin/env python3
"""
Unit tests for configuration system.
"""

import unittest
import tempfile
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add root directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import DockingConfig

class TestDockingConfig(unittest.TestCase):
    """Test DockingConfig functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.config = DockingConfig()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_config_paths(self):
        config = DockingConfig()
        # These should be None or a string path
        self.assertTrue(hasattr(config, 'mgltools_path'))
        self.assertTrue(hasattr(config, 'vina_path'))
        self.assertTrue(hasattr(config, 'gnina_path'))
        self.assertTrue(hasattr(config, 'diffdock_path'))

    def test_validate_mgltools(self):
        config = DockingConfig()
        # Should return True or False, not raise
        result = config.validate_mgltools()
        self.assertIn(result, [True, False])

    def test_get_config_dict(self):
        config = DockingConfig()
        config_dict = config.get_config_dict()
        self.assertIn('mgltools_path', config_dict)
        self.assertIn('vina_path', config_dict)
        self.assertIn('gnina_path', config_dict)
        self.assertIn('diffdock_path', config_dict)
        self.assertIn('output_dir', config_dict)
        self.assertIn('mgltools_available', config_dict)

if __name__ == '__main__':
    unittest.main() 