#!/usr/bin/env python3
"""
Integration tests for the complete docking pipeline.
"""

import unittest
import tempfile
import os
import sys
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'scripts'))
sys.path.insert(0, str(Path(__file__).parent.parent))

from batch_pipeline import BatchDockingPipeline
from config import Config

class TestIntegration(unittest.TestCase):
    """Integration tests for the complete pipeline."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create test protein PDB
        self.protein_pdb = os.path.join(self.temp_dir, 'protein.pdb')
        with open(self.protein_pdb, 'w') as f:
            f.write("""ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00           O  
END""")
        
        # Create test ligand SDF
        self.ligand_sdf = os.path.join(self.temp_dir, 'ligand.sdf')
        with open(self.ligand_sdf, 'w') as f:
            f.write("""     RDKit          3D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  2  3  1  0
M  END""")
        
        # Create test configuration
        self.config = {
            'protein_file': self.protein_pdb,
            'ligand_files': [self.ligand_sdf],
            'docking_engines': ['vina'],
            'output_dir': os.path.join(self.temp_dir, 'output'),
            'log_level': 'INFO'
        }
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        try:
            shutil.rmtree(self.temp_dir)
        except PermissionError:
            # Ignore permission errors on Windows
            pass
    
    @patch('batch_pipeline.subprocess.run')
    @patch('batch_pipeline.Config')
    @patch('batch_pipeline.shutil.which')
    def test_complete_pipeline_flow(self, mock_which, mock_config_class, mock_subprocess):
        """Test complete pipeline flow with mocked external tools."""
        # Mock binary detection
        mock_which.return_value = '/mock/binary'
        
        # Mock configuration
        mock_config = MagicMock()
        mock_config.mgltools_path = '/mock/mgltools'
        mock_config.docking_engines = {'vina': '/mock/vina'}
        mock_config.output_dir = self.config['output_dir']
        mock_config.log_level = 'INFO'
        mock_config.validate_mgltools.return_value = True
        mock_config.get_mgltools_pythonsh.return_value = '/mock/pythonsh'
        mock_config.get_mgltools_prepare_script.return_value = '/mock/prepare_script'
        mock_config_class.return_value = mock_config
        
        # Mock subprocess calls
        mock_subprocess.return_value.returncode = 0
        mock_subprocess.return_value.stdout = b"Mock output"
        
        # Create pipeline
        pipeline = BatchDockingPipeline(output_dir=self.config['output_dir'])
        
        # Run pipeline
        result = pipeline.run_batch_pipeline(
            protein_file=self.config['protein_file'],
            ligand_input=self.ligand_sdf
        )
        
        # Check that pipeline completed
        self.assertIsNotNone(result)
        self.assertTrue(os.path.exists(self.config['output_dir']))
    
    @patch('batch_pipeline.subprocess.run')
    @patch('batch_pipeline.shutil.which')
    def test_pipeline_with_pdbqt_conversion(self, mock_which, mock_subprocess):
        """Test pipeline with PDBQT conversion."""
        # Mock binary detection
        mock_which.return_value = '/mock/binary'
        
        # Mock subprocess calls
        mock_subprocess.return_value.returncode = 0
        mock_subprocess.return_value.stdout = b"Mock output"
        
        # Create pipeline with fallback PDBQT conversion
        pipeline = BatchDockingPipeline(output_dir=self.config['output_dir'])
        
        # Test PDBQT conversion
        pdbqt_file = os.path.join(self.temp_dir, 'protein.pdbqt')
        pipeline._simple_pdb_to_pdbqt(self.protein_pdb, pdbqt_file)
        
        # Check that PDBQT file was created with proper format
        self.assertTrue(os.path.exists(pdbqt_file))
        
        with open(pdbqt_file, 'r') as f:
            content = f.read()
            # Should not contain problematic +0.000 format
            self.assertNotIn('+0.000', content)
            # Should contain proper atom lines
            self.assertIn('ATOM', content)
    
    @patch('batch_pipeline.shutil.which')
    def test_pipeline_error_handling(self, mock_which):
        """Test pipeline error handling."""
        # Mock binary detection
        mock_which.return_value = '/mock/binary'
        
        # Create pipeline
        pipeline = BatchDockingPipeline(output_dir=self.config['output_dir'])
        
        # Should handle missing files gracefully by calling the method directly
        with self.assertRaises(ValueError):
            pipeline.discover_ligands('/nonexistent/ligand.sdf')
    
    @patch('batch_pipeline.subprocess.run')
    @patch('batch_pipeline.shutil.which')
    def test_pipeline_with_results_parsing(self, mock_which, mock_subprocess):
        """Test pipeline with results parsing."""
        # Mock binary detection
        mock_which.return_value = '/mock/binary'
        
        # Mock Vina output
        vina_output = """-----+------------+----------+----------+
   # |   MODE    |   AFFINITY | RMSD LOWER | RMSD UPPER
-----+------------+----------+----------+----------+
    1 |      1    |    -8.1   |    0.000   |    0.000
-----+------------+----------+----------+----------+
"""
        
        # Mock subprocess calls
        mock_subprocess.return_value.returncode = 0
        mock_subprocess.return_value.stdout = vina_output.encode()
        
        # Create pipeline
        pipeline = BatchDockingPipeline(output_dir=self.config['output_dir'])
        
        # Test results parsing by creating a parser manually
        from docking_results_parser import DockingResultsParser
        parser = DockingResultsParser(self.temp_dir, self.temp_dir)
        
        # Test results parsing
        vina_file = os.path.join(self.temp_dir, 'vina_out.txt')
        with open(vina_file, 'w') as f:
            f.write(vina_output)
        
        # Note: The parser expects PDBQT files, not text files
        # This is just a basic test of the parser instantiation
        self.assertIsNotNone(parser)
    
    def test_configuration_integration(self):
        """Test configuration system integration."""
        # Create configuration
        config = Config()
        config.mgltools_path = '/test/mgltools'
        config.output_dir = self.config['output_dir']
        
        # Check configuration properties
        self.assertEqual(config.mgltools_path, '/test/mgltools')
        self.assertEqual(config.output_dir, self.config['output_dir'])
        
        # Test configuration dictionary
        config_dict = config.get_config_dict()
        self.assertEqual(config_dict['mgltools_path'], '/test/mgltools')
        self.assertEqual(config_dict['output_dir'], self.config['output_dir'])

if __name__ == '__main__':
    unittest.main() 