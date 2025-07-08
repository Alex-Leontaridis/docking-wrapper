#!/usr/bin/env python3
"""
Unit tests for PDBQT conversion fixes.
"""

import unittest
import tempfile
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'scripts'))

from prep_structures import _simple_pdb_to_pdbqt
from batch_pipeline import BatchDockingPipeline

class TestPDBQTConversion(unittest.TestCase):
    """Test PDBQT conversion functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a simple test PDB file
        self.test_pdb_content = """ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00           O  
HETATM    5  C   LIG A   2      20.000  20.000  20.000  1.00 20.00           C  
HETATM    6  O   LIG A   2      21.000  21.000  21.000  1.00 20.00           O  
END"""
        
        self.test_pdb_file = os.path.join(self.temp_dir, 'test.pdb')
        with open(self.test_pdb_file, 'w') as f:
            f.write(self.test_pdb_content)
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        try:
            shutil.rmtree(self.temp_dir)
        except PermissionError:
            # Ignore permission errors on Windows
            pass
    
    def test_simple_pdb_to_pdbqt_format(self):
        """Test that PDBQT conversion uses proper AutoDock format."""
        output_file = os.path.join(self.temp_dir, 'test.pdbqt')
        
        # Run conversion
        _simple_pdb_to_pdbqt(self.test_pdb_file, output_file)
        
        # Check that file was created
        self.assertTrue(os.path.exists(output_file))
        
        # Read the output and check format
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        # Check that we have atom lines (RDKit may or may not add hydrogens)
        atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
        self.assertGreaterEqual(len(atom_lines), 6)  # At least the original 6 atoms
        
        # Check that charge format is correct (not +0.000)
        for line in atom_lines:
            # Extract charge part (columns 66-76)
            if len(line) >= 76:
                charge_part = line[66:76].strip()
                # Should not contain the problematic +0.000 format
                self.assertNotIn('+0.000', charge_part)
                # Should be a proper number format
                self.assertTrue(charge_part.replace('.', '').replace('-', '').isdigit())
    
    @unittest.mock.patch('batch_pipeline.shutil.which')
    def test_batch_pipeline_pdbqt_conversion(self, mock_which):
        """Test PDBQT conversion in batch pipeline."""
        # Mock binary detection
        mock_which.return_value = '/mock/binary'
        
        pipeline = BatchDockingPipeline(output_dir=self.temp_dir)
        output_file = os.path.join(self.temp_dir, 'test_batch.pdbqt')
        
        # Run conversion using batch pipeline method
        pipeline._simple_pdb_to_pdbqt(self.test_pdb_file, output_file)
        
        # Check that file was created
        self.assertTrue(os.path.exists(output_file))
        
        # Read and verify format
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
        self.assertGreaterEqual(len(atom_lines), 6)  # At least the original 6 atoms
        
        # Check charge format
        for line in atom_lines:
            if len(line) >= 76:
                charge_part = line[66:76].strip()
                self.assertNotIn('+0.000', charge_part)
    
    def test_pdbqt_conversion_with_different_atom_types(self):
        """Test PDBQT conversion handles different atom types correctly."""
        # Create PDB with various atom types
        pdb_content = """ATOM      1  C   ALA A   1      27.462  24.337   5.045  1.00 20.00           C  
ATOM      2  N   ALA A   1      26.336  25.234   5.234  1.00 20.00           N  
ATOM      3  O   ALA A   1      25.085  24.456   5.567  1.00 20.00           O  
ATOM      4  S   ALA A   1      24.000  25.000   5.000  1.00 20.00           S  
ATOM      5  P   ALA A   1      23.000  26.000   6.000  1.00 20.00           P  
HETATM    6  H   LIG A   2      20.000  20.000  20.000  1.00 20.00           H  
END"""
        
        test_file = os.path.join(self.temp_dir, 'atoms.pdb')
        with open(test_file, 'w') as f:
            f.write(pdb_content)
        
        output_file = os.path.join(self.temp_dir, 'atoms.pdbqt')
        _simple_pdb_to_pdbqt(test_file, output_file)
        
        # Check that all atom types are properly converted
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
        self.assertGreaterEqual(len(atom_lines), 6)  # At least the original 6 atoms
        
        # Check that each line has proper AutoDock atom type at the end
        for line in atom_lines:
            if len(line) >= 78:
                autodock_type = line[77:79].strip()
                self.assertIn(autodock_type, ['C', 'N', 'O', 'S', 'P', 'H'])
    
    def test_pdbqt_conversion_handles_malformed_input(self):
        """Test PDBQT conversion handles malformed input gracefully."""
        # Create malformed PDB
        malformed_pdb = """ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
# This is a comment line
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00           O  
END"""
        
        test_file = os.path.join(self.temp_dir, 'malformed.pdb')
        with open(test_file, 'w') as f:
            f.write(malformed_pdb)
        
        output_file = os.path.join(self.temp_dir, 'malformed.pdbqt')
        
        # Should not raise an exception
        _simple_pdb_to_pdbqt(test_file, output_file)
        
        # Should create output file
        self.assertTrue(os.path.exists(output_file))
        
        # Should have processed valid atom lines (RDKit may or may not add hydrogens)
        with open(output_file, 'r') as f:
            lines = f.readlines()
        
        atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
        self.assertGreaterEqual(len(atom_lines), 4)  # At least the original 4 atoms

if __name__ == '__main__':
    unittest.main() 