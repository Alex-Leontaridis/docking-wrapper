"""
Pytest configuration and fixtures for docking wrapper tests.
"""

import pytest
import tempfile
import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'scripts'))

@pytest.fixture(scope="session")
def test_data_dir():
    """Provide a temporary directory for test data."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield tmp_dir

@pytest.fixture
def sample_protein_pdb(test_data_dir):
    """Create a sample protein PDB file for testing."""
    protein_file = os.path.join(test_data_dir, 'sample_protein.pdb')
    with open(protein_file, 'w') as f:
        f.write("""ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
ATOM      3  C   ALA A   1      25.085  24.456   5.567  1.00 20.00           C  
ATOM      4  O   ALA A   1      24.000  25.000   5.000  1.00 20.00           O  
TER       5      ALA A   1                                                      
END""")
    return protein_file

@pytest.fixture
def sample_ligand_sdf(test_data_dir):
    """Create a sample ligand SDF file for testing."""
    ligand_file = os.path.join(test_data_dir, 'sample_ligand.sdf')
    with open(ligand_file, 'w') as f:
        f.write("""     RDKit          3D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  2  3  1  0
M  END""")
    return ligand_file

@pytest.fixture
def mock_vina_output():
    """Mock Vina docking output."""
    return """-----+------------+----------+----------+
   # |   MODE    |   AFFINITY | RMSD LOWER | RMSD UPPER
-----+------------+----------+----------+----------+
    1 |      1    |    -8.1   |    0.000   |    0.000
    2 |      2    |    -7.5   |    1.200   |    1.800
    3 |      3    |    -7.2   |    2.100   |    2.900
-----+------------+----------+----------+----------+
"""

@pytest.fixture
def mock_gnina_output():
    """Mock GNINA docking output."""
    return """-----+------------+----------+----------+
   # |   MODE    |   AFFINITY | RMSD LOWER | RMSD UPPER
-----+------------+----------+----------+----------+
    1 |      1    |    -7.8   |    0.000   |    0.000
    2 |      2    |    -7.3   |    1.500   |    2.100
-----+------------+----------+----------+----------+
"""

@pytest.fixture
def mock_diffdock_output():
    """Mock DiffDock confidence output."""
    return """pose_1.sdf,0.85
pose_2.sdf,0.72
pose_3.sdf,0.68"""

@pytest.fixture
def sample_config():
    """Sample configuration for testing."""
    return {
        'protein_file': 'test_protein.pdb',
        'ligand_files': ['test_ligand.sdf'],
        'docking_engines': ['vina', 'gnina'],
        'output_dir': 'test_output',
        'log_level': 'INFO',
        'mgltools_path': '/mock/mgltools',
        'docking_engines': {
            'vina': '/mock/vina',
            'gnina': '/mock/gnina',
            'diffdock': '/mock/diffdock'
        }
    }

@pytest.fixture
def mock_docking_tools(test_data_dir):
    """Create mock docking tools for testing."""
    bin_dir = os.path.join(test_data_dir, 'bin')
    ensure_dir(bin_dir)
    
    # Mock Vina
    vina_path = os.path.join(bin_dir, 'vina')
    with open(vina_path, 'w') as f:
        f.write("""#!/bin/bash
echo "Mock Vina output"
echo "-----+------------+----------+----------+"
echo "   # |   MODE    |   AFFINITY | RMSD LOWER | RMSD UPPER"
echo "-----+------------+----------+----------+----------+"
echo "    1 |      1    |    -8.1   |    0.000   |    0.000"
echo "-----+------------+----------+----------+----------+"
""")
    os.chmod(vina_path, 0o755)
    
    # Mock GNINA
    gnina_path = os.path.join(bin_dir, 'gnina')
    with open(gnina_path, 'w') as f:
        f.write("""#!/bin/bash
echo "Mock GNINA output"
echo "-----+------------+----------+----------+"
echo "   # |   MODE    |   AFFINITY | RMSD LOWER | RMSD UPPER"
echo "-----+------------+----------+----------+----------+"
echo "    1 |      1    |    -7.5   |    0.000   |    0.000"
echo "-----+------------+----------+----------+----------+"
""")
    os.chmod(gnina_path, 0o755)
    
    return bin_dir

@pytest.fixture
def mock_subprocess():
    """Mock subprocess for testing external tool calls."""
    with patch('subprocess.run') as mock_run:
        mock_run.return_value.returncode = 0
        mock_run.return_value.stdout = b"Mock output"
        mock_run.return_value.stderr = b""
        yield mock_run

@pytest.fixture(autouse=True)
def setup_test_environment():
    """Set up test environment variables."""
    # Set test environment variables
    os.environ['TESTING'] = 'true'
    os.environ['LOG_LEVEL'] = 'DEBUG'
    
    # Mock external tool paths
    os.environ['VINA_PATH'] = '/mock/vina'
    os.environ['GNINA_PATH'] = '/mock/gnina'
    os.environ['DIFFDOCK_PATH'] = '/mock/diffdock'
    
    yield
    
    # Clean up
    for key in ['TESTING', 'LOG_LEVEL', 'VINA_PATH', 'GNINA_PATH', 'DIFFDOCK_PATH']:
        os.environ.pop(key, None)

@pytest.fixture
def sample_pdbqt_file(test_data_dir):
    """Create a sample PDBQT file for testing."""
    pdbqt_file = os.path.join(test_data_dir, 'sample_protein.pdbqt')
    with open(pdbqt_file, 'w') as f:
        f.write("""ROOT
ATOM      1  C   ALA A   1      27.462  24.337   5.045  0.00  0.00    +0.000 C 
ATOM      2  C   ALA A   1      26.336  25.234   5.234  0.00  0.00    +0.000 C 
ATOM      3  C   ALA A   1      25.085  24.456   5.567  0.00  0.00    +0.000 C 
ATOM      4  O   ALA A   1      24.000  25.000   5.000  0.00  0.00    -0.385 O 
ENDROOT
TORSDOF 0
""")
    return pdbqt_file

# Markers for different test types
def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line("markers", "slow: marks tests as slow")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "unit: marks tests as unit tests")
    config.addinivalue_line("markers", "windows: marks tests as windows-specific")
    config.addinivalue_line("markers", "linux: marks tests as linux-specific")
    config.addinivalue_line("markers", "macos: marks tests as macos-specific")
    config.addinivalue_line("markers", "docking: marks tests as docking-specific")
    config.addinivalue_line("markers", "parser: marks tests as parser-specific") 