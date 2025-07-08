# Docking Wrapper: Installation and Setup Guide

## Overview

This guide provides step-by-step instructions for installing and configuring the docking wrapper with all the recent bug fixes and improvements.

## Prerequisites

### System Requirements
- **Python**: 3.10 (recommended for best compatibility)
- **Operating System**: Windows, macOS, or Linux
- **Memory**: Minimum 4GB RAM (8GB+ recommended)
- **Storage**: At least 2GB free space

### Required External Tools
- **AutoDock Vina**: For molecular docking
- **GNINA**: For CNN-based docking (optional)
- **DiffDock**: For diffusion-based docking (optional)
- **MGLTools**: For structure preparation (optional, has fallback)

## Installation Steps

### 1. Clone the Repository
```bash
git clone <repository-url>
cd docking-wrapper
```

### 2. Install Python Dependencies

#### Option A: Using pip (Recommended)
```bash
# Install all required packages
pip install -r requirements_enhanced.txt

# Or install individually
pip install rdkit-pypi numpy pandas scikit-learn meeko colorama tabulate
```

#### Option B: Using conda
```bash
# Create new environment
conda create -n docking-wrapper python=3.10
conda activate docking-wrapper

# Install packages
conda install -c conda-forge rdkit numpy pandas scikit-learn
pip install meeko colorama tabulate
```

### 3. Install External Docking Tools

#### AutoDock Vina
```bash
# Ubuntu/Debian
sudo apt-get install autodock-vina

# macOS
brew install autodock-vina

# Windows
# Download from: http://vina.scripps.edu/download.html
# Add to PATH environment variable

# Verify installation
vina --version
```

#### GNINA (Optional)
```bash
# Download from: https://github.com/gnina/gnina
# Follow installation instructions in the repository

# Verify installation
gnina --version
```

#### DiffDock (Optional)
```bash
# Clone DiffDock repository
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock

# Install DiffDock dependencies
pip install -r requirements.txt

# Download pre-trained models
# Follow instructions in DiffDock README

# Verify installation
python inference.py --help
```

#### MGLTools (Optional)
```bash
# Download from: http://mgltools.scripps.edu/downloads/downloads/tools/downloads
# Follow platform-specific installation instructions

# Set environment variable (optional)
export MGLTOOLS_PATH="/path/to/mgltools"
```

### 4. Verify Installation

Run the dependency checker to verify all components are properly installed:

```bash
python scripts/check_dependencies.py
```

Expected output:
```
✅ Python Dependencies:
  - rdkit: OK
  - numpy: OK
  - pandas: OK
  - scikit-learn: OK
  - meeko: OK

✅ External Tools:
  - vina: OK (/usr/local/bin/vina)
  - gnina: OK (/usr/local/bin/gnina)
  - mgltools: OK (/path/to/mgltools)

✅ All dependencies satisfied!
```

## Configuration

### 1. Basic Configuration

The system automatically detects most tools, but you can customize the configuration:

```python
from config import DockingConfig

# Create configuration
config = DockingConfig()

# Set custom paths (if needed)
config.mgltools_path = "/custom/path/to/mgltools"
config.output_dir = "/custom/output/directory"

# Save configuration
config.save("my_config.json")
```

### 2. Environment Variables

You can also set environment variables for tool paths:

```bash
# Set MGLTools path
export MGLTOOLS_PATH="/path/to/mgltools"

# Set docking engine paths
export VINA_PATH="/path/to/vina"
export GNINA_PATH="/path/to/gnina"

# Set output directory
export DOCKING_OUTPUT_DIR="/path/to/output"
```

### 3. Configuration File

Create a configuration file `config.json`:

```json
{
    "mgltools_path": "/path/to/mgltools",
    "output_dir": "/path/to/output",
    "log_level": "INFO",
    "docking_engines": {
        "vina": {"enabled": true, "exhaustiveness": 8},
        "gnina": {"enabled": true, "use_gpu": false},
        "diffdock": {"enabled": false}
    }
}
```

## Quick Start

### 1. Test Installation

Run a simple test to verify everything works:

```bash
# Create test files
mkdir test_run
cd test_run

# Create simple protein PDB
echo "ATOM      1  N   ALA A   1      27.462  24.337   5.045  1.00 20.00           N  
ATOM      2  CA  ALA A   1      26.336  25.234   5.234  1.00 20.00           C  
END" > protein.pdb

# Create simple ligand SDF
echo "     RDKit          3D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5000    0.8660    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  2  3  1  0
M  END" > ligand.sdf

# Run test
python ../scripts/batch_pipeline.py \
    --protein protein.pdb \
    --ligands ligand.sdf \
    --engines vina \
    --output test_results
```

### 2. Basic Usage

```python
from scripts.batch_pipeline import BatchDockingPipeline

# Create pipeline
pipeline = BatchDockingPipeline(
    protein_file="protein.pdb",
    ligand_files=["ligand1.sdf", "ligand2.sdf"],
    docking_engines=["vina", "gnina"],
    output_dir="results"
)

# Run docking
results = pipeline.run()

# Check results
print(f"Successfully processed {results['successful']} ligands")
print(f"Failed: {results['failed']} ligands")
```

## Troubleshooting

### Common Issues

#### 1. MGLTools Not Found
```
Warning: MGLTools not found. Protein preparation may be limited.
```

**Solution**: 
- Install MGLTools or set `MGLTOOLS_PATH` environment variable
- The system will use fallback PDBQT conversion if MGLTools is not available

#### 2. Docking Engine Not Found
```
Missing required external binaries: vina
```

**Solution**:
- Install the missing docking engine
- Ensure it's in your PATH
- Or set the appropriate environment variable

#### 3. Python Dependencies Missing
```
ModuleNotFoundError: No module named 'rdkit'
```

**Solution**:
```bash
pip install rdkit-pypi
# or
conda install -c conda-forge rdkit
```

#### 4. PDBQT Conversion Issues
```
Error: Invalid atom type in PDBQT file
```

**Solution**: 
- This should be fixed with the recent updates
- Check that you're using the latest version
- Verify PDB file format is correct

#### 5. Permission Errors
```
Permission denied: /path/to/output
```

**Solution**:
- Check file permissions
- Use a different output directory
- Run with appropriate permissions

### Debug Mode

Enable debug logging for detailed troubleshooting:

```python
from utils.logging import setup_logging

logger = setup_logging(
    name="DebugRun",
    log_file="debug.log",
    level="DEBUG"
)
```

### Running Tests

Run the test suite to verify your installation:

```bash
# Run all tests
python run_tests.py

# Run specific test
python run_tests.py test_pdbqt_conversion
```

## Platform-Specific Notes

### Windows
- Use Windows Subsystem for Linux (WSL) for better compatibility
- Ensure all tools are in PATH environment variable
- Use forward slashes in paths or escape backslashes

### macOS
- Use Homebrew for installing external tools
- MGLTools may need special installation on newer macOS versions
- Consider using conda for Python package management

### Linux
- Use package manager for system dependencies
- Ensure proper permissions for output directories
- Consider using virtual environment for Python packages

## Performance Optimization

### 1. Parallel Processing
```python
# Enable parallel processing
pipeline = BatchDockingPipeline(
    # ... other parameters
    parallel=True,
    max_workers=4  # Adjust based on your system
)
```

### 2. GPU Acceleration
```python
# Enable GPU for GNINA
config = {
    "engines": {
        "gnina": {"enabled": True, "use_gpu": True}
    }
}
```

### 3. Memory Management
- Use smaller ligand batches for large datasets
- Monitor memory usage during long runs
- Consider using SSD storage for better I/O performance

## Support

### Getting Help
1. Check the troubleshooting section above
2. Review the test suite for examples
3. Check log files for detailed error information
4. Run dependency checker to verify setup

### Reporting Issues
When reporting issues, please include:
- Operating system and version
- Python version
- Complete error message
- Log file contents
- Steps to reproduce the issue

### Contributing
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## Updates and Maintenance

### Updating Dependencies
```bash
# Update Python packages
pip install --upgrade -r requirements_enhanced.txt

# Check for updates
python scripts/check_dependencies.py
```

### Backup Configuration
```bash
# Backup your configuration
cp config.json config_backup.json
```

### Clean Installation
```bash
# Remove old installation
rm -rf outputs/
rm -rf logs/
rm *.log

# Reinstall dependencies
pip install --force-reinstall -r requirements_enhanced.txt
``` 