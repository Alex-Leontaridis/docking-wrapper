# EquiBind Production-Ready Script

This document describes the production-ready implementation of the EquiBind molecular docking script, based on the [EquiBind GitHub repository](https://github.com/HannesStark/EquiBind).

## 🚀 Production-Ready Features

### ✅ **Critical Issues Fixed**

1. **✅ Missing Configuration File**
   - Created `tools_config.json` with comprehensive settings
   - Fallback defaults if config file is missing
   - Configurable retry logic, timeouts, and resource limits

2. **✅ Input File Validation**
   - PDB file format validation (ATOM/HETATM records, TITLE)
   - SDF file format validation (molecule count, $$$$ delimiter)
   - MOL2, PDBQT, and PDB ligand format support
   - File size limits (configurable, default 100MB)
   - Content integrity checks

3. **✅ Security Hardening**
   - Path validation (restricts file access to allowed directories)
   - Input sanitization and validation
   - File size limits to prevent DoS attacks
   - Secure temporary directory usage

4. **✅ Error Recovery & Retry Logic**
   - Configurable retry attempts (default: 3)
   - Timeout handling (default: 300 seconds)
   - Graceful degradation on failures
   - Comprehensive error messages and logging

5. **✅ Resource Management**
   - Memory usage monitoring (configurable limit: 8GB)
   - CPU usage monitoring (configurable limit: 80%)
   - Process cleanup on termination
   - Resource availability checks before execution

6. **✅ Comprehensive Logging**
   - Structured logging with configurable levels
   - File and console output
   - Detailed error tracking
   - Performance metrics logging

7. **✅ Health Checks**
   - System resource monitoring
   - Dependency validation
   - Installation verification
   - Periodic health status reporting

8. **✅ Graceful Shutdown**
   - Signal handling (SIGINT, SIGTERM)
   - Process termination with cleanup
   - Resource cleanup on interruption
   - State preservation where possible

## 📋 Configuration

### `tools_config.json`

```json
{
  "tools": {
    "equibind": {
      "repo_path": "EquiBind",
      "conda_env": "equibind",
      "max_retries": 3,
      "timeout_seconds": 300,
      "memory_limit_gb": 8,
      "cpu_limit_percent": 80,
      "log_level": "INFO",
      "health_check_interval": 30,
      "output_formats": [".sdf", ".pdb"],
      "supported_ligand_formats": [".mol2", ".sdf", ".pdbqt", ".pdb"],
      "supported_protein_formats": [".pdb"]
    }
  },
  "logging": {
    "level": "INFO",
    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    "file": "logs/equibind.log"
  },
  "security": {
    "max_file_size_mb": 100,
    "allowed_paths": ["inputs/", "outputs/", "temp/"],
    "validate_file_content": true
  }
}
```

## 🛠️ Usage

### Basic Usage

```bash
python scripts/run_equibind.py \
  --protein inputs/protein.pdb \
  --ligand inputs/ligands/aspirin.sdf \
  --output outputs/equibind_results/
```

### Advanced Options

```bash
python scripts/run_equibind.py \
  --protein inputs/protein.pdb \
  --ligand inputs/ligands/ligand.mol2 \
  --output outputs/equibind_results/ \
  --conda_env equibind \
  --skip-validation
```

### Command Line Arguments

- `--protein`: Protein PDB file (required)
- `--ligand`: Ligand file in SDF, MOL2, PDBQT, or PDB format (required)
- `--output`: Output directory (required)
- `--conda_env`: Conda environment name (default: equibind)
- `--skip-validation`: Skip input file validation (use with caution)

## 🔧 Installation Requirements

### Dependencies

Based on the [EquiBind documentation](https://github.com/HannesStark/EquiBind), the following dependencies are required:

```bash
# Core dependencies
python=3.7
pytorch>=1.10
torchvision
cudatoolkit=10.2
torchaudio
dgl-cuda10.2
rdkit
openbabel
biopython
rdkit
biopandas
pot
dgllife
joblib
pyaml
icecream
matplotlib
tensorboard

# Production dependencies
psutil>=5.9.0
```

### Environment Setup

```bash
# Clone EquiBind repository
git clone https://github.com/HannesStark/EquiBind

# Create conda environment
conda env create -f EquiBind/environment.yml

# Activate environment
conda activate equibind

# Install production dependencies
pip install psutil>=5.9.0
```

## 🧪 Testing

### Run Production Tests

```bash
python test_equibind_production.py
```

### Test Coverage

The test suite validates:

- ✅ Configuration file validation
- ✅ Input file validation (PDB, SDF, MOL2, PDBQT)
- ✅ Security features (path validation, file size limits)
- ✅ Logging setup and functionality
- ✅ Error handling and graceful failure
- ✅ Skip validation flag functionality
- ✅ Help output completeness

## 📊 Monitoring & Logging

### Log Files

- **Location**: `logs/equibind.log`
- **Format**: Structured logging with timestamps
- **Levels**: DEBUG, INFO, WARNING, ERROR, CRITICAL

### Example Log Output

```
2024-01-15 10:30:15 - EquiBind - INFO - Starting EquiBind pose prediction
2024-01-15 10:30:15 - EquiBind - INFO - Configuration loaded successfully
2024-01-15 10:30:16 - EquiBind - INFO - Input files validated successfully
2024-01-15 10:30:17 - EquiBind - INFO - EquiBind installation validated successfully
2024-01-15 10:30:18 - EquiBind - INFO - System resources OK - Memory: 12.5GB, CPU: 45.2%
2024-01-15 10:30:19 - EquiBind - INFO - EquiBind inference attempt 1/3
2024-01-15 10:32:45 - EquiBind - INFO - EquiBind completed successfully with 2 output files
2024-01-15 10:32:45 - EquiBind - INFO - EquiBind completed successfully in 150.23s
```

## 🔒 Security Features

### Path Validation

The script restricts file access to predefined allowed paths:
- `inputs/`
- `outputs/`
- `temp/`

### File Size Limits

- **Default**: 100MB per file
- **Configurable**: Via `tools_config.json`
- **Purpose**: Prevent DoS attacks and resource exhaustion

### Input Sanitization

- File format validation
- Content integrity checks
- Malicious file detection
- Safe temporary directory usage

## ⚡ Performance Features

### Resource Monitoring

- **Memory**: Monitors available RAM (default limit: 8GB)
- **CPU**: Monitors CPU usage (default limit: 80%)
- **Disk**: Validates available disk space
- **Network**: Timeout handling for network operations

### Retry Logic

- **Default**: 3 retry attempts
- **Configurable**: Via `tools_config.json`
- **Backoff**: Exponential backoff between retries
- **Timeout**: Configurable per attempt (default: 300s)

## 🚨 Error Handling

### Graceful Degradation

- Comprehensive error messages
- Detailed logging of failures
- Resource cleanup on errors
- State preservation where possible

### Signal Handling

- **SIGINT**: Graceful shutdown on Ctrl+C
- **SIGTERM**: Clean termination on system signals
- **Process Cleanup**: Automatic cleanup of child processes

## 📈 Production Deployment

### Recommended Settings

For production deployment, consider these settings:

```json
{
  "tools": {
    "equibind": {
      "max_retries": 5,
      "timeout_seconds": 600,
      "memory_limit_gb": 16,
      "cpu_limit_percent": 90,
      "log_level": "WARNING"
    }
  },
  "security": {
    "max_file_size_mb": 500,
    "allowed_paths": ["/data/inputs/", "/data/outputs/", "/tmp/"]
  }
}
```

### Monitoring Integration

The script provides structured output suitable for monitoring systems:

- Exit codes for automation
- Structured logging for log aggregation
- Performance metrics for resource monitoring
- Health check endpoints for load balancers

## 🔄 Integration with Existing Pipeline

The production-ready script integrates seamlessly with the existing docking pipeline:

```python
# Example integration
from scripts.run_equibind import run_equibind_inference, load_config

config = load_config()
success, error = run_equibind_inference(
    protein_path="protein.pdb",
    ligand_path="ligand.sdf", 
    output_path="results/",
    config=config
)

if success:
    print("EquiBind completed successfully")
else:
    print(f"EquiBind failed: {error}")
```

## 📚 References

- [EquiBind GitHub Repository](https://github.com/HannesStark/EquiBind)
- [EquiBind Paper](https://arxiv.org/abs/2202.05146)
- [EquiBind Documentation](https://github.com/HannesStark/EquiBind#readme)

## 🤝 Contributing

When contributing to the production-ready EquiBind script:

1. Follow the existing code style and patterns
2. Add comprehensive tests for new features
3. Update documentation for any changes
4. Ensure backward compatibility
5. Test with various input formats and edge cases

## 📄 License

This production-ready implementation follows the same MIT license as the original EquiBind project. 