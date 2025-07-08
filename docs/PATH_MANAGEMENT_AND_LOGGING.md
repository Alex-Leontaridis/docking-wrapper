# Path Management and Logging System

## Overview

This document describes the centralized path management and consistent logging system implemented to eliminate hardcoded paths and standardize logging across the docking wrapper project.

## Path Management System

### Features

- **Platform-agnostic**: Works on Windows, macOS, and Linux
- **Environment variable support**: Configurable via environment variables
- **Automatic detection**: Finds tools and directories automatically
- **Configuration file**: Persistent path configuration
- **Fallback mechanisms**: Graceful handling when tools are not found

### Core Components

#### PathManager Class

The `PathManager` class provides centralized path management:

```python
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir

# Get the global path manager
path_manager = get_path_manager()

# Get a specific path
vina_path = get_path("vina")
output_dir = get_path("outputs")

# Get absolute path
abs_path = get_absolute_path("scripts")

# Ensure directory exists
log_dir = ensure_dir("logs")
```

#### Environment Variables

The system supports the following environment variables:

- `DOCKING_OUTPUT_DIR`: Output directory for results
- `DOCKING_LOG_DIR`: Directory for log files
- `DOCKING_TEMP_DIR`: Temporary files directory
- `DOCKING_CACHE_DIR`: Cache directory
- `VINA_PATH`: Path to Vina executable
- `GNINA_PATH`: Path to GNINA executable
- `DIFFDOCK_PATH`: Path to DiffDock installation
- `MGLTOOLS_PATH`: Path to MGLTools installation

#### Configuration File

Paths are stored in `config/paths.json`:

```json
{
  "project_root": "/path/to/project",
  "scripts": "/path/to/project/scripts",
  "outputs": "/path/to/outputs",
  "vina": "/usr/local/bin/vina",
  "gnina": "/usr/local/bin/gnina",
  "models": {
    "equibind": "/path/to/EquiBind",
    "neuralplexer": "/path/to/NeuralPLexer"
  }
}
```

### Usage Examples

#### Basic Path Access

```python
from utils.path_manager import get_path, get_absolute_path

# Get tool paths
vina_path = get_path("vina")
gnina_path = get_path("gnina")

# Get model paths
equibind_path = get_path("models", "equibind")

# Get absolute paths
scripts_dir = get_absolute_path("scripts")
```

#### Directory Management

```python
from utils.path_manager import ensure_dir, get_output_path

# Ensure directories exist
log_dir = ensure_dir("logs")
temp_dir = ensure_dir("temp")

# Get output paths
output_file = get_output_path("docking_results", "ligand1", "vina_out.pdbqt")
log_file = get_log_path("batch_pipeline.log")
```

#### Platform Detection

```python
from utils.path_manager import get_path_manager

path_manager = get_path_manager()
platform_info = path_manager.get_platform_info()
tools = path_manager.list_available_tools()

print(f"Platform: {platform_info['system']}")
print(f"Available tools: {tools}")
```

## Logging System

### Features

- **Consistent formatting**: Standardized log formats across all modules
- **Error codes**: Categorized error reporting with numeric codes
- **Performance tracking**: Built-in timing and performance metrics
- **Multiple handlers**: Console and file logging with rotation
- **Third-party noise reduction**: Suppresses noisy library logs

### Core Components

#### DockingLogger Class

The `DockingLogger` class provides enhanced logging capabilities:

```python
from utils.logging import setup_logging, get_logger

# Setup logger with automatic file creation
logger = setup_logging('MyModule')

# Get existing logger
logger = get_logger('MyModule')

# Log messages with error codes
logger.info("Processing started")
logger.warning("Configuration file not found")
logger.error("Docking failed", error_code='VINA_FAILED')
logger.critical("System error", error_code='SYSTEM_ERROR')
```

#### Error Codes

The system defines standardized error codes:

```python
# File/IO errors (1000-1999)
FILE_NOT_FOUND = 1001
FILE_PERMISSION_DENIED = 1002

# Docking engine errors (2000-2999)
VINA_FAILED = 2001
GNINA_FAILED = 2002
ENGINE_NOT_FOUND = 2004

# Structure preparation errors (3000-3999)
PROTEIN_PREP_FAILED = 3001
LIGAND_PREP_FAILED = 3002

# System errors (9000-9999)
SYSTEM_ERROR = 9001
MEMORY_ERROR = 9002
```

#### Performance Tracking

```python
# Time operations
with logger.timer("docking_operation"):
    result = run_docking(protein, ligand)

# Log performance metrics
logger.performance("ligand_preparation", 2.5)

# Get performance statistics
stats = logger.get_performance_stats("docking_operation")
print(f"Average time: {stats['mean']:.3f}s")
```

#### Convenience Functions

```python
from utils.logging import log_startup, log_shutdown, log_error_with_context

# Log module startup
log_startup('BatchPipeline', '1.0.0')

# Log module shutdown with performance summary
log_shutdown('BatchPipeline')

# Log errors with context
try:
    result = risky_operation()
except Exception as e:
    log_error_with_context(e, "During docking", error_code='DOCKING_FAILED')
```

### Log Formats

#### Console Output (Simple)
```
INFO: Processing ligand aspirin
WARNING: MGLTools not found, using fallback
ERROR: [VINA_FAILED] Docking engine failed
```

#### File Output (Detailed)
```
2024-01-15 10:30:45 | INFO     | BatchPipeline        | batch_pipeline.py:123 | Processing ligand aspirin
2024-01-15 10:30:46 | WARNING  | BatchPipeline        | batch_pipeline.py:124 | MGLTools not found, using fallback
2024-01-15 10:30:47 | ERROR    | BatchPipeline        | batch_pipeline.py:125 | [VINA_FAILED] Docking engine failed
```

## Migration Guide

### Updating Existing Code

#### Before (Old Style)
```python
import os
import logging

# Hardcoded paths
output_dir = os.path.join(os.getcwd(), 'outputs')
log_file = os.path.join(output_dir, 'log.txt')

# Basic logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.info("Processing started")
```

#### After (New Style)
```python
from utils.path_manager import get_path, get_log_path
from utils.logging import setup_logging

# Path manager
output_dir = get_path("outputs")
log_file = get_log_path("processing.log")

# Enhanced logging
logger = setup_logging(__name__, log_file)
logger.info("Processing started")
```

### Automatic Migration

Use the provided migration script:

```bash
python scripts/update_paths_and_logging.py
```

This script will automatically:
- Add necessary imports
- Replace hardcoded paths with path manager calls
- Update logging calls to use the new system
- Preserve existing functionality

## Configuration

### Environment Variables

Set these environment variables to customize behavior:

```bash
# Output directories
export DOCKING_OUTPUT_DIR="/custom/outputs"
export DOCKING_LOG_DIR="/custom/logs"

# Tool paths
export VINA_PATH="/usr/local/bin/vina"
export GNINA_PATH="/usr/local/bin/gnina"
export MGLTOOLS_PATH="/opt/mgltools"
```

### Configuration File

Edit `config/paths.json` to customize paths:

```json
{
  "outputs": "/custom/outputs",
  "logs": "/custom/logs",
  "vina": "/custom/path/to/vina",
  "models": {
    "equibind": "/custom/path/to/EquiBind"
  }
}
```

## Best Practices

### Path Management

1. **Always use path manager**: Never hardcode paths in your code
2. **Use environment variables**: For deployment-specific configurations
3. **Check path existence**: Use `ensure_dir()` for directories you need
4. **Handle missing tools gracefully**: Check if tools are available before using them

### Logging

1. **Use appropriate log levels**: INFO for normal operations, WARNING for issues, ERROR for failures
2. **Include error codes**: Use standardized error codes for better debugging
3. **Add context**: Provide meaningful context in log messages
4. **Time operations**: Use the timer context manager for performance tracking
5. **Avoid sensitive data**: Never log passwords, API keys, or personal information

### Error Handling

```python
try:
    result = operation()
except FileNotFoundError as e:
    logger.error(f"File not found: {e}", error_code='FILE_NOT_FOUND')
except PermissionError as e:
    logger.error(f"Permission denied: {e}", error_code='PERMISSION_ERROR')
except Exception as e:
    log_error_with_context(e, "Unexpected error", error_code='SYSTEM_ERROR')
```

## Troubleshooting

### Common Issues

1. **Path not found**: Check if the path exists and is accessible
2. **Tool not found**: Verify tool installation and PATH configuration
3. **Permission errors**: Check file and directory permissions
4. **Log file issues**: Ensure log directory is writable

### Debug Mode

Enable debug logging for troubleshooting:

```python
logger = setup_logging(__name__, level='DEBUG')
```

### Platform-Specific Notes

- **Windows**: Use forward slashes or raw strings for paths
- **macOS**: Check for Homebrew installations in `/opt/homebrew`
- **Linux**: Verify tool installations in `/usr/local/bin` or `/usr/bin`

## Integration with Existing Code

The new system is designed to be backward compatible. Existing code will continue to work, but you can gradually migrate to use the new features:

1. **Phase 1**: Add imports and basic path manager usage
2. **Phase 2**: Replace hardcoded paths with path manager calls
3. **Phase 3**: Update logging to use the new system
4. **Phase 4**: Add error codes and performance tracking

This approach ensures a smooth transition without breaking existing functionality. 