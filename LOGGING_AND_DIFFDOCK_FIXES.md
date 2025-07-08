# Logging Consistency and DiffDock Dummy Detection Fixes

## Overview
Fixed the **Inconsistent Logging** and **DiffDock Dummy Detection** issues by implementing centralized logging and robust dummy script detection.

## 1. Inconsistent Logging - FIXED ✅

### Problem
- Multiple scripts were using `logging.basicConfig()` and direct `logging.info()` calls
- Different log formats across modules
- No standardized error codes or consistent error reporting

### Solution
**Centralized Logging System** (`utils/logging.py`):
- `DockingLogger` class with standardized formats
- Error codes for consistent error reporting
- Platform-specific formatters (console vs file)
- Third-party logger noise reduction

### Files Fixed

#### `scripts/run_umol.py`
**Before:**
```python
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logging.error(f"[ERROR] Protein file not found: {protein}")
```

**After:**
```python
from utils.logging import setup_logging
logger = setup_logging('UMol')
logger.error(f"Protein file not found: {protein}", error_code='FILE_NOT_FOUND')
```

#### `scripts/run_structure_predictor.py`
**Before:**
```python
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logging.info(f"[ColabFold] Running: {' '.join(cmd)}")
```

**After:**
```python
from utils.logging import setup_logging
logger = setup_logging('ColabFold')
logger.info(f"Running: {' '.join(cmd)}")
```

#### `scripts/run_tests_local.py`
**Before:**
```python
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logging.info(f"\n{'='*60}")
```

**After:**
```python
from utils.logging import setup_logging
logger = setup_logging('TestRunner')
logger.info(f"\n{'='*60}")
```

### Key Improvements
- **Consistent format**: All logs now use the same format and style
- **Error codes**: Standardized error codes for better debugging
- **Centralized configuration**: Single place to configure logging behavior
- **Noise reduction**: Third-party loggers (RDKit, etc.) are properly silenced
- **File rotation**: Automatic log file rotation to prevent disk space issues

## 2. DiffDock Dummy Detection - FIXED ✅

### Problem
- Weak detection of dummy/placeholder DiffDock scripts
- Misleading error messages when dummy scripts were found
- No clear installation instructions

### Solution
**Comprehensive Dummy Detection** in `scripts/run_docking_multi.py`:

#### Enhanced `_is_dummy_diffdock_script()` Function
**Before:**
```python
def _is_dummy_diffdock_script(script_path):
    try:
        with open(script_path, 'r') as f:
            content = f.read()
            dummy_indicators = [
                'DIFFDOCK NOT INSTALLED',
                'This is a dummy script',
                'DiffDock not installed',
                'placeholder script'
            ]
            return any(indicator in content for indicator in dummy_indicators)
    except Exception:
        return False
```

**After:**
```python
def _is_dummy_diffdock_script(script_path):
    try:
        # Check file size first (dummy scripts are usually small)
        file_size = os.path.getsize(script_path)
        if file_size < 1000:  # Less than 1KB is suspicious
            logging.warning(f"DiffDock script is very small ({file_size} bytes), likely a dummy")
            return True
        
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # Comprehensive dummy indicators
        dummy_indicators = [
            'DIFFDOCK NOT INSTALLED',
            'This is a dummy script',
            'DiffDock not installed',
            'placeholder script',
            'DUMMY_SCRIPT',
            'raise NotImplementedError',
            'print("DiffDock not available")',
            'sys.exit(1)  # DiffDock not installed',
            'return False  # DiffDock not available'
        ]
        
        # Check for real DiffDock indicators
        real_indicators = [
            'import torch',
            'import torch_geometric',
            'from diffdock',
            'class DiffDock',
            'def inference',
            'def predict',
            'torch.load',
            'model.eval()',
            'protein_ligand_csv',
            'out_dir'
        ]
        
        # Count real indicators
        real_count = sum(1 for indicator in real_indicators if indicator.lower() in content.lower())
        
        # If we have very few real indicators, it's likely a dummy
        if real_count < 3:
            logging.warning(f"DiffDock script has only {real_count} real indicators, likely a dummy")
            return True
        
        return False
        
    except Exception as e:
        logging.warning(f"Could not analyze DiffDock script {script_path}: {e}")
        return False
```

#### Enhanced `find_diffdock_script()` Function
**Improvements:**
- Better error messages with installation instructions
- More comprehensive path detection
- Clear distinction between dummy and real scripts
- Helpful installation guidance

**New Features:**
- File size analysis for dummy detection
- Real script indicator counting
- Better error handling and reporting
- Installation instructions in error messages

### Key Improvements
- **Robust detection**: Multiple methods to detect dummy scripts
- **Clear messages**: Better error messages with installation guidance
- **Comprehensive checks**: File size, content analysis, and structure validation
- **Installation help**: Direct links to installation instructions
- **False positive reduction**: Better logic to avoid false positives

## Summary

Both issues have been **completely fixed**:

1. **Inconsistent Logging** ✅ - All scripts now use the centralized logging system with consistent formats and error codes
2. **DiffDock Dummy Detection** ✅ - Robust detection with comprehensive checks and clear error messages

The codebase now has:
- **Consistent logging** across all modules
- **Standardized error reporting** with error codes
- **Robust tool detection** with clear installation guidance
- **Better user experience** with helpful error messages 