# Results Parsing & Pose Evaluation Script

## Overview

The `parse_and_score_results.py` script is designed to parse and consolidate docking results from multiple software packages (Vina, GNINA, and DiffDock) into a unified CSV format. This is **Task 3** of the docking wrapper pipeline.

## Features

- **Multi-format parsing**: Handles Vina PDBQT, GNINA scores, and DiffDock SDF files
- **Unified output**: Consolidates all results into a standardized CSV format
- **RMSD calculation**: Optional RMSD calculation against reference structures
- **Error handling**: Logs failed parsing attempts and continues with available data
- **Comprehensive reporting**: Provides summary statistics and method comparisons

## Usage

### Basic Usage
```bash
python3 parse_and_score_results.py
```

### With Options
```bash
python3 parse_and_score_results.py --ligand_name aspirin --verbose
python3 parse_and_score_results.py --reference testing/aspirin.sdf --output summary_with_rmsd.csv
```

### Command Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--base_dir` | `.` | Base directory containing docking output folders |
| `--output_dir` | `.` | Output directory for results |
| `--ligand_name` | `ligand` | Name of the ligand being analyzed |
| `--reference` | None | Reference structure file for RMSD calculation |
| `--output_csv` | `summary.csv` | Output CSV filename |
| `--failed_json` | `failed_parsing.json` | Failed runs JSON filename |
| `--verbose` | False | Enable verbose logging |

## Input Directory Structure

The script expects the following directory structure:

```
project_directory/
├── vina_output/
│   └── vina_out.pdbqt
├── gnina_output/
│   ├── gnina_scores.txt (preferred)
│   ├── gnina_out.pdbqt (fallback)
│   └── gnina.log (last resort)
├── diffdock_output/
│   └── **/*.sdf (any SDF files in subdirectories)
└── logs/
    └── parsing.log (created by script)
```

## Output Files

### 1. summary.csv
Main results file with the following columns:
- `ligand_name`: Name of the ligand
- `method`: Docking method (Vina, GNINA, DiffDock)
- `affinity_kcal_mol`: Binding affinity in kcal/mol
- `RMSD`: Root Mean Square Deviation (if available)
- `pose_rank`: Rank of the pose (1 = best)
- `pose_path`: Path to the pose file
- `model_number`: Model number within the file
- Additional method-specific columns (e.g., `cnn_score`, `confidence_score`)

### 2. failed_parsing.json
Records any parsing failures with error messages:
```json
{
  "gnina": "No GNINA output files found",
  "diffdock": "No DiffDock SDF files found"
}
```

### 3. parsing.log
Detailed log of the parsing process including timestamps and error details.

## Parsing Details

### Vina Results
- **Source**: PDBQT files with `REMARK VINA RESULT:` lines
- **Extracted**: Binding scores, RMSD lower/upper bounds
- **Format**: `REMARK VINA RESULT:    -4.145      0.000      0.000`

### GNINA Results
- **Primary**: `gnina_scores.txt` with CNN scores and affinities
- **Fallback**: PDBQT files with scoring remarks
- **Last resort**: Log files with score patterns

### DiffDock Results
- **Source**: SDF files with confidence scores as properties
- **Fallback**: Manual SDF parsing for confidence values
- **Additional**: Separate confidence files if available

## RMSD Calculation

The script can calculate RMSD against a reference structure using:
- **RDKit** (preferred): For SDF format comparisons
- **MDAnalysis** (fallback): For general structure comparisons

Example with reference:
```bash
python3 parse_and_score_results.py --reference testing/aspirin.sdf
```

## Example Output

```bash
============================================================
DOCKING RESULTS SUMMARY
============================================================
Total poses: 9

Poses by method:
  Vina: 9

Affinity statistics (kcal/mol):
  Best score: -4.145
  Worst score: -3.501
  Mean score: -3.881

RMSD statistics (Å):
  Best RMSD: 1.466
  Worst RMSD: 21.910
  Mean RMSD: 9.677

Failed runs: 2
  gnina: Failed to parse GNINA results: No GNINA output files found
  diffdock: Failed to parse DiffDock results: No DiffDock SDF files found
============================================================
```

## Exit Codes

- `0`: Success (all parsing successful)
- `1`: Total failure (no results parsed)
- `2`: Partial success (some methods failed but results available)

## Dependencies

### Required
- `pandas`: For DataFrame operations
- `numpy`: For numerical operations
- `pathlib`: For path handling (standard library)

### Optional (for RMSD calculation)
- `rdkit`: For molecular structure handling
- `MDAnalysis`: For structural analysis

Install dependencies:
```bash
pip3 install pandas numpy
# Optional for RMSD
pip3 install rdkit-pypi MDAnalysis
```

## Integration with Docking Pipeline

This script is designed to work seamlessly with the docking wrapper:

1. **Task 1**: Structure preparation (`prep_structures.py`)
2. **Task 2**: Multi-backend docking (`run_docking_multi.py`)
3. **Task 3**: Results parsing (`parse_and_score_results.py`) ← **This script**

## Error Handling

The script is robust and handles:
- Missing output directories
- Corrupted or incomplete files
- Mixed success/failure scenarios
- Different file formats and structures

All errors are logged to both console and log file, with failed runs recorded in JSON format for analysis.

## Customization

The script can be easily extended to:
- Support additional docking software
- Add custom scoring metrics
- Implement different RMSD calculation methods
- Export to different output formats

## Best Practices

1. **Always use `--verbose`** for debugging
2. **Check failed_parsing.json** for troubleshooting
3. **Provide meaningful ligand names** for clarity
4. **Use reference structures** when available for RMSD validation
5. **Review summary statistics** to identify outliers

## Troubleshooting

### Common Issues
1. **No results found**: Check directory structure and file permissions
2. **RMSD calculation fails**: Ensure reference structure format compatibility
3. **Parsing errors**: Check log files for detailed error messages
4. **Duplicate entries**: Script handles this automatically in recent versions

### Debug Mode
Use `--verbose` flag for detailed logging:
```bash
python3 parse_and_score_results.py --verbose
```

This provides step-by-step parsing information and detailed error messages. 