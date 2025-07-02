# Molecular Docking Pipeline

A comprehensive, production-ready molecular docking pipeline supporting AutoDock Vina, GNINA, and DiffDock with batch processing capabilities.

## 📁 Clean Directory Structure

```
/docking_wrapper/
├── inputs/                    # Input data
│   ├── protein.pdb           # Target protein structure
│   └── ligands/              # Ligand molecules (.smi, .sdf, .mol2)
├── outputs/                  # All pipeline results
│   ├── ligandX/              # Individual ligand results
│   │   ├── vina/            # AutoDock Vina outputs
│   │   ├── gnina/           # GNINA outputs
│   │   └── diffdock/        # DiffDock outputs
│   └── final_summary.csv    # Aggregated results
├── logs/                     # Centralized logging
│   ├── batch_log.txt        # Main processing log
│   └── failed_runs.json     # Failed run tracking
├── scripts/                  # Core pipeline components
│   ├── prep_structures.py   # Task 1: Structure preparation
│   ├── run_docking_multi.py # Task 2: Multi-engine docking
│   ├── parse_and_score_results.py # Task 3: Results parsing
│   └── batch_pipeline.py    # Task 4: Batch orchestration
├── docs/                     # Documentation and examples
├── Dockerfile               # Container definition
├── pipeline_config.json     # Pipeline configuration
└── run_batch_pipeline.py    # Main entry point
```

## 🚀 Quick Start

### Basic Usage
```bash
# Run batch docking with sample data
python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/

# Enable all docking engines
python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock

# Use custom configuration
python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --config pipeline_config.json
```

### Docker Usage
```bash
# Build container
docker build -t docking-pipeline .

# Run with volume mounts
docker run -v $(pwd)/inputs:/app/inputs \
           -v $(pwd)/outputs:/app/outputs \
           -v $(pwd)/logs:/app/logs \
           docking-pipeline \
           --protein inputs/protein.pdb --ligands inputs/ligands/
```

## 📋 Pipeline Tasks

### Task 1: Structure Preparation
```bash
python3 scripts/prep_structures.py inputs/protein.pdb inputs/ligands/ligand1.smi
```

### Task 2: Multi-Engine Docking
```bash
python3 scripts/run_docking_multi.py
```

### Task 3: Results Parsing
```bash
python3 scripts/parse_and_score_results.py
```

### Task 4: Batch Processing
```bash
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

## ⚙️ Configuration

Edit `pipeline_config.json` to customize:

```json
{
    "engines": {
        "vina": {"enabled": true, "exhaustiveness": 8},
        "gnina": {"enabled": false, "use_gpu": false},
        "diffdock": {"enabled": false}
    },
    "box": {
        "auto_detect": true,
        "default_size": [25.0, 25.0, 25.0]
    },
    "parallel": {
        "max_workers": 4
    }
}
```

## 📊 Output Analysis

After running the pipeline:

1. **`outputs/final_summary.csv`** - Combined results from all ligands
2. **`logs/batch_log.txt`** - Detailed processing logs
3. **`outputs/ligandX/`** - Individual detailed results per ligand

## 🧪 Testing

```bash
# Run validation tests
python3 test_batch_pipeline.py
```

## 📚 Documentation

Complete documentation available in `docs/`:
- `docs/README_directory_structure.md` - Directory layout details
- `docs/README_batch_pipeline.md` - Batch processing guide
- `docs/README_parse_results.md` - Results analysis guide
- `docs/TASK4_COMPLETION_SUMMARY.md` - Implementation summary

## 🔧 Requirements

- Python 3.8+
- AutoDock Vina
- MGLTools (for structure preparation)
- Optional: GNINA, DiffDock, Docker

## 📈 Features

- **Multi-Engine Support**: Vina, GNINA, DiffDock
- **Batch Processing**: Handle multiple ligands automatically
- **Parallel Execution**: Configurable worker processes
- **Error Handling**: Comprehensive logging and recovery
- **Docker Ready**: Complete containerized environment
- **Flexible Input**: SMILES, SDF, MOL2, PDB formats
- **Auto Box Detection**: Intelligent binding site detection

