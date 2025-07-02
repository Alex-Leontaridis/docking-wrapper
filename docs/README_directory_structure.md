# Molecular Docking Pipeline - Directory Structure

## 📁 Directory Organization

This project follows a clean, organized structure that separates inputs, outputs, scripts, and logs:

```
/docking_wrapper/
├── inputs/                          # Input data
│   ├── protein.pdb                 # Target protein structure
│   └── ligands/                    # Ligand molecules
│       ├── ligand1.smi            # SMILES format ligand
│       ├── ligand2.smi            # Another SMILES ligand
│       └── aspirin.sdf            # SDF format ligand
├── outputs/                        # All pipeline outputs
│   ├── ligand1/                   # Results for ligand1
│   │   ├── vina/                  # AutoDock Vina results
│   │   ├── gnina/                 # GNINA results
│   │   └── diffdock/              # DiffDock results
│   ├── ligand2/                   # Results for ligand2
│   └── final_summary.csv          # Aggregated results
├── logs/                          # All log files
│   ├── preprocessing_log.txt      # Structure preparation logs
│   ├── batch_log.txt             # Batch processing logs
│   └── failed_runs.json          # Failed run details
├── scripts/                       # All pipeline scripts
│   ├── prep_structures.py        # Structure preparation
│   ├── run_docking_multi.py       # Docking execution
│   ├── parse_and_score_results.py # Results parsing
│   └── batch_pipeline.py          # Main batch orchestrator
├── Dockerfile                     # Container definition
├── pipeline_config.json           # Pipeline configuration
└── run_batch_pipeline.py          # Wrapper script
```

## 🚀 Usage

### Quick Start
```bash
# Run batch docking with the organized structure
python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/

# Or use the script directly
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

### Individual Task Scripts
```bash
# Structure preparation (Task 1)
python3 scripts/prep_structures.py inputs/protein.pdb inputs/ligands/ligand1.smi

# Docking (Task 2)
python3 scripts/run_docking_multi.py

# Results parsing (Task 3)
python3 scripts/parse_and_score_results.py
```

## 📂 Directory Details

### `inputs/`
- **protein.pdb**: Your target protein structure
- **ligands/**: Directory containing all ligand files
  - Supports `.smi`, `.sdf`, `.mol2`, `.pdb` formats
  - Can mix different file formats in the same directory

### `outputs/`
- **ligandX/**: Individual ligand results
  - **vina/**: AutoDock Vina output files (.pdbqt, .log)
  - **gnina/**: GNINA output files (.sdf, .log)
  - **diffdock/**: DiffDock output files (.sdf, .log)
- **final_summary.csv**: Aggregated results from all ligands and engines

### `logs/`
- **batch_log.txt**: Main batch processing log with timestamps
- **preprocessing_log.txt**: Detailed structure preparation logs
- **failed_runs.json**: JSON file tracking any failed ligand processing

### `scripts/`
Contains all the main pipeline components:
- **batch_pipeline.py**: Main orchestrator for batch processing
- **prep_structures.py**: Structure preparation and validation
- **run_docking_multi.py**: Multi-engine docking execution
- **parse_and_score_results.py**: Results parsing and analysis

## 🔧 Configuration

The `pipeline_config.json` file controls all pipeline behavior:

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

## 🎯 Benefits of This Structure

1. **Clear Separation**: Inputs, outputs, scripts, and logs are clearly separated
2. **Scalability**: Easy to add new ligands by dropping files in `inputs/ligands/`
3. **Organization**: Results are automatically organized by ligand name
4. **Maintainability**: Scripts are centralized and easily accessible
5. **Docker Ready**: Structure works seamlessly with containerization

## 🐳 Docker Usage

The directory structure is optimized for Docker deployment:

```bash
# Build the container
docker build -t docking-pipeline .

# Run with volume mounts
docker run -v $(pwd)/inputs:/app/inputs \
           -v $(pwd)/outputs:/app/outputs \
           -v $(pwd)/logs:/app/logs \
           docking-pipeline \
           --protein inputs/protein.pdb \
           --ligands inputs/ligands/
```

## 📊 Output Analysis

After running the pipeline, check:

1. **`outputs/final_summary.csv`**: Combined results from all ligands
2. **`logs/batch_log.txt`**: Processing details and any issues
3. **`outputs/ligandX/`**: Individual detailed results for each ligand

This structure makes it easy to scale from single ligands to hundreds of compounds while maintaining organization and traceability. 