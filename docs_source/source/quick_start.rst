Quick Start Guide
=================

Installation
------------

1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/docking-wrapper.git
   cd docking-wrapper
   ```

2. **Install Python dependencies:**
   ```bash
   pip install -r requirements_enhanced.txt
   ```

3. **Install external tools (optional):**
   ```bash
   python scripts/install_backends.py
   ```

Basic Usage
-----------

**Run a simple docking job:**
```bash
python scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

**Enable specific engines:**
```bash
python scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock
```

**Run with ML models:**
```bash
python scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-equibind --enable-neuralplexer
```

**Full enhanced pipeline:**
```bash
python scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-all-ml --enable-all-analysis --enable-consensus --enable-confidence
```

Docker Usage
------------

**Quick start with Docker:**
```bash
docker run -it --rm -v $(pwd):/workspace yourusername/docking-wrapper:latest
```

**Build and run locally:**
```bash
cd docker
./build_docker.sh
docker run --rm -v $(pwd)/../:/workspace docking-pipeline:latest python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

Configuration
-------------

Create a custom configuration file `pipeline_config.json`:
```json
{
  "docking": {
    "exhaustiveness": 8,
    "num_poses": 10,
    "energy_range": 3
  },
  "ml_models": {
    "equibind": {
      "enabled": true,
      "num_poses": 5
    },
    "neuralplexer": {
      "enabled": true,
      "confidence_threshold": 0.7
    }
  },
  "analysis": {
    "interactions": {
      "enabled": true
    },
    "druggability": {
      "enabled": true
    }
  }
}
```

Output Structure
----------------

The pipeline generates organized outputs:
```
outputs/
├── docking_results/
│   └── ligand_name/
│       ├── vina_output/
│       ├── gnina_output/
│       ├── diffdock_output/
│       └── logs/
├── ml_results/
│   └── ligand_name/
│       ├── equibind_output/
│       ├── neuralplexer_output/
│       └── umol_output/
├── analysis_results/
│   └── ligand_name/
│       ├── interactions/
│       ├── druggability/
│       └── consensus/
└── summary_reports/
    ├── final_summary.csv
    └── confidence_scores.json
```

Next Steps
----------

- Read the :doc:`user_guide` for detailed usage instructions
- Check the :doc:`api_reference` for programmatic access
- Explore the :doc:`developer_guide` for contributing
- Review the :doc:`architecture` for system design 