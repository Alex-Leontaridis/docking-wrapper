# Batch Molecular Docking Pipeline & Docker Environment

A comprehensive, production-ready molecular docking system supporting **AutoDock Vina**, **GNINA**, and **DiffDock** with automated batch processing, intelligent parameter detection, and containerized deployment.

## 🚀 Features

### Pipeline Capabilities
- **Multi-Engine Support**: AutoDock Vina, GNINA, and DiffDock
- **Batch Processing**: Process hundreds of ligands automatically
- **Intelligent Box Detection**: 4-tier automatic parameter detection
- **Parallel Processing**: Multi-core ligand processing
- **Comprehensive Logging**: Detailed execution tracking
- **Error Recovery**: Robust error handling and continuation
- **Results Aggregation**: Automatic summary generation

### Docker Environment
- **Zero Setup**: Complete containerized environment
- **All Dependencies**: RDKit, OpenBabel, Meeko, pandas, and more
- **Production Ready**: Optimized for high-throughput processing
- **Jupyter Integration**: Optional notebook environment for analysis
- **Resource Management**: Configurable CPU/memory limits

## 📦 Quick Start

### Using Docker (Recommended)

1. **Build the container:**
```bash
docker build -t docking-pipeline .
```

2. **Run batch docking:**
```bash
# Basic usage - Vina only
docker run -v $(pwd):/workspace docking-pipeline \
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/

# With all engines enabled
docker run -v $(pwd):/workspace docking-pipeline \
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --enable-gnina --enable-diffdock
```

3. **Using Docker Compose:**
```bash
# Start the pipeline service
docker-compose up -d docking-pipeline

# Execute docking
docker-compose exec docking-pipeline \
  python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/

# Start Jupyter for analysis (optional)
docker-compose --profile jupyter up -d
```

### Native Installation

1. **Install dependencies:**
```bash
# Core Python packages
pip3 install numpy pandas rdkit-pypi biopython scikit-learn meeko

# Docking software
# - Install AutoDock Vina from https://github.com/ccsb-scripps/AutoDock-Vina
# - Install GNINA from https://github.com/gnina/gnina
# - Install DiffDock from https://github.com/gcorso/DiffDock
```

2. **Run pipeline:**
```bash
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/
```

## 📖 Usage Guide

### Basic Usage

```bash
# Process single ligand
python3 batch_pipeline.py --protein receptor.pdb --ligands aspirin.sdf

# Process directory of ligands
python3 batch_pipeline.py --protein receptor.pdb --ligands ligand_library/

# Enable all docking engines
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --enable-gnina --enable-diffdock

# Use custom configuration
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --config custom_config.json

# Serial processing (no parallelization)
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --serial
```

### Advanced Options

```bash
# Control parallel processing
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --max-workers 8

# Verbose logging
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --verbose

# Custom output directory
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --output-dir results_2024
```

## 📁 Output Structure

The pipeline creates a comprehensive output structure:

```
outputs/
├── logs/
│   ├── batch_log.txt              # Main execution log
│   └── batch_log.json             # Detailed execution data
├── prepared_structures/
│   ├── receptor_prepared.pdbqt     # Prepared protein
│   ├── ligand1_prepared.pdbqt      # Prepared ligands
│   └── ligand2_prepared.pdbqt
├── docking_results/
│   ├── ligand1/
│   │   ├── vina_output/
│   │   │   └── vina_out.pdbqt      # Vina poses
│   │   ├── gnina_output/
│   │   │   ├── gnina_out.pdbqt     # GNINA poses
│   │   │   └── gnina_scores.txt    # GNINA scores
│   │   └── diffdock_output/
│   │       └── *.sdf               # DiffDock poses
│   └── ligand2/
│       └── ...
├── parsed_results/
│   ├── ligand1/
│   │   ├── ligand1_summary.csv     # Per-ligand results
│   │   └── ligand1_failed.json     # Failed engines
│   └── ligand2/
│       └── ...
└── final_summary.csv               # Aggregated results
```

## ⚙️ Configuration

The pipeline uses `pipeline_config.json` for advanced configuration:

### Engine Settings
```json
{
  "engines": {
    "vina": {
      "enabled": true,
      "exhaustiveness": 8,
      "num_modes": 9
    },
    "gnina": {
      "enabled": false,
      "use_gpu": false,
      "cnn_scoring": "rescore"
    },
    "diffdock": {
      "enabled": false,
      "inference_steps": 20,
      "samples_per_complex": 10
    }
  }
}
```

### Box Detection
```json
{
  "box": {
    "auto_detect": true,
    "default_size": [25.0, 25.0, 25.0],
    "padding": 8.0
  }
}
```

### Performance Tuning
```json
{
  "parallel": {
    "max_workers": 4,
    "chunk_size": 1
  },
  "timeouts": {
    "vina": 1800,
    "gnina": 3600,
    "diffdock": 7200
  }
}
```

## 🧪 Supported Input Formats

### Proteins
- `.pdb` - Protein Data Bank format
- `.pdbqt` - AutoDock PDBQT format (pre-prepared)

### Ligands
- `.smi` - SMILES strings
- `.sdf` - Structure Data Format
- `.mol2` - Tripos MOL2 format

## 📊 Results Analysis

### CSV Output Format

The `final_summary.csv` contains:

| Column | Description |
|--------|-------------|
| `ligand_name` | Name of the ligand |
| `method` | Docking engine (Vina, GNINA, DiffDock) |
| `affinity_kcal_mol` | Binding affinity (kcal/mol) |
| `RMSD` | Root mean square deviation |
| `pose_rank` | Rank of the pose |
| `pose_path` | Path to the pose file |
| `global_rank` | Overall ranking across all results |

### Analysis Examples

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
df = pd.read_csv('outputs/final_summary.csv')

# Best poses per ligand
best_poses = df.groupby('ligand_name')['affinity_kcal_mol'].min()

# Method comparison
method_stats = df.groupby('method')['affinity_kcal_mol'].describe()

# Plot affinity distribution
plt.figure(figsize=(10, 6))
df.boxplot(column='affinity_kcal_mol', by='method')
plt.title('Binding Affinity Distribution by Method')
plt.show()
```

## 🔧 Docker Environment Details

### Included Software

| Software | Version | Purpose |
|----------|---------|---------|
| AutoDock Vina | Latest | Classical docking |
| GNINA | Latest | CNN-enhanced docking |
| DiffDock | Latest | Diffusion-based docking |
| RDKit | Latest | Molecular manipulation |
| OpenBabel | Latest | Format conversion |
| Meeko | Latest | Structure preparation |
| Python | 3.10 | Runtime environment |

### Container Resources

```yaml
# Default resource limits
deploy:
  resources:
    limits:
      cpus: '8.0'
      memory: 16G
    reservations:
      cpus: '2.0'
      memory: 4G
```

### Volume Mounts

- `/workspace` - Your working directory
- Persistent caches for RDKit and conda packages

## 🚀 Performance Optimization

### Parallel Processing

```bash
# Use all available cores
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --max-workers $(nproc)

# Conservative approach for memory-limited systems
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --max-workers 2
```

### Large-Scale Processing

For thousands of ligands:

1. **Enable only fast engines** (Vina only)
2. **Use appropriate timeouts** in config
3. **Monitor disk space** for outputs
4. **Process in batches** if memory is limited

### Memory Management

```json
{
  "performance": {
    "memory_limit_gb": 8,
    "disk_space_check": true,
    "min_free_space_gb": 5
  }
}
```

## 🐛 Troubleshooting

### Common Issues

1. **"No ligands found"**
   - Check file formats are supported (.smi, .sdf, .mol2)
   - Verify directory path is correct

2. **"Protein preparation failed"**
   - Ensure PDB file is valid
   - Check for missing atoms or corrupted structure

3. **"Docker build fails"**
   - Ensure sufficient disk space (>20GB)
   - Check internet connection for downloads

4. **"Out of memory"**
   - Reduce `max_workers` in config
   - Process ligands in smaller batches

### Debug Mode

```bash
# Enable verbose logging
python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ \
  --verbose

# Check detailed logs
tail -f outputs/logs/batch_log.txt
```

### Validation

```bash
# Test with single ligand first
python3 batch_pipeline.py --protein receptor.pdb --ligands test_ligand.sdf

# Verify Docker environment
docker run docking-pipeline python3 -c "import rdkit, pandas; print('OK')"
```

## 📈 Benchmarking

### Performance Expectations

| Engine | Ligands/hour | Memory Usage | Notes |
|--------|--------------|--------------|-------|
| Vina | 100-200 | ~500MB | Fast, reliable |
| GNINA | 20-50 | ~2GB | Slower, higher quality |
| DiffDock | 5-15 | ~4GB | Slowest, novel approach |

### Scaling Guidelines

- **Small scale** (< 100 ligands): Use all engines
- **Medium scale** (100-1000 ligands): Vina + GNINA
- **Large scale** (> 1000 ligands): Vina only

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🔗 References

- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
- [GNINA](https://github.com/gnina/gnina)
- [DiffDock](https://github.com/gcorso/DiffDock)
- [RDKit](https://github.com/rdkit/rdkit)
- [Meeko](https://github.com/forlilab/Meeko)

## 📞 Support

For issues and questions:
1. Check the troubleshooting section
2. Review logs in `outputs/logs/`
3. Open an issue with detailed error information
4. Include system specifications and input files

---

**Happy Docking! 🧬** 