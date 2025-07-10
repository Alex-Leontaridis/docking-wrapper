# Convexia Docking Wrapper Documentation

## 🚀 Project Overview

This project provides a robust, production-ready molecular docking pipeline supporting AutoDock Vina, GNINA, and DiffDock. It features batch processing, intelligent parameter detection, and containerized deployment for reproducible research.

---

## 📁 Codebase Structure

```
convexia-docking-wrapper/
├── inputs/                # Input data (proteins, ligands)
├── outputs/               # All pipeline outputs (results, logs, summaries)
├── logs/                  # Log files (preprocessing, batch, errors)
├── scripts/               # Main pipeline scripts
│   ├── batch_pipeline.py
│   ├── prep_structures.py
│   ├── run_docking_multi.py
│   └── parse_and_score_results.py
├── pipeline_config.json   # Pipeline configuration
├── docker/                # Dockerfiles, build scripts, compose
├── docs/                  # Documentation (this folder)
└── ...
```

---

## 🖥️ How to Run the Pipeline

### Native (Local) Execution

Install dependencies (see [detailed install guide](./BACKEND_INSTALLATION.md)):
```bash
pip3 install -r requirements.txt  # or see docs for conda setup
```

Run the batch pipeline:
```bash
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

Prepare structures:
```bash
python3 scripts/prep_structures.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.smi
```

---
## Quick Start (Docker)

1. **Pull the Docker image:**
   ```bash
   docker pull yourdockerhubusername/convexia-docking-wrapper:latest
   ```
2. **Run the container:**
   ```bash
   docker run -it --rm -v $(pwd):/workspace yourdockerhubusername/convexia-docking-wrapper:latest
   ```
   - Replace `$(pwd)` with your project directory if needed.
   - The code and data will be available in `/workspace` inside the container.

## Build the Docker Image (if you want to build locally)

```bash
docker build -t yourdockerhubusername/convexia-docking-wrapper:latest -f Dockerfile .
```

## Run from Source (Advanced)

1. Clone the repository:
   ```bash
   git clone https://github.com/yourgithubusername/convexia-docking-wrapper.git
   cd convexia-docking-wrapper
   ```
2. Install dependencies (see `environment.yml` or `requirements.txt`):
   ```bash
   pip3 install -r requirements.txt
   # or use conda
   conda env create -f environment.yml
   conda activate convexia-docking-wrapper
   ```

---
## 📝 CLI Commands

### 1. Batch Pipeline (`scripts/batch_pipeline.py`)

**Basic usage:**
```bash
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

**All options:**
| Option                | Description                                                        |
|-----------------------|--------------------------------------------------------------------|
| `--protein`           | Path to protein file (.pdb or .pdbqt) **[required]**               |
| `--ligands`           | Path to ligand file or directory **[required]**                    |
| `--config`            | Path to pipeline configuration JSON file                           |
| `--output-dir`        | Output directory (default: outputs)                                |
| `--enable-gnina`      | Enable GNINA docking                                               |
| `--enable-diffdock`   | Enable DiffDock docking                                            |
| `--serial`            | Process ligands serially (default: parallel)                       |
| `--max-workers`       | Maximum number of parallel workers                                 |
| `--verbose`           | Enable verbose logging                                             |

**Examples:**
```bash
# With custom config
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --config pipeline_config.json

# Enable all engines
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock

# Serial processing
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --serial
```

---

### 2. Structure Preparation (`scripts/prep_structures.py`)

**Basic usage:**
```bash
python3 scripts/prep_structures.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.smi
```

**All options:**
| Option             | Description                                                        |
|--------------------|--------------------------------------------------------------------|
| `--protein`        | Input protein file (.pdb or .pdbqt) **[required]**                 |
| `--ligand`         | Input ligand file (.smi, .sdf, or .mol2) **[required unless using --batch_ligands]** |
| `--batch_ligands`  | Multiple ligand files for batch processing                         |
| `--output_dir`     | Output directory for batch processing (default: .)                 |

**Examples:**
```bash
# Single ligand
python3 scripts/prep_structures.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.smi

# Batch mode
python3 scripts/prep_structures.py --protein inputs/protein.pdb --batch_ligands inputs/ligands/*.smi --output_dir outputs/prepared_structures/
```

---

### 3. Results Parsing & Scoring (`scripts/parse_and_score_results.py`)

**Basic usage:**
```bash
python3 scripts/parse_and_score_results.py --base_dir outputs/docking_results/ligand1 --ligand_name ligand1
```

**All options:**
| Option           | Description                                                        |
|------------------|--------------------------------------------------------------------|
| `--base_dir`     | Base directory containing docking output folders (default: .)      |
| `--batch_dirs`   | Multiple base directories for batch processing                     |
| `--output_dir`   | Output directory for results (default: .)                          |
| `--ligand_name`  | Name of the ligand being analyzed (default: ligand)                |
| `--reference`    | Reference structure file for RMSD calculation (SDF, PDB, etc.)     |
| `--output_csv`   | Output CSV filename (default: summary.csv)                         |
| `--failed_json`  | Failed runs JSON filename (default: failed_parsing.json)           |
| `--verbose`, `-v`| Enable verbose logging                                             |

**Examples:**
```bash
# Single directory
python3 scripts/parse_and_score_results.py --base_dir outputs/docking_results/ligand1 --ligand_name ligand1

# With reference structure
python3 scripts/parse_and_score_results.py --base_dir outputs/docking_results/ligand1 --ligand_name ligand1 --reference inputs/ligands/aspirin.sdf

# Batch processing
python3 scripts/parse_and_score_results.py --batch_dirs outputs/docking_results/ligand1 outputs/docking_results/ligand2 --output_dir outputs/parsed_results/
```

---

### 4. Docker Usage

**Minimal image:**
```bash
cd docker
./build_docker_minimal.sh
docker run --rm -v $(pwd)/../:/workspace docking-pipeline-minimal:latest \
    python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/
```

**Full image:**
```bash
cd docker
./build_docker.sh
docker run --rm -v $(pwd)/../:/workspace docking-pipeline:latest \
    python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock
```

**Docker Compose:**
```bash
cd docker
docker-compose up --build
```

---

## 🐳 Dockerfiles & Build Scripts

- **docker/Dockerfile**: Full environment (Vina, GNINA, DiffDock, all dependencies)
- **docker/Dockerfile.minimal**: Minimal (Vina only, smaller image)
- **docker/Dockerfile.simple**: (Simplified, see file for details)
- **docker/docker-compose.yml**: Multi-service orchestration (pipeline, minimal, jupyter, GPU)
- **docker/build_docker.sh**: Build full image
- **docker/build_docker_minimal.sh**: Build minimal image
- **docker/test_docker_setup.sh**: Test Docker images

---

## 🧪 Testing & Examples

Test files and example data are available in [`docs/testing/`](./testing/).

---

## 📄 License

**All Rights Reserved**

This software and its documentation are proprietary and confidential. All rights, title, and interest in and to this software, including all intellectual property rights, are and will remain the exclusive property of the original author.

**No Rights Granted:**
- No permission is granted to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of this software.
- No permission is granted to use this software for any commercial or non-commercial purpose.
- No permission is granted to create derivative works based on this software.

**Unauthorized Use Prohibited:**
Any unauthorized use, reproduction, or distribution of this software is strictly prohibited and may result in legal action.

For permission to use this software, contact the original author directly.

