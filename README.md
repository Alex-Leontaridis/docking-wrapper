# Molecular Docking Wrapper System

A comprehensive multi-backend molecular docking system supporting AutoDock Vina, GNINA, and DiffDock with robust error handling and automatic parameter detection.

## 🚀 **Features**

- **Multi-Backend Support**: AutoDock Vina, GNINA, and DiffDock
- **Intelligent Box Detection**: **4-tier automatic parameter detection system**
  1. **Bound Ligand Detection**: Uses existing ligands for precise targeting
  2. **Cavity Detection**: AI-powered cavity identification using geometric analysis
  3. **Protein Center**: Geometric center with conservative box sizing
  4. **Default Fallback**: Robust default parameters for any edge cases
- **Robust Error Handling**: Continues execution even if individual backends fail
- **Comprehensive Logging**: Detailed logs with timestamps and error tracking
- **Input Validation**: Checks file formats and required dependencies
- **PDBQT Formatting**: Automatic cleanup of malformed PDBQT files from MGLTools
- **Production Scale**: Designed for processing thousands of proteins automatically

## 📦 **Installation & Setup**

### Prerequisites
```bash
# Install AutoDock Vina
brew install autodock-vina

# Install MGLTools for structure preparation
# Download from: http://mgltools.scripps.edu/downloads
```

### Dependencies
```bash
pip3 install numpy biopython scikit-learn
```

## 🔧 **Structure Preparation**

The system includes automatic structure preparation:

```bash
# Prepare protein and ligand structures
python3 prep_structures.py testing/1ubq.pdb testing/aspirin.sdf
```

This will generate:
- `protein_prepped.pdbqt` (cleaned and formatted)
- `ligand_prepped.pdbqt` (prepared ligand)

## 🎯 **Usage**

### Basic Usage (Vina only)
```bash
python3 run_docking_multi.py \
    --protein protein_prepped.pdbqt \
    --ligand ligand_prepped.pdbqt \
    --center_x 20 --center_y 20 --center_z 20 \
    --size_x 20 --size_y 20 --size_z 20
```

### Multi-Backend Usage
```bash
python3 run_docking_multi.py \
    --protein protein_prepped.pdbqt \
    --ligand ligand_prepped.pdbqt \
    --center_x 20 --center_y 20 --center_z 20 \
    --size_x 20 --size_y 20 --size_z 20 \
    --use_gnina \
    --use_diffdock
```

### **Automatic Box Detection (Recommended for Production)**
For thousands of proteins, simply omit box parameters and let the system detect them automatically:
```bash
python3 run_docking_multi.py \
    --protein protein_prepped.pdbqt \
    --ligand ligand_prepped.pdbqt
```

### Automatic Box Detection with Multiple Backends
```bash
python3 run_docking_multi.py \
    --protein protein_prepped.pdbqt \
    --ligand ligand_prepped.pdbqt \
    --use_gnina \
    --use_diffdock
```

## 📁 **Output Structure**

```
./
├── vina_output/
│   └── vina_out.pdbqt          # Docked poses with scores
├── gnina_output/
│   ├── gnina_out.pdbqt         # GNINA docked poses
│   ├── gnina.log               # Detailed GNINA log
│   └── gnina_scores.txt        # Scoring results
├── diffdock_output/
│   └── [complex structures]     # DiffDock results
└── logs/
    ├── docking_run.log         # Complete execution log
    └── failed_runs.json        # Error details for debugging
```

## 🎛️ **Command Line Options**

| Option | Description | Required |
|--------|-------------|----------|
| `--protein` | Path to receptor PDBQT file | ✅ |
| `--ligand` | Path to ligand PDBQT file | ✅ |
| `--center_x/y/z` | Grid center coordinates | ⚠️ (if no bound ligand) |
| `--size_x/y/z` | Grid dimensions (Å) | ⚠️ (if no bound ligand) |
| `--use_gnina` | Enable GNINA backend | ❌ |
| `--use_diffdock` | Enable DiffDock backend | ❌ |
| `--output_dir` | Output directory (default: current) | ❌ |

## 📊 **Backend Status**

| Backend | Status | Requirements | Notes |
|---------|--------|--------------|-------|
| **AutoDock Vina** | ✅ **Working** | `vina` in PATH | Fast, reliable docking |
| **GNINA** | ⚠️ **Partial** | CUDA + Build required | CPU-only build possible |
| **DiffDock** | ⚠️ **Partial** | CUDA + PyTorch | Deep learning approach |

## 🔧 **GNINA Setup** (Optional)

GNINA requires compilation from source:
```bash
cd gnina
mkdir build && cd build
cmake .. -DCUDA=OFF  # For CPU-only
make -j$(nproc)
```

## 🌊 **DiffDock Setup** (Optional)

DiffDock requires PyTorch and CUDA:
```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
pip3 install -r requirements.txt
```

## 🎯 **Automatic Box Parameter Detection**

The system uses a sophisticated 4-tier strategy to automatically determine optimal docking parameters:

### Strategy 1: Bound Ligand Detection
- Scans for HETATM records with ≤50 heavy atoms
- Calculates geometric center and bounding box
- Adds 8Å padding for optimal search space
- **Best for**: Proteins with co-crystallized ligands

### Strategy 2: Cavity Detection
- Uses 3D grid-based analysis to identify binding cavities
- Employs DBSCAN clustering to find distinct pockets
- Selects largest cavity by volume
- Adds 10Å padding for comprehensive coverage
- **Best for**: Apo structures without bound ligands

### Strategy 3: Protein Geometric Center
- Calculates protein center of mass
- Uses 40% of protein span (20-30Å per dimension)
- Conservative approach for unknown binding sites
- **Best for**: Large proteins or unusual cases

### Strategy 4: Default Parameters
- Last resort: center (0,0,0) with 25Å³ box
- Ensures system never fails completely
- **Best for**: Edge cases or corrupted files

### Example Output
INFO: Strategy 2: Using largest detected cavity (volume: 145600.0 Ų)
INFO: Extracted box center: (30.82, 29.37, 16.50), size: (50.00, 52.00, 56.00)

## 📝 **Example Output**

```
=== Docking Summary ===
VINA: SUCCESS (time: 2.28s)
GNINA: SKIPPED/FAILED (time: 0.0s)
  Reason: GNINA binary not found. Please install GNINA or build it from source.
DIFFDOCK: SKIPPED/FAILED (time: 3.38s)
  Reason: CUDA dependencies not available on macOS
======================
```

## 🏆 **Key Features Implemented**

1. **Robust Error Handling**: Each backend runs independently
2. **Automatic Cleanup**: PDBQT formatting fixes for MGLTools output
3. **Smart Detection**: Automatic binding site parameter extraction
4. **Comprehensive Logging**: Full execution tracking and debugging
5. **Flexible Backend Selection**: Mix and match docking methods
6. **Cross-Platform**: Works on macOS, Linux (with appropriate dependencies)

## 🚧 **Known Limitations**

- GNINA and DiffDock require CUDA for optimal performance
- macOS lacks CUDA support (affects GNINA/DiffDock)
- Large proteins may require more memory for DiffDock
- Box parameters required if no bound ligand is present

## 📞 **Support**

- Check `logs/docking_run.log` for detailed execution info
- Review `logs/failed_runs.json` for specific error details
- Ensure all dependencies are properly installed
- Verify input file formats (PDBQT for proteins/ligands)

---

**Status**: Production ready for AutoDock Vina with optional GNINA/DiffDock support 