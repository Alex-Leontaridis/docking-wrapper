# Docking Wrapper

## Overview
This project provides a batch molecular docking pipeline supporting AutoDock Vina, GNINA, and DiffDock. It includes structure preparation, docking, and results parsing, with robust logging and error handling.

---

## Installation

1. **Python dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```
2. **External binaries:**
   - **AutoDock Vina**: Must be installed and available in your PATH.
   - **GNINA**: Must be installed and available in your PATH (or use the provided dummy for testing).
   - **DiffDock**: Must be installed and available in your PATH (or use the provided dummy for testing).
   - **MGLTools** (for protein prep, optional): Place in `~/mgltools_1.5.7_MacOS-X/` or adjust script paths.

Scripts will check for these binaries at startup and exit with a clear error if not found.


### Making GNINA Executable

After downloading GNINA, you may need to make it executable:

```bash
# Make gnina executable (Linux/macOS)
chmod +x /path/to/gnina

# Or if using the project's bin directory
chmod +x bin/gnina
```

### Fixing Output Directory Permissions

If you encounter permission errors with DiffDock or other output directories:

```bash
# Navigate to your project directory
cd ~/Desktop/docking-wrapper

# Fix permissions for the output directory
chmod -R 755 test_bug_fix_outputs/

# Make sure the user has write permissions
chmod -R u+w test_bug_fix_outputs/
```
---

## CLI Usage

### Enhanced Batch Pipeline (Recommended)
The enhanced pipeline now integrates all ML/Analysis tools into a single unified workflow with 7 comprehensive stages.

#### **Traditional Docking Only**
Basic molecular docking with Vina, GNINA, and DiffDock:

```bash
# Apply to test_bug_fix scenario
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --output-dir test_bug_fix_outputs/
# Run
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock --output-dir test_bug_fix_outputs/
```

#### **Full Enhanced Pipeline**
Complete decision-ready insight engine with ML models, analysis tools, consensus, and confidence scoring:

```bash
# Apply to test_bug_fix scenario
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-all-ml --enable-all-analysis --output-dir test_bug_fix_outputs/
# Run
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-all-ml --enable-all-analysis --enable-consensus --enable-confidence --output-dir test_bug_fix_outputs/
```

#### **Custom Configuration**
Use a custom configuration file for specific requirements:

```bash
# Apply to test_bug_fix scenario
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --config pipeline_config.json --output-dir test_bug_fix_outputs/
# Run
python3 scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --config pipeline_config.json --enable-equibind --enable-neuralplexer --enable-boltz2 --enable-interactions --output-dir test_bug_fix_outputs/
```

### Pipeline Stages
The enhanced pipeline includes 7 comprehensive stages:

1. **Structure Preparation**: Convert protein and ligands to PDBQT format
2. **Traditional Docking**: Run Vina, GNINA, DiffDock (if enabled)
3. **ML Model Docking**: Run EquiBind, NeuralPLexer, UMol (if enabled)
4. **Analysis Tools**: Run Boltz2, Interactions, Druggability (if enabled)
5. **Consensus Analysis**: RMSD-based pose clustering (if enabled)
6. **Confidence Scoring**: Composite confidence score 0-100 (if enabled)
7. **Results Parsing**: Generate comprehensive summaries

### Command-Line Options
- **Traditional engines**: `--enable-gnina --enable-diffdock`
- **ML models**: `--enable-equibind --enable-neuralplexer --enable-umol --enable-structure-predictor`
- **Analysis tools**: `--enable-boltz2 --enable-interactions --enable-druggability`
- **Convenience flags**: `--enable-all-ml --enable-all-analysis`
- **Processing**: `--serial --max-workers 4 --verbose`

- **Config file options:** See `pipeline_config.json` or the script docstring for all options. The config is validated at startup.
- **Where GNINA and DiffDock are called:**
  - GNINA and DiffDock are invoked in `scripts/run_docking_multi.py` via the functions `run_gnina` and `run_diffdock`.
  - The batch pipeline enables these via `--enable-gnina` and `--enable-diffdock`.

### 2. Structure Preparation
Prepare proteins and ligands for docking.

```bash
python3 scripts/prep_structures.py --protein protein.pdb --ligand ligand.smi
# Batch mode
python3 scripts/prep_structures.py --protein protein.pdb --batch_ligands ligand1.smi ligand2.sdf --output_dir prepped/
```

### 3. Docking (Single or Batch)
Run docking for prepared structures.

```bash
python3 scripts/run_docking_multi.py --protein protein_prepped.pdbqt --ligand ligand_prepped.pdbqt
# Batch mode
python3 scripts/run_docking_multi.py --batch_proteins p1.pdbqt p2.pdbqt --batch_ligands l1.pdbqt l2.pdbqt --use_gnina --use_diffdock
```

- **GNINA and DiffDock location:**
  - GNINA: `run_gnina` in `scripts/run_docking_multi.py`
  - DiffDock: `run_diffdock` in `scripts/run_docking_multi.py`

### 4. Results Parsing
Parse and score docking results.

```bash
python3 scripts/parse_and_score_results.py --base_dir outputs/
# Batch mode
python3 scripts/parse_and_score_results.py --batch_dirs out1/ out2/ --output_dir summary/
```

---

## Individual Script Usage (Advanced)

The following commands can be run individually for specific tasks or custom workflows. However, the **recommended approach** is to use the enhanced batch pipeline above, which integrates all these tools automatically.

### 1. EquiBind Pose Prediction
```bash
python3 scripts/run_equibind.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.sdf --output outputs/equibind/
```

### 2. NeuralPLexer-2 Pose Prediction
```bash
python3 scripts/run_neuralplexer.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.sdf --output outputs/neuralplexer/
```

### 3. UMol Structure-Free Pose Prediction
```bash
python3 scripts/run_umol.py --protein inputs/protein.fasta --ligand inputs/ligands/aspirin.sdf --output outputs/umol/
```

### 4. Structure Prediction (AlphaFold2/OpenFold/ESMFold)
```bash
python3 scripts/run_structure_predictor.py --protein inputs/protein.fasta --output_dir outputs/structures/
```

### 5. Binding Affinity Prediction (Boltz2)
```bash
python3 scripts/run_boltz2.py --protein inputs/protein.fasta --ligand inputs/ligands/aspirin.smi --output outputs/affinity.json
```

### 6. Protein-Ligand Interaction Analysis (PLIP/RDKit)
```bash
python3 scripts/extract_interactions.py --pdb outputs/equibind/equibind_pose.pdb --ligand LIG --protein PROT --output_dir outputs/interactions/
```

### 7. Druggability Scoring (fpocket)
```bash
python3 scripts/run_druggability.py --protein outputs/equibind/equibind_pose.pdb --output outputs/druggability.json
```

### 8. Consensus Analysis (RMSD-based)
```bash
python3 scripts/model_consensus.py --poses outputs/equibind/equibind_pose.pdb outputs/neuralplexer/neuralplexer_pose.pdb outputs/umol/umol_pose.pdb --output outputs/consensus.json --ligand_id LIG1
```

### 9. Composite Confidence Scoring
```bash
python3 scripts/compute_confidence.py --consensus_json outputs/consensus.json --druggability_json outputs/druggability.json --affinity_json outputs/affinity.json --interaction_json outputs/interactions/LIG_PROT.json --output outputs/confidence.json --ligand_id LIG1
```

---

## Enhanced Output Structure

The enhanced pipeline generates comprehensive, organized outputs:

```
test_bug_fix_outputs/
├── docking_results/
│   └── ligand_name/
│       ├── vina_output/          # Traditional docking results
│       ├── gnina_output/         # GNINA docking results
│       ├── diffdock_output/      # DiffDock results
│       ├── ml_models/            # ML model poses
│       │   ├── equibind/
│       │   ├── neuralplexer/
│       │   └── umol/
│       ├── analysis/             # Analysis results
│       │   ├── boltz2/
│       │   ├── interactions/
│       │   └── druggability/
│       ├── consensus/            # Pose consensus analysis
│       └── confidence/           # Confidence scores
├── parsed_results/               # Parsed summaries
├── prepared_structures/          # Prepared protein/ligands
└── logs/                         # Comprehensive logs
```

### Key Output Files
- **Traditional**: CSV summaries with poses and affinity scores
- **ML Models**: PDB/SDF files with predicted poses
- **Analysis**: JSON files with detailed metrics
- **Consensus**: JSON with pose agreement analysis
- **Confidence**: JSON with composite confidence scores (0-100)

---

## Dummy Testing for CI/Review

This project includes **dummy scripts** for GNINA and DiffDock to allow reviewers and CI systems to test the pipeline without the real binaries:

- **Dummy GNINA:** `bin/gnina` (bash script)
- **Dummy DiffDock:** `DiffDock/inference.py` (Python script)

To use the dummies:
```bash
export PATH="$PWD/bin:$PATH"
```
- The pipeline will call the dummy scripts, which simulate successful runs and create expected output files.
- This allows full workflow testing, including error handling and summary generation, without requiring GPU or special licenses.

**To use real binaries:**
- Install GNINA and DiffDock, and ensure their real executables are in your PATH before the dummies.
- Remove or rename the dummy scripts if needed.

---

## Configuration
- All config options are documented in `pipeline_config.json` and in the docstrings of `scripts/batch_pipeline.py`.
- The pipeline validates config files and CLI arguments at startup, exiting with an error if invalid.

---

## Error Handling
- All scripts provide robust error handling and will print clear error messages for missing files, invalid arguments, or missing dependencies.
- If a required external binary (vina, gnina, diffdock) is not found, the script will exit with an error.

---

## External Binaries Check
At startup, the scripts check for required binaries using `shutil.which()` or similar. If a binary is missing, a clear error is printed and the script exits.

---

## Project Status
- **Enhanced pipeline ready for production use.**
- **Complete ML/Analysis integration** with decision-ready insights.
- **Unified workflow** with 7 comprehensive stages.
- **Backward compatibility** with existing traditional workflows.
- **Dummy scripts included** for GNINA and DiffDock to enable full workflow testing in any environment.
- **Comprehensive error handling** and resource management.

---

## License

**All Rights Reserved**

This software and its documentation are proprietary and confidential. All rights, title, and interest in and to this software, including all intellectual property rights, are and will remain the exclusive property of the original author.

**No Rights Granted:**
- No permission is granted to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of this software.
- No permission is granted to use this software for any commercial or non-commercial purpose.
- No permission is granted to create derivative works based on this software.

**Unauthorized Use Prohibited:**
Any unauthorized use, reproduction, or distribution of this software is strictly prohibited and may result in legal action.

For permission to use this software, contact the original author directly.
