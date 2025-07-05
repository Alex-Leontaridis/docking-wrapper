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

---

## CLI Usage

### 1. Batch Pipeline
Orchestrates the full workflow: structure prep, docking, and parsing.

```bash
python3 scripts/batch_pipeline.py --protein receptor.pdb --ligands ligand_dir/
# With custom config
python3 scripts/batch_pipeline.py --protein receptor.pdb --ligands ligands/ --config pipeline_config.json
# Enable all engines
python3 scripts/batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-gnina --enable-diffdock
# Serial processing
python3 scripts/batch_pipeline.py --protein receptor.pdb --ligands ligands/ --serial
```

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

## Enhanced Model & Analysis Commands

The following commands use the new ML-based docking, scoring, and analysis scripts. These are not yet integrated into the main batch pipeline, but can be run individually or orchestrated via custom scripts.

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
- **Ready for review.**
- Dummy scripts included for GNINA and DiffDock to enable full workflow testing in any environment.

---

