# Enhanced Docking Wrapper - Decision-Ready Insight Engine

This enhanced version of the docking-wrapper integrates multiple state-of-the-art ML models for comprehensive protein-ligand docking analysis, transforming the system into a decision-ready insight engine.

## 🚀 New Features

### 1. **Multiple Docking Methods**
- **DiffDock**: Already implemented (baseline)
- **EquiBind**: Structure-based pose prediction
- **NeuralPLexer-2**: High-accuracy neural network-based docking
- **UMol**: Structure-free pose prediction (when no structure available)

### 2. **Binding Affinity Estimation**
- **Boltz2**: State-of-the-art binding energy prediction
- Fallback mechanisms for robustness

### 3. **Structure Prediction**
- **ColabFold**: AlphaFold2-based structure prediction
- **OpenFold**: Open-source AlphaFold2 implementation
- **ESMFold**: Meta's protein structure prediction
- Intelligent caching system

### 4. **Interaction Analysis**
- **PLIP**: Protein-ligand interaction profiler
- **RDKit**: Chemical informatics fallback
- H-bonds, π-π interactions, van der Waals contacts

### 5. **Druggability Scoring**
- **fpocket**: Binding site druggability analysis
- Multi-parameter scoring algorithm (0-1 scale)

### 6. **Consensus & Confidence Engine**
- RMSD-based clustering for pose agreement
- Composite scoring with weighted confidence metrics
- Ligand efficiency calculations (LE, SILE, LLE)

### 7. **Enhanced Output Format**
- CSV summary with comprehensive metrics
- JSON details for downstream processing
- Standardized format across all outputs

## 📁 Enhanced Scripts

### New Docking Scripts
- `scripts/run_equibind.py` - EquiBind pose prediction
- `scripts/run_neuralplexer.py` - NeuralPLexer-2 docking
- `scripts/run_umol.py` - UMol structure-free docking
- `scripts/run_structure_predictor.py` - Protein structure prediction

### Analysis Scripts
- `scripts/run_boltz2.py` - Binding affinity prediction
- `scripts/extract_interactions.py` - Protein-ligand interaction analysis
- `scripts/run_druggability.py` - Binding site druggability scoring

### Consensus & Confidence
- `scripts/model_consensus.py` - Compare poses across methods
- `scripts/compute_confidence.py` - Generate confidence scores (0-100)

### Enhanced Output
- `scripts/parse_and_score_results.py` - Updated with new metrics

## 🛠️ Installation

### 1. Install Python Dependencies
```bash
pip install -r requirements.txt
```

### 2. Install External Tools (Optional)
```bash
# Ubuntu/Debian
sudo apt-get install fpocket

# macOS
brew install fpocket

# Conda
conda install -c conda-forge fpocket
```

### 3. Model Setup
Each ML model may require additional setup:
- **EquiBind**: Follow official installation guide
- **NeuralPLexer**: Install from source
- **UMol**: Download pre-trained models
- **Boltz2**: Install via pip or conda

## 🚀 Usage

### Individual Scripts

#### 1. Pose Prediction
```bash
# EquiBind
python scripts/run_equibind.py --protein inputs/protein.pdb --ligand inputs/ligands/ligand.sdf --output outputs/equibind_pose.pdb

# NeuralPLexer
python scripts/run_neuralplexer.py --protein inputs/protein.pdb --ligand inputs/ligands/ligand.sdf --output outputs/neuralplexer_pose.pdb

# UMol (structure-free)
python scripts/run_umol.py --protein inputs/protein.fasta --ligand inputs/ligands/ligand.sdf --output outputs/umol_pose.pdb
```

#### 2. Structure Prediction
```bash
python scripts/run_structure_predictor.py --protein inputs/protein.fasta --output_dir outputs/structures/
```

#### 3. Binding Affinity
```bash
python scripts/run_boltz2.py --protein inputs/protein.fasta --ligand inputs/ligands/ligand.smi --output outputs/affinity.json
```

#### 4. Interaction Analysis
```bash
python scripts/extract_interactions.py --pdb outputs/pose.pdb --ligand LIG --protein PROT --output_dir outputs/interactions/
```

#### 5. Druggability Scoring
```bash
python scripts/run_druggability.py --protein outputs/pose.pdb --output outputs/druggability.json
```

#### 6. Consensus Analysis
```bash
python scripts/model_consensus.py --poses outputs/equibind_pose.pdb outputs/neuralplexer_pose.pdb outputs/umol_pose.pdb --output outputs/consensus.json --ligand_id LIG1
```

#### 7. Confidence Scoring
```bash
python scripts/compute_confidence.py --consensus_json outputs/consensus.json --druggability_json outputs/druggability.json --affinity_json outputs/affinity.json --interaction_json outputs/interactions/LIG_PROT.json --output outputs/confidence.json --ligand_id LIG1
```

### Complete Enhanced Pipeline

```bash
# 1. Structure prediction (if needed)
python scripts/run_structure_predictor.py --protein inputs/protein.fasta --output_dir outputs/structures/

# 2. Multiple pose predictions
python scripts/run_equibind.py --protein outputs/structures/protein.pdb --ligand inputs/ligands/ligand.sdf --output outputs/equibind_pose.pdb
python scripts/run_neuralplexer.py --protein outputs/structures/protein.pdb --ligand inputs/ligands/ligand.sdf --output outputs/neuralplexer_pose.pdb
python scripts/run_umol.py --protein outputs/structures/protein.pdb --ligand inputs/ligands/ligand.sdf --output outputs/umol_pose.pdb

# 3. Binding affinity
python scripts/run_boltz2.py --protein inputs/protein.fasta --ligand inputs/ligands/ligand.smi --output outputs/affinity.json

# 4. Interaction analysis
python scripts/extract_interactions.py --pdb outputs/equibind_pose.pdb --ligand LIG --protein PROT --output_dir outputs/interactions/

# 5. Druggability scoring
python scripts/run_druggability.py --protein outputs/equibind_pose.pdb --output outputs/druggability.json

# 6. Consensus analysis
python scripts/model_consensus.py --poses outputs/equibind_pose.pdb outputs/neuralplexer_pose.pdb outputs/umol_pose.pdb --output outputs/consensus.json --ligand_id LIG1

# 7. Confidence scoring
python scripts/compute_confidence.py --consensus_json outputs/consensus.json --druggability_json outputs/druggability.json --affinity_json outputs/affinity.json --interaction_json outputs/interactions/LIG_PROT.json --output outputs/confidence.json --ligand_id LIG1

# 8. Final summary
python scripts/parse_and_score_results.py --input_dir outputs/ --output_csv outputs/final_summary.csv --output_json outputs/summary.json
```

## 📊 Enhanced Output Formats

### CSV Summary (`final_summary.csv`)
```csv
ligand_id,confidence_score,consensus_score,druggability_score,affinity_kcal_per_mol,LE,model_voting,mean_rmsd,num_clusters,h_bonds,pi_pi,vdW,structure_source
LIG1,85,0.8,0.75,-8.5,0.34,0.8,1.2,2,"ASP155;ARG198","TRP156","LEU83;VAL91",inputs/protein.fasta
```

### JSON Summary (`summary.json`)
```json
{
  "LIG1": {
    "ligand_id": "LIG1",
    "confidence_score": 85,
    "consensus_score": 0.8,
    "druggability_score": 0.75,
    "affinity_kcal_per_mol": -8.5,
    "LE": 0.34,
    "model_voting": 0.8,
    "mean_rmsd": 1.2,
    "num_clusters": 2,
    "h_bonds": ["ASP155", "ARG198"],
    "pi_pi": ["TRP156"],
    "vdW": ["LEU83", "VAL91"],
    "structure_source": "inputs/protein.fasta",
    "summary": {
      "confidence": {...},
      "consensus": {...},
      "druggability": {...},
      "affinity": {...},
      "interaction": {...}
    }
  }
}
```

## 🔧 Configuration

The system uses `config_enhanced.py` for configuration:
- Model paths and thresholds
- Output directories
- RMSD and confidence thresholds

## 🧪 Testing

Test the enhanced system with the provided example files:
```bash
# Test individual components
python scripts/run_equibind.py --protein inputs/protein.pdb --ligand inputs/ligands/aspirin.sdf --output test_outputs/equibind_test.pdb

# Test complete pipeline
python scripts/batch_pipeline.py --config pipeline_config.json
```

## 📈 Key Metrics

### Confidence Score (0-100)
- **Pose convergence**: RMSD ≤ 2Å
- **Druggability**: ≥ 0.75
- **Ligand efficiency**: LE, SILE, LLE thresholds
- **Model agreement**: % of models agreeing on top pose

### Consensus Score (0-1)
- RMSD-based clustering
- Pose agreement analysis
- Weighted voting across methods

### Druggability Score (0-1)
- Pocket volume and quality
- Multi-parameter algorithm
- fpocket-based analysis

## 🐛 Troubleshooting

### Common Issues
1. **External tool not found**: Scripts include fallback mechanisms
2. **Memory issues**: Reduce batch sizes, use CPU-only mode
3. **File format errors**: Ensure correct input formats

### Getting Help
- Check script help: `python scripts/script_name.py --help`
- Review log files in `logs/` directory
- Ensure all dependencies are installed

## 🔄 Migration from Original System

The enhanced system is backward compatible with your existing:
- `batch_pipeline.py` - Enhanced with new metrics
- `parse_and_score_results.py` - Updated output format
- Configuration files - Extended with new options

## 📚 Acknowledgments

- **EquiBind**: Authors of the EquiBind model
- **NeuralPLexer**: NeuralPLexer development team
- **UMol**: UMol research group
- **Boltz2**: Boltz2 developers
- **ColabFold**: ColabFold team
- **OpenFold**: OpenFold contributors
- **ESMFold**: Meta AI Research
- **fpocket**: fpocket development team
- **PLIP**: PLIP authors
- **RDKit**: RDKit contributors

## 📄 License

This project is licensed under the MIT License. 