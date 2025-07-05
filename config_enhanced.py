# config.py

import os

# Shared paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODELS_DIR = os.path.join(BASE_DIR, 'models')
OUTPUT_DIR = os.path.join(BASE_DIR, 'outputs')
DATA_DIR = os.path.join(BASE_DIR, 'data')

# Model paths (placeholders)
DIFFDOCK_PATH = os.path.join(MODELS_DIR, 'diffdock')
EQUIBIND_PATH = os.path.join(MODELS_DIR, 'equibind')
NEURALPLEXER2_PATH = os.path.join(MODELS_DIR, 'neuralplexer2')
UMOL_PATH = os.path.join(MODELS_DIR, 'umol')
BOLTZ2_PATH = os.path.join(MODELS_DIR, 'boltz2')

# Thresholds and constants
RMSD_THRESHOLD = 2.0  # Angstroms
DRUGGABILITY_THRESHOLD = 0.5
CONFIDENCE_THRESHOLD = 0.7

# Output files
FINAL_SUMMARY_CSV = os.path.join(OUTPUT_DIR, 'final_summary.csv')
SUMMARY_JSON_DIR = os.path.join(OUTPUT_DIR, 'summaries')

# Add more as needed for new models or outputs 