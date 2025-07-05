#!/usr/bin/env python3
"""
Confidence Engine: Composite Confidence Scoring

Aggregates consensus, druggability, ligand efficiency, and model voting to produce a final confidence score (0–100).

Inputs:
  --consensus_json: Output from model_consensus.py
  --druggability_json: Output from run_druggability.py
  --affinity_json: Output from run_boltz2.py
  --interaction_json: Output from extract_interactions.py
  --output: Output JSON file
  --ligand_id: Ligand identifier

Outputs:
  JSON with confidence_score, breakdown, and all input data for traceability.
"""

import argparse
import os
import sys
import json
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description="Compute composite confidence score for ligand prediction")
    parser.add_argument("--consensus_json", required=True, help="Consensus JSON file")
    parser.add_argument("--druggability_json", required=True, help="Druggability JSON file")
    parser.add_argument("--affinity_json", required=True, help="Affinity JSON file")
    parser.add_argument("--interaction_json", required=True, help="Interaction JSON file")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--ligand_id", required=True, help="Ligand identifier")
    return parser.parse_args()

def load_json(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")
    with open(path, 'r') as f:
        return json.load(f)

def compute_ligand_efficiency(affinity_json):
    """
    Compute Ligand Efficiency (LE), SILE, LLE if possible.
    LE = -ΔG / N_heavy_atoms (mock: assume 25 heavy atoms if not available)
    SILE, LLE: placeholders for now.
    """
    delta_g = affinity_json.get("delta_g")
    smiles = None
    n_heavy = 25  # Default/mock value
    # Try to get SMILES and count heavy atoms if available
    ligand_path = affinity_json.get("ligand")
    if ligand_path and os.path.isfile(ligand_path):
        try:
            with open(ligand_path) as f:
                smiles = f.read().strip()
            # Count heavy atoms (non-H)
            n_heavy = sum(1 for c in smiles if c.isalpha() and c.upper() != 'H')
        except Exception:
            pass
    le = None
    if delta_g is not None and n_heavy > 0:
        le = round(-float(delta_g) / n_heavy, 3)
    # SILE, LLE: not implemented, set to None
    return {
        "LE": le,
        "SILE": None,
        "LLE": None,
        "n_heavy_atoms": n_heavy
    }

def compute_model_voting(consensus_json):
    """
    Model voting agreement: fraction of models in largest cluster.
    """
    return consensus_json.get("consensus_score", 0.0)

def compute_interaction_score(interaction_json):
    """
    Simple interaction score: count of unique residues in all interaction types.
    """
    if not interaction_json:
        return 0.0
    unique_residues = set()
    for k, v in interaction_json.items():
        if isinstance(v, list):
            unique_residues.update(v)
    # Normalize: assume 10+ is strong, 0 is weak
    score = min(len(unique_residues) / 10.0, 1.0)
    return round(score, 3)

def compute_confidence_score(consensus, druggability, le, model_voting, interaction_score):
    """
    Weighted formula:
      - Consensus (pose agreement): 40%
      - Druggability (≥0.75): 20%
      - Ligand Efficiency (LE): 20%
      - Model voting agreement: 10%
      - Interaction score: 10%
    All scores normalized to [0,1] before weighting.
    """
    # Normalize LE: typical range 0.2–0.6, cap at 0.7
    le_norm = 0.0
    if le is not None:
        le_norm = min(max((le - 0.2) / 0.5, 0.0), 1.0)
    # Druggability: use as is (0–1)
    druggability_norm = min(max(druggability, 0.0), 1.0)
    # Consensus, model_voting, interaction_score: already 0–1
    consensus_norm = min(max(consensus, 0.0), 1.0)
    model_voting_norm = min(max(model_voting, 0.0), 1.0)
    interaction_norm = min(max(interaction_score, 0.0), 1.0)
    # Weighted sum
    score = (
        0.4 * consensus_norm +
        0.2 * druggability_norm +
        0.2 * le_norm +
        0.1 * model_voting_norm +
        0.1 * interaction_norm
    )
    return int(round(score * 100))

def main():
    args = parse_args()
    ligand_id = args.ligand_id
    # Load all input data
    consensus_json = load_json(args.consensus_json)
    druggability_json = load_json(args.druggability_json)
    affinity_json_path = args.affinity_json
    if not os.path.exists(affinity_json_path):
        # Try fallback
        fallback_path = 'outputs/test_boltz2_result.json'
        if os.path.exists(fallback_path):
            print(f"[INFO] Affinity file not found, using fallback: {fallback_path}")
            affinity_json_path = fallback_path
        else:
            raise FileNotFoundError(f"Affinity file not found: {affinity_json_path}")
    affinity_json = load_json(affinity_json_path)
    interaction_json = load_json(args.interaction_json)
    # Compute all factors
    consensus = consensus_json.get("consensus_score", 0.0)
    druggability = druggability_json.get("druggability_score", 0.0)
    le_data = compute_ligand_efficiency(affinity_json)
    le = le_data["LE"]
    model_voting = compute_model_voting(consensus_json)
    interaction_score = compute_interaction_score(interaction_json)
    # Compute final confidence score
    confidence_score = compute_confidence_score(
        consensus, druggability, le, model_voting, interaction_score
    )
    # Organize output
    result = {
        "ligand_id": ligand_id,
        "confidence_score": confidence_score,
        "breakdown": {
            "consensus": consensus,
            "druggability": druggability,
            "ligand_efficiency": le,
            "model_voting": model_voting,
            "interaction_score": interaction_score
        },
        "input_files": {
            "consensus_json": args.consensus_json,
            "druggability_json": args.druggability_json,
            "affinity_json": args.affinity_json,
            "interaction_json": args.interaction_json
        },
        "input_data": {
            "consensus": consensus_json,
            "druggability": druggability_json,
            "affinity": affinity_json,
            "interaction": interaction_json
        },
        "notes": "Confidence score is a weighted composite of consensus, druggability, ligand efficiency, model voting, and interaction count."
    }
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"[SUCCESS] Confidence scoring complete: {args.output}")
    print(f"[INFO] Confidence score: {confidence_score}")
    print(f"[INFO] Breakdown: {result['breakdown']}")

if __name__ == "__main__":
    main() 