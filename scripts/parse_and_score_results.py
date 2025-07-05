#!/usr/bin/env python3
"""
Output Formatter: Consolidate Model Outputs

Merges all model outputs (affinity, druggability, consensus, confidence, interactions, etc.)
into a structured CSV and JSON summary for each ligand.

Inputs:
  --input_dir: Directory containing all per-ligand JSON outputs
  --output_csv: Path to final summary CSV
  --output_json: Path to final summary JSON
  --ligand_ids: Optional, comma-separated list of ligand IDs to include

Outputs:
  - final_summary.csv: One row per ligand, all key metrics
  - summary.json: Per-ligand full details
"""

import argparse
import os
import sys
import json
import glob
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Consolidate model outputs into summary CSV and JSON")
    parser.add_argument("--input_dir", required=True, help="Directory with per-ligand JSON outputs")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_json", required=True, help="Output JSON file")
    parser.add_argument("--ligand_ids", default=None, help="Comma-separated list of ligand IDs to include (optional)")
    return parser.parse_args()

def find_json_files(input_dir, pattern):
    return glob.glob(os.path.join(input_dir, pattern))

def load_json(path):
    if not os.path.isfile(path):
        return None
    with open(path, 'r') as f:
        return json.load(f)

def safe_get(d, *keys, default="N/A"):
    for k in keys:
        if d is None or k not in d:
            return default
        d = d[k]
    return d

def collect_ligand_data(input_dir, ligand_ids=None):
    """
    For each ligand, collect all relevant JSON outputs and merge into a summary dict.
    """
    # Find all confidence JSONs as the anchor
    confidence_files = find_json_files(input_dir, "*confidence.json")
    summary = {}
    for conf_file in confidence_files:
        conf = load_json(conf_file)
        ligand_id = safe_get(conf, "ligand_id", default=None)
        if not ligand_id:
            continue
        if ligand_ids and ligand_id not in ligand_ids:
            continue
        # Try to find other files by ligand_id
        base = ligand_id
        consensus = load_json(os.path.join(input_dir, f"{base}_consensus.json"))
        druggability = load_json(os.path.join(input_dir, f"{base}_druggability.json"))
        affinity = load_json(os.path.join(input_dir, f"{base}_affinity.json"))
        interaction = load_json(os.path.join(input_dir, f"{base}_interactions.json"))
        # Compose summary
        summary[ligand_id] = {
            "ligand_id": ligand_id,
            "confidence_score": safe_get(conf, "confidence_score"),
            "consensus_score": safe_get(consensus, "consensus_score"),
            "druggability_score": safe_get(druggability, "druggability_score"),
            "affinity_kcal_per_mol": safe_get(affinity, "delta_g"),
            "LE": safe_get(conf, "breakdown", "ligand_efficiency"),
            "model_voting": safe_get(conf, "breakdown", "model_voting"),
            "mean_rmsd": safe_get(consensus, "mean_rmsd"),
            "num_clusters": safe_get(consensus, "num_clusters"),
            "h_bonds": safe_get(interaction, "H-bond"),
            "pi_pi": safe_get(interaction, "pi-pi"),
            "vdW": safe_get(interaction, "vdW"),
            "structure_source": safe_get(affinity, "protein"),
            "summary": {
                "confidence": conf,
                "consensus": consensus,
                "druggability": druggability,
                "affinity": affinity,
                "interaction": interaction
            }
        }
    return summary

def write_csv(summary, output_csv):
    # Define CSV columns
    columns = [
        "ligand_id", "confidence_score", "consensus_score", "druggability_score", "affinity_kcal_per_mol",
        "LE", "model_voting", "mean_rmsd", "num_clusters", "h_bonds", "pi_pi", "vdW", "structure_source"
    ]
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for ligand_id, row in summary.items():
            # Flatten LE if it's a dict
            if isinstance(row["LE"], dict):
                row["LE"] = row["LE"].get("LE", "N/A")
            # Convert lists to string
            for k in ["h_bonds", "pi_pi", "vdW"]:
                if isinstance(row[k], list):
                    row[k] = ";".join(row[k])
            writer.writerow({k: row.get(k, "N/A") for k in columns})

def write_json(summary, output_json):
    with open(output_json, 'w') as f:
        json.dump(summary, f, indent=2)

def main():
    args = parse_args()
    ligand_ids = None
    if args.ligand_ids:
        ligand_ids = [x.strip() for x in args.ligand_ids.split(",") if x.strip()]
    summary = collect_ligand_data(args.input_dir, ligand_ids)
    if not summary:
        print("[ERROR] No ligand data found.")
        sys.exit(1)
    write_csv(summary, args.output_csv)
    write_json(summary, args.output_json)
    print(f"[SUCCESS] Final summary written to {args.output_csv} and {args.output_json}")
    print(f"[INFO] Ligands summarized: {list(summary.keys())}")

if __name__ == "__main__":
    main() 