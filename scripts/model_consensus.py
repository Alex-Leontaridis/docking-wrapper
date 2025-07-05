#!/usr/bin/env python3
"""
Consensus Engine: Pose Clustering and Consensus Scoring

Aligns and clusters predicted poses for a ligand, computes consensus score based on RMSD clustering.

Inputs:
  --poses: List of PDB files (predicted poses from different models)
  --output: Output JSON file for consensus results
  --ligand_id: Ligand identifier
  --rmsd_threshold: RMSD threshold for clustering (default: 2.0 Å)

Outputs:
  JSON with consensus_score, cluster assignments, mean RMSD, and breakdown.
"""

import argparse
import os
import sys
import json
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from sklearn.cluster import DBSCAN
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Pose consensus scoring via RMSD clustering")
    parser.add_argument("--poses", nargs='+', required=True, help="List of PDB files (poses)")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--ligand_id", required=True, help="Ligand identifier")
    parser.add_argument("--rmsd_threshold", type=float, default=2.0, help="RMSD threshold for clustering (Å)")
    return parser.parse_args()

def get_ligand_atoms(structure):
    """Extract ligand atoms (HETATM, not water) from a PDB structure."""
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == 'H_' and residue.resname not in ('HOH', 'WAT'):
                    for atom in residue:
                        atoms.append(atom)
    return atoms

def load_ligand_coords(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('lig', pdb_file)
    atoms = get_ligand_atoms(structure)
    if not atoms:
        raise ValueError(f"No ligand atoms found in {pdb_file}")
    coords = np.array([atom.get_coord() for atom in atoms])
    return coords

def compute_pairwise_rmsd(coord_list):
    """Compute pairwise RMSD matrix for a list of coordinate arrays."""
    n = len(coord_list)
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            # Align and compute RMSD
            si = Superimposer()
            # Use min length (in case of atom mismatch)
            min_len = min(len(coord_list[i]), len(coord_list[j]))
            si.set_atoms(coord_list[i][:min_len], coord_list[j][:min_len])
            rmsd = si.rms
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd
    return rmsd_matrix

def cluster_poses(rmsd_matrix, threshold):
    """Cluster poses using DBSCAN with RMSD threshold."""
    # DBSCAN expects a distance matrix; set eps=threshold
    clustering = DBSCAN(eps=threshold, min_samples=1, metric='precomputed')
    labels = clustering.fit_predict(rmsd_matrix)
    return labels

def consensus_score_from_labels(labels):
    """Consensus score: fraction of poses in largest cluster."""
    counts = np.bincount(labels)
    largest = counts.max()
    return round(largest / len(labels), 3)

def main():
    args = parse_args()
    pose_files = args.poses
    ligand_id = args.ligand_id
    rmsd_threshold = args.rmsd_threshold
    output_json = args.output

    # Validate input files
    for f in pose_files:
        if not os.path.isfile(f):
            print(f"[ERROR] Pose file not found: {f}")
            sys.exit(1)

    # Load ligand coordinates
    coords_list = []
    failed_files = []
    for f in pose_files:
        try:
            coords = load_ligand_coords(f)
            coords_list.append(coords)
        except Exception as e:
            print(f"[WARNING] Failed to load ligand from {f}: {e}")
            failed_files.append(f)

    if len(coords_list) < 2:
        print("[ERROR] Need at least 2 valid poses for consensus scoring.")
        sys.exit(1)

    # Compute pairwise RMSD
    rmsd_matrix = compute_pairwise_rmsd(coords_list)
    mean_rmsd = float(np.mean(rmsd_matrix[np.triu_indices(len(coords_list), 1)]))

    # Cluster poses
    labels = cluster_poses(rmsd_matrix, rmsd_threshold)
    consensus_score = consensus_score_from_labels(labels)

    # Organize output
    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[int(label)].append(os.path.basename(pose_files[idx]))

    result = {
        "ligand_id": ligand_id,
        "pose_files": [os.path.basename(f) for f in pose_files],
        "failed_files": failed_files,
        "rmsd_threshold": rmsd_threshold,
        "consensus_score": consensus_score,
        "mean_rmsd": round(mean_rmsd, 3),
        "num_clusters": int(len(set(labels))),
        "cluster_assignments": {str(k): v for k, v in clusters.items()},
        "labels": labels.tolist(),
        "notes": "Consensus score = fraction of poses in largest cluster (RMSD ≤ threshold)"
    }

    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"[SUCCESS] Consensus scoring complete: {output_json}")
    print(f"[INFO] Consensus score: {consensus_score}")
    print(f"[INFO] Mean RMSD: {mean_rmsd:.3f} Å")
    print(f"[INFO] Clusters: {result['cluster_assignments']}")

if __name__ == "__main__":
    main() 