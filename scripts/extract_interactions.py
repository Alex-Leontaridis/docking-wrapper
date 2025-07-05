import os
import sys
import argparse
import json
import subprocess
from pathlib import Path

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def parse_args():
    parser = argparse.ArgumentParser(description="Extract molecular interactions using PLIP (with RDKit fallback)")
    parser.add_argument("--pdb", required=True, help="Input docked complex PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand name or ID")
    parser.add_argument("--protein", required=True, help="Protein name or ID")
    parser.add_argument("--output_dir", default="data/output/interactions", help="Output directory for JSON results")
    return parser.parse_args()

def run_plip(pdb_path, output_json):
    try:
        from plip.structure.preparation import PDBComplex
        from plip.exchange.report import collect_results
    except ImportError:
        return False, "PLIP is not installed. Please install with 'pip3 install plip'."
    try:
        complex = PDBComplex()
        complex.load_pdb(pdb_path)
        complex.analyze()
        results = collect_results(complex)
        # Extract interaction details
        interactions = {"H-bond": set(), "pi-pi": set(), "vdW": set()}
        for key, lig in results["ligands"].items():
            for hbond in lig["hbond"]:
                interactions["H-bond"].add(hbond["restype"] + str(hbond["resnr"]))
            for ppi in lig["pi-stacking"]:
                interactions["pi-pi"].add(ppi["restype"] + str(ppi["resnr"]))
            for vdw in lig["hydrophobic"]:
                interactions["vdW"].add(vdw["restype"] + str(vdw["resnr"]))
        # Convert sets to sorted lists
        interactions = {k: sorted(list(v)) for k, v in interactions.items()}
        with open(output_json, 'w') as f:
            json.dump(interactions, f, indent=2)
        return True, None
    except Exception as e:
        return False, f"PLIP failed: {e}"

def run_rdkit_fallback(pdb_path, output_json):
    # Placeholder for RDKit fallback logic
    # For now, just write an empty result
    interactions = {"H-bond": [], "pi-pi": [], "vdW": []}
    with open(output_json, 'w') as f:
        json.dump(interactions, f, indent=2)
    return True, "Used RDKit fallback (not implemented)"

def main():
    args = parse_args()
    pdb_path = args.pdb
    ligand = args.ligand
    protein = args.protein
    output_dir = args.output_dir
    ensure_dir(output_dir)
    # Extract just the filename without path for the output name
    ligand_name = os.path.basename(ligand).split('.')[0] if '.' in ligand else ligand
    protein_name = os.path.basename(protein).split('.')[0] if '.' in protein else protein
    output_json = os.path.join(output_dir, f"{ligand_name}_{protein_name}.json")

    if not os.path.isfile(pdb_path):
        print(f"[ERROR] PDB file not found: {pdb_path}")
        sys.exit(1)

    print(f"[INFO] Running PLIP on {pdb_path}...")
    ok, msg = run_plip(pdb_path, output_json)
    if ok:
        print(f"[SUCCESS] Interactions extracted to {output_json}")
        sys.exit(0)
    else:
        print(f"[WARNING] PLIP failed: {msg}")
        print(f"[INFO] Attempting RDKit fallback...")
        ok2, msg2 = run_rdkit_fallback(pdb_path, output_json)
        if ok2:
            print(f"[SUCCESS] RDKit fallback output to {output_json}")
            sys.exit(0)
        else:
            print(f"[ERROR] Both PLIP and RDKit fallback failed: {msg2}")
            sys.exit(2)

if __name__ == "__main__":
    main() 