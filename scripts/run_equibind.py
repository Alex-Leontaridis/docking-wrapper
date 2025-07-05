import argparse
import os
import subprocess
import sys
import time
import shutil
import json
import tempfile
from pathlib import Path

# Load configuration
def load_config():
    config_path = "tools_config.json"
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return json.load(f)
    return {}

EQUIBIND_REPO = "EquiBind"
INFERENCE_SCRIPT = "inference.py"
CONFIG_PATH = "configs_clean/inference.yml"

def check_equibind_installed():
    """Check if EquiBind is properly installed"""
    config = load_config()
    equibind_config = config.get("tools", {}).get("equibind", {})
    repo_path = equibind_config.get("repo_path", EQUIBIND_REPO)
    
    # Check if repository exists
    if not os.path.exists(repo_path):
        return False, f"EquiBind repository not found at {repo_path}"
    
    # Check if inference script exists
    inference_path = os.path.join(repo_path, INFERENCE_SCRIPT)
    if not os.path.exists(inference_path):
        return False, f"EquiBind inference script not found at {inference_path}"
    
    # Check if conda environment exists
    conda_env = equibind_config.get("conda_env", "equibind")
    try:
        result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
        if conda_env not in result.stdout:
            return False, f"EquiBind conda environment '{conda_env}' not found"
    except:
        return False, "Conda not available"
    
    return True, inference_path

def prepare_equibind_input(protein_path, ligand_path, output_dir):
    """Prepare input files for EquiBind"""
    # EquiBind expects specific input format
    # Copy files to output directory with proper names
    protein_dest = os.path.join(output_dir, "protein.pdb")
    ligand_dest = os.path.join(output_dir, "ligand.sdf")
    
    shutil.copy2(protein_path, protein_dest)
    shutil.copy2(ligand_path, ligand_dest)
    
    return protein_dest, ligand_dest

def run_equibind_inference(protein_path, ligand_path, output_path, conda_env="equibind"):
    """Run EquiBind inference using conda environment"""
    config = load_config()
    equibind_config = config.get("tools", {}).get("equibind", {})
    repo_path = equibind_config.get("repo_path", EQUIBIND_REPO)
    conda_env = equibind_config.get("conda_env", conda_env)
    
    # Create temporary directory for EquiBind input/output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Prepare input files
        protein_dest, ligand_dest = prepare_equibind_input(protein_path, ligand_path, tmpdir)
        
        # Run EquiBind inference
        cmd = [
            "conda", "run", "-n", conda_env,
            "python", os.path.join(repo_path, INFERENCE_SCRIPT),
            "--protein", protein_dest,
            "--ligand", ligand_dest,
            "--out", tmpdir
        ]
        
        print(f"[EquiBind] Running: {' '.join(cmd)}")
        
        try:
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Find output files
            output_files = []
            for root, _, files in os.walk(tmpdir):
                for file in files:
                    if file.endswith('.pdb') or file.endswith('.sdf'):
                        output_files.append(os.path.join(root, file))
            
            if output_files:
                # Copy output to destination
                for output_file in output_files:
                    dest_file = os.path.join(output_path, os.path.basename(output_file))
                    shutil.copy2(output_file, dest_file)
                
                # Save log
                log_file = os.path.join(output_path, "equibind.log")
                with open(log_file, 'w') as f:
                    f.write(proc.stdout + "\n" + proc.stderr)
                
                return True, None
            else:
                return False, "EquiBind did not produce output files"
                
        except subprocess.CalledProcessError as e:
            return False, f"EquiBind failed: {e.stderr}"
        except Exception as e:
            return False, f"EquiBind error: {e}"

def main():
    parser = argparse.ArgumentParser(description="EquiBind pose prediction")
    parser.add_argument("--protein", required=True, help="Protein PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand SDF file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--conda_env", default="equibind", help="Conda environment name")
    args = parser.parse_args()
    
    # Check inputs
    if not os.path.isfile(args.protein):
        print(f"[ERROR] Protein file not found: {args.protein}")
        sys.exit(1)
    if not os.path.isfile(args.ligand):
        print(f"[ERROR] Ligand file not found: {args.ligand}")
        sys.exit(1)
    
    # Check EquiBind installation
    ok, result = check_equibind_installed()
    if not ok:
        print(f"[ERROR] {result}")
        print("[INFO] Install EquiBind with: python install_real_tools.py")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Run EquiBind
    start_time = time.time()
    success, error = run_equibind_inference(args.protein, args.ligand, args.output, args.conda_env)
    
    if success:
        elapsed = time.time() - start_time
        print(f"[SUCCESS] EquiBind completed in {elapsed:.2f}s")
        print(f"[SUCCESS] Output written to {args.output}")
    else:
        print(f"[ERROR] EquiBind failed: {error}")
        sys.exit(1)

if __name__ == "__main__":
    main() 