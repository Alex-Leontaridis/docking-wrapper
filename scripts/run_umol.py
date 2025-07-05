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

UMOL_REPO = "Umol"
PREDICT_SCRIPT = "predict.sh"

def check_umol_installed():
    """Check if UMol is properly installed"""
    config = load_config()
    umol_config = config.get("tools", {}).get("umol", {})
    repo_path = umol_config.get("repo_path", UMOL_REPO)
    
    # Check if repository exists
    if not os.path.exists(repo_path):
        return False, f"UMol repository not found at {repo_path}"
    
    # Check for predict script
    script_path = os.path.join(repo_path, PREDICT_SCRIPT)
    if os.path.isfile(script_path):
        return True, script_path
    
    # Check for Python inference script
    inference_script = os.path.join(repo_path, "inference.py")
    if os.path.isfile(inference_script):
        return True, inference_script
    
    # Check for main script
    main_script = os.path.join(repo_path, "main.py")
    if os.path.isfile(main_script):
        return True, main_script
    
    return False, "UMol not properly installed"

def prepare_umol_input(protein_path, ligand_path, output_dir):
    """Prepare input files for UMol"""
    # UMol expects specific input format
    # Copy files to output directory with proper names
    protein_dest = os.path.join(output_dir, "protein.pdb")
    ligand_dest = os.path.join(output_dir, "ligand.sdf")
    
    shutil.copy2(protein_path, protein_dest)
    shutil.copy2(ligand_path, ligand_dest)
    
    return protein_dest, ligand_dest

def run_umol_inference(protein_path, ligand_path, output_path):
    """Run UMol inference"""
    config = load_config()
    umol_config = config.get("tools", {}).get("umol", {})
    repo_path = umol_config.get("repo_path", UMOL_REPO)
    
    # Check installation
    ok, result = check_umol_installed()
    if not ok:
        return False, result
    
    # Create temporary directory for UMol input/output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Prepare input files
        protein_dest, ligand_dest = prepare_umol_input(protein_path, ligand_path, tmpdir)
        
        # Run UMol inference
        if result.endswith('.sh'):
            # Use bash script
            cmd = ["bash", result, protein_dest, ligand_dest, tmpdir]
        elif result.endswith('.py'):
            # Use Python script
            cmd = [
                sys.executable, result,
                "--protein", protein_dest,
                "--ligand", ligand_dest,
                "--output", tmpdir
            ]
        else:
            return False, f"Unknown UMol script type: {result}"
        
        print(f"[UMol] Running: {' '.join(cmd)}")
        
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
                log_file = os.path.join(output_path, "umol.log")
                with open(log_file, 'w') as f:
                    f.write(proc.stdout + "\n" + proc.stderr)
                
                return True, None
            else:
                return False, "UMol did not produce output files"
                
        except subprocess.CalledProcessError as e:
            return False, f"UMol failed: {e.stderr}"
        except Exception as e:
            return False, f"UMol error: {e}"

def main():
    parser = argparse.ArgumentParser(description="UMol pose prediction")
    parser.add_argument("--protein", required=True, help="Protein PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand SDF file")
    parser.add_argument("--output", required=True, help="Output directory")
    args = parser.parse_args()
    
    # Check inputs
    if not os.path.isfile(args.protein):
        print(f"[ERROR] Protein file not found: {args.protein}")
        sys.exit(1)
    if not os.path.isfile(args.ligand):
        print(f"[ERROR] Ligand file not found: {args.ligand}")
        sys.exit(1)
    
    # Check UMol installation
    ok, result = check_umol_installed()
    if not ok:
        print(f"[ERROR] {result}")
        print("[INFO] Install UMol with: python install_real_tools.py")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Run UMol
    start_time = time.time()
    success, error = run_umol_inference(args.protein, args.ligand, args.output)
    
    if success:
        elapsed = time.time() - start_time
        print(f"[SUCCESS] UMol completed in {elapsed:.2f}s")
        print(f"[SUCCESS] Output written to {args.output}")
    else:
        print(f"[ERROR] UMol failed: {error}")
        sys.exit(1)

if __name__ == "__main__":
    main() 