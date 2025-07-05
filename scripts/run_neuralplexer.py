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

NEURALPLEXER_REPO = "NeuralPLexer"
INFERENCE_CLI = "neuralplexer-inference"
MODEL_CKPT = "models/complex_structure_prediction.ckpt"

def check_neuralplexer_installed():
    """Check if NeuralPLexer is properly installed"""
    config = load_config()
    neuralplexer_config = config.get("tools", {}).get("neuralplexer", {})
    repo_path = neuralplexer_config.get("repo_path", NEURALPLEXER_REPO)
    
    # Check if repository exists
    if not os.path.exists(repo_path):
        return False, f"NeuralPLexer repository not found at {repo_path}"
    
    # Check for installed CLI
    cli_path = shutil.which(INFERENCE_CLI)
    if cli_path:
        return True, cli_path
    
    # Check for inference script in repository
    inference_script = os.path.join(repo_path, "inference.py")
    if os.path.exists(inference_script):
        return True, inference_script
    
    # Check for makefile (indicates proper installation)
    makefile_path = os.path.join(repo_path, "Makefile")
    if os.path.exists(makefile_path):
        return True, "makefile"
    
    return False, "NeuralPLexer not properly installed"

def prepare_neuralplexer_input(protein_path, ligand_path, output_dir):
    """Prepare input files for NeuralPLexer"""
    # NeuralPLexer expects specific input format
    # Copy files to output directory with proper names
    protein_dest = os.path.join(output_dir, "protein.pdb")
    ligand_dest = os.path.join(output_dir, "ligand.sdf")
    
    shutil.copy2(protein_path, protein_dest)
    shutil.copy2(ligand_path, ligand_dest)
    
    return protein_dest, ligand_dest

def run_neuralplexer_inference(protein_path, ligand_path, output_path):
    """Run NeuralPLexer inference"""
    config = load_config()
    neuralplexer_config = config.get("tools", {}).get("neuralplexer", {})
    repo_path = neuralplexer_config.get("repo_path", NEURALPLEXER_REPO)
    
    # Check installation
    ok, result = check_neuralplexer_installed()
    if not ok:
        return False, result
    
    # Create temporary directory for NeuralPLexer input/output
    with tempfile.TemporaryDirectory() as tmpdir:
        # Prepare input files
        protein_dest, ligand_dest = prepare_neuralplexer_input(protein_path, ligand_path, tmpdir)
        
        # Run NeuralPLexer inference
        if result == "makefile":
            # Use make command
            cmd = [
                "make", "-C", repo_path,
                "inference",
                f"PROTEIN={protein_dest}",
                f"LIGAND={ligand_dest}",
                f"OUTPUT={tmpdir}"
            ]
        elif result.endswith('.py'):
            # Use Python script
            cmd = [
                sys.executable, result,
                "--protein", protein_dest,
                "--ligand", ligand_dest,
                "--output", tmpdir
            ]
        else:
            # Use CLI
            cmd = [
                result,
                "--protein", protein_dest,
                "--ligand", ligand_dest,
                "--output", tmpdir
            ]
        
        print(f"[NeuralPLexer] Running: {' '.join(cmd)}")
        
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
                log_file = os.path.join(output_path, "neuralplexer.log")
                with open(log_file, 'w') as f:
                    f.write(proc.stdout + "\n" + proc.stderr)
                
                return True, None
            else:
                return False, "NeuralPLexer did not produce output files"
                
        except subprocess.CalledProcessError as e:
            return False, f"NeuralPLexer failed: {e.stderr}"
        except Exception as e:
            return False, f"NeuralPLexer error: {e}"

def main():
    parser = argparse.ArgumentParser(description="NeuralPLexer pose prediction")
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
    
    # Check NeuralPLexer installation
    ok, result = check_neuralplexer_installed()
    if not ok:
        print(f"[ERROR] {result}")
        print("[INFO] Install NeuralPLexer with: python install_real_tools.py")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Run NeuralPLexer
    start_time = time.time()
    success, error = run_neuralplexer_inference(args.protein, args.ligand, args.output)
    
    if success:
        elapsed = time.time() - start_time
        print(f"[SUCCESS] NeuralPLexer completed in {elapsed:.2f}s")
        print(f"[SUCCESS] Output written to {args.output}")
    else:
        print(f"[ERROR] NeuralPLexer failed: {error}")
        sys.exit(1)

if __name__ == "__main__":
    main() 