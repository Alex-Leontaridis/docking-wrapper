import argparse
import os
import sys
import json
import time
import tempfile
import subprocess
from pathlib import Path

# Load configuration
def load_config():
    config_path = "tools_config.json"
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return json.load(f)
    return {}

# Try to import real Boltz2
try:
    from boltz import predict
    BOLTZ_AVAILABLE = True
except ImportError:
    BOLTZ_AVAILABLE = False
    print("[INFO] Boltz2 not available, will try to install or use fallback")

def check_boltz2_installed():
    """Check if Boltz2 is properly installed"""
    config = load_config()
    boltz_config = config.get("boltz2", {})
    repo_path = boltz_config.get("repo_path", "boltz")
    
    # Check if Boltz2 is importable
    if BOLTZ_AVAILABLE:
        return True, "imported"
    
    # Check if repository exists (official path)
    if os.path.exists(repo_path):
        return True, repo_path
    
    # Check if old repository path exists
    if os.path.exists("Boltz2"):
        return True, "Boltz2"
    
    return False, "Boltz2 not installed"

def install_boltz2():
    """Try to install Boltz2 using official methods"""
    print("[INFO] Attempting to install Boltz2 using official methods...")
    
    # Method 1: PyPI installation (recommended)
    print("[INFO] Trying PyPI installation...")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "boltz", "-U"], 
                      check=True, capture_output=True)
        print("[SUCCESS] Boltz2 installed via PyPI")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] PyPI installation failed: {e}")
    
    # Method 2: GitHub installation
    print("[INFO] Trying GitHub installation...")
    try:
        subprocess.run(["git", "clone", "https://github.com/jwohlwend/boltz.git", "boltz"], 
                      check=True, capture_output=True)
        os.chdir("boltz")
        subprocess.run([sys.executable, "-m", "pip", "install", "-e", "."], 
                      check=True, capture_output=True)
        os.chdir("..")
        print("[SUCCESS] Boltz2 installed from GitHub")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[WARNING] GitHub installation failed: {e}")
        if os.path.exists("boltz"):
            os.chdir("..")  # Make sure we're back in the right directory
    
    print("[ERROR] All installation methods failed")
    print("[INFO] Consider using Docker: https://hub.docker.com/r/boltz/boltz")
    return False

def mock_predict(yaml_path):
    """Mock Boltz2 prediction for testing"""
    import random
    # Simulate realistic binding affinity predictions
    # Typical ΔG values range from -15 to +5 kcal/mol
    delta_g = random.uniform(-12.0, -2.0)
    confidence = random.uniform(0.6, 0.95)
    
    return {
        "affinity_pred_value": round(delta_g, 2),
        "affinity_probability_binary": round(confidence, 3)
    }

def check_inputs(protein, ligand):
    if not os.path.isfile(protein):
        print(f"[ERROR] Protein file not found: {protein}")
        sys.exit(1)
    if not os.path.isfile(ligand):
        print(f"[ERROR] Ligand file not found: {ligand}")
        sys.exit(1)
    if not (protein.endswith('.fasta') or protein.endswith('.fa')):
        print(f"[ERROR] Protein must be a FASTA file (.fasta or .fa): {protein}")
        sys.exit(1)
    if not ligand.endswith('.smi'):
        print(f"[ERROR] Ligand must be a SMILES file (.smi): {ligand}")
        sys.exit(1)

def make_yaml(protein, ligand, yaml_path):
    yaml_content = f"""protein:
  fasta_path: {os.path.abspath(protein)}
ligand:
  smiles_path: {os.path.abspath(ligand)}
properties:
  - affinity
"""
    with open(yaml_path, "w") as f:
        f.write(yaml_content)

def run_boltz2_prediction(yaml_path):
    """Run Boltz2 prediction"""
    # Check if Boltz2 is available
    if BOLTZ_AVAILABLE:
        try:
            return predict(yaml_path), False
        except Exception as e:
            print(f"[WARNING] Boltz2 prediction failed: {e}")
            return mock_predict(yaml_path), True
    else:
        # Try to install Boltz2
        if install_boltz2():
            try:
                from boltz import predict
                return predict(yaml_path), False
            except ImportError:
                pass
        
        # Use mock as fallback
        return mock_predict(yaml_path), True

def main():
    parser = argparse.ArgumentParser(description="Binding affinity estimation using Boltz2")
    parser.add_argument("--protein", required=True, help="Protein FASTA file (.fasta)")
    parser.add_argument("--ligand", required=True, help="Ligand SMILES file (.smi)")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--ligand_id", default="LIG", help="Ligand ID")
    parser.add_argument("--protein_id", default="PROT", help="Protein ID")
    parser.add_argument("--force_mock", action="store_true", help="Force use of mock prediction")
    args = parser.parse_args()

    check_inputs(args.protein, args.ligand)
    
    # Check Boltz2 installation
    ok, result = check_boltz2_installed()
    if not ok and not args.force_mock:
        print(f"[ERROR] {result}")
        print("[INFO] Install Boltz2 with: python install_real_tools.py")
        print("[INFO] Or use --force_mock for testing")
        sys.exit(1)
    
    start = time.time()
    result = {
        "protein": args.protein,
        "ligand": args.ligand,
        "protein_id": args.protein_id,
        "ligand_id": args.ligand_id,
        "model": "Boltz2",
        "fallback": False,
    }
    
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            yaml_path = os.path.join(tmpdir, "input.yaml")
            make_yaml(args.protein, args.ligand, yaml_path)
            
            if args.force_mock:
                pred = mock_predict(yaml_path)
                result["fallback"] = True
                result["mock_prediction"] = True
            else:
                pred, is_fallback = run_boltz2_prediction(yaml_path)
                if is_fallback:
                    result["fallback"] = True
                    result["mock_prediction"] = True
            
            # Parse output (see Boltz2 docs for exact keys)
            affinity = pred.get("affinity_pred_value", None)
            confidence = pred.get("affinity_probability_binary", None)
            result["delta_g"] = affinity
            result["confidence"] = confidence
            
    except Exception as e:
        result["delta_g"] = None
        result["confidence"] = None
        result["fallback"] = True
        result["error"] = str(e)
    
    elapsed = time.time() - start
    result["runtime_sec"] = round(elapsed, 2)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    with open(args.output, "w") as f:
        json.dump(result, f, indent=2)
    
    if result.get("fallback", False):
        print(f"[Boltz2] Mock prediction completed in {elapsed:.2f}s")
        print(f"[Boltz2] Output written to {args.output} (MOCK)")
    else:
        print(f"[Boltz2] Real prediction completed in {elapsed:.2f}s")
        print(f"[Boltz2] Output written to {args.output}")

if __name__ == "__main__":
    main() 