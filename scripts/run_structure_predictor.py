import argparse
import os
import sys
import hashlib
import subprocess
import time
import json
import shutil
import tempfile
from pathlib import Path

# Load configuration
def load_config():
    config_path = "tools_config.json"
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            return json.load(f)
    return {}

# Helper: hash FASTA content for unique cache key
def fasta_hash(fasta_path):
    with open(fasta_path, 'rb') as f:
        content = f.read()
    return hashlib.sha256(content).hexdigest()[:16]

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def check_cached_pdb(cache_dir, protein_id, hash_id):
    pdb_path = os.path.join(cache_dir, f"{protein_id}_{hash_id}.pdb")
    log_path = os.path.join(cache_dir, f"{protein_id}_{hash_id}.log.json")
    if os.path.isfile(pdb_path):
        return pdb_path, log_path
    return None, None

def log_result(log_path, model, runtime, error=None):
    log = {
        "model_used": model,
        "runtime_sec": round(runtime, 2),
        "error": error
    }
    with open(log_path, 'w') as f:
        json.dump(log, f, indent=2)

def find_pdb_file(search_dir):
    """Recursively find the first .pdb file in search_dir."""
    for root, _, files in os.walk(search_dir):
        for fname in files:
            if fname.endswith('.pdb'):
                return os.path.join(root, fname)
    return None

def run_colabfold(fasta, out_pdb):
    """Run ColabFold structure prediction"""
    config = load_config()
    colabfold_config = config.get("tools", {}).get("colabfold", {})
    
    # Check if ColabFold is available
    try:
        import colabfold
        colabfold_available = True
    except ImportError:
        colabfold_available = False
    
    if not colabfold_available:
        return False, "ColabFold not installed. Run: pip install colabfold"
    
    tmp_outdir = os.path.dirname(out_pdb)
    ensure_dir(tmp_outdir)
    
    try:
        # Use colabfold_batch command
        cmd = [
            "colabfold_batch",
            "--stop-at-score", "85",  # Stop when confident
            "--num-models", "1",      # Generate 1 model for speed
            fasta,
            tmp_outdir
        ]
        
        print(f"[ColabFold] Running: {' '.join(cmd)}")
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Find the generated PDB file
        pdb_found = find_pdb_file(tmp_outdir)
        if pdb_found:
            shutil.move(pdb_found, out_pdb)
            # Save log
            with open(out_pdb + ".colabfold.log", 'w') as f:
                f.write(proc.stdout + "\n" + proc.stderr)
            return True, None
        else:
            return False, "ColabFold did not produce a .pdb file"
            
    except subprocess.CalledProcessError as e:
        return False, f"ColabFold failed: {e.stderr}"
    except Exception as e:
        return False, f"ColabFold error: {e}"

def run_openfold(fasta, out_pdb):
    """Run OpenFold structure prediction"""
    config = load_config()
    openfold_config = config.get("tools", {}).get("openfold", {})
    repo_path = openfold_config.get("repo_path", "OpenFold")
    
    if not os.path.exists(repo_path):
        return False, f"OpenFold repository not found at {repo_path}"
    
    tmp_outdir = os.path.dirname(out_pdb)
    ensure_dir(tmp_outdir)
    
    try:
        # Use OpenFold inference script
        cmd = [
            sys.executable,
            os.path.join(repo_path, "run_pretrained_openfold.py"),
            fasta,
            tmp_outdir,
            "--model_device", "cpu",  # Use CPU for compatibility
            "--config_preset", "model_1_ptm",
            "--openfold_checkpoint_path", os.path.join(repo_path, "openfold/resources/openfold_params/finetuning_ptm_1.pt")
        ]
        
        print(f"[OpenFold] Running: {' '.join(cmd)}")
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Find the generated PDB file
        pdb_found = find_pdb_file(tmp_outdir)
        if pdb_found:
            shutil.move(pdb_found, out_pdb)
            # Save log
            with open(out_pdb + ".openfold.log", 'w') as f:
                f.write(proc.stdout + "\n" + proc.stderr)
            return True, None
        else:
            return False, "OpenFold did not produce a .pdb file"
            
    except subprocess.CalledProcessError as e:
        return False, f"OpenFold failed: {e.stderr}"
    except Exception as e:
        return False, f"OpenFold error: {e}"

def run_esmfold(fasta, out_pdb):
    """Run ESMFold structure prediction"""
    config = load_config()
    esmfold_config = config.get("tools", {}).get("esmfold", {})
    
    # Check if ESM is available
    try:
        import esm
        esm_available = True
    except ImportError:
        esm_available = False
    
    if not esm_available:
        return False, "ESM not installed. Run: pip install fair-esm"
    
    tmp_outdir = os.path.dirname(out_pdb)
    ensure_dir(tmp_outdir)
    
    try:
        # Use ESMFold command
        cmd = [
            sys.executable, "-m", "esm.pretrained.esmfold_v1",
            fasta,
            tmp_outdir,
            "--num-recycles", "1"  # Reduce for speed
        ]
        
        print(f"[ESMFold] Running: {' '.join(cmd)}")
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Find the generated PDB file
        pdb_found = find_pdb_file(tmp_outdir)
        if pdb_found:
            shutil.move(pdb_found, out_pdb)
            # Save log
            with open(out_pdb + ".esmfold.log", 'w') as f:
                f.write(proc.stdout + "\n" + proc.stderr)
            return True, None
        else:
            return False, "ESMFold did not produce a .pdb file"
            
    except subprocess.CalledProcessError as e:
        return False, f"ESMFold failed: {e.stderr}"
    except Exception as e:
        return False, f"ESMFold error: {e}"

def main():
    parser = argparse.ArgumentParser(description="Predict 3D structure from FASTA using ColabFold/OpenFold/ESMFold with caching.")
    parser.add_argument("--protein", required=True, help="Input protein FASTA file (.fasta)")
    parser.add_argument("--output_dir", default="data/cache/structures/", help="Directory to save predicted PDB")
    parser.add_argument("--protein_id", default=None, help="Protein ID for naming (default: basename of FASTA)")
    parser.add_argument("--prefer", choices=["colabfold", "openfold", "esmfold"], 
                       help="Prefer a specific tool (default: try all in order)")
    args = parser.parse_args()

    fasta = args.protein
    if not os.path.isfile(fasta):
        print(f"[ERROR] FASTA file not found: {fasta}")
        sys.exit(1)
    ensure_dir(args.output_dir)
    hash_id = fasta_hash(fasta)
    protein_id = args.protein_id or os.path.splitext(os.path.basename(fasta))[0]
    pdb_path = os.path.join(args.output_dir, f"{protein_id}_{hash_id}.pdb")
    log_path = os.path.join(args.output_dir, f"{protein_id}_{hash_id}.log.json")

    # Caching
    if os.path.isfile(pdb_path):
        print(f"[CACHE] Structure already exists: {pdb_path}")
        print(f"[CACHE] Log: {log_path}")
        sys.exit(0)

    # Define tool order based on preference
    if args.prefer:
        tools = [(args.prefer, globals()[f"run_{args.prefer}"])]
        # Add other tools as fallback
        all_tools = [("colabfold", run_colabfold), ("openfold", run_openfold), ("esmfold", run_esmfold)]
        for tool_name, runner in all_tools:
            if tool_name != args.prefer:
                tools.append((tool_name, runner))
    else:
        tools = [
            ("colabfold", run_colabfold),
            ("openfold", run_openfold),
            ("esmfold", run_esmfold)
        ]

    # Try models in order
    start = time.time()
    for model, runner in tools:
        print(f"[INFO] Trying {model}...")
        ok, err = runner(fasta, pdb_path)
        if ok:
            elapsed = time.time() - start
            print(f"[SUCCESS] {model} prediction complete. Output: {pdb_path}")
            log_result(log_path, model, elapsed)
            sys.exit(0)
        else:
            print(f"[WARNING] {model} failed: {err}")
    
    # All failed
    elapsed = time.time() - start
    log_result(log_path, None, elapsed, error="All models failed.")
    print(f"[ERROR] All structure predictors failed. See log: {log_path}")
    print("[INFO] Try installing tools with: python install_real_tools.py")
    sys.exit(2)

if __name__ == "__main__":
    main() 