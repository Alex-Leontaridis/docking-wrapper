#!/usr/bin/env python3
"""
Druggability Scoring using fpocket

Analyzes protein-ligand complexes to determine binding site druggability.
Uses fpocket CLI tool to identify and score potential binding pockets.
"""

import argparse
import os
import sys
import json
import subprocess
import tempfile
import time
from pathlib import Path

def ensure_dir(path):
    """Ensure directory exists"""
    if not os.path.exists(path):
        os.makedirs(path)

def check_fpocket_installed():
    """Check if fpocket is available in PATH"""
    try:
        result = subprocess.run(["fpocket", "-h"], 
                              capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False

def install_fpocket():
    """Attempt to install fpocket"""
    print("[INFO] Attempting to install fpocket...")
    
    # Try conda installation first
    try:
        subprocess.run(["conda", "install", "-c", "conda-forge", "fpocket", "-y"], 
                      check=True, capture_output=True)
        print("[SUCCESS] fpocket installed via conda")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[WARNING] Conda installation failed")
    
    # Try apt-get (Ubuntu/Debian)
    try:
        subprocess.run(["sudo", "apt-get", "update"], check=True, capture_output=True)
        subprocess.run(["sudo", "apt-get", "install", "-y", "fpocket"], 
                      check=True, capture_output=True)
        print("[SUCCESS] fpocket installed via apt-get")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[WARNING] apt-get installation failed")
    
    # Try brew (macOS)
    try:
        subprocess.run(["brew", "install", "fpocket"], check=True, capture_output=True)
        print("[SUCCESS] fpocket installed via brew")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("[WARNING] brew installation failed")
    
    print("[ERROR] Could not install fpocket automatically")
    print("[INFO] Please install manually:")
    print("  - Ubuntu/Debian: sudo apt-get install fpocket")
    print("  - macOS: brew install fpocket")
    print("  - Conda: conda install -c conda-forge fpocket")
    return False

def parse_fpocket_output(fpocket_dir):
    """Parse fpocket output to extract druggability scores"""
    scores = []
    
    # Look for fpocket output files
    pdb_files = list(Path(fpocket_dir).glob("*.pdb"))
    info_files = list(Path(fpocket_dir).glob("*.txt"))
    
    if not pdb_files and not info_files:
        return None, "No fpocket output files found"
    
    # Parse info files for pocket statistics
    for info_file in info_files:
        if "info" in info_file.name:
            try:
                with open(info_file, 'r') as f:
                    lines = f.readlines()
                
                # Extract pocket information
                pocket_data = {}
                for line in lines:
                    if line.startswith("Pocket"):
                        parts = line.strip().split()
                        if len(parts) >= 8:
                            pocket_num = int(parts[1])
                            # Extract key metrics
                            try:
                                pocket_data[pocket_num] = {
                                    "volume": float(parts[2]),
                                    "area": float(parts[3]),
                                    "opening": float(parts[4]),
                                    "apol_area": float(parts[5]),
                                    "hydrophobicity": float(parts[6]),
                                    "polarity": float(parts[7])
                                }
                            except (ValueError, IndexError):
                                continue
                
                # Calculate druggability score based on pocket properties
                for pocket_num, data in pocket_data.items():
                    score = calculate_druggability_score(data)
                    scores.append({
                        "pocket_id": pocket_num,
                        "druggability_score": score,
                        "volume": data["volume"],
                        "area": data["area"],
                        "opening": data["opening"],
                        "apol_area": data["apol_area"],
                        "hydrophobicity": data["hydrophobicity"],
                        "polarity": data["polarity"]
                    })
                
            except Exception as e:
                return None, f"Error parsing fpocket info: {e}"
    
    if not scores:
        return None, "No valid pocket data found"
    
    # Return the best scoring pocket
    best_pocket = max(scores, key=lambda x: x["druggability_score"])
    return best_pocket, None

def calculate_druggability_score(pocket_data):
    """
    Calculate druggability score based on pocket properties.
    Returns a score between 0 and 1.
    
    Scoring criteria based on typical druggable pocket characteristics:
    - Volume: 100-1000 Å³ (optimal range)
    - Area: 200-2000 Å² (optimal range)
    - Hydrophobicity: 0.5-1.0 (optimal range)
    - Polarity: 0.2-0.8 (optimal range)
    """
    volume = pocket_data["volume"]
    area = pocket_data["area"]
    hydrophobicity = pocket_data["hydrophobicity"]
    polarity = pocket_data["polarity"]
    
    # Volume score (0-1)
    if 100 <= volume <= 1000:
        volume_score = 1.0
    elif volume < 50 or volume > 2000:
        volume_score = 0.0
    else:
        # Linear interpolation
        if volume < 100:
            volume_score = volume / 100
        else:
            volume_score = max(0, (2000 - volume) / 1000)
    
    # Area score (0-1)
    if 200 <= area <= 2000:
        area_score = 1.0
    elif area < 100 or area > 3000:
        area_score = 0.0
    else:
        # Linear interpolation
        if area < 200:
            area_score = area / 200
        else:
            area_score = max(0, (3000 - area) / 1000)
    
    # Hydrophobicity score (0-1)
    if 0.5 <= hydrophobicity <= 1.0:
        hydrophobicity_score = 1.0
    elif hydrophobicity < 0.0 or hydrophobicity > 1.5:
        hydrophobicity_score = 0.0
    else:
        # Linear interpolation
        if hydrophobicity < 0.5:
            hydrophobicity_score = hydrophobicity / 0.5
        else:
            hydrophobicity_score = max(0, (1.5 - hydrophobicity) / 0.5)
    
    # Polarity score (0-1)
    if 0.2 <= polarity <= 0.8:
        polarity_score = 1.0
    elif polarity < 0.0 or polarity > 1.0:
        polarity_score = 0.0
    else:
        # Linear interpolation
        if polarity < 0.2:
            polarity_score = polarity / 0.2
        else:
            polarity_score = max(0, (1.0 - polarity) / 0.2)
    
    # Weighted average (volume and area are most important)
    final_score = (
        0.4 * volume_score +
        0.3 * area_score +
        0.2 * hydrophobicity_score +
        0.1 * polarity_score
    )
    
    return round(final_score, 3)

def run_fpocket_analysis(pdb_path, output_dir):
    """Run fpocket analysis on PDB file"""
    try:
        # Create temporary directory for fpocket output
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run fpocket
            cmd = ["fpocket", "-f", pdb_path, "-o", tmpdir]
            print(f"[fpocket] Running: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                return None, f"fpocket failed: {result.stderr}"
            
            # Parse results
            pocket_data, error = parse_fpocket_output(tmpdir)
            if error:
                return None, error
            
            return pocket_data, None
            
    except subprocess.TimeoutExpired:
        return None, "fpocket analysis timed out"
    except Exception as e:
        return None, f"fpocket error: {e}"

def mock_druggability_analysis(pdb_path):
    """Mock druggability analysis for testing"""
    import random
    
    # Simulate realistic druggability scores
    # Most binding sites have scores between 0.3 and 0.9
    score = random.uniform(0.3, 0.9)
    
    return {
        "pocket_id": 1,
        "druggability_score": round(score, 3),
        "volume": random.uniform(200, 800),
        "area": random.uniform(300, 1500),
        "opening": random.uniform(0.5, 2.0),
        "apol_area": random.uniform(100, 600),
        "hydrophobicity": random.uniform(0.4, 0.9),
        "polarity": random.uniform(0.2, 0.7),
        "mock_prediction": True
    }

def main():
    parser = argparse.ArgumentParser(description="Druggability scoring using fpocket")
    parser.add_argument("--protein", required=True, help="Input protein PDB file (with bound ligand)")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--protein_id", default=None, help="Protein ID for naming")
    parser.add_argument("--ligand_id", default=None, help="Ligand ID for naming")
    parser.add_argument("--force_mock", action="store_true", help="Force use of mock prediction")
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isfile(args.protein):
        print(f"[ERROR] Protein file not found: {args.protein}")
        sys.exit(1)
    
    if not args.protein.endswith('.pdb'):
        print(f"[ERROR] Protein must be a PDB file: {args.protein}")
        sys.exit(1)
    
    # Set default IDs if not provided
    protein_id = args.protein_id or os.path.splitext(os.path.basename(args.protein))[0]
    ligand_id = args.ligand_id or "LIG"
    
    # Initialize result structure
    start_time = time.time()
    result = {
        "protein": args.protein,
        "protein_id": protein_id,
        "ligand_id": ligand_id,
        "tool": "fpocket",
        "fallback": False,
        "mock_prediction": False
    }
    
    try:
        # Check if fpocket is available
        if not args.force_mock and check_fpocket_installed():
            print("[INFO] fpocket found, running analysis...")
            pocket_data, error = run_fpocket_analysis(args.protein, "data/output/druggability")
            
            if pocket_data:
                result.update(pocket_data)
                result["success"] = True
            else:
                print(f"[WARNING] fpocket analysis failed: {error}")
                print("[INFO] Using mock prediction as fallback...")
                result.update(mock_druggability_analysis(args.protein))
                result["fallback"] = True
                result["mock_prediction"] = True
                result["error"] = error
        else:
            if not args.force_mock:
                print("[WARNING] fpocket not found")
                if not install_fpocket():
                    print("[INFO] Using mock prediction as fallback...")
                    result.update(mock_druggability_analysis(args.protein))
                    result["fallback"] = True
                    result["mock_prediction"] = True
                    result["error"] = "fpocket not installed"
                else:
                    # Try again after installation
                    pocket_data, error = run_fpocket_analysis(args.protein, "data/output/druggability")
                    if pocket_data:
                        result.update(pocket_data)
                        result["success"] = True
                    else:
                        result.update(mock_druggability_analysis(args.protein))
                        result["fallback"] = True
                        result["mock_prediction"] = True
                        result["error"] = error
            else:
                print("[INFO] Using mock prediction (forced)")
                result.update(mock_druggability_analysis(args.protein))
                result["mock_prediction"] = True
        
        # Add timing information
        result["runtime_sec"] = round(time.time() - start_time, 2)
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"[SUCCESS] Druggability analysis complete: {args.output}")
        print(f"[INFO] Druggability score: {result.get('druggability_score', 'N/A')}")
        
        if result.get("fallback"):
            print("[INFO] Used fallback prediction - consider installing fpocket for real analysis")
        
        return 0
        
    except Exception as e:
        result["error"] = str(e)
        result["runtime_sec"] = round(time.time() - start_time, 2)
        
        # Save error result
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"[ERROR] Druggability analysis failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 