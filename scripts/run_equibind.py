#!/usr/bin/env python3
"""
EquiBind Integration Script

This script integrates EquiBind docking engine into the docking wrapper pipeline.
EquiBind is a machine learning-based method for predicting protein-ligand binding.

Author: Development Team
Date: 2024
"""

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Optional, Dict, Any

# Add the project root to the path
sys.path.append(str(Path(__file__).parent.parent))

from utils.logging import setup_logger
from utils.path_manager import PathManager

class EquiBindRunner:
    """EquiBind docking engine integration."""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """Initialize EquiBind runner.
        
        Args:
            config: Configuration dictionary for EquiBind
        """
        self.config = config or {}
        self.logger = setup_logger(__name__)
        self.path_manager = PathManager()
        
        # EquiBind specific paths
        self.equibind_dir = Path("EquiBind")
        self.model_path = self.equibind_dir / "runs" / "flexible_self_docking" / "best_checkpoint.pt"
        
    def check_dependencies(self) -> bool:
        """Check if EquiBind dependencies are available.
        
        Returns:
            True if all dependencies are available, False otherwise
        """
        self.logger.info("Checking EquiBind dependencies...")
        
        # Check if EquiBind directory exists
        if not self.equibind_dir.exists():
            self.logger.error(f"EquiBind directory not found: {self.equibind_dir}")
            return False
            
        # Check if model checkpoint exists
        if not self.model_path.exists():
            self.logger.error(f"EquiBind model not found: {self.model_path}")
            return False
            
        # Check if inference script exists
        inference_script = self.equibind_dir / "inference.py"
        if not inference_script.exists():
            self.logger.error(f"EquiBind inference script not found: {inference_script}")
            return False
            
        self.logger.info("All EquiBind dependencies are available")
        return True
        
    def prepare_inputs(self, protein_path: str, ligand_path: str, output_dir: str) -> Dict[str, str]:
        """Prepare inputs for EquiBind.
        
        Args:
            protein_path: Path to protein structure file
            ligand_path: Path to ligand structure file
            output_dir: Output directory for results
            
        Returns:
            Dictionary with prepared input paths
        """
        self.logger.info("Preparing inputs for EquiBind...")
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Prepare protein input
        protein_input = output_path / "protein_equibind.pdb"
        self._prepare_protein(protein_path, protein_input)
        
        # Prepare ligand input
        ligand_input = output_path / "ligand_equibind.sdf"
        self._prepare_ligand(ligand_path, ligand_input)
        
        return {
            "protein": str(protein_input),
            "ligand": str(ligand_input),
            "output": str(output_path)
        }
        
    def _prepare_protein(self, input_path: str, output_path: Path) -> None:
        """Prepare protein structure for EquiBind.
        
        Args:
            input_path: Input protein file path
            output_path: Output protein file path
        """
        # For now, just copy the file
        # In a real implementation, you might need to convert formats
        import shutil
        shutil.copy2(input_path, output_path)
        self.logger.debug(f"Prepared protein: {output_path}")
        
    def _prepare_ligand(self, input_path: str, output_path: Path) -> None:
        """Prepare ligand structure for EquiBind.
        
        Args:
            input_path: Input ligand file path
            output_path: Output ligand file path
        """
        # For now, just copy the file
        # In a real implementation, you might need to convert formats
        import shutil
        shutil.copy2(input_path, output_path)
        self.logger.debug(f"Prepared ligand: {output_path}")
        
    def run_docking(self, protein_path: str, ligand_path: str, output_dir: str) -> Dict[str, Any]:
        """Run EquiBind docking.
        
        Args:
            protein_path: Path to protein structure file
            ligand_path: Path to ligand structure file
            output_dir: Output directory for results
            
        Returns:
            Dictionary with docking results
        """
        self.logger.info("Starting EquiBind docking...")
        
        # Check dependencies
        if not self.check_dependencies():
            raise RuntimeError("EquiBind dependencies not available")
            
        # Prepare inputs
        inputs = self.prepare_inputs(protein_path, ligand_path, output_dir)
        
        # Run EquiBind inference
        try:
            result = self._run_inference(inputs)
            self.logger.info("EquiBind docking completed successfully")
            return result
        except Exception as e:
            self.logger.error(f"EquiBind docking failed: {e}")
            raise
            
    def _run_inference(self, inputs: Dict[str, str]) -> Dict[str, Any]:
        """Run EquiBind inference.
        
        Args:
            inputs: Dictionary with input paths
            
        Returns:
            Dictionary with inference results
        """
        # Change to EquiBind directory
        original_dir = os.getcwd()
        os.chdir(self.equibind_dir)
        
        try:
            # Build command for EquiBind inference
            cmd = [
                sys.executable, "inference.py",
                "--protein", inputs["protein"],
                "--ligand", inputs["ligand"],
                "--output", inputs["output"],
                "--model_path", str(self.model_path)
            ]
            
            # Add configuration options
            if self.config.get("num_poses"):
                cmd.extend(["--num_poses", str(self.config["num_poses"])])
            if self.config.get("confidence_threshold"):
                cmd.extend(["--confidence_threshold", str(self.config["confidence_threshold"])])
                
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            
            # Run the command
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse results
            return self._parse_results(inputs["output"])
            
        finally:
            # Restore original directory
            os.chdir(original_dir)
            
    def _parse_results(self, output_dir: str) -> Dict[str, Any]:
        """Parse EquiBind results.
        
        Args:
            output_dir: Directory containing results
            
        Returns:
            Dictionary with parsed results
        """
        output_path = Path(output_dir)
        
        # Look for output files
        results = {
            "poses": [],
            "confidence_scores": [],
            "output_files": []
        }
        
        # Find pose files
        for pose_file in output_path.glob("pose_*.sdf"):
            results["poses"].append(str(pose_file))
            results["output_files"].append(str(pose_file))
            
        # Find confidence file
        confidence_file = output_path / "equibind_confidence.txt"
        if confidence_file.exists():
            results["confidence_scores"] = self._read_confidence_scores(confidence_file)
            results["output_files"].append(str(confidence_file))
            
        # Find main output file
        main_output = output_path / "equibind_out.sdf"
        if main_output.exists():
            results["main_output"] = str(main_output)
            results["output_files"].append(str(main_output))
            
        return results
        
    def _read_confidence_scores(self, confidence_file: Path) -> list:
        """Read confidence scores from file.
        
        Args:
            confidence_file: Path to confidence scores file
            
        Returns:
            List of confidence scores
        """
        scores = []
        try:
            with open(confidence_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and line.replace('.', '').replace('-', '').isdigit():
                        scores.append(float(line))
        except Exception as e:
            self.logger.warning(f"Could not read confidence scores: {e}")
            
        return scores

def main():
    """Main function for EquiBind integration."""
    import argparse
    
    parser = argparse.ArgumentParser(description="EquiBind Docking Integration")
    parser.add_argument("--protein", required=True, help="Path to protein structure file")
    parser.add_argument("--ligand", required=True, help="Path to ligand structure file")
    parser.add_argument("--output", required=True, help="Output directory for results")
    parser.add_argument("--config", help="Path to configuration file")
    parser.add_argument("--num_poses", type=int, default=3, help="Number of poses to generate")
    parser.add_argument("--confidence_threshold", type=float, default=0.5, help="Confidence threshold")
    
    args = parser.parse_args()
    
    # Load configuration
    config = {
        "num_poses": args.num_poses,
        "confidence_threshold": args.confidence_threshold
    }
    
    # Create EquiBind runner
    runner = EquiBindRunner(config)
    
    # Run docking
    try:
        results = runner.run_docking(args.protein, args.ligand, args.output)
        print(f"EquiBind docking completed successfully")
        print(f"Generated {len(results['poses'])} poses")
        print(f"Output files: {results['output_files']}")
    except Exception as e:
        print(f"EquiBind docking failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 