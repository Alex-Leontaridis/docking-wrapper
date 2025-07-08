#!/usr/bin/env python3
"""
Docking Results Parser

Parses output files from different docking engines (Vina, GNINA, DiffDock)
and generates summary DataFrames with pose information and scores.
"""

import os
import re
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem

class DockingResultsParser:
    """Parser for docking engine output files."""
    
    def __init__(self, base_dir: str, output_dir: str):
        """
        Initialize the parser.
        
        Args:
            base_dir: Directory containing docking results
            output_dir: Directory to save parsed results
        """
        self.base_dir = Path(base_dir)
        self.output_dir = Path(output_dir)
        self.failed_runs = {}
        self.logger = logging.getLogger(__name__)
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_summary(self, ligand_name: str) -> pd.DataFrame:
        """
        Generate summary DataFrame for a ligand.
        
        Args:
            ligand_name: Name of the ligand
            
        Returns:
            DataFrame with pose information and scores
        """
        self.logger.info("Starting results parsing...")
        
        all_poses = []
        
        # Parse Vina results
        vina_poses = self._parse_vina_results(ligand_name)
        all_poses.extend(vina_poses)
        
        # Parse GNINA results
        gnina_poses = self._parse_gnina_results(ligand_name)
        all_poses.extend(gnina_poses)
        
        # Parse DiffDock results
        diffdock_poses = self._parse_diffdock_results(ligand_name)
        all_poses.extend(diffdock_poses)
        
        # Create summary DataFrame
        if all_poses:
            df = pd.DataFrame(all_poses)
            self.logger.info(f"Generated summary with {len(df)} total poses")
            
            # Save summary to CSV
            summary_file = self.output_dir / f"{ligand_name}_summary.csv"
            df.to_csv(summary_file, index=False)
            self.logger.info(f"Results saved to: {summary_file}")
            
            return df
        else:
            self.logger.warning("No poses found in any docking results")
            return pd.DataFrame()
    
    def _parse_vina_results(self, ligand_name: str) -> List[Dict]:
        """Parse Vina PDBQT output file."""
        vina_file = self.base_dir / "vina_output" / "vina_out.pdbqt"
        
        if not vina_file.exists():
            self.logger.error("Failed to parse Vina results: No Vina output files found")
            return []
        
        poses = []
        try:
            with open(vina_file, 'r') as f:
                content = f.read()
            
            # Extract affinity scores from REMARK lines
            affinity_pattern = r'REMARK VINA RESULT:\s+([-\d.]+)'
            scores = re.findall(affinity_pattern, content)
            
            # Count models (poses)
            model_pattern = r'MODEL\s+(\d+)'
            models = re.findall(model_pattern, content)
            
            for i, (score, model_num) in enumerate(zip(scores, models)):
                poses.append({
                    'ligand_name': ligand_name,
                    'engine': 'vina',
                    'pose_id': int(model_num),
                    'affinity_kcal_per_mol': float(score),
                    'confidence_score': None,  # Vina doesn't provide confidence
                    'pose_file': str(vina_file),
                    'status': 'success'
                })
            
            self.logger.info(f"Parsed {len(poses)} Vina poses from {vina_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to parse Vina results: {e}")
            self.failed_runs['vina'] = str(e)
        
        return poses
    
    def _parse_gnina_results(self, ligand_name: str) -> List[Dict]:
        """Parse GNINA output files."""
        gnina_dir = self.base_dir / "gnina_output"
        
        if not gnina_dir.exists():
            self.logger.error("Failed to parse GNINA results: No GNINA output files found")
            return []
        
        poses = []
        try:
            # Look for GNINA output files (typically .sdf or .pdbqt)
            gnina_files = list(gnina_dir.glob("*.sdf")) + list(gnina_dir.glob("*.pdbqt"))
            
            if not gnina_files:
                self.logger.error("Failed to parse GNINA results: No GNINA output files found")
                return []
            
            # For now, we'll create a basic entry since GNINA parsing is complex
            # In a full implementation, you'd parse the actual GNINA output format
            poses.append({
                'ligand_name': ligand_name,
                'engine': 'gnina',
                'pose_id': 1,
                'affinity_kcal_per_mol': None,  # Would extract from GNINA output
                'confidence_score': None,  # Would extract from GNINA output
                'pose_file': str(gnina_files[0]),
                'status': 'success'
            })
            
            self.logger.info(f"Parsed {len(poses)} GNINA poses")
            
        except Exception as e:
            self.logger.error(f"Failed to parse GNINA results: {e}")
            self.failed_runs['gnina'] = str(e)
        
        return poses
    
    def _parse_diffdock_results(self, ligand_name: str) -> List[Dict]:
        """Parse DiffDock output files."""
        diffdock_dir = self.base_dir / "diffdock_output"
        
        if not diffdock_dir.exists():
            self.logger.error("Failed to parse DiffDock results: No DiffDock SDF files found")
            return []
        
        poses = []
        try:
            # Parse confidence file
            confidence_file = diffdock_dir / "diffdock_confidence.txt"
            confidence_scores = []
            
            if confidence_file.exists():
                with open(confidence_file, 'r') as f:
                    lines = f.readlines()
                
                # Skip header line and parse confidence scores
                for line in lines[1:]:  # Skip header
                    line = line.strip()
                    if line and ',' in line:
                        try:
                            parts = line.split(',')
                            if len(parts) >= 2:
                                confidence = float(parts[1])
                                confidence_scores.append(confidence)
                        except (ValueError, IndexError):
                            continue
            
            # Look for SDF files with poses
            sdf_files = list(diffdock_dir.glob("*.sdf"))
            
            if not sdf_files:
                self.logger.error("Failed to parse DiffDock results: No DiffDock SDF files found")
                return []
            
            # Create pose entries
            for i, sdf_file in enumerate(sdf_files):
                confidence = confidence_scores[i] if i < len(confidence_scores) else None
                
                poses.append({
                    'ligand_name': ligand_name,
                    'engine': 'diffdock',
                    'pose_id': i + 1,
                    'affinity_kcal_per_mol': None,  # DiffDock doesn't provide affinity
                    'confidence_score': confidence,
                    'pose_file': str(sdf_file),
                    'status': 'success'
                })
            
            self.logger.info(f"Parsed {len(poses)} DiffDock poses")
            
        except Exception as e:
            self.logger.error(f"Failed to parse DiffDock results: {e}")
            self.failed_runs['diffdock'] = str(e)
        
        return poses
    
    def get_failed_runs(self) -> Dict[str, str]:
        """Get dictionary of failed runs and their error messages."""
        return self.failed_runs.copy() 