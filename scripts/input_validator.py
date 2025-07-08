#!/usr/bin/env python3
"""
Input validation module for molecular docking pipeline.
Provides comprehensive validation for protein and ligand files.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import logging
from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context
from pathlib import Path
from typing import Tuple, List, Dict, Optional
import tempfile
import shutil

# RDKit imports for molecular validation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not available - molecular validation will be limited")

# BioPython imports for protein validation
try:
    from Bio import PDB
    from Bio.PDB.PDBParser import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logger.warning("BioPython not available - protein validation will be limited")


class InputValidator:
    """Comprehensive input validation for molecular docking files."""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the validator.
        
        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        
        # Supported file extensions
        self.SUPPORTED_PROTEIN_EXT = {'.pdb', '.pdbqt'}
        self.SUPPORTED_LIGAND_EXT = {'.smi', '.sdf', '.mol2', '.pdb', '.pdbqt'}
        
        # File size limits (in MB)
        self.MAX_PROTEIN_SIZE_MB = 100
        self.MAX_LIGAND_SIZE_MB = 50
        
        # Validation results cache
        self._validation_cache = {}
    
    def validate_protein_file(self, file_path: str) -> Tuple[bool, str]:
        """
        Validate a protein file for docking.
        
        Args:
            file_path: Path to the protein file
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        cache_key = f"protein_{file_path}"
        if cache_key in self._validation_cache:
            return self._validation_cache[cache_key]
        
        try:
            # Basic file checks
            if not os.path.exists(file_path):
                result = (False, f"Protein file does not exist: {file_path}")
                self._validation_cache[cache_key] = result
                return result
            
            if not os.access(file_path, os.R_OK):
                result = (False, f"Protein file is not readable: {file_path}")
                self._validation_cache[cache_key] = result
                return result
            
            # File extension check
            ext = Path(file_path).suffix.lower()
            if ext not in self.SUPPORTED_PROTEIN_EXT:
                result = (False, f"Unsupported protein file format: {ext}. Supported: {', '.join(self.SUPPORTED_PROTEIN_EXT)}")
                self._validation_cache[cache_key] = result
                return result
            
            # File size check
            file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
            if file_size_mb > self.MAX_PROTEIN_SIZE_MB:
                result = (False, f"Protein file too large: {file_size_mb:.1f}MB (max: {self.MAX_PROTEIN_SIZE_MB}MB)")
                self._validation_cache[cache_key] = result
                return result
            
            # Content validation based on file type
            if ext == '.pdb':
                valid, error = self._validate_pdb_content(file_path)
            elif ext == '.pdbqt':
                valid, error = self._validate_pdbqt_content(file_path)
            else:
                valid, error = True, ""
            
            result = (valid, error)
            self._validation_cache[cache_key] = result
            return result
            
        except Exception as e:
            result = (False, f"Unexpected error validating protein file: {e}")
            self._validation_cache[cache_key] = result
            return result
    
    def validate_ligand_file(self, file_path: str) -> Tuple[bool, str]:
        """
        Validate a ligand file for docking.
        
        Args:
            file_path: Path to the ligand file
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        cache_key = f"ligand_{file_path}"
        if cache_key in self._validation_cache:
            return self._validation_cache[cache_key]
        
        try:
            # Basic file checks
            if not os.path.exists(file_path):
                result = (False, f"Ligand file does not exist: {file_path}")
                self._validation_cache[cache_key] = result
                return result
            
            if not os.access(file_path, os.R_OK):
                result = (False, f"Ligand file is not readable: {file_path}")
                self._validation_cache[cache_key] = result
                return result
            
            # File extension check
            ext = Path(file_path).suffix.lower()
            if ext not in self.SUPPORTED_LIGAND_EXT:
                result = (False, f"Unsupported ligand file format: {ext}. Supported: {', '.join(self.SUPPORTED_LIGAND_EXT)}")
                self._validation_cache[cache_key] = result
                return result
            
            # File size check
            file_size_mb = os.path.getsize(file_path) / (1024 * 1024)
            if file_size_mb > self.MAX_LIGAND_SIZE_MB:
                result = (False, f"Ligand file too large: {file_size_mb:.1f}MB (max: {self.MAX_LIGAND_SIZE_MB}MB)")
                self._validation_cache[cache_key] = result
                return result
            
            # Content validation based on file type
            if ext == '.sdf':
                valid, error = self._validate_sdf_content(file_path)
            elif ext == '.smi':
                valid, error = self._validate_smi_content(file_path)
            elif ext == '.mol2':
                valid, error = self._validate_mol2_content(file_path)
            elif ext in ['.pdb', '.pdbqt']:
                valid, error = self._validate_pdb_content(file_path)
            else:
                valid, error = True, ""
            
            result = (valid, error)
            self._validation_cache[cache_key] = result
            return result
            
        except Exception as e:
            result = (False, f"Unexpected error validating ligand file: {e}")
            self._validation_cache[cache_key] = result
            return result
    
    def _validate_pdb_content(self, file_path: str) -> Tuple[bool, str]:
        """Validate PDB file content."""
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                return False, "PDB file is empty"
            
            # Check for basic PDB structure
            has_atom = any(line.startswith('ATOM') for line in lines)
            has_hetatm = any(line.startswith('HETATM') for line in lines)
            
            if not (has_atom or has_hetatm):
                return False, "PDB file contains no ATOM or HETATM records"
            
            # Validate with BioPython if available
            if BIOPYTHON_AVAILABLE:
                try:
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure('test', file_path)
                    if len(list(structure.get_atoms())) == 0:
                        return False, "PDB file contains no valid atoms"
                except Exception as e:
                    return False, f"BioPython PDB validation failed: {e}"
            
            return True, ""
            
        except Exception as e:
            return False, f"Error reading PDB file: {e}"
    
    def _validate_pdbqt_content(self, file_path: str) -> Tuple[bool, str]:
        """Validate PDBQT file content."""
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                return False, "PDBQT file is empty"
            
            # Check for basic PDBQT structure
            has_atom = any(line.startswith('ATOM') for line in lines)
            has_hetatm = any(line.startswith('HETATM') for line in lines)
            
            if not (has_atom or has_hetatm):
                return False, "PDBQT file contains no ATOM or HETATM records"
            
            # Check for AutoDock atom types (should be present in PDBQT)
            has_autodock_types = any(len(line) > 78 and line[77:79].strip() for line in lines if line.startswith(('ATOM', 'HETATM')))
            if not has_autodock_types:
                self.logger.warning("PDBQT file may not have proper AutoDock atom types")
            
            return True, ""
            
        except Exception as e:
            return False, f"Error reading PDBQT file: {e}"
    
    def _validate_sdf_content(self, file_path: str) -> Tuple[bool, str]:
        """Validate SDF file content."""
        if not RDKIT_AVAILABLE:
            return True, "RDKit not available for SDF validation"
        
        try:
            suppl = Chem.SDMolSupplier(file_path, removeHs=False)
            mol_count = 0
            valid_mols = 0
            
            for mol in suppl:
                mol_count += 1
                if mol is not None:
                    valid_mols += 1
                    # Basic molecular validation
                    if mol.GetNumAtoms() == 0:
                        return False, f"SDF file contains molecule with no atoms at position {mol_count}"
                    if mol.GetNumBonds() == 0 and mol.GetNumAtoms() > 1:
                        return False, f"SDF file contains disconnected atoms at position {mol_count}"
            
            if mol_count == 0:
                return False, "SDF file contains no molecules"
            
            if valid_mols == 0:
                return False, "SDF file contains no valid molecules"
            
            if valid_mols < mol_count:
                self.logger.warning(f"SDF file contains {mol_count - valid_mols} invalid molecules out of {mol_count}")
            
            return True, ""
            
        except Exception as e:
            return False, f"Error reading SDF file: {e}"
    
    def _validate_smi_content(self, file_path: str) -> Tuple[bool, str]:
        """Validate SMILES file content."""
        if not RDKIT_AVAILABLE:
            return True, "RDKit not available for SMILES validation"
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                return False, "SMILES file is empty"
            
            valid_smiles = 0
            for line_num, line in enumerate(lines, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse SMILES (format: SMILES [name])
                parts = line.split()
                if not parts:
                    continue
                
                smiles = parts[0]
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        valid_smiles += 1
                        if mol.GetNumAtoms() == 0:
                            return False, f"SMILES file contains molecule with no atoms at line {line_num}"
                    else:
                        return False, f"SMILES file contains invalid SMILES at line {line_num}: {smiles}"
                except Exception as e:
                    return False, f"Error parsing SMILES at line {line_num}: {e}"
            
            if valid_smiles == 0:
                return False, "SMILES file contains no valid molecules"
            
            return True, ""
            
        except Exception as e:
            return False, f"Error reading SMILES file: {e}"
    
    def _validate_mol2_content(self, file_path: str) -> Tuple[bool, str]:
        """Validate MOL2 file content."""
        if not RDKIT_AVAILABLE:
            return True, "RDKit not available for MOL2 validation"
        
        try:
            mol = Chem.MolFromMol2File(file_path, removeHs=False)
            if mol is None:
                return False, "MOL2 file contains no valid molecules"
            
            if mol.GetNumAtoms() == 0:
                return False, "MOL2 file contains molecule with no atoms"
            
            return True, ""
            
        except Exception as e:
            return False, f"Error reading MOL2 file: {e}"
    
    def validate_directory(self, dir_path: str, file_type: str = "any") -> Tuple[bool, str]:
        """
        Validate a directory containing input files.
        
        Args:
            dir_path: Path to the directory
            file_type: Type of files to look for ("protein", "ligand", or "any")
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            if not os.path.exists(dir_path):
                return False, f"Directory does not exist: {dir_path}"
            
            if not os.path.isdir(dir_path):
                return False, f"Path is not a directory: {dir_path}"
            
            if not os.access(dir_path, os.R_OK):
                return False, f"Directory is not readable: {dir_path}"
            
            # Find files based on type
            if file_type == "protein":
                extensions = self.SUPPORTED_PROTEIN_EXT
            elif file_type == "ligand":
                extensions = self.SUPPORTED_LIGAND_EXT
            else:
                extensions = self.SUPPORTED_PROTEIN_EXT | self.SUPPORTED_LIGAND_EXT
            
            files = []
            for ext in extensions:
                files.extend(Path(dir_path).glob(f"*{ext}"))
                files.extend(Path(dir_path).glob(f"*{ext.upper()}"))
            
            if not files:
                return False, f"No supported files found in directory: {dir_path}"
            
            return True, f"Found {len(files)} supported files"
            
        except Exception as e:
            return False, f"Error validating directory: {e}"
    
    def clear_cache(self):
        """Clear the validation cache."""
        self._validation_cache.clear()


# Convenience functions for backward compatibility
def validate_protein_file(file_path: str) -> Tuple[bool, str]:
    """Validate a protein file."""
    validator = InputValidator()
    return validator.validate_protein_file(file_path)


def validate_ligand_file(file_path: str) -> Tuple[bool, str]:
    """Validate a ligand file."""
    validator = InputValidator()
    return validator.validate_ligand_file(file_path)


def validate_directory(dir_path: str, file_type: str = "any") -> Tuple[bool, str]:
    """Validate a directory."""
    validator = InputValidator()
    return validator.validate_directory(dir_path, file_type) 