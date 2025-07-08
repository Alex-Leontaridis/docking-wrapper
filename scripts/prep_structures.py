#!/usr/bin/env python3
import argparse
import os
import sys
import logging
import subprocess
import shutil
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Setup logging first
from utils.logging import setup_logging as setup_docking_logging
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir

# Create logs directory if it doesn't exist
ensure_dir('logs')

# Configure logging
logger = setup_docking_logging(__name__)

# Import config after logging is set up
from config import config

# Import resource management utilities
try:
    from resource_manager import temp_file, temp_directory, register_temp_file, cleanup_all
    RESOURCE_MANAGEMENT_AVAILABLE = True
except ImportError:
    RESOURCE_MANAGEMENT_AVAILABLE = False
    logger.warning("Resource management utilities not available - using basic cleanup")

try:
    # Try to import meeko with proper error handling
    try:
        from meeko import MoleculePreparation
        MEEKO_AVAILABLE = True
        logger.info("Meeko successfully imported")
    except ImportError as e:
        if "rdkit.six" in str(e):
            logger.warning("Meeko import failed due to rdkit.six module issue. This is a known compatibility problem.")
            logger.warning("Falling back to RDKit-based ligand preparation.")
            MEEKO_AVAILABLE = False
        else:
            logger.warning(f"Meeko import failed: {e}")
            MEEKO_AVAILABLE = False
    except Exception as e:
        logger.warning(f"Unexpected error importing Meeko: {e}")
        MEEKO_AVAILABLE = False
except Exception as e:
    logger.warning(f"Failed to check Meeko availability: {e}")
    MEEKO_AVAILABLE = False

SUPPORTED_PROTEIN_EXT = {'.pdb', '.pdbqt'}
SUPPORTED_LIGAND_EXT = {'.smi', '.sdf', '.mol2'}


def validate_file(file_path, valid_exts, file_type):
    if not os.path.isfile(file_path):
        logger.error(f"{file_type} file '{file_path}' does not exist.")
        sys.exit(1)
    if not os.access(file_path, os.R_OK):
        logger.error(f"{file_type} file '{file_path}' is not readable.")
        sys.exit(1)
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in valid_exts:
        logger.error(f"{file_type} file '{file_path}' must have one of the following extensions: {', '.join(valid_exts)}.")
        sys.exit(1)
    return file_path


def validate_output_dir(output_file):
    output_dir = os.path.dirname(os.path.abspath(output_file)) or '.'
    if not os.access(output_dir, os.W_OK):
        logger.error(f"Output directory '{output_dir}' is not writable.")
        sys.exit(1)


def clean_pdbqt_formatting(pdbqt_file):
    """
    Clean up PDBQT formatting issues that can cause Vina parsing errors.
    - Fix atom names that are too long (> 4 characters)
    - Remove problematic alternate conformations
    - Ensure proper coordinate formatting
    """
    if RESOURCE_MANAGEMENT_AVAILABLE:
        # Use resource management utilities
        with temp_file(suffix='.pdbqt') as temp_file_path:
            _clean_pdbqt_content(pdbqt_file, temp_file_path)
            # Move cleaned file to original location
            shutil.move(temp_file_path, pdbqt_file)
            logger.info(f"Cleaned PDBQT formatting in {pdbqt_file}")
    else:
        # Fallback to original implementation
        _clean_pdbqt_formatting_fallback(pdbqt_file)


def _clean_pdbqt_content(input_file: str, output_file: str):
    """Clean PDBQT content and write to output file."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            if line.startswith(('ATOM', 'HETATM')):
                # Extract fields according to PDB/PDBQT format
                record_type = line[0:6].strip()
                atom_num = line[6:11].strip()
                atom_name = line[12:16].strip()
                alt_loc = line[16:17].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                res_num = line[22:26].strip()
                insert_code = line[26:27].strip()
                
                # Skip alternate conformations (A, B, etc.) except the first one
                if alt_loc and alt_loc not in ['', ' ', 'A']:
                    continue
                
                # Fix atom names that are too long
                if len(atom_name) > 4:
                    # Common fixes for hydrogen names
                    if atom_name.startswith('HG1A'):
                        atom_name = 'HG11'
                    elif atom_name.startswith('HG1B'):
                        atom_name = 'HG12'
                    elif atom_name.startswith('HG2A'):
                        atom_name = 'HG21'
                    elif atom_name.startswith('HG2B'):
                        atom_name = 'HG22'
                    elif 'HD1' in atom_name:
                        atom_name = 'HD1'
                    elif 'HD2' in atom_name:
                        atom_name = 'HD2'
                    else:
                        # Generic truncation
                        atom_name = atom_name[:4]
                
                try:
                    # Parse coordinates
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occupancy = line[54:60].strip()
                    temp_factor = line[60:66].strip()
                    charge = line[66:76].strip() if len(line) > 66 else ""
                    element = line[76:78].strip() if len(line) > 76 else ""
                    
                    # Reconstruct the line with proper formatting
                    new_line = f"{record_type:<6}{atom_num:>5} {atom_name:<4}{alt_loc:1}{res_name:>3} {chain_id:1}{res_num:>4}{insert_code:1}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:>6}{temp_factor:>6}{charge:>10}{element:>2}\n"
                    outfile.write(new_line)
                    
                except (ValueError, IndexError) as e:
                    logger.warning(f"Skipping malformed line {line_num}: {line.strip()}")
                    continue
            else:
                # Keep non-ATOM/HETATM lines as-is
                outfile.write(line)


def _clean_pdbqt_formatting_fallback(pdbqt_file):
    """Fallback implementation for PDBQT cleaning without resource management."""
    import tempfile
    import shutil
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pdbqt')
    
    try:
        _clean_pdbqt_content(pdbqt_file, temp_file.name)
        temp_file.close()
        
        # Replace original file with cleaned version
        shutil.move(temp_file.name, pdbqt_file)
        logger.info(f"Cleaned PDBQT formatting in {pdbqt_file}")
        
    except Exception as e:
        logger.error(f"Error cleaning PDBQT file: {e}")
        if os.path.exists(temp_file.name):
            os.unlink(temp_file.name)
        raise


def _simple_pdb_to_pdbqt(pdb_file, pdbqt_file):
    """Improved PDB to PDBQT conversion using RDKit for atom typing and Gasteiger charge calculation."""
    try:
        # Load the PDB file as an RDKit molecule
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            logger.error(f"RDKit failed to load protein from {pdb_file}")
            raise ValueError("Invalid PDB file for fallback conversion")
        
        # Add hydrogens if not present
        mol = Chem.AddHs(mol)
        # Compute Gasteiger charges
        AllChem.ComputeGasteigerCharges(mol)
        
        # Map atomic numbers to basic AutoDock types
        autodock_type_map = {
            6: "C",   # Carbon
            7: "N",   # Nitrogen
            8: "O",   # Oxygen
            16: "S",  # Sulfur
            15: "P",  # Phosphorus
            1: "H",   # Hydrogen
        }
        
        pdbqt_lines = []
        atom_idx = 1
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atom_name = atom.GetSymbol()
            res_name = atom.GetPDBResidueInfo().GetResidueName() if atom.GetPDBResidueInfo() else "UNK"
            chain_id = atom.GetPDBResidueInfo().GetChainId() if atom.GetPDBResidueInfo() else "A"
            res_num = atom.GetPDBResidueInfo().GetResidueNumber() if atom.GetPDBResidueInfo() else 1
            # Get Gasteiger charge
            try:
                charge = float(atom.GetProp('_GasteigerCharge'))
                if np.isnan(charge) or np.isinf(charge):
                    charge = 0.0
            except Exception:
                charge = 0.0
            # Map to AutoDock type
            autodock_type = autodock_type_map.get(atom.GetAtomicNum(), "C")
            # Compose PDBQT line (ATOM/HETATM format)
            pdbqt_line = (
                f"ATOM  {atom_idx:5d} {atom_name:<4} {res_name:>3} {chain_id:1}{res_num:4d}    "
                f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    "
                f"{charge:>6.3f} {autodock_type}\n"
            )
            pdbqt_lines.append(pdbqt_line)
            atom_idx += 1
        # Write out the PDBQT file
        with open(pdbqt_file, 'w') as f:
            f.writelines(pdbqt_lines)
        logger.info(f"Improved PDB to PDBQT conversion completed: {pdbqt_file}")
    except Exception as e:
        logger.error(f"Improved PDB to PDBQT conversion failed: {e}")
        raise


def prepare_protein(protein_file):
    ext = os.path.splitext(protein_file)[1].lower()
    output_file = 'protein_prepped.pdbqt'
    validate_output_dir(output_file)
    if ext == '.pdb':
        try:
            logger.info(f"Preparing protein: cleaning, adding hydrogens, assigning Gasteiger charges, converting to PDBQT...")
            mgltools_pythonsh = config.get_mgltools_pythonsh()
            prepare_script = config.get_mgltools_prepare_script()
            
            # Check if MGLTools is available
            if mgltools_pythonsh and prepare_script and os.path.exists(mgltools_pythonsh) and os.path.exists(prepare_script):
                cmd = [
                    mgltools_pythonsh, prepare_script,
                    '-r', protein_file,
                    '-o', output_file,
                    '-A', 'hydrogens',  # Add hydrogens
                    '-U', 'waters',     # Remove waters
                ]
                logger.info(f"Running MGLTools: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logger.error(f"MGLTools failed: {result.stderr}\n{result.stdout}")
                    if 'IndexError' in result.stderr or 'Unable to assign HAD type' in result.stderr:
                        logger.error("Protein structure appears to be missing atoms or is corrupted. Please repair the PDB using PDBFixer (https://pdbfixer.openmm.org/) and try again.")
                    sys.exit(1)
            else:
                # Fallback for Docker/Linux environments - use meeko for basic conversion
                logger.info("MGLTools not found. Using meeko for basic protein preparation...")
                try:
                    from Bio import PDB
                    from rdkit import Chem
                    from meeko import PDBQTWriterLegacy
                    
                    # Load PDB file and convert to PDBQT
                    mol = Chem.MolFromPDBFile(protein_file, removeHs=False)
                    if mol is None:
                        logger.error(f"Failed to load protein from {protein_file}")
                        sys.exit(1)
                    
                    # Add hydrogens if not present
                    mol = Chem.AddHs(mol)
                    
                    # Write PDBQT using meeko
                    writer = PDBQTWriterLegacy()
                    pdbqt_string = writer.write_string(mol)
                    
                    with open(output_file, 'w') as f:
                        f.write(pdbqt_string)
                    
                    logger.info(f"Protein prepared using meeko fallback method")
                    
                except ImportError as e:
                    logger.error(f"Required packages not available for protein preparation: {e}")
                    logger.error("Please install MGLTools or ensure BioPython and meeko are available")
                    sys.exit(1)
                except Exception as e:
                    logger.error(f"Meeko protein preparation failed: {e}")
                    # Try simple PDB to PDBQT conversion as last resort
                    logger.info("Attempting basic PDB to PDBQT conversion...")
                    try:
                        _simple_pdb_to_pdbqt(protein_file, output_file)
                    except Exception as e2:
                        logger.error(f"All protein preparation methods failed: {e2}")
                        sys.exit(1)
            
            # Clean up PDBQT formatting issues
            logger.info("Cleaning PDBQT formatting...")
            clean_pdbqt_formatting(output_file)
            
            logger.info(f"Protein cleaned (waters removed), hydrogens added, Gasteiger charges assigned, and saved as: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Error during protein preparation: {e}")
            sys.exit(1)
    elif ext == '.pdbqt':
        try:
            if os.path.getsize(protein_file) == 0:
                logger.error(f"PDBQT file '{protein_file}' is empty.")
                sys.exit(1)
            logger.info(f"PDBQT file '{protein_file}' is valid.")
            return protein_file
        except Exception as e:
            logger.error(f"Error validating PDBQT file: {e}")
            sys.exit(1)
    else:
        logger.error(f"Unsupported protein file format: {ext}")
        sys.exit(1)


def batch_prepare_ligands(ligand_files, output_dir='.'):
    """
    Prepare multiple ligands in batch mode for scalability.
    
    Args:
        ligand_files: List of ligand file paths
        output_dir: Directory to save prepared ligands
    
    Returns:
        List of prepared ligand file paths
    """
    prepared_ligands = []
    ensure_dir(output_dir)
    
    for i, ligand_file in enumerate(ligand_files):
        try:
            logger.info(f"Processing ligand {i+1}/{len(ligand_files)}: {ligand_file}")
            
            # Generate unique output filename
            base_name = os.path.splitext(os.path.basename(ligand_file))[0]
            output_file = os.path.join(output_dir, f"{base_name}_prepped.pdbqt")
            
            # Prepare individual ligand
            prepared_ligand = prepare_ligand_single(ligand_file, output_file)
            prepared_ligands.append(prepared_ligand)
            
        except Exception as e:
            logger.error(f"Failed to prepare ligand {ligand_file}: {e}")
            continue
    
    logger.info(f"Successfully prepared {len(prepared_ligands)} out of {len(ligand_files)} ligands")
    return prepared_ligands


def prepare_ligand(ligand_file):
    """Prepare a single ligand using default filename."""
    return prepare_ligand_single(ligand_file, 'ligand_prepped.pdbqt')


def prepare_ligand_single(ligand_file, output_file):
    """Prepare a single ligand file - refactored for batch processing."""
    ext = os.path.splitext(ligand_file)[1].lower()
    validate_output_dir(output_file)
    
    try:
        if ext == '.smi':
            logger.info(f"Ligand input detected as SMILES. Converting to RDKit molecule.")
            with open(ligand_file) as f:
                smiles = f.readline().strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Failed to parse SMILES from {ligand_file}")
                raise ValueError("Invalid SMILES")
            mol = Chem.AddHs(mol)
            logger.info("Generating 3D conformer...")
            conformer_success = (AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == 0)
            if not conformer_success:
                # Try generating tautomers as fallback
                logger.warning("3D conformer generation failed. Attempting to generate tautomers and retry...")
                try:
                    from rdkit.Chem import rdMolStandardize
                    enumerator = rdMolStandardize.TautomerEnumerator()
                    tautomers = enumerator.Enumerate(mol)
                    for t in tautomers:
                        t = Chem.AddHs(t)
                        if AllChem.EmbedMolecule(t, AllChem.ETKDG()) == 0:
                            mol = t
                            conformer_success = True
                            logger.info("Successfully generated conformer for a tautomer.")
                            break
                except Exception as e:
                    logger.warning(f"Tautomer generation failed: {e}")
            if not conformer_success:
                logger.error("3D conformer generation failed for all tautomers.")
                raise ValueError("Failed conformer generation")
            logger.info("Performing energy minimization...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logger.warning("Energy minimization may not have fully converged.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        elif ext == '.sdf':
            logger.info(f"Ligand input detected as SDF. Loading molecule(s).")
            suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
            mols = [m for m in suppl if m is not None]
            if not mols:
                logger.error(f"No valid molecules found in {ligand_file}")
                raise ValueError("No valid molecules")
            mol = mols[0]
            if mol.GetNumConformers() == 0:
                logger.info("Generating 3D conformer for SDF molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logger.error("3D conformer generation failed for SDF molecule.")
                    raise ValueError("Failed conformer generation")
            logger.info("Performing energy minimization for SDF molecule...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logger.warning("Energy minimization may not have fully converged for SDF molecule.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        elif ext == '.mol2':
            logger.info(f"Ligand input detected as MOL2. Loading molecule.")
            mol = Chem.MolFromMol2File(ligand_file, removeHs=False)
            if mol is None:
                logger.error(f"Failed to load molecule from {ligand_file}")
                raise ValueError("Failed to load MOL2")
            if mol.GetNumConformers() == 0:
                logger.info("Generating 3D conformer for MOL2 molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logger.error("3D conformer generation failed for MOL2 molecule.")
                    raise ValueError("Failed conformer generation")
            logger.info("Performing energy minimization for MOL2 molecule...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logger.warning("Energy minimization may not have fully converged for MOL2 molecule.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        else:
            logger.error(f"Unsupported ligand file format: {ext}")
            raise ValueError("Unsupported format")
    except Exception as e:
        logger.error(f"Error during ligand preparation: {e}")
        # For single ligand mode, exit; for batch mode, continue
        if output_file == 'ligand_prepped.pdbqt':
            sys.exit(1)
        else:
            raise


def convert_ligand_to_pdbqt(mol, pdbqt_file):
    if MEEKO_AVAILABLE:
        try:
            # Updated Meeko API - use the new interface only
            prep = MoleculePreparation()
            pdbqt_str = prep.write_string(mol)
            with open(pdbqt_file, 'w') as f:
                f.write(pdbqt_str)
            logger.info(f"Ligand converted and saved as: {pdbqt_file}")
            return
        except Exception as e:
            logger.warning(f"Failed to convert ligand to PDBQT using Meeko: {e}")
            logger.info("Falling back to RDKit-based conversion...")
    
    # Fallback: Use RDKit-based conversion
    try:
        # Add hydrogens if not present
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformer if not present
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        
        # Optimize geometry
        AllChem.UFFOptimizeMolecule(mol)
        
        # Write as PDB first, then convert to PDBQT format
        temp_pdb = pdbqt_file.replace('.pdbqt', '_temp.pdb')
        Chem.MolToPDBFile(mol, temp_pdb)
        
        # Simple PDB to PDBQT conversion
        with open(temp_pdb, 'r') as f_in, open(pdbqt_file, 'w') as f_out:
            for line in f_in:
                if line.startswith(('ATOM', 'HETATM')):
                    # Basic PDBQT formatting
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip()
                    
                    # Determine AutoDock atom type
                    if element == 'C':
                        autodock_type = 'C'
                    elif element == 'N':
                        autodock_type = 'N'
                    elif element == 'O':
                        autodock_type = 'O'
                    elif element == 'S':
                        autodock_type = 'S'
                    elif element == 'P':
                        autodock_type = 'P'
                    elif element == 'H':
                        autodock_type = 'H'
                    else:
                        autodock_type = 'C'
                    
                    # Format as PDBQT
                    pdb_part = line[:76].rstrip()
                    charge = 0.000
                    pdbqt_line = f"{pdb_part}  {charge:>6.3f} {autodock_type}\n"
                    f_out.write(pdbqt_line)
                else:
                    # Keep non-ATOM lines
                    f_out.write(line)
        
        # Clean up temp file
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)
        
        logger.info(f"Ligand converted using RDKit fallback and saved as: {pdbqt_file}")
        
    except Exception as e:
        logger.error(f"Failed to convert ligand to PDBQT using RDKit fallback: {e}")
        raise


def main():
    # Check for MGLTools (optional, but warn if missing)
    if not config.validate_mgltools():
        print("Warning: MGLTools not found. Protein preparation may be limited.")
        print("To install MGLTools, visit: http://mgltools.scripps.edu/downloads/downloads/tools/downloads")
        print("Or set MGLTOOLS_PATH environment variable to point to your installation.")

    parser = argparse.ArgumentParser(
        description="Preprocess protein and ligand files for docking."
    )
    parser.add_argument(
        '--protein', required=True, help="Input protein file (.pdb or .pdbqt)"
    )
    parser.add_argument(
        '--ligand', required=True, help="Input ligand file (.smi, .sdf, or .mol2)"
    )
    parser.add_argument(
        '--batch_ligands', nargs='+', help="Multiple ligand files for batch processing"
    )
    parser.add_argument(
        '--output_dir', default='.', help="Output directory for batch processing"
    )
    args = parser.parse_args()

    logger.info("Starting structure preparation...")

    # Validate input files
    protein_file = validate_file(args.protein, SUPPORTED_PROTEIN_EXT, "Protein")
    
    # Protein preparation
    prepared_protein = prepare_protein(protein_file)
    logger.info(f"Prepared protein file: {prepared_protein}")

    # Ligand preparation - support both single and batch mode
    if args.batch_ligands:
        # Batch mode - prepare multiple ligands
        logger.info(f"Batch mode: preparing {len(args.batch_ligands)} ligands")
        for ligand_file in args.batch_ligands:
            validate_file(ligand_file, SUPPORTED_LIGAND_EXT, "Ligand")
        
        prepared_ligands = batch_prepare_ligands(args.batch_ligands, args.output_dir)
        logger.info(f"Batch preparation complete: {len(prepared_ligands)} ligands prepared")
        
    else:
        # Single ligand mode
        ligand_file = validate_file(args.ligand, SUPPORTED_LIGAND_EXT, "Ligand")
        prepared_ligand = prepare_ligand(ligand_file)
        logger.info(f"Ligand prepared: {prepared_ligand}")

    logger.info("Structure preparation completed successfully")


if __name__ == "__main__":
    main()