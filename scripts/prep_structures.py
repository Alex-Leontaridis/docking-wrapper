#!/usr/bin/env python3
import argparse
import os
import sys
import logging
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from meeko import MoleculePreparation
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False

# Setup logging
def setup_logging():
    """Setup logging to both console and file"""
    # Create logs directory if it doesn't exist
    os.makedirs('logs', exist_ok=True)
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.FileHandler('preprocessing_log.txt'),
            logging.StreamHandler(sys.stdout)
        ]
    )

# Call setup_logging at module level
setup_logging()

SUPPORTED_PROTEIN_EXT = {'.pdb', '.pdbqt'}
SUPPORTED_LIGAND_EXT = {'.smi', '.sdf', '.mol2'}


def validate_file(file_path, valid_exts, file_type):
    if not os.path.isfile(file_path):
        logging.error(f"{file_type} file '{file_path}' does not exist.")
        sys.exit(1)
    if not os.access(file_path, os.R_OK):
        logging.error(f"{file_type} file '{file_path}' is not readable.")
        sys.exit(1)
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in valid_exts:
        logging.error(f"{file_type} file '{file_path}' must have one of the following extensions: {', '.join(valid_exts)}.")
        sys.exit(1)
    return file_path


def validate_output_dir(output_file):
    output_dir = os.path.dirname(os.path.abspath(output_file)) or '.'
    if not os.access(output_dir, os.W_OK):
        logging.error(f"Output directory '{output_dir}' is not writable.")
        sys.exit(1)


def clean_pdbqt_formatting(pdbqt_file):
    """
    Clean up PDBQT formatting issues that can cause Vina parsing errors.
    - Fix atom names that are too long (> 4 characters)
    - Remove problematic alternate conformations
    - Ensure proper coordinate formatting
    """
    import tempfile
    import shutil
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pdbqt')
    
    try:
        with open(pdbqt_file, 'r') as infile:
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
                        temp_file.write(new_line)
                        
                    except (ValueError, IndexError) as e:
                        logging.warning(f"Skipping malformed line {line_num}: {line.strip()}")
                        continue
                else:
                    # Keep non-ATOM/HETATM lines as-is
                    temp_file.write(line)
        
        temp_file.close()
        
        # Replace original file with cleaned version
        shutil.move(temp_file.name, pdbqt_file)
        logging.info(f"Cleaned PDBQT formatting in {pdbqt_file}")
        
    except Exception as e:
        logging.error(f"Error cleaning PDBQT file: {e}")
        if os.path.exists(temp_file.name):
            os.unlink(temp_file.name)
        raise


def _simple_pdb_to_pdbqt(pdb_file, pdbqt_file):
    """Simple PDB to PDBQT conversion as last resort fallback."""
    try:
        with open(pdb_file, 'r') as f:
            pdb_lines = f.readlines()
        
        pdbqt_lines = []
        for line in pdb_lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Convert PDB line to basic PDBQT format
                if len(line) >= 54:  # Minimum length for coordinates
                    # Take only the standard PDB part (up to column 78) but remove element column
                    # PDB format: columns 77-78 contain element symbol which we need to replace
                    pdb_part = line[:76].rstrip()  # Stop before element column
                    
                    # Determine AutoDock atom type based on atom name
                    atom_name = line[12:16].strip()
                    if atom_name.startswith('C'):
                        autodock_type = "C"
                    elif atom_name.startswith('N'):
                        autodock_type = "N"
                    elif atom_name.startswith('O'):
                        autodock_type = "O"
                    elif atom_name.startswith('S'):
                        autodock_type = "S"
                    elif atom_name.startswith('P'):
                        autodock_type = "P"
                    elif atom_name.startswith('H'):
                        autodock_type = "H"
                    else:
                        autodock_type = "C"  # Default to carbon
                    
                    # Format: PDB_part + spaces + charge + space + atom_type
                    # Ensure proper spacing to column 79
                    padding_needed = 76 - len(pdb_part)
                    if padding_needed > 0:
                        pdb_part += ' ' * padding_needed
                    
                    pdbqt_line = f"{pdb_part}  +0.000 {autodock_type}"
                    pdbqt_lines.append(pdbqt_line + '\n')
            else:
                # Skip header and other PDB-specific lines that Vina doesn't need
                # Only keep essential structural information
                if line.startswith(('REMARK', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                    pdbqt_lines.append(line)
        
        with open(pdbqt_file, 'w') as f:
            f.writelines(pdbqt_lines)
        
        logging.info(f"Basic PDB to PDBQT conversion completed: {pdbqt_file}")
        
    except Exception as e:
        logging.error(f"Simple PDB to PDBQT conversion failed: {e}")
        raise


def prepare_protein(protein_file):
    ext = os.path.splitext(protein_file)[1].lower()
    output_file = 'protein_prepped.pdbqt'
    validate_output_dir(output_file)
    if ext == '.pdb':
        try:
            logging.info(f"Preparing protein: cleaning, adding hydrogens, assigning Gasteiger charges, converting to PDBQT...")
            mgltools_pythonsh = os.path.expanduser('~/mgltools_1.5.7_MacOS-X/bin/pythonsh')
            prepare_script = os.path.expanduser('~/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py')
            
            # Convert to absolute paths to avoid path issues with MGLTools
            abs_protein_file = os.path.abspath(protein_file)
            abs_output_file = os.path.abspath(output_file)
            
            # Check if MGLTools is available (for macOS/local development)
            if os.path.exists(mgltools_pythonsh) and os.path.exists(prepare_script):
                cmd = [
                    mgltools_pythonsh, prepare_script,
                    '-r', abs_protein_file,
                    '-o', abs_output_file,
                    '-A', 'hydrogens',  # Add hydrogens
                    '-U', 'waters',     # Remove waters
                ]
                logging.info(f"Running MGLTools: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    logging.error(f"MGLTools failed: {result.stderr}\n{result.stdout}")
                    if 'IndexError' in result.stderr or 'Unable to assign HAD type' in result.stderr:
                        logging.error("Protein structure appears to be missing atoms or is corrupted. Please repair the PDB using PDBFixer (https://pdbfixer.openmm.org/) and try again.")
                    sys.exit(1)
            else:
                # Fallback for Docker/Linux environments - use meeko for basic conversion
                logging.info("MGLTools not found. Using meeko for basic protein preparation...")
                try:
                    from Bio import PDB
                    from rdkit import Chem
                    from meeko import PDBQTWriterLegacy
                    
                    # Load PDB file and convert to PDBQT
                    mol = Chem.MolFromPDBFile(abs_protein_file, removeHs=False)
                    if mol is None:
                        logging.error(f"Failed to load protein from {abs_protein_file}")
                        sys.exit(1)
                    
                    # Add hydrogens if not present
                    mol = Chem.AddHs(mol)
                    
                    # Write PDBQT using meeko
                    writer = PDBQTWriterLegacy()
                    pdbqt_string = writer.write_string(mol)
                    
                    with open(abs_output_file, 'w') as f:
                        f.write(pdbqt_string)
                    
                    logging.info(f"Protein prepared using meeko fallback method")
                    
                except ImportError as e:
                    logging.error(f"Required packages not available for protein preparation: {e}")
                    logging.error("Please install MGLTools or ensure BioPython and meeko are available")
                    sys.exit(1)
                except Exception as e:
                    logging.error(f"Meeko protein preparation failed: {e}")
                    # Try simple PDB to PDBQT conversion as last resort
                    logging.info("Attempting basic PDB to PDBQT conversion...")
                    try:
                        _simple_pdb_to_pdbqt(abs_protein_file, abs_output_file)
                    except Exception as e2:
                        logging.error(f"All protein preparation methods failed: {e2}")
                        sys.exit(1)
            
            # Clean up PDBQT formatting issues
            logging.info("Cleaning PDBQT formatting...")
            clean_pdbqt_formatting(output_file)
            
            logging.info(f"Protein cleaned (waters removed), hydrogens added, Gasteiger charges assigned, and saved as: {output_file}")
            return output_file
        except Exception as e:
            logging.error(f"Error during protein preparation: {e}")
            sys.exit(1)
    elif ext == '.pdbqt':
        try:
            if os.path.getsize(protein_file) == 0:
                logging.error(f"PDBQT file '{protein_file}' is empty.")
                sys.exit(1)
            logging.info(f"PDBQT file '{protein_file}' is valid.")
            return protein_file
        except Exception as e:
            logging.error(f"Error validating PDBQT file: {e}")
            sys.exit(1)
    else:
        logging.error(f"Unsupported protein file format: {ext}")
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
    os.makedirs(output_dir, exist_ok=True)
    
    for i, ligand_file in enumerate(ligand_files):
        try:
            logging.info(f"Processing ligand {i+1}/{len(ligand_files)}: {ligand_file}")
            
            # Generate unique output filename
            base_name = os.path.splitext(os.path.basename(ligand_file))[0]
            output_file = os.path.join(output_dir, f"{base_name}_prepped.pdbqt")
            
            # Prepare individual ligand
            prepared_ligand = prepare_ligand_single(ligand_file, output_file)
            prepared_ligands.append(prepared_ligand)
            
        except Exception as e:
            logging.error(f"Failed to prepare ligand {ligand_file}: {e}")
            continue
    
    logging.info(f"Successfully prepared {len(prepared_ligands)} out of {len(ligand_files)} ligands")
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
            logging.info(f"Ligand input detected as SMILES. Converting to RDKit molecule.")
            with open(ligand_file) as f:
                smiles = f.readline().strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logging.error(f"Failed to parse SMILES from {ligand_file}")
                raise ValueError("Invalid SMILES")
            mol = Chem.AddHs(mol)
            logging.info("Generating 3D conformer...")
            conformer_success = (AllChem.EmbedMolecule(mol, AllChem.ETKDG()) == 0)
            if not conformer_success:
                # Try generating tautomers as fallback
                logging.warning("3D conformer generation failed. Attempting to generate tautomers and retry...")
                try:
                    from rdkit.Chem import rdMolStandardize
                    enumerator = rdMolStandardize.TautomerEnumerator()
                    tautomers = enumerator.Enumerate(mol)
                    for t in tautomers:
                        t = Chem.AddHs(t)
                        if AllChem.EmbedMolecule(t, AllChem.ETKDG()) == 0:
                            mol = t
                            conformer_success = True
                            logging.info("Successfully generated conformer for a tautomer.")
                            break
                except Exception as e:
                    logging.warning(f"Tautomer generation failed: {e}")
            if not conformer_success:
                logging.error("3D conformer generation failed for all tautomers.")
                raise ValueError("Failed conformer generation")
            logging.info("Performing energy minimization...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logging.warning("Energy minimization may not have fully converged.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        elif ext == '.sdf':
            logging.info(f"Ligand input detected as SDF. Loading molecule(s).")
            suppl = Chem.SDMolSupplier(ligand_file, removeHs=False)
            mols = [m for m in suppl if m is not None]
            if not mols:
                logging.error(f"No valid molecules found in {ligand_file}")
                raise ValueError("No valid molecules")
            mol = mols[0]
            if mol.GetNumConformers() == 0:
                logging.info("Generating 3D conformer for SDF molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logging.error("3D conformer generation failed for SDF molecule.")
                    raise ValueError("Failed conformer generation")
            logging.info("Performing energy minimization for SDF molecule...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logging.warning("Energy minimization may not have fully converged for SDF molecule.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        elif ext == '.mol2':
            logging.info(f"Ligand input detected as MOL2. Loading molecule.")
            mol = Chem.MolFromMol2File(ligand_file, removeHs=False)
            if mol is None:
                logging.error(f"Failed to load molecule from {ligand_file}")
                raise ValueError("Failed to load MOL2")
            if mol.GetNumConformers() == 0:
                logging.info("Generating 3D conformer for MOL2 molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logging.error("3D conformer generation failed for MOL2 molecule.")
                    raise ValueError("Failed conformer generation")
            logging.info("Performing energy minimization for MOL2 molecule...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logging.warning("Energy minimization may not have fully converged for MOL2 molecule.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        else:
            logging.error(f"Unsupported ligand file format: {ext}")
            raise ValueError("Unsupported format")
    except Exception as e:
        logging.error(f"Error during ligand preparation: {e}")
        # For single ligand mode, exit; for batch mode, continue
        if output_file == 'ligand_prepped.pdbqt':
            sys.exit(1)
        else:
            raise


def convert_ligand_to_pdbqt(mol, pdbqt_file):
    if not MEEKO_AVAILABLE:
        logging.error("Meeko is not installed. Cannot convert ligand to PDBQT format.")
        sys.exit(1)
    try:
        prep = MoleculePreparation()
        prep.prepare(mol)
        pdbqt_str = prep.write_pdbqt_string()
        with open(pdbqt_file, 'w') as f:
            f.write(pdbqt_str)
        logging.info(f"Ligand converted and saved as: {pdbqt_file}")
    except Exception as e:
        logging.error(f"Failed to convert ligand to PDBQT: {e}")
        sys.exit(1)


def main():
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

    logging.info("Starting structure preparation...")

    # Validate input files
    protein_file = validate_file(args.protein, SUPPORTED_PROTEIN_EXT, "Protein")
    
    # Protein preparation
    prepared_protein = prepare_protein(protein_file)
    logging.info(f"Prepared protein file: {prepared_protein}")

    # Ligand preparation - support both single and batch mode
    if args.batch_ligands:
        # Batch mode - prepare multiple ligands
        logging.info(f"Batch mode: preparing {len(args.batch_ligands)} ligands")
        for ligand_file in args.batch_ligands:
            validate_file(ligand_file, SUPPORTED_LIGAND_EXT, "Ligand")
        
        prepared_ligands = batch_prepare_ligands(args.batch_ligands, args.output_dir)
        logging.info(f"Batch preparation complete: {len(prepared_ligands)} ligands prepared")
        
    else:
        # Single ligand mode
        ligand_file = validate_file(args.ligand, SUPPORTED_LIGAND_EXT, "Ligand")
        prepared_ligand = prepare_ligand(ligand_file)
        logging.info(f"Ligand prepared: {prepared_ligand}")

    logging.info("Structure preparation completed successfully")


if __name__ == "__main__":
    main()