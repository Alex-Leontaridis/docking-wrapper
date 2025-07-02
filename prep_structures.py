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
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

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


def prepare_protein(protein_file):
    ext = os.path.splitext(protein_file)[1].lower()
    output_file = 'protein_prepped.pdbqt'
    validate_output_dir(output_file)
    if ext == '.pdb':
        try:
            logging.info(f"Preparing protein: cleaning, adding hydrogens, assigning Gasteiger charges, converting to PDBQT...")
            mgltools_pythonsh = os.path.expanduser('~/mgltools_1.5.7_MacOS-X/bin/pythonsh')
            prepare_script = os.path.expanduser('~/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py')
            cmd = [
                mgltools_pythonsh, prepare_script,
                '-r', protein_file,
                '-o', output_file,
                '-A', 'hydrogens',  # Add hydrogens
                '-U', 'waters',     # Remove waters
            ]
            logging.info(f"Running: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logging.error(f"MGLTools failed: {result.stderr}\n{result.stdout}")
                if 'IndexError' in result.stderr or 'Unable to assign HAD type' in result.stderr:
                    logging.error("Protein structure appears to be missing atoms or is corrupted. Please repair the PDB using PDBFixer (https://pdbfixer.openmm.org/) and try again.")
                sys.exit(1)
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


def prepare_ligand(ligand_file):
    ext = os.path.splitext(ligand_file)[1].lower()
    output_file = 'ligand_prepped.pdbqt'
    validate_output_dir(output_file)
    try:
        if ext == '.smi':
            logging.info(f"Ligand input detected as SMILES. Converting to RDKit molecule.")
            with open(ligand_file) as f:
                smiles = f.readline().strip().split()[0]
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logging.error(f"Failed to parse SMILES from {ligand_file}")
                sys.exit(1)
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
                sys.exit(1)
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
                sys.exit(1)
            mol = mols[0]
            if mol.GetNumConformers() == 0:
                logging.info("Generating 3D conformer for SDF molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logging.error("3D conformer generation failed for SDF molecule.")
                    sys.exit(1)
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
                sys.exit(1)
            if mol.GetNumConformers() == 0:
                logging.info("Generating 3D conformer for MOL2 molecule...")
                if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                    logging.error("3D conformer generation failed for MOL2 molecule.")
                    sys.exit(1)
            logging.info("Performing energy minimization for MOL2 molecule...")
            if AllChem.UFFOptimizeMolecule(mol) != 0:
                logging.warning("Energy minimization may not have fully converged for MOL2 molecule.")
            convert_ligand_to_pdbqt(mol, output_file)
            return output_file
        else:
            logging.error(f"Unsupported ligand file format: {ext}")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error during ligand preparation: {e}")
        sys.exit(1)


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
        '-protein', required=True, help="Input protein file (.pdb or .pdbqt)"
    )
    parser.add_argument(
        '-ligand', required=True, help="Input ligand file (.smi, .sdf, or .mol2)"
    )
    args = parser.parse_args()

    # Validate input files
    protein_file = validate_file(args.protein, SUPPORTED_PROTEIN_EXT, "Protein")
    ligand_file = validate_file(args.ligand, SUPPORTED_LIGAND_EXT, "Ligand")

    # Protein preparation
    prepared_protein = prepare_protein(protein_file)
    logging.info(f"Prepared protein file: {prepared_protein}")

    # Ligand preparation
    prepared_ligand = prepare_ligand(ligand_file)
    logging.info(f"Ligand prepared: {prepared_ligand}")

if __name__ == "__main__":
    main() 