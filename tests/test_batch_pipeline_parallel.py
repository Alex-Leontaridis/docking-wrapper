import os
import json
import pytest
from scripts.batch_pipeline import BatchDockingPipeline

@pytest.fixture
def setup_test_inputs(tmp_path):
    # Create dummy protein and ligand files
    protein = tmp_path / "protein.pdbqt"
    protein.write_text("REMARK DUMMY PROTEIN")
    ligands = []
    for i in range(3):
        ligand = tmp_path / f"ligand{i}.pdbqt"
        ligand.write_text("REMARK DUMMY LIGAND")
        ligands.append(ligand)
    # Add one invalid ligand (file does not exist)
    invalid_ligand = tmp_path / "badligand.pdbqt"
    return str(protein), [str(l) for l in ligands] + [str(invalid_ligand)], tmp_path

def write_ligand_list(ligands, path):
    ligand_list_file = path / "ligands.txt"
    with open(ligand_list_file, "w") as f:
        for l in ligands:
            f.write(f"{l}\n")
    return ligand_list_file

def write_config(path, max_workers=2):
    config = {
        "engines": {"vina": {"enabled": False}, "gnina": {"enabled": False}, "diffdock": {"enabled": False}},
        "ml_models": {},
        "parallel": {"max_workers": max_workers},
        "box": {"auto_detect": False, "default_size": [20, 20, 20]},
        "analysis": {}
    }
    config_file = path / "config.json"
    with open(config_file, "w") as f:
        json.dump(config, f)
    return config_file

def test_parallel_batch_pipeline_success(setup_test_inputs):
    protein, ligands, tmp_path = setup_test_inputs
    ligand_list_file = write_ligand_list(ligands, tmp_path)
    config_file = write_config(tmp_path, max_workers=2)
    pipeline = BatchDockingPipeline(config_path=str(config_file), output_dir=str(tmp_path))
    result = pipeline.run_batch_pipeline(protein, str(ligand_list_file), parallel=True)
    assert result["total_ligands"] == 4
    assert result["failed_ligands"] == 1
    assert result["successful_ligands"] == 3
    assert os.path.exists(os.path.join(tmp_path, "final_summary.csv"))

def test_serial_batch_pipeline_success(setup_test_inputs):
    protein, ligands, tmp_path = setup_test_inputs
    ligand_list_file = write_ligand_list(ligands, tmp_path)
    config_file = write_config(tmp_path, max_workers=1)
    pipeline = BatchDockingPipeline(config_path=str(config_file), output_dir=str(tmp_path))
    result = pipeline.run_batch_pipeline(protein, str(ligand_list_file), parallel=False)
    assert result["total_ligands"] == 4
    assert result["failed_ligands"] == 1
    assert result["successful_ligands"] == 3
    assert os.path.exists(os.path.join(tmp_path, "final_summary.csv")) 