import pytest
from pathlib import Path
import tempfile

# 1. Component Registration and Discovery
def test_component_registry_registration():
    from core import component_registry, ComponentType
    engines = component_registry.get_components_by_type(ComponentType.DOCKING_ENGINE)
    models = component_registry.get_components_by_type(ComponentType.ML_MODEL)
    analyzers = component_registry.get_components_by_type(ComponentType.ANALYZER)
    assert any(e.name == "vina" for e in engines)
    assert any(m.name == "equibind" for m in models)
    assert any(a.name == "plip" for a in analyzers)

# 2. Component Availability
def test_component_availability():
    from core import component_registry
    for name in ["vina", "gnina", "diffdock", "equibind", "umol", "plip", "fpocket"]:
        comp = component_registry.get_component(name)
        assert comp is not None
        status = comp.check_availability()
        assert status in comp.get_status().__class__

# 3. Docking Engine Interface
def test_docking_engine_run(tmp_path):
    from core import component_registry, DockingResult
    vina = component_registry.get_component("vina")
    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    receptor.write_text("RECEPTOR DUMMY")
    ligand.write_text("LIGAND DUMMY")
    result = vina.run(receptor, ligand, output_dir=tmp_path, center=(0,0,0), size=(10,10,10))
    assert isinstance(result, DockingResult)
    assert hasattr(result, "success")

# 4. ML Model Interface
def test_ml_model_predict(tmp_path):
    from core import component_registry, PredictionResult
    equibind = component_registry.get_component("equibind")
    pdb = tmp_path / "complex.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N")
    result = equibind.predict(pdb, output_dir=tmp_path)
    assert hasattr(result, "success")

# 5. Analyzer Interface
def test_analyzer_analyze(tmp_path):
    from core import component_registry, AnalysisResult
    plip = component_registry.get_component("plip")
    pose = tmp_path / "pose.pdb"
    pose.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N")
    result = plip.analyze(pose, output_dir=tmp_path)
    assert hasattr(result, "success")

# 6. Pipeline Integration Test
def test_pipeline_integration(tmp_path):
    from scripts.batch_pipeline import BatchDockingPipeline
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "ligand.sdf"
    protein.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N")
    ligand.write_text("$$$$\n")
    pipeline = BatchDockingPipeline(config_path=None, output_dir=str(tmp_path))
    result = pipeline.process_single_ligand(("ligand", ligand), str(protein))
    assert result["stages"]["preparation"] is not False

# 7. Error Handling
def test_pipeline_missing_ligand(tmp_path):
    from scripts.batch_pipeline import BatchDockingPipeline
    protein = tmp_path / "protein.pdb"
    protein.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N")
    pipeline = BatchDockingPipeline(config_path=None, output_dir=str(tmp_path))
    result = pipeline.process_single_ligand(("missing_ligand", tmp_path / "missing.sdf"), str(protein))
    assert result["stages"]["preparation"] is False
    assert "Ligand file not found" in result["errors"]["preparation"]

# 8. Registry Dynamic Selection
def test_pipeline_dynamic_selection(tmp_path):
    from scripts.batch_pipeline import BatchDockingPipeline
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "ligand.sdf"
    protein.write_text("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N")
    ligand.write_text("$$$$\n")
    pipeline = BatchDockingPipeline(config_path=None, output_dir=str(tmp_path))
    pipeline.config["engines"]["vina"]["enabled"] = False
    pipeline.config["engines"]["gnina"]["enabled"] = False
    pipeline.config["engines"]["diffdock"]["enabled"] = False
    result = pipeline.process_single_ligand(("ligand", ligand), str(protein))
    assert "No docking engines or ML models enabled" in result["errors"]["docking"] 