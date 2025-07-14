from pydantic import BaseSettings, Field
from typing import Optional

class AppConfig(BaseSettings):
    vina_path: str = Field(..., description="Path to Vina executable")
    mgltools_path: str = Field(..., description="Path to MGLTools directory")
    gnina_path: str = Field(..., description="Path to GNINA executable")
    diffdock_path: str = Field(..., description="Path to DiffDock executable")
    output_dir: str = Field("./output", description="Default output directory")
    exhaustiveness: int = Field(8, description="Docking exhaustiveness")
    num_poses: int = Field(10, description="Number of docking poses")

    class Config:
        env_prefix = ""
        env_file = ".env"
        case_sensitive = False 