#!/usr/bin/env python3
"""
Core Interfaces for Docking System Components

This module defines abstract base classes for all docking system components,
providing a unified interface for different implementations.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
import logging
from dataclasses import dataclass
from enum import Enum


class ComponentType(Enum):
    """Enumeration of component types."""
    DOCKING_ENGINE = "docking_engine"
    ML_MODEL = "ml_model"
    ANALYZER = "analyzer"


class ComponentStatus(Enum):
    """Enumeration of component status."""
    AVAILABLE = "available"
    UNAVAILABLE = "unavailable"
    ERROR = "error"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class DockingResult:
    """Standardized result structure for docking operations."""
    success: bool
    poses: List[Path]  # Paths to pose files
    scores: List[float]  # Docking scores
    confidence: Optional[float] = None
    metadata: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    execution_time: Optional[float] = None


@dataclass
class PredictionResult:
    """Standardized result structure for ML model predictions."""
    success: bool
    predictions: Dict[str, Any]  # Model-specific predictions
    confidence: Optional[float] = None
    metadata: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    execution_time: Optional[float] = None


@dataclass
class AnalysisResult:
    """Standardized result structure for analysis operations."""
    success: bool
    metrics: Dict[str, Any]  # Analysis metrics
    visualizations: Optional[List[Path]] = None  # Paths to visualization files
    metadata: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
    execution_time: Optional[float] = None


class BaseComponent(ABC):
    """Base class for all docking system components."""
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None):
        """Initialize component with name and configuration."""
        self.name = name
        self.config = config or {}
        self.logger = logging.getLogger(f"{self.__class__.__name__}.{name}")
        self.status = ComponentStatus.UNAVAILABLE
        
    @abstractmethod
    def check_availability(self) -> ComponentStatus:
        """Check if the component is available and ready to use."""
        pass
    
    @abstractmethod
    def get_capabilities(self) -> Dict[str, Any]:
        """Get component capabilities and supported features."""
        pass
    
    def validate_inputs(self, **kwargs) -> bool:
        """Validate input parameters for the component."""
        return True
    
    def get_status(self) -> ComponentStatus:
        """Get current component status."""
        return self.status
    
    def set_status(self, status: ComponentStatus):
        """Set component status."""
        self.status = status
        self.logger.debug(f"Component status changed to: {status.value}")


class DockingEngine(BaseComponent):
    """Abstract base class for docking engines."""
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.component_type = ComponentType.DOCKING_ENGINE
    
    @abstractmethod
    def run(self, receptor: Path, ligand: Path, **kwargs) -> DockingResult:
        """
        Run docking calculation.
        
        Args:
            receptor: Path to receptor file
            ligand: Path to ligand file
            **kwargs: Additional parameters specific to the engine
            
        Returns:
            DockingResult with poses and scores
        """
        pass
    
    @abstractmethod
    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get supported file formats for receptor and ligand."""
        pass
    
    def prepare_inputs(self, receptor: Path, ligand: Path, **kwargs) -> Dict[str, Path]:
        """Prepare input files for docking (can be overridden by subclasses)."""
        return {"receptor": receptor, "ligand": ligand}


class MLModel(BaseComponent):
    """Abstract base class for machine learning models."""
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.component_type = ComponentType.ML_MODEL
    
    @abstractmethod
    def predict(self, complex_pdb: Path, **kwargs) -> PredictionResult:
        """
        Run prediction on protein-ligand complex.
        
        Args:
            complex_pdb: Path to protein-ligand complex PDB file
            **kwargs: Additional parameters specific to the model
            
        Returns:
            PredictionResult with model predictions
        """
        pass
    
    @abstractmethod
    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load the ML model."""
        pass
    
    def get_model_info(self) -> Dict[str, Any]:
        """Get information about the loaded model."""
        return {
            "name": self.name,
            "type": self.component_type.value,
            "loaded": self.status == ComponentStatus.AVAILABLE
        }


class Analyzer(BaseComponent):
    """Abstract base class for analysis tools."""
    
    def __init__(self, name: str, config: Optional[Dict[str, Any]] = None):
        super().__init__(name, config)
        self.component_type = ComponentType.ANALYZER
    
    @abstractmethod
    def analyze(self, pose_file: Path, **kwargs) -> AnalysisResult:
        """
        Analyze docking pose or complex.
        
        Args:
            pose_file: Path to pose file (PDB, SDF, etc.)
            **kwargs: Additional parameters specific to the analyzer
            
        Returns:
            AnalysisResult with analysis metrics
        """
        pass
    
    @abstractmethod
    def get_supported_analyses(self) -> List[str]:
        """Get list of supported analysis types."""
        pass
    
    def validate_pose_file(self, pose_file: Path) -> bool:
        """Validate if the pose file is supported by this analyzer."""
        if not pose_file.exists():
            return False
        
        # Basic validation - can be overridden by subclasses
        return pose_file.suffix.lower() in ['.pdb', '.sdf', '.mol2', '.pdbqt']


class ComponentRegistry:
    """Registry for managing docking system components."""
    
    def __init__(self):
        self.components: Dict[str, BaseComponent] = {}
        self.logger = logging.getLogger(__name__)
    
    def register(self, component: BaseComponent):
        """Register a component in the registry."""
        self.components[component.name] = component
        self.logger.info(f"Registered component: {component.name} ({component.component_type.value})")
    
    def get_component(self, name: str) -> Optional[BaseComponent]:
        """Get a component by name."""
        return self.components.get(name)
    
    def get_components_by_type(self, component_type: ComponentType) -> List[BaseComponent]:
        """Get all components of a specific type."""
        return [comp for comp in self.components.values() if comp.component_type == component_type]
    
    def list_components(self) -> Dict[ComponentType, List[str]]:
        """List all registered components by type."""
        result = {comp_type: [] for comp_type in ComponentType}
        for component in self.components.values():
            result[component.component_type].append(component.name)
        return result
    
    def check_all_availability(self) -> Dict[str, ComponentStatus]:
        """Check availability of all registered components."""
        status = {}
        for name, component in self.components.items():
            try:
                status[name] = component.check_availability()
            except Exception as e:
                self.logger.error(f"Error checking availability of {name}: {e}")
                status[name] = ComponentStatus.ERROR
        return status


# Global component registry
component_registry = ComponentRegistry() 