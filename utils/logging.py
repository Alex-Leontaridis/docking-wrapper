#!/usr/bin/env python3
"""
Centralized Logging Configuration

Provides consistent logging across all modules with standardized formats,
error codes, and configuration options.
"""

import logging
import sys
import os
from pathlib import Path
from typing import Optional, Dict, Any
from logging.handlers import RotatingFileHandler

# Standard error codes for consistent error reporting
ERROR_CODES = {
    # File/IO errors (1000-1999)
    'FILE_NOT_FOUND': 1001,
    'FILE_PERMISSION_DENIED': 1002,
    'FILE_CORRUPTED': 1003,
    'DIRECTORY_NOT_FOUND': 1004,
    
    # Docking engine errors (2000-2999)
    'VINA_FAILED': 2001,
    'GNINA_FAILED': 2002,
    'DIFFDOCK_FAILED': 2003,
    'ENGINE_NOT_FOUND': 2004,
    'ENGINE_TIMEOUT': 2005,
    
    # Structure preparation errors (3000-3999)
    'PROTEIN_PREP_FAILED': 3001,
    'LIGAND_PREP_FAILED': 3002,
    'PDBQT_CONVERSION_FAILED': 3003,
    'MGLTOOLS_NOT_FOUND': 3004,
    
    # Parsing errors (4000-4999)
    'PARSE_FAILED': 4001,
    'OUTPUT_MALFORMED': 4002,
    'NO_POSES_FOUND': 4003,
    
    # Configuration errors (5000-5999)
    'CONFIG_INVALID': 5001,
    'MISSING_PARAMETER': 5002,
    'INVALID_VALUE': 5003,
    
    # System errors (9000-9999)
    'SYSTEM_ERROR': 9001,
    'MEMORY_ERROR': 9002,
    'TIMEOUT_ERROR': 9003,
}

class DockingLogger:
    """Centralized logger for the docking wrapper."""
    
    def __init__(self, name: str = 'DockingWrapper', 
                 log_file: Optional[str] = None,
                 level: str = 'INFO',
                 max_file_size_mb: int = 100,
                 backup_count: int = 5):
        """
        Initialize the logger.
        
        Args:
            name: Logger name
            log_file: Path to log file (optional)
            level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
            max_file_size_mb: Maximum log file size in MB
            backup_count: Number of backup files to keep
        """
        self.name = name
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, level.upper()))
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Create formatters
        self._create_formatters()
        
        # Add handlers
        self._add_console_handler()
        if log_file:
            self._add_file_handler(log_file, max_file_size_mb, backup_count)
        
        # Set specific logger levels for noisy libraries
        self._configure_third_party_loggers()
    
    def _create_formatters(self):
        """Create standardized formatters."""
        # Detailed formatter for file logging
        self.detailed_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)-8s | %(name)-20s | %(filename)s:%(lineno)d | %(message)s'
        )
        
        # Simple formatter for console output
        self.simple_formatter = logging.Formatter(
            '%(levelname)s: %(message)s'
        )
        
        # Error formatter with error codes
        self.error_formatter = logging.Formatter(
            '[ERROR-%(error_code)s] %(asctime)s | %(name)s | %(message)s'
        )
    
    def _add_console_handler(self):
        """Add console handler."""
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(self.simple_formatter)
        console_handler.setLevel(logging.INFO)
        self.logger.addHandler(console_handler)
    
    def _add_file_handler(self, log_file: str, max_file_size_mb: int, backup_count: int):
        """Add file handler with rotation."""
        # Ensure log directory exists
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=max_file_size_mb * 1024 * 1024,
            backupCount=backup_count
        )
        file_handler.setFormatter(self.detailed_formatter)
        file_handler.setLevel(logging.DEBUG)
        self.logger.addHandler(file_handler)
    
    def _configure_third_party_loggers(self):
        """Configure third-party loggers to reduce noise."""
        noisy_loggers = [
            'rdkit',
            'urllib3',
            'requests',
            'matplotlib',
            'PIL',
            'tensorflow',
            'torch'
        ]
        
        for logger_name in noisy_loggers:
            logging.getLogger(logger_name).setLevel(logging.WARNING)
    
    def info(self, message: str, **kwargs):
        """Log info message."""
        self.logger.info(message, **kwargs)
    
    def warning(self, message: str, **kwargs):
        """Log warning message."""
        self.logger.warning(message, **kwargs)
    
    def error(self, message: str, error_code: Optional[str] = None, **kwargs):
        """Log error message with optional error code."""
        if error_code:
            extra = {'error_code': ERROR_CODES.get(error_code, 'UNKNOWN')}
            kwargs['extra'] = extra
            self.logger.error(f"[{error_code}] {message}", **kwargs)
        else:
            self.logger.error(message, **kwargs)
    
    def debug(self, message: str, **kwargs):
        """Log debug message."""
        self.logger.debug(message, **kwargs)
    
    def critical(self, message: str, error_code: Optional[str] = None, **kwargs):
        """Log critical message with optional error code."""
        if error_code:
            extra = {'error_code': ERROR_CODES.get(error_code, 'UNKNOWN')}
            kwargs['extra'] = extra
            self.logger.critical(f"[{error_code}] {message}", **kwargs)
        else:
            self.logger.critical(message, **kwargs)

def setup_logging(name: str = 'DockingWrapper',
                 log_file: Optional[str] = None,
                 level: str = 'INFO',
                 max_file_size_mb: int = 100,
                 backup_count: int = 5) -> DockingLogger:
    """
    Setup centralized logging.
    
    Args:
        name: Logger name
        log_file: Path to log file (optional)
        level: Logging level
        max_file_size_mb: Maximum log file size in MB
        backup_count: Number of backup files to keep
    
    Returns:
        Configured DockingLogger instance
    """
    return DockingLogger(
        name=name,
        log_file=log_file,
        level=level,
        max_file_size_mb=max_file_size_mb,
        backup_count=backup_count
    )

def get_logger(name: str = None) -> DockingLogger:
    """
    Get an existing logger or create a new one.
    
    Args:
        name: Logger name (optional)
    
    Returns:
        DockingLogger instance
    """
    if name is None:
        name = 'DockingWrapper'
    
    # Check if logger already exists
    existing_logger = logging.getLogger(name)
    if existing_logger.handlers:
        # Return existing logger
        return DockingLogger(name)
    else:
        # Create new logger with defaults
        return setup_logging(name)

# Global logger instance
_global_logger = None

def get_global_logger() -> DockingLogger:
    """Get the global logger instance."""
    global _global_logger
    if _global_logger is None:
        _global_logger = setup_logging()
    return _global_logger

def set_global_logger(logger: DockingLogger):
    """Set the global logger instance."""
    global _global_logger
    _global_logger = logger 