#!/usr/bin/env python3
"""
Centralized Logging Configuration

Provides consistent logging across all modules with standardized formats,
error codes, and configuration options. Integrates with path management.
"""

import structlog
import logging
import sys
import os
import time
from pathlib import Path
from typing import Optional, Dict, Any, Union
from logging.handlers import RotatingFileHandler
import threading
from contextlib import contextmanager

# Import path manager (lazy import to avoid circular dependency)
PATH_MANAGER_AVAILABLE = None

def _get_path_manager():
    """Lazy import of path manager to avoid circular dependency."""
    global PATH_MANAGER_AVAILABLE
    if PATH_MANAGER_AVAILABLE is None:
        try:
            from .path_manager import get_path_manager
            PATH_MANAGER_AVAILABLE = get_path_manager
        except ImportError:
            PATH_MANAGER_AVAILABLE = False
    return PATH_MANAGER_AVAILABLE if PATH_MANAGER_AVAILABLE is not False else None

# Standard error codes for consistent error reporting
ERROR_CODES = {
    # File/IO errors (1000-1999)
    'FILE_NOT_FOUND': 1001,
    'FILE_PERMISSION_DENIED': 1002,
    'FILE_CORRUPTED': 1003,
    'DIRECTORY_NOT_FOUND': 1004,
    'PATH_NOT_FOUND': 1005,
    
    # Docking engine errors (2000-2999)
    'VINA_FAILED': 2001,
    'GNINA_FAILED': 2002,
    'DIFFDOCK_FAILED': 2003,
    'ENGINE_NOT_FOUND': 2004,
    'ENGINE_TIMEOUT': 2005,
    'ENGINE_CONFIG_ERROR': 2006,
    
    # Structure preparation errors (3000-3999)
    'PROTEIN_PREP_FAILED': 3001,
    'LIGAND_PREP_FAILED': 3002,
    'PDBQT_CONVERSION_FAILED': 3003,
    'MGLTOOLS_NOT_FOUND': 3004,
    'MOLECULE_INVALID': 3005,
    
    # Parsing errors (4000-4999)
    'PARSE_FAILED': 4001,
    'OUTPUT_MALFORMED': 4002,
    'NO_POSES_FOUND': 4003,
    'FORMAT_UNSUPPORTED': 4004,
    
    # Configuration errors (5000-5999)
    'CONFIG_INVALID': 5001,
    'MISSING_PARAMETER': 5002,
    'INVALID_VALUE': 5003,
    'CONFIG_NOT_FOUND': 5004,
    
    # ML/Analysis errors (6000-6999)
    'ML_MODEL_FAILED': 6001,
    'ANALYSIS_FAILED': 6002,
    'MODEL_NOT_FOUND': 6003,
    'PREDICTION_FAILED': 6004,
    
    # System errors (9000-9999)
    'SYSTEM_ERROR': 9001,
    'MEMORY_ERROR': 9002,
    'TIMEOUT_ERROR': 9003,
    'PERMISSION_ERROR': 9004,
    'NETWORK_ERROR': 9005,
}

class DockingLogger:
    """Centralized logger for the docking wrapper."""
    
    def __init__(self, name: str = 'DockingWrapper', 
                 log_file: Optional[str] = None,
                 level: str = 'INFO',
                 max_file_size_mb: int = 100,
                 backup_count: int = 5,
                 console_output: bool = True,
                 detailed_format: bool = True):
        """
        Initialize the logger.
        
        Args:
            name: Logger name
            log_file: Path to log file (optional)
            level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
            max_file_size_mb: Maximum log file size in MB
            backup_count: Number of backup files to keep
            console_output: Whether to output to console
            detailed_format: Whether to use detailed formatting
        """
        self.name = name
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, level.upper()))
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Create formatters
        self._create_formatters(detailed_format)
        
        # Add handlers
        if console_output:
            self._add_console_handler()
        
        if log_file:
            self._add_file_handler(log_file, max_file_size_mb, backup_count)
        
        # Set specific logger levels for noisy libraries
        self._configure_third_party_loggers()
        
        # Performance tracking
        self._performance_data = {}
        self._lock = threading.Lock()
    
    def _create_formatters(self, detailed_format: bool):
        """Create standardized formatters."""
        if detailed_format:
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
        else:
            # Simple formatters
            self.detailed_formatter = logging.Formatter(
                '%(asctime)s | %(levelname)s | %(message)s'
            )
            self.simple_formatter = logging.Formatter(
                '%(levelname)s: %(message)s'
            )
            self.error_formatter = logging.Formatter(
                '[ERROR-%(error_code)s] %(message)s'
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
            'torch',
            'openmm',
            'mdtraj',
            'numpy',
            'scipy',
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
    
    def performance(self, operation: str, duration: float, **kwargs):
        """Log performance metrics."""
        with self._lock:
            if operation not in self._performance_data:
                self._performance_data[operation] = []
            self._performance_data[operation].append(duration)
        
        self.logger.info(f"Performance: {operation} took {duration:.3f}s", **kwargs)
    
    def get_performance_stats(self, operation: str) -> Dict[str, float]:
        """Get performance statistics for an operation."""
        with self._lock:
            if operation not in self._performance_data:
                return {}
            
            durations = self._performance_data[operation]
            if not durations:
                return {}
            
            return {
                'count': len(durations),
                'total': sum(durations),
                'mean': sum(durations) / len(durations),
                'min': min(durations),
                'max': max(durations),
            }
    
    @contextmanager
    def timer(self, operation: str):
        """Context manager for timing operations."""
        start_time = time.time()
        try:
            yield
        finally:
            duration = time.time() - start_time
            self.performance(operation, duration)
    
    def section(self, title: str, level: str = 'info'):
        """Log a section header."""
        separator = '=' * 60
        message = f"\n{separator}\n{title}\n{separator}"
        getattr(self, level)(message)
    
    def subsection(self, title: str, level: str = 'info'):
        """Log a subsection header."""
        separator = '-' * 40
        message = f"\n{separator}\n{title}\n{separator}"
        getattr(self, level)(message)


def setup_logging():
    logging.basicConfig(
        format="%(message)s",
        stream=sys.stdout,
        level=logging.INFO,
    )
    structlog.configure(
        processors=[
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.stdlib.add_log_level,
            structlog.stdlib.add_logger_name,
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.JSONRenderer(),
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )
    return structlog.get_logger()


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


def log_function_call(func):
    """Decorator to log function calls with timing."""
    def wrapper(*args, **kwargs):
        logger = get_logger('FunctionCall')
        func_name = func.__name__
        
        with logger.timer(f"function_call:{func_name}"):
            logger.debug(f"Calling {func_name} with args={args}, kwargs={kwargs}")
            try:
                result = func(*args, **kwargs)
                logger.debug(f"{func_name} returned successfully")
                return result
            except Exception as e:
                logger.error(f"{func_name} failed: {e}", error_code='SYSTEM_ERROR')
                raise
    
    return wrapper


def log_class_methods(cls):
    """Class decorator to add logging to all methods."""
    for attr_name in dir(cls):
        attr = getattr(cls, attr_name)
        if callable(attr) and not attr_name.startswith('_'):
            setattr(cls, attr_name, log_function_call(attr))
    return cls


# Convenience functions for common logging patterns
def log_startup(module_name: str, version: str = "unknown"):
    """Log module startup information."""
    logger = get_logger(module_name)
    logger.section(f"Starting {module_name} (version: {version})")
    
    try:
        path_manager_func = _get_path_manager()
        if path_manager_func:
            path_manager = path_manager_func()
            platform_info = path_manager.get_platform_info()
            logger.info(f"Platform: {platform_info['system']} {platform_info['release']}")
            logger.info(f"Python: {platform_info['python_version']}")
            
            tools = path_manager.list_available_tools()
            logger.info("Available tools:")
            for tool, available in tools.items():
                status = "✓" if available else "✗"
                logger.info(f"  {status} {tool}")
    except Exception as e:
        logger.warning(f"Could not get platform info: {e}")


def log_shutdown(module_name: str):
    """Log module shutdown information."""
    logger = get_logger(module_name)
    logger.section(f"Shutting down {module_name}")
    
    # Log performance statistics if available
    if hasattr(logger, 'get_performance_stats'):
        logger.info("Performance summary:")
        for operation in ['function_call', 'docking', 'preparation', 'parsing']:
            stats = logger.get_performance_stats(operation)
            if stats:
                logger.info(f"  {operation}: {stats['count']} calls, "
                          f"avg {stats['mean']:.3f}s, "
                          f"total {stats['total']:.3f}s")


def log_error_with_context(error: Exception, context: str = "", error_code: str = "SYSTEM_ERROR"):
    """Log an error with context information."""
    logger = get_logger('ErrorHandler')
    message = f"{context}: {str(error)}" if context else str(error)
    logger.error(message, error_code=error_code)
    
    # Log traceback in debug mode
    import traceback
    logger.debug(f"Traceback: {traceback.format_exc()}")


# Global error handler
def setup_global_error_handler():
    """Setup global error handling for uncaught exceptions."""
    def handle_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            # Don't log keyboard interrupts
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        
        logger = get_logger('GlobalErrorHandler')
        logger.critical("Uncaught exception", error_code='SYSTEM_ERROR')
        logger.debug(f"Exception type: {exc_type}")
        logger.debug(f"Exception value: {exc_value}")
        logger.debug(f"Traceback: {exc_traceback}")
        
        # Call the original exception hook
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
    
    sys.excepthook = handle_exception 