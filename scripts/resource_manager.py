#!/usr/bin/env python3
"""
Resource management utilities for molecular docking pipeline.
Handles temporary file cleanup and resource tracking.
"""

import os
import tempfile
import shutil
import atexit
import logging
import threading
from pathlib import Path
from typing import List, Set, Optional, Callable
from contextlib import contextmanager


class ResourceManager:
    """Manages temporary files and resources with automatic cleanup."""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the resource manager.
        
        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self._temp_files: Set[str] = set()
        self._temp_dirs: Set[str] = set()
        self._cleanup_callbacks: List[Callable] = []
        self._lock = threading.Lock()
        
        # Register cleanup on exit
        atexit.register(self.cleanup_all)
    
    def register_temp_file(self, file_path: str) -> None:
        """
        Register a temporary file for cleanup.
        
        Args:
            file_path: Path to the temporary file
        """
        with self._lock:
            self._temp_files.add(file_path)
            self.logger.debug(f"Registered temp file for cleanup: {file_path}")
    
    def register_temp_dir(self, dir_path: str) -> None:
        """
        Register a temporary directory for cleanup.
        
        Args:
            dir_path: Path to the temporary directory
        """
        with self._lock:
            self._temp_dirs.add(dir_path)
            self.logger.debug(f"Registered temp directory for cleanup: {dir_path}")
    
    def register_cleanup_callback(self, callback: Callable) -> None:
        """
        Register a cleanup callback function.
        
        Args:
            callback: Function to call during cleanup
        """
        with self._lock:
            self._cleanup_callbacks.append(callback)
            self.logger.debug(f"Registered cleanup callback: {callback.__name__}")
    
    def cleanup_file(self, file_path: str) -> bool:
        """
        Clean up a specific temporary file.
        
        Args:
            file_path: Path to the file to clean up
            
        Returns:
            True if cleanup was successful, False otherwise
        """
        try:
            if os.path.exists(file_path):
                os.unlink(file_path)
                self.logger.debug(f"Cleaned up temp file: {file_path}")
                return True
        except Exception as e:
            self.logger.warning(f"Failed to clean up temp file {file_path}: {e}")
            return False
        finally:
            with self._lock:
                self._temp_files.discard(file_path)
        return True
    
    def cleanup_directory(self, dir_path: str) -> bool:
        """
        Clean up a specific temporary directory.
        
        Args:
            dir_path: Path to the directory to clean up
            
        Returns:
            True if cleanup was successful, False otherwise
        """
        try:
            if os.path.exists(dir_path):
                shutil.rmtree(dir_path)
                self.logger.debug(f"Cleaned up temp directory: {dir_path}")
                return True
        except Exception as e:
            self.logger.warning(f"Failed to clean up temp directory {dir_path}: {e}")
            return False
        finally:
            with self._lock:
                self._temp_dirs.discard(dir_path)
        return True
    
    def cleanup_all(self) -> None:
        """Clean up all registered temporary files and directories."""
        self.logger.info("Starting cleanup of all temporary resources")
        
        # Run cleanup callbacks
        with self._lock:
            callbacks = self._cleanup_callbacks.copy()
        
        for callback in callbacks:
            try:
                callback()
                self.logger.debug(f"Executed cleanup callback: {callback.__name__}")
            except Exception as e:
                self.logger.warning(f"Cleanup callback {callback.__name__} failed: {e}")
        
        # Clean up temporary files
        with self._lock:
            temp_files = self._temp_files.copy()
            temp_dirs = self._temp_dirs.copy()
        
        for file_path in temp_files:
            self.cleanup_file(file_path)
        
        for dir_path in temp_dirs:
            self.cleanup_directory(dir_path)
        
        self.logger.info("Completed cleanup of all temporary resources")
    
    def get_registered_resources(self) -> dict:
        """
        Get information about currently registered resources.
        
        Returns:
            Dictionary with resource information
        """
        with self._lock:
            return {
                'temp_files': list(self._temp_files),
                'temp_dirs': list(self._temp_dirs),
                'cleanup_callbacks': len(self._cleanup_callbacks)
            }


# Global resource manager instance
_global_resource_manager = ResourceManager()


@contextmanager
def temp_file(suffix: str = '', prefix: str = 'tmp', delete: bool = True):
    """
    Context manager for temporary files with automatic cleanup.
    
    Args:
        suffix: File suffix
        prefix: File prefix
        delete: Whether to delete the file on exit
        
    Yields:
        Path to the temporary file
    """
    temp_file_path = None
    try:
        temp_file_path = tempfile.NamedTemporaryFile(
            mode='w', 
            suffix=suffix, 
            prefix=prefix, 
            delete=False
        ).name
        _global_resource_manager.register_temp_file(temp_file_path)
        yield temp_file_path
    finally:
        if delete and temp_file_path and os.path.exists(temp_file_path):
            _global_resource_manager.cleanup_file(temp_file_path)


@contextmanager
def temp_directory():
    """
    Context manager for temporary directories with automatic cleanup.
    
    Yields:
        Path to the temporary directory
    """
    temp_dir_path = None
    try:
        temp_dir_path = tempfile.mkdtemp()
        _global_resource_manager.register_temp_dir(temp_dir_path)
        yield temp_dir_path
    finally:
        if temp_dir_path and os.path.exists(temp_dir_path):
            _global_resource_manager.cleanup_directory(temp_dir_path)


def safe_cleanup(func: Callable) -> Callable:
    """
    Decorator to ensure cleanup happens even if function fails.
    
    Args:
        func: Function to wrap
        
    Returns:
        Wrapped function
    """
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.error(f"Function {func.__name__} failed: {e}")
            raise
        finally:
            # Trigger cleanup on function exit
            _global_resource_manager.cleanup_all()
    
    return wrapper


class FileBackup:
    """Context manager for backing up files before modification."""
    
    def __init__(self, file_path: str, backup_suffix: str = '.backup'):
        """
        Initialize file backup.
        
        Args:
            file_path: Path to the file to backup
            backup_suffix: Suffix for the backup file
        """
        self.file_path = file_path
        self.backup_path = f"{file_path}{backup_suffix}"
        self.backup_created = False
    
    def __enter__(self):
        """Create backup of the file."""
        if os.path.exists(self.file_path):
            shutil.copy2(self.file_path, self.backup_path)
            self.backup_created = True
            logging.debug(f"Created backup: {self.backup_path}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Restore backup if an exception occurred."""
        if exc_type is not None and self.backup_created:
            logging.warning(f"Restoring backup due to exception: {exc_type}")
            shutil.copy2(self.backup_path, self.file_path)
            logging.info(f"Restored file from backup: {self.file_path}")
        
        # Clean up backup file
        if self.backup_created and os.path.exists(self.backup_path):
            os.unlink(self.backup_path)
            logging.debug(f"Removed backup file: {self.backup_path}")


def cleanup_on_exit():
    """Ensure cleanup happens on program exit."""
    _global_resource_manager.cleanup_all()


# Register cleanup on exit
atexit.register(cleanup_on_exit)


# Convenience functions
def register_temp_file(file_path: str) -> None:
    """Register a temporary file for cleanup."""
    _global_resource_manager.register_temp_file(file_path)


def register_temp_dir(dir_path: str) -> None:
    """Register a temporary directory for cleanup."""
    _global_resource_manager.register_temp_dir(dir_path)


def cleanup_all() -> None:
    """Clean up all registered resources."""
    _global_resource_manager.cleanup_all()


def get_resource_info() -> dict:
    """Get information about registered resources."""
    return _global_resource_manager.get_registered_resources() 