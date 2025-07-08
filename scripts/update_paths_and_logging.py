#!/usr/bin/env python3
"""
Script to update all Python files to use the new path manager and consistent logging system.
This script removes hardcoded paths and standardizes logging across the project.
"""

import os
import re
import sys
from pathlib import Path
from typing import List, Dict, Tuple

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

def find_python_files(directory: Path) -> List[Path]:
    """Find all Python files in the directory and subdirectories."""
    python_files = []
    for root, dirs, files in os.walk(directory):
        # Skip certain directories
        dirs[:] = [d for d in dirs if not d.startswith('.') and d not in ['__pycache__', 'venv', 'env', 'node_modules']]
        
        for file in files:
            if file.endswith('.py'):
                python_files.append(Path(root) / file)
    
    return python_files

def update_imports(content: str) -> str:
    """Update import statements to use new path manager and logging."""
    
    # Replace old logging imports
    old_logging_imports = [
        r'import logging',
        r'from logging import .*',
        r'logging\.basicConfig\(.*\)',
        r'logging\.getLogger\(.*\)',
    ]
    
    # Add new imports if not present
    if 'from utils.path_manager import' not in content:
        content = re.sub(
            r'(import os\n)',
            r'\1from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir\n',
            content
        )
    
    if 'from utils.logging import' not in content:
        content = re.sub(
            r'(import logging\n)',
            r'\1from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context\n',
            content
        )
    
    return content

def update_hardcoded_paths(content: str) -> str:
    """Replace hardcoded paths with path manager calls."""
    
    # Common hardcoded path patterns
    path_replacements = [
        # os.path.join patterns
        (r'os\.path\.join\(os\.getcwd\(\), \'([^\']+)\'\)', r'get_path_manager().join("\1")'),
        (r'os\.path\.join\(os\.getcwd\(\), "([^"]+)"\)', r'get_path_manager().join("\1")'),
        
        # Direct path assignments
        (r'(\w+)_dir = ["\']([^"\']+)["\']', r'\1_dir = get_path("\2")'),
        (r'(\w+)_path = ["\']([^"\']+)["\']', r'\1_path = get_path("\2")'),
        
        # os.path.abspath patterns
        (r'os\.path\.abspath\(["\']([^"\']+)["\']\)', r'get_absolute_path("\1")'),
        
        # Directory creation patterns
        (r'os\.makedirs\(([^,]+), exist_ok=True\)', r'ensure_dir(\1)'),
        (r'Path\(([^)]+)\)\.mkdir\(exist_ok=True\)', r'ensure_dir(\1)'),
    ]
    
    for pattern, replacement in path_replacements:
        content = re.sub(pattern, replacement, content)
    
    return content

def update_logging_calls(content: str) -> str:
    """Update logging calls to use the new logging system."""
    
    # Replace basic logging setup
    content = re.sub(
        r'logging\.basicConfig\([^)]+\)',
        'logger = setup_logging(__name__)',
        content
    )
    
    # Replace logger creation
    content = re.sub(
        r'logger = logging\.getLogger\([\'"]([^\'"]+)[\'"]\)',
        r'logger = setup_logging("\1")',
        content
    )
    
    # Replace direct logging calls with logger calls
    logging_calls = [
        (r'logging\.info\(([^)]+)\)', r'logger.info(\1)'),
        (r'logging\.error\(([^)]+)\)', r'logger.error(\1)'),
        (r'logging\.warning\(([^)]+)\)', r'logger.warning(\1)'),
        (r'logging\.debug\(([^)]+)\)', r'logger.debug(\1)'),
        (r'logging\.critical\(([^)]+)\)', r'logger.critical(\1)'),
    ]
    
    for pattern, replacement in logging_calls:
        content = re.sub(pattern, replacement, content)
    
    return content

def update_file(file_path: Path) -> Tuple[bool, List[str]]:
    """Update a single Python file."""
    changes = []
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        
        # Skip files that are already updated
        if 'from utils.path_manager import' in content and 'from utils.logging import' in content:
            return False, []
        
        # Apply updates
        content = update_imports(content)
        content = update_hardcoded_paths(content)
        content = update_logging_calls(content)
        
        # Check if content changed
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            
            changes.append(f"Updated {file_path}")
            return True, changes
        
        return False, changes
        
    except Exception as e:
        changes.append(f"Error updating {file_path}: {e}")
        return False, changes

def main():
    """Main function to update all Python files."""
    print("Updating Python files to use new path manager and logging system...")
    
    # Find all Python files
    python_files = find_python_files(project_root)
    
    # Filter out files we don't want to modify
    exclude_patterns = [
        'utils/path_manager.py',
        'utils/logging.py',
        'scripts/update_paths_and_logging.py',
        '__pycache__',
        'venv',
        'env',
        'node_modules',
        'Umol/',  # External library
        'EquiBind/',  # External library
    ]
    
    filtered_files = []
    for file_path in python_files:
        file_str = str(file_path.relative_to(project_root))
        if not any(pattern in file_str for pattern in exclude_patterns):
            filtered_files.append(file_path)
    
    print(f"Found {len(filtered_files)} Python files to update")
    
    # Update files
    updated_count = 0
    all_changes = []
    
    for file_path in filtered_files:
        updated, changes = update_file(file_path)
        if updated:
            updated_count += 1
        all_changes.extend(changes)
    
    # Print summary
    print(f"\nUpdated {updated_count} files:")
    for change in all_changes:
        print(f"  {change}")
    
    print(f"\nTotal files processed: {len(filtered_files)}")
    print(f"Files updated: {updated_count}")
    print(f"Files unchanged: {len(filtered_files) - updated_count}")

if __name__ == "__main__":
    main() 