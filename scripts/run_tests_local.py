#!/usr/bin/env python3
"""
Local Test Runner

Runs tests and checks locally before committing.
"""

import os
from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
import sys
import subprocess
import importlib
from pathlib import Path
from typing import List, Dict, Any

# Add the project root to the path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from utils.logging import setup_logging

def run_command(cmd: List[str], description: str) -> bool:
    """Run a command and return success status."""
    logger = setup_logging('TestRunner')
    
    logger.info(f"\n{'='*60}")
    logger.info(f"Running: {description}")
    logger.info(f"Command: {' '.join(cmd)}")
    logger.info(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info(f"\n✅ {description} completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"\n❌ {description} failed with exit code {e.returncode}")
        logger.error(f"Error output: {e.stderr}")
        return False
    except FileNotFoundError:
        logger.error(f"\n❌ Command not found: {cmd[0]}")
        logger.error("Please install the required dependencies:")
        logger.error("pip install pytest pytest-cov pytest-mock pytest-asyncio")
        return False

def check_python_packages() -> bool:
    """Check if required Python packages are installed."""
    logger = setup_logging('TestRunner')
    
    required_packages = [
        'pytest',
        'pytest-cov', 
        'pytest-mock',
        'pytest-asyncio',
        'rdkit',
        'numpy',
        'pandas',
        'scikit-learn',
        'meeko',
        'biopython'
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            importlib.import_module(package.replace('-', '_'))
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        logger.error("❌ Missing required packages:")
        for package in missing_packages:
            logger.error(f"  - {package}")
        logger.error("\nInstall them with:")
        logger.error(f"pip install {' '.join(missing_packages)}")
        return False
    
    logger.info("✅ All required Python packages are installed")
    return True

def run_tests():
    """Run all tests and checks."""
    logger = setup_logging('TestRunner')
    
    logger.info("🧪 DOCKING WRAPPER - LOCAL TEST RUNNER")
    logger.info("=" * 60)
    
    # Check if we're in the right directory
    if not os.path.exists('pytest.ini'):
        logger.error("❌ pytest.ini not found. Please run this script from the project root.")
        sys.exit(1)
    
    # Check Python packages
    if not check_python_packages():
        sys.exit(1)
    
    # Define test commands
    test_commands = [
        (['python', '-m', 'pytest', 'tests/', '-v'], 'Unit Tests'),
        (['python', '-m', 'pytest', 'tests/', '--cov=scripts', '--cov=utils', '--cov-report=html'], 'Coverage Report'),
        (['python', '-m', 'pytest', 'tests/', '--cov=scripts', '--cov=utils', '--cov-report=term-missing'], 'Coverage Summary'),
        (['python', 'scripts/check_dependencies.py'], 'Dependency Check'),
        (['python', '-c', 'from config import Config; print("Config loaded successfully")'], 'Config Test'),
    ]
    
    # Run tests
    results = []
    for cmd, description in test_commands:
        success = run_command(cmd, description)
        results.append((description, success))
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("TEST SUMMARY")
    logger.info(f"{'='*60}")
    
    all_passed = all(success for _, success in results)
    
    if all_passed:
        logger.info("✅ All tests and checks passed!")
        logger.info("\n🚀 Your code is ready for commit!")
    else:
        logger.error("❌ Some tests or checks failed!")
        logger.error("\n🔧 Please fix the issues before committing.")
        
        # Show failed tests
        failed_tests = [desc for desc, success in results if not success]
        logger.error("Failed tests:")
        for test in failed_tests:
            logger.error(f"  - {test}")
    
    # Coverage report location
    if os.path.exists('htmlcov/'):
        logger.info(f"\n📊 Coverage report generated in: {project_root}/htmlcov/")
        logger.info("Open htmlcov/index.html in your browser to view it.")
    
    return all_passed

def main():
    """Main function."""
    success = run_tests()
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main() 