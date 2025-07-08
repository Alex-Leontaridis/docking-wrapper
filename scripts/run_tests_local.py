#!/usr/bin/env python3
"""
Local test runner for the docking wrapper project.
Provides an easy interface to run tests with different options.
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path
import logging

# Setup logging if not already present
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')

def run_command(cmd, description):
    """Run a command and handle errors."""
    logging.info(f"\n{'='*60}")
    logging.info(f"Running: {description}")
    logging.info(f"Command: {' '.join(cmd)}")
    logging.info(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        logging.info(f"\n✅ {description} completed successfully!")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"\n❌ {description} failed with exit code {e.returncode}")
        return False
    except FileNotFoundError:
        logging.error(f"\n❌ Command not found: {cmd[0]}")
        logging.error("Please install the required dependencies:")
        logging.error("pip install pytest pytest-cov pytest-mock pytest-asyncio")
        return False

def check_dependencies():
    """Check if required dependencies are installed."""
    required_packages = ['pytest', 'pytest-cov', 'pytest-mock']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package.replace('-', '_'))
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        logging.error("❌ Missing required packages:")
        for package in missing_packages:
            logging.error(f"  - {package}")
        logging.error("\nInstall them with:")
        logging.error(f"pip install {' '.join(missing_packages)}")
        return False
    
    return True

def main():
    """Main test runner."""
    parser = argparse.ArgumentParser(description='Local test runner for docking wrapper')
    parser.add_argument('--type', choices=['unit', 'integration', 'docking', 'all'], 
                       default='all', help='Type of tests to run')
    parser.add_argument('--coverage', action='store_true', 
                       help='Generate coverage report')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Verbose output')
    parser.add_argument('--fast', action='store_true', 
                       help='Skip slow tests')
    parser.add_argument('--platform', choices=['linux', 'windows', 'macos'], 
                       help='Run platform-specific tests only')
    parser.add_argument('--lint', action='store_true', 
                       help='Run linting checks')
    parser.add_argument('--security', action='store_true', 
                       help='Run security checks')
    
    args = parser.parse_args()
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Get project root
    project_root = Path(__file__).parent.parent
    tests_dir = project_root / 'tests'
    
    logging.info("🧪 DOCKING WRAPPER - LOCAL TEST RUNNER")
    logging.info("=" * 60)
    
    success = True
    
    # Build pytest command
    pytest_cmd = ['python', '-m', 'pytest']
    
    if args.verbose:
        pytest_cmd.append('-v')
    else:
        pytest_cmd.append('-q')
    
    # Add test type filters
    if args.type == 'unit':
        pytest_cmd.extend(['-m', 'unit'])
    elif args.type == 'integration':
        pytest_cmd.extend(['-m', 'integration'])
    elif args.type == 'docking':
        pytest_cmd.extend(['-m', 'docking'])
    
    # Add platform filter
    if args.platform:
        pytest_cmd.extend(['-m', args.platform])
    
    # Skip slow tests
    if args.fast:
        pytest_cmd.extend(['-m', 'not slow'])
    
    # Add coverage
    if args.coverage:
        pytest_cmd.extend([
            '--cov=scripts',
            '--cov=utils',
            '--cov-report=term-missing',
            '--cov-report=html:htmlcov'
        ])
    
    # Add test directory
    pytest_cmd.append(str(tests_dir))
    
    # Run tests
    if not run_command(pytest_cmd, f"Running {args.type} tests"):
        success = False
    
    # Run linting
    if args.lint:
        lint_commands = [
            (['black', '--check', '--diff', 'scripts/', 'tests/', 'utils/'], 
             'Black code formatting check'),
            (['flake8', 'scripts/', 'tests/', 'utils/', '--max-line-length=88', '--extend-ignore=E203,W503'], 
             'Flake8 code style check'),
            (['mypy', 'scripts/', 'utils/', '--ignore-missing-imports', '--no-strict-optional'], 
             'MyPy type checking')
        ]
        
        for cmd, description in lint_commands:
            if not run_command(cmd, description):
                success = False
    
    # Run security checks
    if args.security:
        security_commands = [
            (['bandit', '-r', 'scripts/', 'utils/', '-f', 'json', '-o', 'bandit-report.json'], 
             'Bandit security linting'),
            (['safety', 'check', '--json', '--output', 'safety-report.json'], 
             'Safety dependency vulnerability check')
        ]
        
        for cmd, description in security_commands:
            if not run_command(cmd, description):
                success = False
    
    # Summary
    logging.info(f"\n{'='*60}")
    logging.info("TEST SUMMARY")
    logging.info(f"{'='*60}")
    
    if success:
        logging.info("✅ All tests and checks passed!")
        logging.info("\n�� Your code is ready for commit!")
    else:
        logging.error("❌ Some tests or checks failed!")
        logging.error("\n🔧 Please fix the issues before committing.")
    
    if args.coverage:
        logging.info(f"\n📊 Coverage report generated in: {project_root}/htmlcov/")
        logging.info("Open htmlcov/index.html in your browser to view it.")
    
    return 0 if success else 1

if __name__ == '__main__':
    sys.exit(main()) 