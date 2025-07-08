# Testing Guide for Docking Wrapper

This document provides comprehensive information about the testing setup for the docking wrapper project.

## Overview

The project uses a multi-layered testing approach:

- **Unit Tests**: Test individual functions and classes
- **Integration Tests**: Test the complete pipeline with mocked external tools
- **GitHub Actions**: Automated CI/CD testing on multiple platforms
- **Security Tests**: Vulnerability scanning and dependency checks

## Test Structure

```
tests/
├── __init__.py
├── conftest.py              # Pytest configuration and fixtures
├── test_config.py           # Configuration system tests
├── test_docking_results_parser.py  # Results parsing tests
├── test_integration.py      # Integration tests
├── test_pdbqt_conversion.py # PDBQT conversion tests
└── README.md               # This file
```

## Running Tests

### Local Testing

1. **Install test dependencies**:
   ```bash
   pip install -r requirements.txt
   pip install pytest pytest-cov pytest-mock pytest-asyncio
   ```

2. **Run all tests**:
   ```bash
   # Using pytest
   pytest tests/ -v
   
   # Using the test runner
   python run_tests.py
   ```

3. **Run specific test categories**:
   ```bash
   # Unit tests only
   pytest tests/ -v -m "unit"
   
   # Integration tests only
   pytest tests/ -v -m "integration"
   
   # Docking-specific tests
   pytest tests/ -v -m "docking"
   
   # Exclude slow tests
   pytest tests/ -v -m "not slow"
   ```

4. **Run with coverage**:
   ```bash
   pytest tests/ --cov=scripts --cov=utils --cov-report=html
   ```

### GitHub Actions Testing

The project includes several GitHub Actions workflows:

1. **Main Tests** (`.github/workflows/tests.yml`):
   - Runs on multiple Python versions (3.8-3.10)
   - Tests on Ubuntu and Windows
   - Includes linting and security checks
   - Generates coverage reports

2. **Docking Tests** (`.github/workflows/docking-tests.yml`):
   - Specialized for docking functionality
   - Uses mocked external tools
   - Performance testing
   - Integration testing

3. **Dependencies** (`.github/workflows/dependencies.yml`):
   - Weekly dependency updates
   - Security scanning
   - License checking

## Test Categories

### Unit Tests (`-m "unit"`)
- Test individual functions and methods
- Mock external dependencies
- Fast execution
- High coverage

### Integration Tests (`-m "integration"`)
- Test complete workflows
- Use mocked external tools
- Test data flow between components
- Verify end-to-end functionality

### Docking Tests (`-m "docking"`)
- Test docking-specific functionality
- Mock Vina, GNINA, DiffDock outputs
- Test result parsing
- Test PDBQT conversion

### Slow Tests (`-m "slow"`)
- Tests that take longer to run
- May involve file I/O or complex calculations
- Can be excluded with `-m "not slow"`

## Test Fixtures

The `conftest.py` file provides several useful fixtures:

- `sample_protein_pdb`: Sample protein PDB file
- `sample_ligand_sdf`: Sample ligand SDF file
- `mock_vina_output`: Mock Vina docking results
- `mock_gnina_output`: Mock GNINA docking results
- `mock_diffdock_output`: Mock DiffDock confidence scores
- `mock_subprocess`: Mock subprocess calls
- `test_data_dir`: Temporary directory for test data

## Mocking Strategy

Since the project depends on external docking tools (Vina, GNINA, DiffDock), we use extensive mocking:

1. **External Tools**: Mock subprocess calls to external binaries
2. **File I/O**: Use temporary files and directories
3. **Configuration**: Mock configuration files and environment variables
4. **Network Calls**: Mock any API calls or downloads

## Coverage Requirements

- **Minimum Coverage**: 70% (set in `pytest.ini`)
- **Coverage Reports**: Generated as HTML and XML
- **Coverage Upload**: Automatically uploaded to Codecov

## Continuous Integration

### Triggers
- Push to main/master/develop branches
- Pull requests to main/master/develop branches
- Manual workflow dispatch
- Scheduled dependency updates (weekly)

### Platforms
- Ubuntu Latest
- Windows Latest
- Python 3.8, 3.9, 3.10

### Artifacts
- Coverage reports
- Security scan reports
- Dependency trees
- Test logs

## Debugging Tests

### Running Tests in Debug Mode
```bash
pytest tests/ -v -s --pdb
```

### Viewing Test Output
```bash
pytest tests/ -v -s
```

### Running Single Test
```bash
pytest tests/test_integration.py::TestIntegration::test_complete_pipeline_flow -v
```

## Adding New Tests

### Test Naming Convention
- Test files: `test_*.py`
- Test classes: `Test*`
- Test methods: `test_*`

### Test Structure
```python
import pytest
from unittest.mock import patch

@pytest.mark.unit
def test_function_name():
    """Test description."""
    # Arrange
    # Act
    # Assert
    pass

@pytest.mark.integration
class TestFeature:
    def test_feature_workflow(self, sample_protein_pdb, mock_subprocess):
        """Test complete feature workflow."""
        # Test implementation
        pass
```

### Test Markers
Use appropriate markers for test categorization:
- `@pytest.mark.unit`: Unit tests
- `@pytest.mark.integration`: Integration tests
- `@pytest.mark.docking`: Docking-specific tests
- `@pytest.mark.slow`: Slow-running tests

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure project root is in Python path
2. **Missing Dependencies**: Install all requirements
3. **Permission Errors**: Check file permissions on test files
4. **Mock Issues**: Verify mock setup in fixtures

### Environment Setup
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install pytest pytest-cov pytest-mock

# Run tests
pytest tests/ -v
```

## Performance Considerations

- Unit tests should run in < 1 second
- Integration tests should run in < 30 seconds
- Use `@pytest.mark.slow` for longer tests
- Mock expensive operations (file I/O, network calls)

## Security Testing

The CI pipeline includes:
- **Bandit**: Security linting
- **Safety**: Dependency vulnerability scanning
- **License Check**: License compliance verification

## Contributing

When adding new features:
1. Write tests first (TDD approach)
2. Ensure all tests pass
3. Maintain or improve coverage
4. Add appropriate test markers
5. Update this documentation if needed 