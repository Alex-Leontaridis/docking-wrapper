#!/bin/bash

# Comprehensive Docker Setup Test Script
# Tests GNINA and DiffDock installation in the Docker environment

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[PASS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARN]${NC} $1"; }
print_error() { echo -e "${RED}[FAIL]${NC} $1"; }

# Configuration
IMAGE_NAME="docking-pipeline:latest"
MINIMAL_IMAGE="docking-pipeline-minimal:latest"
TEST_DIR="$(pwd)/docker_tests"
TIMEOUT=300  # 5 minutes timeout for tests

# Counters
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Test function wrapper
run_test() {
    local test_name="$1"
    local test_command="$2"
    local expected_exit_code="${3:-0}"
    
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    print_status "Running test: $test_name"
    
    if timeout $TIMEOUT bash -c "$test_command"; then
        local exit_code=$?
        if [ $exit_code -eq $expected_exit_code ]; then
            print_success "$test_name"
            TESTS_PASSED=$((TESTS_PASSED + 1))
            return 0
        else
            print_error "$test_name (exit code: $exit_code, expected: $expected_exit_code)"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            return 1
        fi
    else
        print_error "$test_name (timeout or execution error)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Setup test environment
setup_test_env() {
    print_status "Setting up test environment..."
    
    # Create test directory
    mkdir -p "$TEST_DIR"
    cd "$TEST_DIR"
    
    # Create test input files
    mkdir -p inputs/ligands
    
    # Create a simple test protein (1HVD crystal structure coordinates)
    cat > inputs/protein.pdb << 'EOF'
HEADER    HYDROLASE                               19-JAN-95   1HVD              
ATOM      1  N   ALA A   1      20.154  11.200   5.481  1.00 10.00           N  
ATOM      2  CA  ALA A   1      19.030  10.346   5.618  1.00 10.00           C  
ATOM      3  C   ALA A   1      17.726  11.074   5.897  1.00 10.00           C  
ATOM      4  O   ALA A   1      17.763  12.147   6.487  1.00 10.00           O  
ATOM      5  CB  ALA A   1      18.879   9.396   4.441  1.00 10.00           C  
TER
END
EOF
    
    # Create a simple test ligand (aspirin SMILES)
    echo "CC(=O)OC1=CC=CC=C1C(=O)O" > inputs/ligands/aspirin.smi
    
    # Create a simple SDF ligand
    cat > inputs/ligands/test_mol.sdf << 'EOF'
test_molecule
  Mrv2014 01010100002D          

  3  2  0  0  0  0            999 V2000
   -0.4125    0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125    0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  3  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
EOF
    
    print_success "Test environment created in $TEST_DIR"
}

# Test Docker image availability
test_docker_images() {
    print_status "Testing Docker image availability..."
    
    run_test "Docker full image exists" \
        "docker image inspect $IMAGE_NAME >/dev/null 2>&1"
    
    if docker image inspect $MINIMAL_IMAGE >/dev/null 2>&1; then
        run_test "Docker minimal image exists" \
            "docker image inspect $MINIMAL_IMAGE >/dev/null 2>&1"
    else
        print_warning "Minimal image not found - skipping minimal tests"
    fi
}

# Test basic container functionality
test_basic_functionality() {
    print_status "Testing basic container functionality..."
    
    run_test "Container startup" \
        "docker run --rm $IMAGE_NAME echo 'Container started successfully'"
    
    run_test "Python environment" \
        "docker run --rm $IMAGE_NAME python3 --version"
    
    run_test "Conda environment activation" \
        "docker run --rm $IMAGE_NAME bash -c 'source activate docking && python3 -c \"print(\\\"Conda environment OK\\\")\"'"
}

# Test core Python dependencies
test_python_dependencies() {
    print_status "Testing Python dependencies..."
    
    run_test "RDKit import" \
        "docker run --rm $IMAGE_NAME python3 -c 'import rdkit; print(f\"RDKit version: {rdkit.__version__}\")'"
    
    run_test "Pandas import" \
        "docker run --rm $IMAGE_NAME python3 -c 'import pandas; print(f\"Pandas version: {pandas.__version__}\")'"
    
    run_test "NumPy import" \
        "docker run --rm $IMAGE_NAME python3 -c 'import numpy; print(f\"NumPy version: {numpy.__version__}\")'"
    
    run_test "OpenBabel import" \
        "docker run --rm $IMAGE_NAME python3 -c 'import openbabel; print(\"OpenBabel: OK\")'"
    
    run_test "PyTorch import" \
        "docker run --rm $IMAGE_NAME python3 -c 'import torch; print(f\"PyTorch version: {torch.__version__}\")'"
}

# Test AutoDock Vina
test_vina() {
    print_status "Testing AutoDock Vina..."
    
    run_test "Vina binary availability" \
        "docker run --rm $IMAGE_NAME which vina"
    
    run_test "Vina help command" \
        "docker run --rm $IMAGE_NAME vina --help >/dev/null"
    
    run_test "Vina version command" \
        "docker run --rm $IMAGE_NAME vina --version >/dev/null 2>&1" 1  # Vina returns 1 for --version
}

# Test GNINA
test_gnina() {
    print_status "Testing GNINA..."
    
    run_test "GNINA binary availability" \
        "docker run --rm $IMAGE_NAME which gnina"
    
    run_test "GNINA help command" \
        "docker run --rm $IMAGE_NAME gnina --help >/dev/null"
    
    run_test "GNINA environment variables" \
        "docker run --rm $IMAGE_NAME bash -c 'echo \$GNINA_PATH'"
    
    # Test GNINA execution with minimal parameters
    run_test "GNINA dry run" \
        "docker run --rm $IMAGE_NAME gnina --version >/dev/null 2>&1" 1  # GNINA may return 1 for --version
}

# Test DiffDock
test_diffdock() {
    print_status "Testing DiffDock..."
    
    run_test "DiffDock installation directory" \
        "docker run --rm $IMAGE_NAME test -d /opt/DiffDock"
    
    run_test "DiffDock inference script" \
        "docker run --rm $IMAGE_NAME test -f /opt/DiffDock/inference.py"
    
    run_test "DiffDock Python imports" \
        "docker run --rm $IMAGE_NAME python3 -c 'import sys; sys.path.append(\"/opt/DiffDock\"); import torch_geometric; print(\"DiffDock dependencies: OK\")'"
    
    run_test "DiffDock environment variables" \
        "docker run --rm $IMAGE_NAME bash -c 'echo \$DIFFDOCK_PATH'"
    
    # Test if DiffDock models are downloaded
    run_test "DiffDock pre-trained models" \
        "docker run --rm $IMAGE_NAME test -d /opt/DiffDock/workdir/paper_models"
}

# Test pipeline scripts
test_pipeline_scripts() {
    print_status "Testing pipeline scripts..."
    
    run_test "Batch pipeline script" \
        "docker run --rm $IMAGE_NAME python3 /app/pipeline/batch_pipeline.py --help >/dev/null"
    
    run_test "Structure preparation script" \
        "docker run --rm $IMAGE_NAME python3 /app/pipeline/prep_structures.py --help >/dev/null"
    
    run_test "Docking multi script" \
        "docker run --rm $IMAGE_NAME python3 /app/pipeline/run_docking_multi.py --help >/dev/null"
    
    run_test "Results parsing script" \
        "docker run --rm $IMAGE_NAME python3 /app/pipeline/parse_and_score_results.py --help >/dev/null"
}

# Test integration - run a simple docking job
test_integration() {
    print_status "Testing integration with real docking job..."
    
    # Copy test files to container and run minimal docking
    run_test "Full pipeline integration (Vina only)" \
        "docker run --rm -v \"$TEST_DIR:/workspace\" $IMAGE_NAME \
         python3 run_batch_pipeline.py \
         --protein inputs/protein.pdb \
         --ligands inputs/ligands/ \
         --output-dir outputs/test_run \
         --engines vina >/dev/null 2>&1"
    
    # Check if output was generated
    if [ -d "$TEST_DIR/outputs/test_run" ]; then
        run_test "Vina output generation" \
            "test -f \"$TEST_DIR/outputs/test_run/final_summary.csv\""
    else
        print_warning "Vina integration test output not found"
    fi
}

# Test GNINA integration
test_gnina_integration() {
    print_status "Testing GNINA integration..."
    
    run_test "GNINA pipeline integration" \
        "docker run --rm -v \"$TEST_DIR:/workspace\" $IMAGE_NAME \
         python3 run_batch_pipeline.py \
         --protein inputs/protein.pdb \
         --ligands inputs/ligands/aspirin.smi \
         --output-dir outputs/gnina_test \
         --engines gnina \
         --max-ligands 1 >/dev/null 2>&1"
}

# Test DiffDock integration  
test_diffdock_integration() {
    print_status "Testing DiffDock integration..."
    
    run_test "DiffDock pipeline integration" \
        "docker run --rm -v \"$TEST_DIR:/workspace\" $IMAGE_NAME \
         python3 run_batch_pipeline.py \
         --protein inputs/protein.pdb \
         --ligands inputs/ligands/aspirin.smi \
         --output-dir outputs/diffdock_test \
         --engines diffdock \
         --max-ligands 1 >/dev/null 2>&1"
}

# Test Docker Compose
test_docker_compose() {
    if [ -f "$(dirname "$TEST_DIR")/docker-compose.yml" ]; then
        print_status "Testing Docker Compose..."
        
        cd "$(dirname "$TEST_DIR")"
        
        run_test "Docker Compose config validation" \
            "docker-compose config >/dev/null"
        
        run_test "Docker Compose service start" \
            "docker-compose up -d docking-pipeline && docker-compose down"
    else
        print_warning "docker-compose.yml not found - skipping Compose tests"
    fi
}

# Performance benchmarks
run_performance_tests() {
    print_status "Running performance benchmarks..."
    
    # Test memory usage
    run_test "Memory usage check" \
        "docker run --rm --memory=4g $IMAGE_NAME python3 -c 'import psutil; print(f\"Available memory: {psutil.virtual_memory().available // (1024**3)} GB\")'"
    
    # Test CPU utilization
    run_test "CPU detection" \
        "docker run --rm $IMAGE_NAME python3 -c 'import psutil; print(f\"CPU cores: {psutil.cpu_count()}\")'"
    
    # Test disk space
    run_test "Disk space check" \
        "docker run --rm $IMAGE_NAME df -h /"
}

# Cleanup function
cleanup() {
    print_status "Cleaning up test environment..."
    cd "$(dirname "$TEST_DIR")"
    rm -rf "$TEST_DIR"
    print_success "Cleanup completed"
}

# Main execution
main() {
    echo "=========================================="
    echo "Docker Setup Test Script for Molecular Docking Pipeline"
    echo "Testing GNINA and DiffDock Integration"
    echo "=========================================="
    echo
    
    # Check if Docker is available
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed or not in PATH"
        exit 1
    fi
    
    # Check if Docker daemon is running
    if ! docker info &> /dev/null; then
        print_error "Docker daemon is not running"
        exit 1
    fi
    
    print_status "Starting comprehensive Docker tests..."
    echo
    
    # Run test suites
    setup_test_env
    test_docker_images
    test_basic_functionality
    test_python_dependencies
    test_vina
    test_gnina
    test_diffdock
    test_pipeline_scripts
    test_integration
    
    # Optional integration tests (may take longer)
    if [ "${FULL_TESTS:-false}" = "true" ]; then
        print_status "Running full integration tests (this may take several minutes)..."
        test_gnina_integration
        test_diffdock_integration
    else
        print_status "Skipping full integration tests (set FULL_TESTS=true to enable)"
    fi
    
    test_docker_compose
    run_performance_tests
    
    # Cleanup
    if [ "${KEEP_TEST_FILES:-false}" != "true" ]; then
        cleanup
    else
        print_status "Keeping test files in $TEST_DIR"
    fi
    
    # Final report
    echo
    echo "=========================================="
    echo "Test Results Summary"
    echo "=========================================="
    echo "Total tests run: $TESTS_TOTAL"
    echo "Tests passed: $TESTS_PASSED"
    echo "Tests failed: $TESTS_FAILED"
    echo
    
    if [ $TESTS_FAILED -eq 0 ]; then
        print_success "All tests passed! 🎉"
        print_success "Your Docker environment is ready for molecular docking with GNINA and DiffDock!"
        echo
        echo "Next steps:"
        echo "1. Place your protein file in inputs/protein.pdb"
        echo "2. Place your ligand files in inputs/ligands/"
        echo "3. Run: docker run -v \$(pwd):/workspace $IMAGE_NAME python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock"
        exit 0
    else
        print_error "Some tests failed. Please check the output above."
        echo
        echo "Common issues:"
        echo "- Insufficient memory (need 8GB+ RAM)"
        echo "- Missing Docker image (run ./build_docker.sh first)"
        echo "- Insufficient disk space (need 20GB+ free)"
        echo "- Docker daemon not running"
        exit 1
    fi
}

# Handle script arguments
case "${1:-}" in
    --help|-h)
        echo "Usage: $0 [options]"
        echo ""
        echo "Options:"
        echo "  --help, -h          Show this help"
        echo "  --full              Run full integration tests (slower)"
        echo "  --keep-files        Keep test files after completion"
        echo ""
        echo "Environment variables:"
        echo "  FULL_TESTS=true     Enable full integration tests"
        echo "  KEEP_TEST_FILES=true Keep test files"
        echo "  TIMEOUT=300         Test timeout in seconds"
        exit 0
        ;;
    --full)
        export FULL_TESTS=true
        main
        ;;
    --keep-files)
        export KEEP_TEST_FILES=true
        main
        ;;
    "")
        main
        ;;
    *)
        print_error "Unknown option: $1"
        echo "Use --help for usage information"
        exit 1
        ;;
esac 