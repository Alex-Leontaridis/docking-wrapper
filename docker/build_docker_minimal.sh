#!/bin/bash

# Minimal Molecular Docking Pipeline Docker Build Script
# This script builds a lightweight Docker environment for the batch docking pipeline

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Configuration
IMAGE_NAME="docking-pipeline-minimal"
IMAGE_TAG="latest"
FULL_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"

print_status "Starting minimal Docker build for Molecular Docking Pipeline"
print_status "Image: ${FULL_IMAGE_NAME}"

# Check prerequisites
print_status "Checking prerequisites..."

# Check if Docker is installed and running
if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed. Please install Docker first."
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    print_error "Docker daemon is not running. Please start Docker."
    exit 1
fi

# Check if required files exist
required_files=("Dockerfile.minimal" "scripts/batch_pipeline.py" "scripts/prep_structures.py" "scripts/run_docking_multi.py" "scripts/parse_and_score_results.py" "pipeline_config.json")
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        print_error "Required file not found: $file"
        exit 1
    fi
done

print_success "Prerequisites check passed"

# Build Docker image
print_status "Building minimal Docker image (this should take 5-10 minutes)..."

# Start timer
start_time=$(date +%s)

# Build with progress output
docker build \
    -f Dockerfile.minimal \
    -t "$FULL_IMAGE_NAME" \
    . 2>&1 | tee build_minimal.log

# Check build result
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    # Calculate build time
    end_time=$(date +%s)
    build_time=$((end_time - start_time))
    build_time_min=$((build_time / 60))
    build_time_sec=$((build_time % 60))
    
    print_success "Docker image built successfully in ${build_time_min}m ${build_time_sec}s"
else
    print_error "Docker build failed. Check build_minimal.log for details."
    exit 1
fi

# Test the built image
print_status "Testing the built image..."

# Test basic functionality
if docker run --rm "$FULL_IMAGE_NAME" python3 -c "import rdkit, pandas, numpy; print('Basic imports: OK')" &> /dev/null; then
    print_success "Basic Python imports test passed"
else
    print_error "Basic Python imports test failed"
    exit 1
fi

# Test pipeline help
if docker run --rm -v $(pwd):/workspace "$FULL_IMAGE_NAME" python3 /workspace/scripts/batch_pipeline.py --help &> /dev/null; then
    print_success "Pipeline help test passed"
else
    print_error "Pipeline help test failed"
    exit 1
fi

# Show image information
print_status "Image information:"
docker images "$FULL_IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

# Show usage examples
print_success "Minimal Docker image ready!"
echo
print_status "Usage examples:"
echo
echo "# Basic usage:"
echo "docker run --rm -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/"
echo
echo "# Interactive shell:"
echo "docker run --rm -it -v \$(pwd):/workspace $FULL_IMAGE_NAME bash"
echo
echo "# Test with sample data:"
echo "docker run --rm -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --serial"
echo

# Cleanup build log if successful
if [ -f "build_minimal.log" ]; then
    read -p "Remove build_minimal.log? (Y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        rm build_minimal.log
        print_status "build_minimal.log removed"
    fi
fi

print_success "Minimal build script completed successfully!"

# Show next steps
print_status "Next steps:"
echo "1. Test with your data: docker run --rm -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 /workspace/scripts/batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/"
echo "2. See README_batch_pipeline.md for detailed usage instructions"
echo "3. This minimal image includes: Python 3.10, RDKit, AutoDock Vina, essential scientific libraries" 