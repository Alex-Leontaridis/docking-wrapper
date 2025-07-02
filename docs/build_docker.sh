#!/bin/bash

# Molecular Docking Pipeline Docker Build Script
# This script builds the complete Docker environment for the batch docking pipeline

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
IMAGE_NAME="docking-pipeline"
IMAGE_TAG="latest"
FULL_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"

# Get build metadata
BUILD_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
VCS_REF=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

print_status "Starting Docker build for Molecular Docking Pipeline"
print_status "Image: ${FULL_IMAGE_NAME}"
print_status "Build date: ${BUILD_DATE}"
print_status "VCS ref: ${VCS_REF}"

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

# Check available disk space (require at least 20GB)
available_space=$(df . | tail -1 | awk '{print $4}')
required_space=$((20 * 1024 * 1024))  # 20GB in KB

if [ "$available_space" -lt "$required_space" ]; then
    print_warning "Available disk space is less than 20GB. Build might fail."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_error "Build cancelled by user."
        exit 1
    fi
fi

# Check if required files exist
required_files=("Dockerfile" "batch_pipeline.py" "prep_structures.py" "run_docking_multi.py" "parse_and_score_results.py" "pipeline_config.json")
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        print_error "Required file not found: $file"
        exit 1
    fi
done

print_success "Prerequisites check passed"

# Build Docker image
print_status "Building Docker image (this may take 20-30 minutes)..."

# Start timer
start_time=$(date +%s)

# Build with progress output
docker build \
    --progress=plain \
    --build-arg BUILD_DATE="$BUILD_DATE" \
    --build-arg VCS_REF="$VCS_REF" \
    -t "$FULL_IMAGE_NAME" \
    . 2>&1 | tee build.log

# Check build result
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    # Calculate build time
    end_time=$(date +%s)
    build_time=$((end_time - start_time))
    build_time_min=$((build_time / 60))
    build_time_sec=$((build_time % 60))
    
    print_success "Docker image built successfully in ${build_time_min}m ${build_time_sec}s"
else
    print_error "Docker build failed. Check build.log for details."
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
if docker run --rm "$FULL_IMAGE_NAME" python3 /app/pipeline/batch_pipeline.py --help &> /dev/null; then
    print_success "Pipeline help test passed"
else
    print_error "Pipeline help test failed"
    exit 1
fi

# Test Vina availability
if docker run --rm "$FULL_IMAGE_NAME" bash -c "source activate docking && vina --help" &> /dev/null; then
    print_success "Vina availability test passed"
else
    print_warning "Vina availability test failed (might still work)"
fi

# Show image information
print_status "Image information:"
docker images "$FULL_IMAGE_NAME" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

# Show usage examples
print_success "Docker image ready!"
echo
print_status "Usage examples:"
echo
echo "# Basic usage:"
echo "docker run -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/"
echo
echo "# With all engines:"
echo "docker run -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/ --enable-gnina --enable-diffdock"
echo
echo "# Using Docker Compose:"
echo "docker-compose up -d"
echo "docker-compose exec docking-pipeline python3 batch_pipeline.py --protein receptor.pdb --ligands ligands/"
echo
echo "# Interactive shell:"
echo "docker run -it -v \$(pwd):/workspace $FULL_IMAGE_NAME bash"
echo

# Cleanup build log if successful
if [ -f "build.log" ]; then
    read -p "Remove build.log? (Y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        rm build.log
        print_status "build.log removed"
    fi
fi

print_success "Build script completed successfully!"

# Optional: Tag for different versions
if [ -n "$VCS_REF" ] && [ "$VCS_REF" != "unknown" ]; then
    print_status "Tagging image with git commit: ${IMAGE_NAME}:${VCS_REF}"
    docker tag "$FULL_IMAGE_NAME" "${IMAGE_NAME}:${VCS_REF}"
fi

# Show next steps
print_status "Next steps:"
echo "1. Test with your data: docker run -v \$(pwd):/workspace $FULL_IMAGE_NAME python3 batch_pipeline.py --protein your_protein.pdb --ligands your_ligands/"
echo "2. See README_batch_pipeline.md for detailed usage instructions"
echo "3. Use docker-compose.yml for persistent services" 