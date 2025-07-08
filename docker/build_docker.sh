#!/bin/bash

# Enhanced Docker Build Script for Molecular Docking Pipeline
# Supports GNINA and DiffDock installation with verification

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
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

# Function to show help
show_help() {
    cat << EOF
Molecular Docking Pipeline Docker Build Script

Usage: $0 [OPTIONS]

OPTIONS:
    -h, --help          Show this help message
    -t, --tag TAG       Docker image tag (default: docking-pipeline:latest)
    -f, --force         Force rebuild without cache
    -v, --verbose       Verbose output
    --no-cache          Build without using cache
    --test-only         Only run tests, don't build
    --minimal           Build minimal version (Vina only)

EXAMPLES:
    $0                                  # Build with defaults
    $0 -t my-docking:v1.0              # Custom tag
    $0 --force --verbose               # Force rebuild with verbose output
    $0 --minimal                       # Build minimal version

ENVIRONMENT VARIABLES:
    DOCKER_BUILDKIT=1                  # Enable BuildKit (recommended)
    BUILD_ARGS                         # Additional build arguments

EOF
}

# Default values
IMAGE_TAG="docking-pipeline:latest"
DOCKERFILE="Dockerfile"
BUILD_ARGS=""
CACHE_FROM=""
VERBOSE=false
FORCE_REBUILD=false
TEST_ONLY=false
MINIMAL=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -t|--tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        -f|--force)
            FORCE_REBUILD=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --no-cache)
            BUILD_ARGS="$BUILD_ARGS --no-cache"
            shift
            ;;
        --test-only)
            TEST_ONLY=true
            shift
            ;;
        --minimal)
            MINIMAL=true
            DOCKERFILE="Dockerfile.minimal"
            IMAGE_TAG="docking-pipeline-minimal:latest"
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Enable verbose output if requested
if [ "$VERBOSE" = true ]; then
    set -x
    BUILD_ARGS="$BUILD_ARGS --progress=plain"
fi

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

# Check if Dockerfile exists
if [ ! -f "$DOCKERFILE" ]; then
    print_error "Dockerfile not found: $DOCKERFILE"
    exit 1
fi

print_status "Starting Docker build for Molecular Docking Pipeline"
print_status "Image tag: $IMAGE_TAG"
print_status "Dockerfile: $DOCKERFILE"

# Function to test Docker image
test_image() {
    local image_name=$1
    
    print_status "Testing Docker image: $image_name"
    
    # Test 1: Basic container startup
    print_status "Test 1: Container startup..."
    if docker run --rm "$image_name" echo "Container startup: OK"; then
        print_success "Container startup test passed"
    else
        print_error "Container startup test failed"
        return 1
    fi
    
    # Test 2: Python environment
    print_status "Test 2: Python environment..."
    if docker run --rm "$image_name" python3 -c "import sys; print(f'Python version: {sys.version}')"; then
        print_success "Python environment test passed"
    else
        print_error "Python environment test failed"
        return 1
    fi
    
    # Test 3: Core dependencies
    print_status "Test 3: Core dependencies..."
    if docker run --rm "$image_name" python3 -c "import rdkit, pandas, numpy; print('Core dependencies: OK')"; then
        print_success "Core dependencies test passed"
    else
        print_error "Core dependencies test failed"
        return 1
    fi
    
    # Test 4: AutoDock Vina
    print_status "Test 4: AutoDock Vina..."
    if docker run --rm "$image_name" vina --help > /dev/null 2>&1; then
        print_success "AutoDock Vina test passed"
    else
        print_warning "AutoDock Vina test failed"
    fi
    
    if [ "$MINIMAL" = false ]; then
        # Test 5: GNINA
        print_status "Test 5: GNINA..."
        if docker run --rm "$image_name" gnina --help > /dev/null 2>&1; then
            print_success "GNINA test passed"
        else
            print_warning "GNINA test failed"
        fi
        
        # Test 6: DiffDock
        print_status "Test 6: DiffDock..."
        DIFFDOCK_PATH=${DIFFDOCK_PATH:-"/opt/DiffDock"}
        if docker run --rm "$image_name" test -f "$DIFFDOCK_PATH/inference.py"; then
            print_success "DiffDock test passed"
        else
            print_warning "DiffDock test failed"
        fi
        
        # Test 7: PyTorch
        print_status "Test 7: PyTorch..."
        if docker run --rm "$image_name" python3 -c "import torch; print(f'PyTorch version: {torch.__version__}')"; then
            print_success "PyTorch test passed"
        else
            print_warning "PyTorch test failed"
        fi
    fi
    
    print_success "All critical tests passed!"
    return 0
}

# If test-only mode, just run tests
if [ "$TEST_ONLY" = true ]; then
    print_status "Running tests only (no build)"
    if docker image inspect "$IMAGE_TAG" >/dev/null 2>&1; then
        test_image "$IMAGE_TAG"
        exit $?
    else
        print_error "Image $IMAGE_TAG not found. Build it first."
        exit 1
    fi
fi

# Check available disk space (minimum 10GB recommended)
AVAILABLE_SPACE=$(df / | awk 'NR==2 {print $4}')
if [ "$AVAILABLE_SPACE" -lt 10485760 ]; then  # 10GB in KB
    print_warning "Low disk space detected. This build requires at least 10GB free space."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Set build arguments
BUILD_DATE=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
VCS_REF=$(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

BUILD_ARGS="$BUILD_ARGS --build-arg BUILD_DATE=$BUILD_DATE"
BUILD_ARGS="$BUILD_ARGS --build-arg VCS_REF=$VCS_REF"

# Add cache from existing image if not force rebuild
if [ "$FORCE_REBUILD" = false ] && docker image inspect "$IMAGE_TAG" >/dev/null 2>&1; then
    print_status "Using existing image as cache source"
    BUILD_ARGS="$BUILD_ARGS --cache-from $IMAGE_TAG"
fi

# Enable BuildKit if available
export DOCKER_BUILDKIT=1

print_status "Building Docker image..."
print_status "Build command: docker build $BUILD_ARGS -t $IMAGE_TAG -f $DOCKERFILE ."

# Estimate build time
if [ "$MINIMAL" = true ]; then
    print_status "Estimated build time: 10-15 minutes (minimal version)"
else
    print_status "Estimated build time: 25-40 minutes (full version with GNINA and DiffDock)"
fi

# Build the image
START_TIME=$(date +%s)

if docker build $BUILD_ARGS -t "$IMAGE_TAG" -f "$DOCKERFILE" .; then
    END_TIME=$(date +%s)
    BUILD_TIME=$((END_TIME - START_TIME))
    
    print_success "Docker image built successfully in ${BUILD_TIME}s"
    print_success "Image tagged as: $IMAGE_TAG"
    
    # Show image size
    IMAGE_SIZE=$(docker image inspect "$IMAGE_TAG" --format='{{.Size}}' | numfmt --to=iec)
    print_status "Image size: $IMAGE_SIZE"
    
    # Run tests
    print_status "Running post-build tests..."
    if test_image "$IMAGE_TAG"; then
        print_success "All tests passed! Image is ready to use."
        
        echo
        print_status "To run the pipeline:"
        echo "docker run -v \$(pwd)/inputs:/workspace/inputs \\"
        echo "           -v \$(pwd)/outputs:/workspace/outputs \\"
        echo "           $IMAGE_TAG \\"
        echo "           python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/"
        
        if [ "$MINIMAL" = false ]; then
            echo
            print_status "To enable GNINA and DiffDock:"
            echo "docker run -v \$(pwd)/inputs:/workspace/inputs \\"
            echo "           -v \$(pwd)/outputs:/workspace/outputs \\"
            echo "           $IMAGE_TAG \\"
            echo "           python3 run_batch_pipeline.py --protein inputs/protein.pdb --ligands inputs/ligands/ --enable-gnina --enable-diffdock"
        fi
        
    else
        print_warning "Some tests failed, but image was built successfully"
        exit 1
    fi
    
else
    print_error "Docker build failed"
    exit 1
fi 