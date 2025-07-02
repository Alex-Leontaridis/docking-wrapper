#!/bin/bash

set -e

# Helper function for colored output
echo_info() { echo -e "\033[1;34m[INFO]\033[0m $1"; }
echo_warn() { echo -e "\033[1;33m[WARN]\033[0m $1"; }
echo_err()  { echo -e "\033[1;31m[ERR ]\033[0m $1"; }

# 1. Miniconda installation check
if ! command -v conda &> /dev/null; then
    echo_info "Miniconda not found. Installing via Homebrew..."
    if ! command -v brew &> /dev/null; then
        echo_err "Homebrew is required but not found. Please install Homebrew first."
        exit 1
    fi
    brew install --cask miniconda
    eval "$('/opt/homebrew/Caskroom/miniconda/base/bin/conda' 'shell.bash' 'hook')"
else
    echo_info "Miniconda/conda already installed."
fi

# 2. Conda initialization for bash
if ! grep -q 'conda initialize' ~/.bash_profile 2>/dev/null; then
    echo_info "Initializing conda for bash..."
    conda init bash
    source ~/.bash_profile
else
    echo_info "Conda already initialized in bash profile."
    source ~/.bash_profile
fi

# 3. Create and activate docking environment
if conda info --envs | grep -q '^docking'; then
    echo_info "Conda environment 'docking' already exists."
else
    echo_info "Creating conda environment 'docking' with Python 3.10..."
    conda create -y -n docking python=3.10
fi

echo_info "Activating 'docking' environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate docking

# 4. Install dependencies for GNINA
conda install -y -c conda-forge cmake boost git

# 5. GNINA installation/build
if [ -d "gnina" ] && [ -x "gnina/build/gnina" ]; then
    echo_info "GNINA already built in ./gnina/build/gnina."
else
    if [ ! -d "gnina" ]; then
        echo_info "Cloning GNINA repository..."
        git clone https://github.com/gnina/gnina.git
    fi
    cd gnina
    mkdir -p build
    cd build
    if [ ! -x "gnina" ]; then
        echo_info "Building GNINA..."
        cmake ..
        make -j$(sysctl -n hw.ncpu)
    fi
    cd ../..
fi

# 6. DiffDock installation
if [ -d "DiffDock" ] && [ -x "DiffDock/scripts/run_diffdock.py" ]; then
    echo_info "DiffDock already present."
else
    if [ ! -d "DiffDock" ]; then
        echo_info "Cloning DiffDock repository..."
        git clone https://github.com/gnina/DiffDock.git
    fi
    cd DiffDock
    echo_info "Installing DiffDock Python dependencies..."
    pip install -r requirements.txt
    cd ..
fi

# 7. AutoDock Vina installation
if command -v vina &> /dev/null; then
    echo_info "AutoDock Vina already installed."
else
    echo_info "Attempting to install AutoDock Vina from source..."
    if [ ! -d "vina" ]; then
        git clone https://github.com/ccsb-scripps/AutoDock-Vina.git vina
    fi
    cd vina
    mkdir -p build
    cd build
    cmake ..
    make -j$(sysctl -n hw.ncpu)
    if [ -f "vina" ]; then
        echo_info "AutoDock Vina built successfully."
        sudo cp vina /usr/local/bin/
    else
        echo_err "AutoDock Vina build failed. Please check the output above."
    fi
    cd ../..
fi

echo_info "Docking environment setup complete!" 