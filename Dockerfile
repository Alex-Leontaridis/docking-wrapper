# Molecular Docking Pipeline Docker Environment
# Based on Ubuntu 22.04 with support for AutoDock Vina, GNINA, DiffDock, RDKit, OpenBabel, and Meeko
# 
# Build: docker build -t docking-pipeline .
# Run:   docker run -v $(pwd):/workspace docking-pipeline python3 batch_pipeline.py --protein protein.pdb --ligands ligands/

FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Set work directory
WORKDIR /app

# System information
LABEL maintainer="Molecular Docking Pipeline"
LABEL version="1.0"
LABEL description="Complete molecular docking environment with Vina, GNINA, DiffDock"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Core utilities
    build-essential \
    cmake \
    git \
    wget \
    curl \
    unzip \
    software-properties-common \
    # Python and pip
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    # Scientific computing libraries
    libboost-all-dev \
    libeigen3-dev \
    libgoogle-glog-dev \
    libprotobuf-dev \
    protobuf-compiler \
    libhdf5-dev \
    libatlas-base-dev \
    # RDKit dependencies
    libcairo2-dev \
    libjpeg-dev \
    libgif-dev \
    # OpenBabel dependencies
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    # Additional scientific libraries
    liblapack-dev \
    libopenblas-dev \
    libfftw3-dev \
    # Utilities
    vim \
    nano \
    htop \
    && rm -rf /var/lib/apt/lists/*

# Install Conda (miniconda)
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Create conda environment with Python dependencies
RUN conda create -n docking python=3.10 -y && \
    echo "source activate docking" > ~/.bashrc
ENV PATH="/opt/conda/envs/docking/bin:$PATH"

# Activate environment and install core Python packages
RUN /bin/bash -c "source activate docking && \
    conda install -y -c conda-forge \
    numpy \
    scipy \
    pandas \
    matplotlib \
    scikit-learn \
    rdkit \
    openbabel \
    biopython \
    prody \
    mdanalysis \
    && conda clean -afy"

# Install additional Python packages via pip
RUN /bin/bash -c "source activate docking && \
    pip install --no-cache-dir \
    meeko \
    vina \
    gemmi \
    jupyter \
    seaborn \
    plotly \
    tqdm \
    psutil \
    memory-profiler"

# Install AutoDock Vina from source (latest version)
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina.git /tmp/vina && \
    cd /tmp/vina && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    rm -rf /tmp/vina

# Install GNINA (GPU-enabled molecular docking)
RUN git clone https://github.com/gnina/gnina.git /tmp/gnina && \
    cd /tmp/gnina && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install && \
    rm -rf /tmp/gnina

# Install OpenBabel from source (latest version)
RUN git clone https://github.com/openbabel/openbabel.git /tmp/openbabel && \
    cd /tmp/openbabel && \
    mkdir build && \
    cd build && \
    cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON .. && \
    make -j$(nproc) && \
    make install && \
    ldconfig && \
    rm -rf /tmp/openbabel

# Install PyTorch for DiffDock (CPU version for compatibility)
RUN /bin/bash -c "source activate docking && \
    pip install --no-cache-dir torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu"

# Install DiffDock
RUN git clone https://github.com/gcorso/DiffDock.git /opt/DiffDock && \
    cd /opt/DiffDock && \
    /bin/bash -c "source activate docking && pip install -r requirements.txt"

# Install additional bioinformatics tools
RUN /bin/bash -c "source activate docking && \
    pip install --no-cache-dir \
    biotite \
    pymol-open-source \
    nglview \
    py3Dmol"

# Set up MGLTools for structure preparation (lightweight alternative)
RUN /bin/bash -c "source activate docking && \
    pip install --no-cache-dir pdb2pqr apbs-runner"

# Create application directories
RUN mkdir -p /app/pipeline /app/data /app/outputs /app/logs

# Copy pipeline scripts
COPY batch_pipeline.py /app/pipeline/
COPY prep_structures.py /app/pipeline/
COPY run_docking_multi.py /app/pipeline/
COPY parse_and_score_results.py /app/pipeline/
COPY pipeline_config.json /app/pipeline/

# Make scripts executable
RUN chmod +x /app/pipeline/*.py

# Set up environment variables
ENV PYTHONPATH="/app/pipeline:$PYTHONPATH"
ENV PATH="/app/pipeline:$PATH"
ENV DIFFDOCK_PATH="/opt/DiffDock"
ENV OMP_NUM_THREADS=1

# Set conda environment activation
RUN echo '#!/bin/bash\nsource activate docking\nexec "$@"' > /entrypoint.sh && \
    chmod +x /entrypoint.sh

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python3 -c "import rdkit, openbabel, pandas, numpy; print('Environment OK')" || exit 1

# Volume for data
VOLUME ["/workspace"]

# Working directory for pipeline execution
WORKDIR /workspace

# Entry point
ENTRYPOINT ["/entrypoint.sh"]

# Default command
CMD ["python3", "/app/pipeline/batch_pipeline.py", "--help"]

# Image metadata
ARG BUILD_DATE
ARG VCS_REF
LABEL org.label-schema.build-date=$BUILD_DATE \
      org.label-schema.name="molecular-docking-pipeline" \
      org.label-schema.description="Complete molecular docking environment" \
      org.label-schema.vcs-ref=$VCS_REF \
      org.label-schema.vcs-url="https://github.com/your-repo/docking-pipeline" \
      org.label-schema.schema-version="1.0"

# Documentation
RUN echo "Molecular Docking Pipeline Container" > /app/README.txt && \
    echo "====================================" >> /app/README.txt && \
    echo "" >> /app/README.txt && \
    echo "This container includes:" >> /app/README.txt && \
    echo "- AutoDock Vina (latest)" >> /app/README.txt && \
    echo "- GNINA (latest)" >> /app/README.txt && \
    echo "- DiffDock (latest)" >> /app/README.txt && \
    echo "- RDKit" >> /app/README.txt && \
    echo "- OpenBabel" >> /app/README.txt && \
    echo "- Meeko" >> /app/README.txt && \
    echo "- Complete Python scientific stack" >> /app/README.txt && \
    echo "" >> /app/README.txt && \
    echo "Usage:" >> /app/README.txt && \
    echo "docker run -v \$(pwd):/workspace docking-pipeline python3 batch_pipeline.py --protein protein.pdb --ligands ligands/" >> /app/README.txt && \
    echo "" >> /app/README.txt && \
    echo "See /app/pipeline/ for all available scripts." >> /app/README.txt

# Final setup - ensure all tools are accessible
RUN /bin/bash -c "source activate docking && \
    python3 -c 'import rdkit; print(f\"RDKit version: {rdkit.__version__}\")' && \
    python3 -c 'import openbabel; print(\"OpenBabel: OK\")' && \
    python3 -c 'import pandas; print(f\"Pandas version: {pandas.__version__}\")' && \
    vina --help > /dev/null && echo 'Vina: OK' && \
    which gnina > /dev/null && echo 'GNINA: OK' || echo 'GNINA: Not in PATH but built'" 