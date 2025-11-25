#!/usr/bin/env bash

################################################################################
# RNA-seq Docker Image Build Automation Script
# 
# This script automates the building of a Docker image for RNA-seq analysis
# with all necessary bioinformatics tools
#
# Author: kiennm121103
# Date: November 24, 2025
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
IMAGE_NAME="${IMAGE_NAME:-rnaseq-pipeline}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
DOCKER_REGISTRY="${DOCKER_REGISTRY:-}"
DOCKERFILE_PATH="${DOCKERFILE_PATH:-./Dockerfile}"
BUILD_CONTEXT="${BUILD_CONTEXT:-.}"
NO_CACHE="${NO_CACHE:-false}"
PUSH_IMAGE="${PUSH_IMAGE:-false}"
PLATFORM="${PLATFORM:-linux/amd64}"

# Function to print colored messages
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to print usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Build Docker image for RNA-seq analysis pipeline

OPTIONS:
    -n, --name NAME           Docker image name (default: rnaseq-pipeline)
    -t, --tag TAG             Docker image tag (default: latest)
    -r, --registry REGISTRY   Docker registry to push to (e.g., docker.io/username)
    -f, --file DOCKERFILE     Path to Dockerfile (default: ./Dockerfile)
    -c, --context PATH        Build context path (default: .)
    -p, --platform PLATFORM   Target platform (default: linux/amd64)
    --no-cache                Build without using cache
    --push                    Push image to registry after build
    -h, --help                Show this help message

EXAMPLES:
    # Basic build
    $0

    # Build with custom name and tag
    $0 -n myrnaseq -t v1.0.0

    # Build and push to Docker Hub
    $0 -r docker.io/username --push

    # Build for multiple platforms
    $0 -p linux/amd64,linux/arm64

ENVIRONMENT VARIABLES:
    IMAGE_NAME              Docker image name
    IMAGE_TAG               Docker image tag
    DOCKER_REGISTRY         Docker registry URL
    NO_CACHE                Set to 'true' to disable cache
    PUSH_IMAGE              Set to 'true' to push after build

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--name)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -t|--tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        -r|--registry)
            DOCKER_REGISTRY="$2"
            shift 2
            ;;
        -f|--file)
            DOCKERFILE_PATH="$2"
            shift 2
            ;;
        -c|--context)
            BUILD_CONTEXT="$2"
            shift 2
            ;;
        -p|--platform)
            PLATFORM="$2"
            shift 2
            ;;
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        --push)
            PUSH_IMAGE=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            ;;
    esac
done

# Function to check if Docker is installed
check_docker() {
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed. Please install Docker first."
        exit 1
    fi
    log_success "Docker is installed: $(docker --version)"
}

# Function to check if Dockerfile exists
check_dockerfile() {
    if [[ ! -f "$DOCKERFILE_PATH" ]]; then
        log_error "Dockerfile not found at: $DOCKERFILE_PATH"
        log_info "Creating a default Dockerfile..."
        create_default_dockerfile
    else
        log_success "Dockerfile found at: $DOCKERFILE_PATH"
    fi
}

# Function to create default Dockerfile
create_default_dockerfile() {
    cat > "$DOCKERFILE_PATH" << 'DOCKERFILEEOF'
FROM ubuntu:22.04

LABEL maintainer="kiennm121103"
LABEL description="Docker image for RNA-seq analysis pipeline"
LABEL version="1.0.0"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Set working directory
WORKDIR /workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    autoconf \
    automake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    libncursesw5-dev \
    python3 \
    python3-pip \
    python3-dev \
    default-jre \
    unzip \
    gzip \
    bzip2 \
    pigz \
    libxml2-dev \
    libgomp1 \
    perl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda/mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/conda && \
    rm Mambaforge-Linux-x86_64.sh && \
    /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --set channel_priority strict

ENV PATH="/opt/conda/bin:${PATH}"

# Install bioinformatics tools via conda
RUN mamba install -y -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    multiqc=1.19 \
    trim-galore=0.6.10 \
    cutadapt=4.4 \
    star=2.7.10b \
    samtools=1.18 \
    hisat2=2.2.1 \
    salmon=1.10.2 \
    kallisto=0.50.1 \
    rsem=1.3.3 \
    subread=2.0.6 \
    htseq=2.0.5 \
    bedtools=2.31.0 \
    deeptools=3.5.4 \
    picard=3.1.0 \
    rseqc=5.0.3 \
    preseq=3.2.0 \
    qualimap=2.2.2d \
    dupradar=1.28.0 \
    sortmerna=4.3.6 \
    bbmap=39.01 \
    fastp=0.23.4 \
    stringtie=2.2.1 \
    && mamba clean -a -y

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    scikit-learn \
    biopython \
    pysam \
    pyBigWig

# Install R and R packages
RUN mamba install -y -c conda-forge -c bioconda \
    r-base=4.3.1 \
    r-essentials \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-tximport \
    bioconductor-genomeinfodb \
    bioconductor-genomicfeatures \
    bioconductor-genomicranges \
    r-ggplot2 \
    r-dplyr \
    r-tidyr \
    r-pheatmap \
    r-rcolorbrewer \
    && mamba clean -a -y

# Create directories for data
RUN mkdir -p /data /results /references

# Set permissions
RUN chmod -R 755 /workspace /data /results /references

# Add version information
RUN echo "RNA-seq Pipeline Docker Image v1.0.0" > /VERSION && \
    echo "Build date: $(date)" >> /VERSION && \
    echo "Tools installed:" >> /VERSION && \
    echo "  - FastQC: $(fastqc --version 2>&1)" >> /VERSION && \
    echo "  - STAR: $(STAR --version 2>&1)" >> /VERSION && \
    echo "  - Samtools: $(samtools --version | head -n1)" >> /VERSION && \
    echo "  - Salmon: $(salmon --version 2>&1)" >> /VERSION

# Default command
CMD ["/bin/bash"]
DOCKERFILEEOF
    log_success "Default Dockerfile created at: $DOCKERFILE_PATH"
}

# Function to build full image name
get_full_image_name() {
    if [[ -n "$DOCKER_REGISTRY" ]]; then
        echo "${DOCKER_REGISTRY}/${IMAGE_NAME}:${IMAGE_TAG}"
    else
        echo "${IMAGE_NAME}:${IMAGE_TAG}"
    fi
}

# Function to build Docker image
build_image() {
    local full_image_name=$(get_full_image_name)
    local build_args="--file $DOCKERFILE_PATH --tag $full_image_name"
    
    if [[ "$NO_CACHE" == "true" ]]; then
        build_args="$build_args --no-cache"
        log_info "Building without cache..."
    fi
    
    if [[ -n "$PLATFORM" ]]; then
        build_args="$build_args --platform $PLATFORM"
        log_info "Building for platform: $PLATFORM"
    fi
    
    log_info "Building Docker image: $full_image_name"
    log_info "Build context: $BUILD_CONTEXT"
    log_info "Dockerfile: $DOCKERFILE_PATH"
    
    if docker build $build_args "$BUILD_CONTEXT"; then
        log_success "Docker image built successfully: $full_image_name"
        return 0
    else
        log_error "Failed to build Docker image"
        return 1
    fi
}

# Function to test the built image
test_image() {
    local full_image_name=$(get_full_image_name)
    log_info "Testing Docker image: $full_image_name"
    
    # Test if image exists
    if ! docker image inspect "$full_image_name" &> /dev/null; then
        log_error "Image not found: $full_image_name"
        return 1
    fi
    
    # Test basic tools
    log_info "Testing installed tools..."
    
    local test_commands=(
        "fastqc --version"
        "STAR --version"
        "samtools --version"
        "salmon --version"
        "python3 --version"
        "Rscript --version"
    )
    
    for cmd in "${test_commands[@]}"; do
        if docker run --rm "$full_image_name" bash -c "$cmd" &> /dev/null; then
            log_success "   $cmd"
        else
            log_warning "   $cmd failed"
        fi
    done
    
    log_success "Image testing completed"
}

# Function to push image to registry
push_image() {
    local full_image_name=$(get_full_image_name)
    
    if [[ -z "$DOCKER_REGISTRY" ]]; then
        log_warning "No registry specified. Skipping push."
        return 0
    fi
    
    log_info "Pushing image to registry: $full_image_name"
    
    if docker push "$full_image_name"; then
        log_success "Image pushed successfully: $full_image_name"
        return 0
    else
        log_error "Failed to push image"
        return 1
    fi
}

# Function to display image information
display_info() {
    local full_image_name=$(get_full_image_name)
    
    log_info "Image build complete!"
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Image Name:    $full_image_name"
    echo "  Size:          $(docker image inspect "$full_image_name" --format='{{.Size}}' | numfmt --to=iec-i --suffix=B 2>/dev/null || echo 'N/A')"
    echo "  Created:       $(docker image inspect "$full_image_name" --format='{{.Created}}')"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    echo "To run the container:"
    echo "  docker run -it --rm $full_image_name"
    echo ""
    echo "To run with mounted volumes:"
    echo "  docker run -it --rm \\"
    echo "    -v \$(pwd)/data:/data \\"
    echo "    -v \$(pwd)/results:/results \\"
    echo "    $full_image_name"
    echo ""
    if [[ "$PUSH_IMAGE" == "true" ]] && [[ -n "$DOCKER_REGISTRY" ]]; then
        echo "Image pushed to: $DOCKER_REGISTRY"
        echo ""
    fi
}

# Function to save build log
save_build_log() {
    local log_file="docker_build_$(date +%Y%m%d_%H%M%S).log"
    log_info "Build log saved to: $log_file"
}

# Main execution
main() {
    log_info "Starting RNA-seq Docker image build automation..."
    echo ""
    
    # Pre-build checks
    check_docker
    check_dockerfile
    
    echo ""
    log_info "Build Configuration:"
    echo "  Image Name:     $IMAGE_NAME"
    echo "  Image Tag:      $IMAGE_TAG"
    echo "  Registry:       ${DOCKER_REGISTRY:-'(none - local only)'}"
    echo "  Dockerfile:     $DOCKERFILE_PATH"
    echo "  Build Context:  $BUILD_CONTEXT"
    echo "  Platform:       $PLATFORM"
    echo "  No Cache:       $NO_CACHE"
    echo "  Push to Registry: $PUSH_IMAGE"
    echo ""
    
    # Confirm build
    read -p "Continue with build? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_warning "Build cancelled by user"
        exit 0
    fi
    
    # Build the image
    if ! build_image; then
        log_error "Build failed"
        exit 1
    fi
    
    # Test the image
    test_image
    
    # Push if requested
    if [[ "$PUSH_IMAGE" == "true" ]]; then
        push_image
    fi
    
    # Display final information
    display_info
    
    log_success "All done!"
}

# Run main function
main "$@"
