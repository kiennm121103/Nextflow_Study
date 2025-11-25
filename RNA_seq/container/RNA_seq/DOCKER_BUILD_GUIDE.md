# Docker Build Automation - Quick Reference

## Quick Start

### Basic Build (Local)
```bash
./build_docker_manually.sh
```

### Build with Custom Name and Tag
```bash
./build_docker_manually.sh -n myrnaseq -t v1.0.0
```

### Build and Push to Docker Hub
```bash
./build_docker_manually.sh -r docker.io/yourusername --push
```

### Build Without Cache
```bash
./build_docker_manually.sh --no-cache
```

## Usage Examples

### 1. Local Development Build
```bash
./build_docker_manually.sh \
  -n rnaseq-dev \
  -t latest
```

### 2. Production Build with Registry
```bash
./build_docker_manually.sh \
  -n rnaseq-pipeline \
  -t 1.0.0 \
  -r docker.io/kiennm121103 \
  --push
```

### 3. Multi-platform Build
```bash
./build_docker_manually.sh \
  -n rnaseq-pipeline \
  -t latest \
  -p linux/amd64,linux/arm64
```

### 4. Using Environment Variables
```bash
export IMAGE_NAME=rnaseq-pipeline
export IMAGE_TAG=v2.0.0
export DOCKER_REGISTRY=docker.io/myorg
export PUSH_IMAGE=true

./build_docker_manually.sh
```

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-n, --name` | Docker image name | `rnaseq-pipeline` |
| `-t, --tag` | Docker image tag | `latest` |
| `-r, --registry` | Docker registry URL | `(none)` |
| `-f, --file` | Path to Dockerfile | `./Dockerfile` |
| `-c, --context` | Build context path | `.` |
| `-p, --platform` | Target platform | `linux/amd64` |
| `--no-cache` | Build without cache | `false` |
| `--push` | Push to registry after build | `false` |
| `-h, --help` | Show help message | - |

## Running the Container

### Interactive Shell
```bash
docker run -it --rm rnaseq-pipeline:latest
```

### With Mounted Volumes
```bash
docker run -it --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/references:/references \
  rnaseq-pipeline:latest
```

### Run Specific Command
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  rnaseq-pipeline:latest \
  fastqc /data/sample.fastq.gz -o /data/
```

### Run with Resource Limits
```bash
docker run -it --rm \
  --cpus="8" \
  --memory="32g" \
  -v $(pwd)/data:/data \
  rnaseq-pipeline:latest
```

## Installed Tools

The Docker image includes:

### Quality Control
- FastQC 0.12.1
- MultiQC 1.19
- RSeQC 5.0.3
- Qualimap 2.2.2d

### Trimming & Filtering
- Trim Galore 0.6.10
- Cutadapt 4.4
- fastp 0.23.4
- SortMeRNA 4.3.6

### Alignment
- STAR 2.7.10b
- HISAT2 2.2.1
- Salmon 1.10.2
- Kallisto 0.50.1

### Quantification
- featureCounts (Subread 2.0.6)
- HTSeq 2.0.5
- RSEM 1.3.3
- StringTie 2.2.1

### Utilities
- SAMtools 1.18
- BEDtools 2.31.0
- deepTools 3.5.4
- Picard 3.1.0

### R/Bioconductor
- DESeq2
- edgeR
- limma
- tximport

### Python Packages
- NumPy, Pandas
- Matplotlib, Seaborn
- Biopython, pysam

## Docker Registry Setup

### Login to Docker Hub
```bash
docker login
```

### Login to GitHub Container Registry
```bash
echo $GITHUB_TOKEN | docker login ghcr.io -u USERNAME --password-stdin
```

### Login to Custom Registry
```bash
docker login myregistry.com
```

## Troubleshooting

### Build Fails with "No Space Left"
```bash
# Clean up Docker
docker system prune -a

# Or specify cleanup in build
docker build --rm --force-rm ...
```

### Permission Denied
```bash
# Add user to docker group
sudo usermod -aG docker $USER
# Re-login to apply changes
```

### Can't Find Dockerfile
```bash
# Specify Dockerfile location
./build_docker_manually.sh -f /path/to/Dockerfile
```

### Push Failed - Not Authenticated
```bash
# Login first
docker login docker.io

# Then build and push
./build_docker_manually.sh -r docker.io/username --push
```

## Advanced Usage

### Build with Custom Dockerfile
```bash
./build_docker_manually.sh \
  -f custom.Dockerfile \
  -n custom-rnaseq \
  -t latest
```

### Build from Different Context
```bash
./build_docker_manually.sh \
  -c ../.. \
  -f dockerfiles/rnaseq.Dockerfile
```

### Skip Interactive Confirmation
```bash
yes | ./build_docker_manually.sh -n rnaseq -t latest
```

## CI/CD Integration

### GitHub Actions
```yaml
name: Build Docker Image

on:
  push:
    branches: [ main ]

jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Build image
        run: |
          cd RNA_seq/container/RNA_seq
          ./build_docker_manually.sh -n rnaseq -t ${{ github.sha }}
```

### GitLab CI
```yaml
build-docker:
  stage: build
  script:
    - cd RNA_seq/container/RNA_seq
    - ./build_docker_manually.sh -n rnaseq -t $CI_COMMIT_SHORT_SHA
```

## Maintenance

### Update Image
```bash
# Build new version
./build_docker_manually.sh -n rnaseq -t v1.1.0

# Tag as latest
docker tag rnaseq:v1.1.0 rnaseq:latest
```

### Clean Old Images
```bash
# Remove dangling images
docker image prune

# Remove all unused images
docker image prune -a
```

### Check Image Size
```bash
docker images rnaseq-pipeline:latest
```

### Inspect Image
```bash
docker inspect rnaseq-pipeline:latest
```

## Support

For issues or questions, please check:
- Docker documentation: https://docs.docker.com/
- Dockerfile best practices: https://docs.docker.com/develop/dev-best-practices/
