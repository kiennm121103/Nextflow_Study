# UV Package Manager Guide

## What is UV?

**UV** is an extremely fast Python package installer and resolver written in Rust. It's 10-100x faster than `pip` and provides:
- ⚡ Lightning-fast package installation
- 🔒 Deterministic dependency resolution with `uv.lock`
- 🎯 Compatible with `pip` and `pyproject.toml`
- 🔄 Drop-in replacement for `pip`, `pip-tools`, and `virtualenv`

## Installation

### Install UV

```bash
# Install UV
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or using pip
pip install uv

# Verify installation
uv --version
```

## Quick Start

### 1. Generate uv.lock (Single-Cell Pipeline)

```bash
cd /home/kiennm/Nextflow_Study/scRNAseq

# Initialize and create lock file
uv sync

# This will:
# - Create .venv/ directory
# - Generate uv.lock with exact package versions
# - Install all dependencies
```

### 2. Generate uv.lock (Bulk RNA-seq Pipeline)

```bash
cd /home/kiennm/Nextflow_Study

# Initialize and create lock file
uv sync

# This will:
# - Create .venv/ directory
# - Generate uv.lock with exact package versions
# - Install all dependencies
```

## UV Commands Cheat Sheet

### Basic Operations

```bash
# Create virtual environment
uv venv

# Sync environment from pyproject.toml (creates uv.lock)
uv sync

# Install from existing uv.lock (frozen - exact versions)
uv sync --frozen

# Activate virtual environment
source .venv/bin/activate

# Deactivate
deactivate
```

### Package Management

```bash
# Add a package (updates pyproject.toml and uv.lock)
uv add numpy

# Add a development dependency
uv add --dev pytest

# Add an optional dependency group
uv add --optional batch harmonypy

# Remove a package
uv remove numpy

# Install a package without adding to pyproject.toml
uv pip install package-name
```

### Updating Dependencies

```bash
# Update all packages to latest compatible versions
uv lock --upgrade

# Update a specific package
uv lock --upgrade-package numpy

# Update and sync
uv lock --upgrade && uv sync
```

### Installing Optional Features

```bash
# Install with development tools
uv sync --extra dev

# Install with batch correction tools
uv sync --extra batch

# Install multiple extras
uv sync --extra dev --extra batch --extra trajectory

# Install everything
uv sync --extra full
```

### Working with Lock Files

```bash
# Generate/update lock file without installing
uv lock

# Show what would be installed
uv pip compile pyproject.toml

# Export to requirements.txt format
uv pip compile pyproject.toml -o requirements.txt

# Verify lock file is up to date
uv lock --check
```

## Project Structure

### Single-Cell RNA-seq Pipeline

```
scRNAseq/
├── pyproject.toml          # Project configuration
├── uv.lock                 # Generated lock file (commit to git!)
├── .venv/                  # Virtual environment (don't commit)
├── requirements.txt        # Legacy format (optional)
└── scRNAseq.nf            # Pipeline script
```

### Bulk RNA-seq Pipeline

```
Nextflow_Study/
├── pyproject.toml          # Project configuration
├── uv.lock                 # Generated lock file (commit to git!)
├── .venv/                  # Virtual environment (don't commit)
├── requirements.txt        # Legacy format (optional)
└── RNAseq.nf              # Pipeline script
```

## Workflow Examples

### Setting Up a New Environment

```bash
# 1. Navigate to project directory
cd /home/kiennm/Nextflow_Study/scRNAseq

# 2. Create lock file and install dependencies
uv sync

# 3. Activate environment
source .venv/bin/activate

# 4. Verify installation
python -c "import scanpy; print(scanpy.__version__)"
python -c "import kb_python; print('kb-python OK')"

# 5. Run pipeline
nextflow run scRNAseq.nf
```

### Reproducing Environment on Another Machine

```bash
# 1. Clone repository with pyproject.toml and uv.lock
git clone <repository>
cd scRNAseq

# 2. Install UV if needed
curl -LsSf https://astral.sh/uv/install.sh | sh

# 3. Install exact versions from lock file
uv sync --frozen

# 4. Activate and use
source .venv/bin/activate
```

### Adding New Dependencies

```bash
# 1. Add package
uv add gseapy

# This automatically:
# - Adds to pyproject.toml
# - Updates uv.lock
# - Installs the package

# 2. Commit changes
git add pyproject.toml uv.lock
git commit -m "Add gseapy for enrichment analysis"
```

### Updating Packages

```bash
# Update all packages
uv lock --upgrade
uv sync

# Update specific package
uv lock --upgrade-package scanpy
uv sync

# Check what changed
git diff uv.lock
```

## Advantages over Conda/Pip

| Feature | UV | Conda | Pip |
|---------|----|----|-----|
| **Speed** | ⚡⚡⚡⚡⚡ | ⚡⚡ | ⚡⚡⚡ |
| **Lock File** | ✅ Built-in | ❌ Requires conda-lock | ❌ Requires pip-tools |
| **Reproducibility** | ✅ Perfect | ⚠️ Good | ❌ Poor without lock |
| **Disk Space** | ✅ Efficient | ❌ Large | ✅ Efficient |
| **Cross-platform** | ✅ Yes | ✅ Yes | ✅ Yes |
| **Binary packages** | ✅ PyPI wheels | ✅ Conda packages | ✅ PyPI wheels |

## Troubleshooting

### UV not found

```bash
# Add to PATH
export PATH="$HOME/.cargo/bin:$PATH"

# Or add to ~/.bashrc
echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Lock file conflicts

```bash
# Regenerate lock file
rm uv.lock
uv lock

# Force update
uv lock --upgrade
```

### Python version mismatch

```bash
# Specify Python version
uv venv --python 3.9

# Or
uv venv --python python3.10
```

### Package conflicts

```bash
# Show resolution details
uv lock --verbose

# Try with different resolver
uv lock --upgrade-package problematic-package
```

## Best Practices

### 1. Always Commit uv.lock

```bash
# .gitignore
.venv/
__pycache__/
*.pyc

# But NOT uv.lock - commit it!
```

### 2. Use --frozen in CI/CD

```yaml
# GitHub Actions example
- name: Install dependencies
  run: |
    pip install uv
    uv sync --frozen
```

### 3. Regular Updates

```bash
# Weekly/monthly
uv lock --upgrade
uv sync
git add uv.lock
git commit -m "Update dependencies"
```

### 4. Separate Optional Dependencies

```toml
[project.optional-dependencies]
dev = ["pytest", "black", "ruff"]
batch = ["harmonypy", "scvi-tools"]
```

### 5. Pin Critical Packages

```toml
dependencies = [
    "scanpy>=1.9.0,<2.0.0",  # Pin major version
    "numpy>=1.21.0,<2.0.0",  # Avoid numpy 2.0 issues
]
```

## Migration from Conda

### Export Conda Environment

```bash
# From conda environment
conda list --export > conda-requirements.txt

# Convert to pyproject.toml format
# (manual process - review each package)
```

### Hybrid Approach

```bash
# Use conda for system dependencies
conda install -c conda-forge nextflow

# Use UV for Python packages
uv sync
```

## Performance Comparison

Example: Installing scanpy + dependencies

| Tool | Time | Disk Usage |
|------|------|-----------|
| **uv** | 15 seconds | 500 MB |
| pip | 2 minutes | 500 MB |
| conda | 5 minutes | 2 GB |

## Additional Resources

- **UV Documentation**: https://github.com/astral-sh/uv
- **pyproject.toml Spec**: https://peps.python.org/pep-0621/
- **Python Packaging**: https://packaging.python.org/

## Quick Reference

```bash
# Essential commands
uv sync                      # Install from pyproject.toml + create lock
uv sync --frozen            # Install from lock file (exact versions)
uv add <package>            # Add package
uv remove <package>         # Remove package
uv lock --upgrade           # Update all packages
uv pip list                 # List installed packages
uv pip freeze              # Export installed packages

# With extras
uv sync --extra dev
uv sync --extra batch --extra trajectory
uv sync --all-extras       # Install all optional dependencies
```
