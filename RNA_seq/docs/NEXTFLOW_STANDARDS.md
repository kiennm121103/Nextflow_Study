# Standard Nextflow Pipeline File Organization

## Overview

This guide shows the standard file organization for professional Nextflow pipelines, following nf-core best practices.

## Standard Directory Structure

```
my-pipeline/                      # Root directory
│
├── .github/                      # GitHub-specific files
│   ├── workflows/               # CI/CD workflows
│   │   ├── ci.yml              # Continuous integration
│   │   └── linting.yml         # Code linting
│   ├── ISSUE_TEMPLATE/         # Issue templates
│   └── PULL_REQUEST_TEMPLATE.md
│
├── assets/                       # Static assets
│   ├── multiqc_config.yaml     # MultiQC configuration
│   ├── schema_input.json       # Input validation schema
│   └── email_template.html     # Email notification template
│
├── bin/                          # Executable scripts
│   ├── scrape_software_versions.py
│   ├── markdown_to_html.py
│   └── check_samplesheet.py
│
├── conf/                         # Configuration files
│   ├── base.config             # Base configuration
│   ├── test.config             # Test configuration
│   ├── test_full.config        # Full test configuration
│   └── modules.config          # Module-specific configs
│
├── docs/                         # Documentation
│   ├── output.md               # Output description
│   ├── usage.md                # Usage guide
│   └── images/                 # Documentation images
│
├── lib/                          # Groovy library scripts
│   ├── nfcore_external_java_deps.jar
│   ├── NfcoreSchema.groovy
│   └── WorkflowMain.groovy
│
├── modules/                      # Process modules
│   ├── local/                  # Local modules
│   │   ├── custom_process.nf
│   │   └── another_process.nf
│   └── nf-core/                # nf-core modules
│       ├── fastqc/
│       │   ├── main.nf
│       │   └── meta.yml
│       └── multiqc/
│           ├── main.nf
│           └── meta.yml
│
├── subworkflows/                 # Reusable subworkflows
│   ├── local/                  # Local subworkflows
│   │   └── input_check.nf
│   └── nf-core/               # nf-core subworkflows
│       └── fastq_qc/
│
├── workflows/                    # Main workflow files
│   └── pipeline.nf             # Main pipeline workflow
│
├── .editorconfig                # Editor configuration
├── .gitattributes              # Git attributes
├── .gitignore                  # Git ignore rules
├── .nf-core.yml                # nf-core configuration
├── .prettierignore             # Prettier ignore
├── .prettierrc.yml             # Prettier configuration
│
├── CHANGELOG.md                 # Version history
├── CITATIONS.md                # Citations and references
├── CODE_OF_CONDUCT.md          # Code of conduct
├── LICENSE                      # License file
├── README.md                    # Main documentation
│
├── main.nf                      # Main pipeline entry point
├── nextflow.config             # Main configuration
├── nextflow_schema.json        # Parameter schema
│
└── pyproject.toml              # Python dependencies (if needed)
```

## Detailed Explanation

### 1. Root Level Files

#### `main.nf` (Required)
The main pipeline entry point.

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    Pipeline Name
========================================================================================
    Github : https://github.com/yourname/pipeline
    Author : Your Name
----------------------------------------------------------------------------------------
*/

// Print help message
def helpMessage() {
    log.info"""
    Usage:
      nextflow run main.nf --input samplesheet.csv --genome GRCh38

    Required arguments:
      --input               Path to input samplesheet
      --genome              Reference genome
    
    Optional arguments:
      --outdir              Output directory [default: ./results]
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Import workflows
include { PIPELINE } from './workflows/pipeline'

// Main workflow
workflow {
    PIPELINE ()
}
```

#### `nextflow.config` (Required)
Main configuration file.

```groovy
// Pipeline parameters
params {
    // Input/output options
    input                      = null
    outdir                     = './results'
    
    // Reference genome options
    genome                     = null
    fasta                      = null
    
    // Other options
    help                       = false
    publish_dir_mode           = 'copy'
    
    // Max resource options
    max_cpus                   = 16
    max_memory                 = '128.GB'
    max_time                   = '240.h'
}

// Load base config
includeConfig 'conf/base.config'

// Profiles
profiles {
    debug {
        cleanup = false
    }
    
    conda {
        conda.enabled          = true
        docker.enabled         = false
        singularity.enabled    = false
    }
    
    docker {
        docker.enabled         = true
        docker.userEmulation   = true
        conda.enabled          = false
        singularity.enabled    = false
    }
    
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        conda.enabled          = false
        docker.enabled         = false
    }
    
    test {
        includeConfig 'conf/test.config'
    }
}

// Export these variables to prevent local Python/R libraries
env {
    PYTHONNOUSERSITE = 1
    R_PROFILE_USER   = "/.Rprofile"
    R_ENVIRON_USER   = "/.Renviron"
}

// Capture exit codes
process.errorStrategy = 'finish'

// Function to ensure that resource requirements don't go beyond maximum
def check_max(obj, type) {
    if (type == 'memory') {
        try {
            if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                return params.max_memory as nextflow.util.MemoryUnit
            else
                return obj
        } catch (all) {
            println "ERROR - Max memory '${params.max_memory}' not valid"
            return obj
        }
    } else if (type == 'time') {
        try {
            if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
            else
                return obj
        } catch (all) {
            println "ERROR - Max time '${params.max_time}' not valid"
            return obj
        }
    } else if (type == 'cpus') {
        try {
            return Math.min(obj, params.max_cpus as int)
        } catch (all) {
            println "ERROR - Max cpus '${params.max_cpus}' not valid"
            return obj
        }
    }
}
```

### 2. Configuration Directory (`conf/`)

#### `conf/base.config`
Base configuration for all profiles.

```groovy
process {
    cpus   = { check_max( 1    * task.attempt, 'cpus'   ) }
    memory = { check_max( 6.GB * task.attempt, 'memory' ) }
    time   = { check_max( 4.h  * task.attempt, 'time'   ) }

    errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries    = 1
    maxErrors     = '-1'

    // Process-specific resource requirements
    withLabel:process_low {
        cpus   = { check_max( 2     * task.attempt, 'cpus'    ) }
        memory = { check_max( 12.GB * task.attempt, 'memory'  ) }
        time   = { check_max( 4.h   * task.attempt, 'time'    ) }
    }
    withLabel:process_medium {
        cpus   = { check_max( 6     * task.attempt, 'cpus'    ) }
        memory = { check_max( 36.GB * task.attempt, 'memory'  ) }
        time   = { check_max( 8.h   * task.attempt, 'time'    ) }
    }
    withLabel:process_high {
        cpus   = { check_max( 12    * task.attempt, 'cpus'    ) }
        memory = { check_max( 72.GB * task.attempt, 'memory'  ) }
        time   = { check_max( 16.h  * task.attempt, 'time'    ) }
    }
    withLabel:process_long {
        time   = { check_max( 20.h  * task.attempt, 'time'    ) }
    }
}
```

#### `conf/test.config`
Test configuration with minimal data.

```groovy
params {
    config_profile_name        = 'Test profile'
    config_profile_description = 'Minimal test dataset'

    // Limit resources
    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data
    input  = 'https://raw.githubusercontent.com/yourname/testdata/samplesheet.csv'
    genome = 'GRCh38'
}
```

### 3. Modules Directory (`modules/`)

#### `modules/local/custom_process.nf`
Example custom module.

```groovy
process CUSTOM_PROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::tool=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tool:1.0' :
        'quay.io/biocontainers/tool:1.0' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool \\
        $args \\
        -i $reads \\
        -r $reference \\
        -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1 | sed 's/^.*tool //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool: \$(tool --version 2>&1 | sed 's/^.*tool //; s/ .*\$//')
    END_VERSIONS
    """
}
```

#### `modules/local/meta.yml`
Module metadata.

```yaml
name: custom_process
description: Process description
keywords:
  - alignment
  - mapping
  - genomics
tools:
  - tool:
      description: Tool description
      homepage: https://tool.homepage.com
      documentation: https://tool.docs.com
      tool_dev_url: https://github.com/tool/tool
      doi: "10.1000/doi"
      licence: ['MIT']

input:
  - meta:
      type: map
      description: Groovy Map containing sample information
  - reads:
      type: file
      description: Input FASTQ files
      pattern: "*.{fastq,fq}.gz"
  - reference:
      type: file
      description: Reference genome
      pattern: "*.{fa,fasta}"

output:
  - meta:
      type: map
      description: Groovy Map containing sample information
  - bam:
      type: file
      description: Output BAM file
      pattern: "*.bam"
  - versions:
      type: file
      description: File containing software versions
      pattern: "versions.yml"

authors:
  - "@yourname"
```

### 4. Workflows Directory (`workflows/`)

#### `workflows/pipeline.nf`
Main workflow.

```groovy
/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

// Local modules
include { CUSTOM_PROCESS } from '../modules/local/custom_process'

// nf-core modules
include { FASTQC } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'

// Subworkflows
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow PIPELINE {
    
    // Create channels
    ch_versions = Channel.empty()
    
    // Input check
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    // FastQC
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    
    // Custom process
    CUSTOM_PROCESS (
        INPUT_CHECK.out.reads,
        params.reference
    )
    ch_versions = ch_versions.mix(CUSTOM_PROCESS.out.versions)
    
    // MultiQC
    MULTIQC (
        FASTQC.out.zip.collect(),
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
```

### 5. Subworkflows Directory (`subworkflows/`)

#### `subworkflows/local/input_check.nf`

```groovy
include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
```

### 6. Binary Scripts (`bin/`)

#### `bin/check_samplesheet.py`

```python
#!/usr/bin/env python3

import sys
import csv
from pathlib import Path

def check_samplesheet(file_in, file_out):
    """
    Check that the samplesheet follows the format:
    sample,fastq_1,fastq_2,single_end
    """
    
    required_columns = {"sample", "fastq_1", "fastq_2", "single_end"}
    
    # Check header
    with open(file_in, "r") as fin:
        reader = csv.DictReader(fin)
        
        # Validate header
        if not required_columns.issubset(reader.fieldnames):
            print(f"ERROR: Missing required columns: {required_columns - set(reader.fieldnames)}")
            sys.exit(1)
        
        # Check rows
        sample_dict = {}
        for i, row in enumerate(reader):
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            single_end = row["single_end"]
            
            # Check sample name
            if not sample:
                print(f"ERROR: Line {i+2}: Sample name is required")
                sys.exit(1)
            
            # Check files exist
            if not Path(fastq_1).exists():
                print(f"ERROR: Line {i+2}: FASTQ file does not exist: {fastq_1}")
                sys.exit(1)
            
            if single_end == "false" and not Path(fastq_2).exists():
                print(f"ERROR: Line {i+2}: FASTQ file does not exist: {fastq_2}")
                sys.exit(1)
            
            # Store samples
            if sample in sample_dict:
                print(f"ERROR: Duplicate sample name: {sample}")
                sys.exit(1)
            sample_dict[sample] = row
    
    # Write validated samplesheet
    with open(file_in, "r") as fin, open(file_out, "w") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            writer.writerow(row)

if __name__ == "__main__":
    check_samplesheet(sys.argv[1], sys.argv[2])
```

### 7. Documentation (`docs/`)

#### `docs/usage.md`
Usage documentation.

```markdown
# Pipeline Usage

## Introduction

This pipeline processes...

## Running the pipeline

### Quick start

```bash
nextflow run main.nf \\
  --input samplesheet.csv \\
  --genome GRCh38 \\
  --outdir results \\
  -profile docker
```

### Samplesheet format

```csv
sample,fastq_1,fastq_2,single_end
SAMPLE1,sample1_R1.fastq.gz,sample1_R2.fastq.gz,false
SAMPLE2,sample2_R1.fastq.gz,sample2_R2.fastq.gz,false
```
```

#### `docs/output.md`
Output documentation.

```markdown
# Pipeline Output

## Introduction

This document describes the output...

## Pipeline results

### FastQC

- `fastqc/`: FastQC HTML reports
- `fastqc/*_fastqc.html`: FastQC report
- `fastqc/*_fastqc.zip`: FastQC data

### Alignment

- `alignment/`: BAM files
- `alignment/*.bam`: Aligned reads
```

### 8. Assets (`assets/`)

#### `assets/multiqc_config.yaml`

```yaml
report_comment: >
    This report has been generated by the pipeline
report_header_info:
    - Application: 'Pipeline Name'
    - Version: '1.0.0'

module_order:
    - fastqc:
        name: 'FastQC (raw)'
        path_filters:
            - '*_fastqc.zip'
```

### 9. GitHub Workflows (`.github/workflows/`)

#### `.github/workflows/ci.yml`

```yaml
name: CI

on:
  push:
    branches: [ main, dev ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Install Nextflow
        run: |
          wget -qO- https://get.nextflow.io | bash
          sudo mv nextflow /usr/local/bin/
      
      - name: Run pipeline with test data
        run: |
          nextflow run main.nf -profile test,docker
```

### 10. Additional Files

#### `.gitignore`

```
work/
results/
.nextflow/
.nextflow.log*
*.pyc
__pycache__/
.DS_Store
```

#### `CHANGELOG.md`

```markdown
# Changelog

## [1.0.0] - 2025-01-01

### Added
- Initial release
- FastQC module
- Alignment module

### Changed
- Updated configuration

### Fixed
- Bug fixes
```

## Best Practices

### 1. Naming Conventions

- **Processes**: UPPERCASE (e.g., `FASTQC`, `STAR_ALIGN`)
- **Workflows**: PascalCase (e.g., `InputCheck`, `AlignReads`)
- **Parameters**: snake_case (e.g., `input_file`, `max_memory`)
- **Channels**: ch_ prefix (e.g., `ch_reads`, `ch_versions`)

### 2. File Organization

- Keep modules atomic (one tool = one module)
- Use subworkflows for repeated patterns
- Separate configuration from code
- Document everything

### 3. Version Control

- Track software versions
- Use semantic versioning (MAJOR.MINOR.PATCH)
- Maintain CHANGELOG.md
- Tag releases

### 4. Testing

- Provide test profile with small dataset
- Include stub runs for faster testing
- Test on different executors

### 5. Documentation

- README: Overview and quick start
- docs/usage.md: Detailed usage
- docs/output.md: Output description
- Module meta.yml: Tool documentation

## Resources

- [nf-core tools](https://nf-co.re/tools)
- [Nextflow patterns](https://nextflow-io.github.io/patterns/)
- [nf-core template](https://github.com/nf-core/tools)
