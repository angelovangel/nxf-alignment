# nxf-alignment

A Nextflow workflow for basecalling (ONT only), aligning, and variant calling for long-read sequencing data (ONT + HiFi).

## Features

- **Basecalling**: Uses Dorado for basecalling (and demultiplexing) with optional adaptive sampling support (ONT)
- **Alignment**: Aligns reads to a reference genome using Dorado aligner (modifications are preserved, ONT or HiFi data)
- **Coverage Analysis**: Calculates per-region coverage statistics with thresholds (1x, 10x, 20x, 30x)
- **Variant Calling**: Uses Clair3 for variant calling (ONT or HiFi data)
- **Variant Annotation**: Uses snpEff for variant annotation (ONT or HiFi data)
- **Interactive HTML Report**: Generates an interactive report with read statistics, coverage, variants and annotations metrics

## Requirements

- **Nextflow** >= 23.04
- **Docker** 
- **NVIDIA GPU** (for basecalling and variants)

## Quick Start

See also [workflow diagram](https://angelovangel.github.io/nxf-alignment/assets/diagram.html) for how parameter use determine the workflow path

#### Basic Workflow (Basecalling + Alignment + Variant Calling)
For an adaptive sampling run, basecalling is done for the accepted reads based on the decision file produced by MinKNOW.

```bash
nextflow run angelovangel/nxf-alignment \
  --pod5 /path/to/pod5/dir \
  --asfile /path/to/AS_decisions.csv # optional
  --model hac \
  --bed /path/to/regions.bed # optional, if provided the report contains coverage analysis per region from bed file 
  --ref /path/to/ref.fasta \
  --variants \
  --annotate
```
#### Barcoded run (Basecalling + Alignment)
For a barcoded run, provide a [samplesheet](#sample-sheet-barcoded-runs) and kit name
```bash
nextflow run angelovangel/nxf-alignment \
  --pod5 /path/to/pod5/dir \
  --model hac,5mC_5hmC \
  --bed /path/to/regions.bed
  --ref /path/to/ref.fasta
  --kit SQK-RBK114-96
  --samplesheet /path/to/samplesheet.csv
```

#### Skip Basecalling (Align Existing BAM/FASTQ + Variant Calling)
If the basecalling has been performed before, the pipeline can be run with the `--reads` parameter. The reads can be in any HTS format, a directory of reads can also be given.

```bash
nextflow run angelovangel/nxf-alignment \
  --reads /path/to/reads.bam \ # can be also a directory with reads
  --ref /path/to/ref.fasta \
  --bedfile /path/to/regions.bed \
  --variants \
  --annotate
```
#### Skip alignment (basecalling only)
Basecalling (for single sample and barcoded runs) can also be performed without alignment, using the `-entry` parameter.
`-entry basecall` will do basecalling (and evt demultiplexing), `-entry report` will do basecalling + report
```bash
nextflow run angelovangel/nxf-alignment \
    --pod5 /path/to/pod5/dir \
    --model hac \
    --kit SQK-RBK114-96 # for barcoded runs only
    --samplesheet /path/to/samplesheet.csv # for barcoded runs only
    -entry report
```


## Parameters

#### Core Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pod5` | path | - | Directory containing POD5 files (required if not using `--reads`) |
| `reads` | path | null | Path to input BAM/FASTQ file(s) or directory (skips basecalling) |
| `ref` | path | - | Reference genome in FASTA format (required) |
| `model` | string | `fast` | Dorado basecall model, see [available models](https://software-docs.nanoporetech.com/dorado/latest/models/list/)|
| `outdir` | string | `results` | Output directory for results |

#### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `asfile` | path | null | Adaptive sampling decisions CSV (if using AS filtering) |
| `kit` | string | null | Barcoding kit name (e.g., `SQK-NBD111-96`). Required for barcoded runs |
| `samplesheet` | path | null | Sample sheet CSV with columns: `sample`, `barcode`. Required for barcoded runs |
| `bed` | path | null | BED file with target regions (auto-generated from reference if not provided) |
| `variants` | boolean | false | Perform variant calling using Clair3 or DeepVariant |
| `variant_caller` | string | `clair3` | Variant caller to use (`clair3` or `deepvariant`, use only with `--variants`) |
| `clair3_platform` | string | `ont` | Platform to use for Clair3 (`ont` or `hifi`) |
| `clair3_model` | string | `r1041_e82_400bps_hac_v500` | Model to use for Clair3 |
| `deepvariant_model` | string | `ONT_R104` | Model to use for DeepVariant |
| `annotate` | boolean | false | Annotate variants using snpEff (use only with `--variants`) |
| `anno_db` | string | `hg38` | Database to use for annotation |

#### Profiles
Predefined set of parameters for common use cases, use with `-profile`:

| Profile | Description |
|---------|-------------|
| `standard` | Standard workflow with Docker GPU support |
| `singularity` | Workflow with Singularity GPU support |
| `revio` | Workflow optimized for HiFi Revio (use with `--reads` to skip basecalling)|

## Output Structure

```
results/
├── 00-basecall/
│   ├── reads.bam              # Basecalled reads
│   └── processed/             # Per-sample BAMs (if barcoded)
├── 01-align/
│   ├── reads.align.bam        # Aligned reads
│   └── reads.align.bam.bai    # BAM index
├── 02-coverage/
│   └── reads.hist.tsv         # Coverage histogram
├── 03-variants/
│   └── reads.variants.vcf     # Variants
└── nxf-alignment-report.html  # Workflow report
```

## Input Files


#### BED File (Optional)
Tab-separated file defining target regions:
```
chr1    1000    5000    GENE_A
chr1    8000    12000   GENE_B
```
If not provided, the workflow auto-generates a BED file covering the entire reference.

#### Sample Sheet (Barcoded Runs)
CSV with minimum columns `sample` and `barcode`:
```
sample,barcode
sample_1,barcode01
sample_2,barcode02
```

#### Adaptive Sampling Decisions File (Optional)
This file is generated by MinKNOW during an adaptive sampling run, and can be found under `runfolder/adaptive_sampling/AS_decisions.csv`

## Enabling docker GPU Support
If you observe the error "could not select device driver with capabilities gpu", additional docker setup is required. The `nvidia-container-toolkit` has to be installed and running on your system. See [here](https://epi2me.nanoporetech.com/epi2me-docs/installation/) and [here](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) for details. 


## Citation

If you use this workflow, please cite:
- **Dorado**: https://github.com/nanoporetech/dorado
- **Bedtools**: Quinlan & Hall, 2010
- **Nextflow**: Di Tommaso et al., 2017
