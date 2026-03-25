# Genotype Imputation Pipeline

Automated pipeline for genotype array imputation using 1000 Genomes or HRC reference panels. Includes genome build detection, liftover, QC, ancestry analysis, phasing, and imputation.

## Features

- ✅ **Multiple input formats**: VCF/VCF.gz and **23andMe .txt/.txt.gz/.zip** (auto-detected)
- ✅ **Automated genome build detection** (hg16-hg38, GRCh37-38, NCBI builds)
- ✅ **Chromosome support**: Autosomes (1-22) + **X chromosome**
- ✅ **Dual reference panels**: 1000 Genomes Phase 3 (81M variants) or HRC (40M variants)
- ✅ **Ancestry-specific processing**: AFR, AMR, EAS, EUR, SAS, and admixed populations
- ✅ **Docker containerization** with all dependencies included
- ✅ **Single chromosome mode** for faster testing/debugging
- ✅ **Partial run support**: Start/stop at any pipeline step (0-6)

## Quick Start (Docker - Recommended)

### 1. Build Image
```bash
docker build -t imputation-pipeline .
```

### 2. Prepare Reference Panel
```bash
# Download 1000G Phase 3 VCF files (autosomes + X)
for chr in {1..22}; do
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz

# Prepare reference (converts to BCF and M3VCF formats)
docker run --rm --privileged \
  -v $(pwd)/vcf:/vcf_input \
  -v $(pwd)/reference:/reference_output \
  --memory=64g --cpus=16 \
  imputation-pipeline \
  ./prepare_reference_panel.sh \
    --vcf-dir /vcf_input \
    --output-dir /reference_output \
    --threads 16
```

**Note**: `--privileged` is required for bcftools/minimac4 multi-threading due to Docker security restrictions.

### 3. Run Pipeline
```bash
# Process all chromosomes (1-22 and X if present)
docker run --rm --privileged \
  -v $(pwd)/data:/data \
  -v $(pwd)/reference:/reference \
  -v $(pwd)/output:/output \
  --memory=64g --cpus=16 \
  imputation-pipeline \
  ./imputation_pipeline.sh \
    --vcf /data/input.vcf.gz \
    --out /output \
    --ref 1KG \
    --ref-path /reference/vcf \
    --ref-bcf /reference/bcf \
    --ref-m3vcf /reference/m3vcf \
    --threads 16

# Or process single chromosome for testing
docker run --rm --privileged \
  -v $(pwd)/data:/data \
  -v $(pwd)/reference:/reference \
  -v $(pwd)/output:/output \
  imputation-pipeline \
  ./imputation_pipeline.sh \
    --vcf /data/input.vcf.gz \
    --out /output \
    --ref 1KG \
    --chr 22 \
    --threads 8
```

## Pipeline Overview

| Step | Description | Output |
|------|-------------|--------|
| 0 | Detect genome build | Build report |
| 1 | Lift to GRCh37 | Lifted genotypes |
| 2 | QC & harmonization | Harmonized VCFs |
| 3 | Ancestry analysis | Split by ancestry (AFR/AMR/EAS/EUR/SAS/mixed) |
| 4 | Missingness/HWE filter | QC'd genotypes |
| 5 | Phasing (Eagle) | Phased haplotypes |
| 6 | Imputation (Minimac4) | **Imputed genotypes (~40M-81M variants)** |

Final outputs: `output/6_impute_{HRC|1KG}/{prefix}/{ancestry}/imputed_{1-22,X}.dose.vcf.gz`

## Options

```bash
Required:
  --vcf PATH       Input file (VCF/VCF.gz or 23andMe .txt/.txt.gz/.zip)
  --out PATH       Output directory

Optional:
  --ref TYPE       Reference panel: HRC or 1KG (default: 1KG)
  --threads N      Number of threads (default: 8)
  --chr N          Process single chromosome (1-22 or X, default: all)
  --ref-path PATH  Custom reference VCF directory
  --ref-bcf PATH   Custom reference BCF directory
  --ref-m3vcf PATH Custom reference M3VCF directory
  --start N        First step (0-6, default: 0)
  --end N          Last step (0-6, default: 6)
  --hwe            Enable HWE filtering (for non-admixed populations)
  --wgs            Enable WGS mode (down-sample variants in ancestry analysis)

Chromosome Processing:
  Pipeline automatically detects and processes chromosomes 1-22 and X (if present).
  Use --chr to process a single chromosome: --chr 22 or --chr X
```

## Examples

```bash
# Full pipeline with default 1KG reference
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --threads 16

# Process 23andMe file (automatically converted to VCF)
./imputation_pipeline.sh --vcf genome_data.txt --out results --threads 16

# Use HRC reference panel instead
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --ref HRC --threads 16

# Test with single chromosome
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --chr 22

# Process only X chromosome
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --ref 1KG --chr X

# Run specific steps (e.g., only imputation after phasing)
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --ref 1KG --start 6 --end 6

# With HWE filtering for homogeneous populations
./imputation_pipeline.sh --vcf sample.vcf.gz --out results --ref 1KG --hwe --threads 16
```

## Requirements

- **RAM**: 64GB+ recommended
- **CPU**: 16+ cores recommended
- **Disk**: ~200GB (50GB reference + 150GB data/output)
- **Docker**: 20.10+

## Output Directory Structure

```
output/
├── 0_convert_23andme/                                     # (if 23andMe input)
│   └── {prefix}.vcf.gz                                    # Converted VCF
├── 0_check_vcf_build/
│   └── {prefix}.BuildChecked                              # Build detection report
├── 1_lift/
│   └── {prefix}.lifted_{build}_to_GRCh37.chr{1-22,X}.bed  # Lifted genotypes
├── 2_GH/
│   └── {prefix}.{lifted}.GH.chr{1-22,X}.vcf.gz            # Harmonized VCFs
├── 3_ancestry/{prefix}/
│   ├── {prefix}.{lifted}.GH.pruned.intersect1KG.5.Q       # Ancestry proportions
│   └── {prefix}.{lifted}.GH.ancestry-{1-5,mixed}.chr*.bed # Split by ancestry
├── 4_split_QC2/{prefix}/
│   └── {prefix}.{lifted}.GH.ancestry-{1-5,mixed}.chr*.bed # QC filtered
├── 5_phase/{prefix}/
│   └── {prefix}.{lifted}.GH.ancestry-*.chr*.phased.vcf.gz # Phased haplotypes
└── 6_impute_{HRC|1KG}/{prefix}/{1-5,mixed}/
    ├── imputed_{1-22,X}.dose.vcf.gz                       # ⭐ Final imputed genotypes
    ├── imputed_{1-22,X}.dose.vcf.gz.tbi                   # Index files
    └── imputed_{1-22,X}.info                              # Quality metrics (R²)
```

**Key outputs**: `6_impute_*/` contains ~40M variants (HRC) or ~81M variants (1KG)

**Ancestry codes**: 1=AFR, 2=AMR, 3=EAS, 4=EUR, 5=SAS, mixed=Admixed

## Troubleshooting

- **Threading errors in Docker**: Use `--privileged` flag (required for bcftools/minimac4)
- **Out of memory**: Increase `--memory=128g` or reduce `--threads`
- **Reference not found**: Prepare reference panel first (see Quick Start step 2)
- **X chromosome missing**: Ensure reference panel includes chrX files (see step 2)

## Examples

```bash
 bash imputation_pipeline.sh --vcf id_600.input_data.zip --out /output_dir --threads 1 --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-m3vcf /reference_prepared/m3vcf --chr 2
```

```bash
docker run --privileged -it --rm -v $PWD:/output_dir -v $PWD:/input_data/ -v $PWD/../1000GRef:/vcf_ref -v $PWD/../1000G_prepared:/reference_prepared imputation-pipeline:latest bash
```

### Run directly from the command line with Docker

```bash
export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared --memory=32g --cpus=16 imputation-pipeline:latest bash imputation_pipeline.sh --vcf /input_data/id_600.input_data.zip --out /output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-m3vcf /reference_prepared/m3vcf --threads 16
```

Single chromosome:

```bash
export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared  imputation-pipeline:latest bash imputation_pipeline.sh --vcf /input_data/id_600.input_data.zip --out /output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-m3vcf /reference_prepared/m3vcf --chr 2 --threads 8
```

Detached mode with log file:

```bash
export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run -d --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared --memory=32g --cpus=16 imputation-pipeline:latest bash -c 'bash imputation_pipeline.sh --vcf /input_data/id_600.input_data.zip --out output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-msav /reference_prepared/msav --threads 16 2>&1 | tee /output_dir/pipeline.log'
```

#### imputation_pipeline_no_ancestry.sh

```bash
export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run -d --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared --memory=32g --cpus=16 imputation-pipeline:latest bash -c 'bash imputation_pipeline_no_ancestry.sh --vcf /input_data/id_600.input_data.zip --out output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-msav /reference_prepared/msav --threads 16 2>&1 | tee /output_dir/pipeline.log'

export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared --memory=32g --cpus=16 imputation-pipeline:latest bash -c 'bash imputation_pipeline_no_ancestry.sh --vcf /input_data/id_600.input_data.zip --out output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-msav /reference_prepared/msav --threads 16 2>&1 | tee /output_dir/pipeline.log'
```

#### preprocess.sh

```bash
export OUTDIR=output_$(date +%Y%m%d_%H%M%S) && mkdir -p $OUTDIR && docker run --rm --privileged -v $PWD/new_imp/input_test_data:/input_data -v $PWD/$OUTDIR:/output_dir -v $PWD/new_imp/1000GRef:/vcf_ref -v $PWD/new_imp/1000G_prepared:/reference_prepared --memory=32g --cpus=16 imputation-pipeline:latest bash -c 'bash preprocess.sh --vcf /input_data/id_600.input_data.zip --out /output_dir --ref-path /vcf_ref --ref-bcf /reference_prepared/bcf --ref-msav /reference_prepared/msav 2>&1 | tee /output_dir/pipeline.log'
```