#!/usr/bin/env bash

################################################################################
# Reference Panel Preparation Script
#
# This script prepares the 1000 Genomes Phase 3 reference panel for use
# with the imputation pipeline. It creates the necessary directory structure
# and converts VCF files to BCF and MSAV formats.
#
# Usage: bash prepare_reference_panel.sh --vcf-dir <path> --output-dir <path>
################################################################################

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

usage() {
    cat << EOF
Usage: bash prepare_reference_panel.sh --vcf-dir <path> --output-dir <path> [OPTIONS]

Required Arguments:
  --vcf-dir PATH          Directory containing downloaded 1000G VCF files
  --output-dir PATH       Output directory for prepared reference panel

Optional Arguments:
  --threads NUM           Number of threads to use (default: 4)
  --skip-bcf              Skip BCF conversion (for phasing)
  --skip-msav             Skip MSAV conversion (for imputation)
  --help                  Show this help message

Example:
  bash prepare_reference_panel.sh \\
    --vcf-dir /path/to/downloaded/vcfs \\
    --output-dir /path/to/reference_panel \\
    --threads 8

EOF
}

# Default values
VCF_DIR=""
OUTPUT_DIR=""
THREADS=4
SKIP_BCF=false
SKIP_MSAV=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf-dir)
            VCF_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --skip-bcf)
            SKIP_BCF=true
            shift
            ;;
        --skip-msav)
            SKIP_MSAV=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ -z "$VCF_DIR" ]]; then
    print_error "VCF directory is required (--vcf-dir)"
    exit 1
fi

if [[ ! -d "$VCF_DIR" ]]; then
    print_error "VCF directory does not exist: $VCF_DIR"
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    print_error "Output directory is required (--output-dir)"
    exit 1
fi

print_info "Starting reference panel preparation..."
print_info "VCF directory: $VCF_DIR"
print_info "Output directory: $OUTPUT_DIR"
print_info "Threads: $THREADS"

# Create directory structure
mkdir -p "$OUTPUT_DIR"/{vcf,bcf,msav,fasta}

# Download reference genome fasta file
if [ ! -f "$OUTPUT_DIR/fasta/human_g1k_v37.fasta.gz" ]; then
    print_info "Downloading reference genome fasta..."
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

    print_info "Fixing reference genome compression format..."
    gunzip human_g1k_v37.fasta.gz
    bgzip human_g1k_v37.fasta
    samtools faidx human_g1k_v37.fasta.gz
    mv human_g1k_v37.fasta.gz "$OUTPUT_DIR/fasta/"
    mv human_g1k_v37.fasta.fai "$OUTPUT_DIR/fasta/"
    print_success "Reference genome fasta downloaded and indexed"
fi

################################################################################
# Step 1: Copy/Link VCF files
################################################################################

print_info "Step 1: Setting up VCF files..."

for chr in {1..22} X; do
    VCF_FILE="$VCF_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

    # For chrX, the filename is different
    if [[ "$chr" == "X" ]]; then
        VCF_FILE="$VCF_DIR/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    fi

    if [[ ! -f "$VCF_FILE" ]]; then
        print_error "Missing VCF file: $VCF_FILE"
        exit 1
    fi

    # Create symlink to save space
    ln -sf "$(realpath "$VCF_FILE")" "$OUTPUT_DIR/vcf/"
    if [[ -f "${VCF_FILE}.tbi" ]]; then
        ln -sf "$(realpath "${VCF_FILE}.tbi")" "$OUTPUT_DIR/vcf/"
    fi
done

print_success "VCF files linked successfully"

################################################################################
# Step 2: Convert to BCF format (for phasing with Eagle)
################################################################################

if [[ "$SKIP_BCF" == false ]]; then
    print_info "Step 2: Converting VCF to BCF format for phasing..."

    for chr in {1..22} X; do
        print_info "Converting chromosome $chr to BCF..."

        INPUT_VCF="$OUTPUT_DIR/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        OUTPUT_BCF="$OUTPUT_DIR/bcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf"

        if [[ "$chr" == "X" ]]; then
            INPUT_VCF="$OUTPUT_DIR/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
            OUTPUT_BCF="$OUTPUT_DIR/bcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.bcf"
        fi

        if [[ ! -f "$OUTPUT_BCF" ]]; then
            print_info "Converting $INPUT_VCF to BCF..."
            bcftools view -Ob -o "$OUTPUT_BCF" "$INPUT_VCF" --threads "$THREADS"
            bcftools index "$OUTPUT_BCF"
        else
            print_info "BCF for chr$chr already exists, skipping..."
        fi
    done

    print_success "BCF conversion completed"
else
    print_info "Skipping BCF conversion (--skip-bcf)"
fi

################################################################################
# Step 3: Convert to msav format (for imputation with Minimac4)
################################################################################

if [[ "$SKIP_MSAV" == false ]]; then
    print_info "Step 3: Converting VCF to msav format for imputation..."
    print_info "This may take a while..."

    for chr in {1..22} X; do
        print_info "Converting chromosome $chr to msav..."
        INPUT_VCF="$OUTPUT_DIR/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        OUTPUT_MSAV="$OUTPUT_DIR/msav/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.msav"

        if [[ "$chr" == "X" ]]; then
            INPUT_VCF="$OUTPUT_DIR/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
            OUTPUT_MSAV="$OUTPUT_DIR/msav/ALL.chrX_PAR1.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav"
        fi

        if [[ ! -f "$OUTPUT_MSAV" ]]; then
            if [[ "$chr" == "X" ]]; then
                # PAR1
                bcftools view -r X:60001-2699520 $INPUT_VCF -Oz -o X_PAR1.vcf.gz
                tabix -p vcf X_PAR1.vcf.gz

                # nonPAR
                bcftools view -r X:2699521-154931043 $INPUT_VCF -Oz -o X_nonPAR.vcf.gz
                tabix -p vcf X_nonPAR.vcf.gz

                # PAR2
                bcftools view -r X:154931044-155260560 $INPUT_VCF -Oz -o X_PAR2.vcf.gz
                tabix -p vcf X_PAR2.vcf.gz

                echo "> Checking male ploidy in nonPAR (just to avoid cursed VCF errors)"
                bcftools +fixploidy X_nonPAR.vcf.gz -- --check || true

                # Optional auto-fix ploidy if needed
                # Uncomment this block if your data screams
                #
                # echo "> Auto-fixing ploidy issues in nonPAR"
                # bcftools +fixploidy X_nonPAR.vcf.gz -- -f > X_nonPAR.fixed.vcf
                # bgzip X_nonPAR.fixed.vcf
                # mv X_nonPAR.fixed.vcf.gz X_nonPAR.vcf.gz
                # tabix -p vcf X_nonPAR.vcf.gz

                minimac4 --compress-reference X_PAR1.vcf.gz > $OUTPUT_DIR/msav/ALL.chrX_PAR1.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav
                minimac4 --compress-reference X_nonPAR.vcf.gz > $OUTPUT_DIR/msav/ALL.chrX_nonPAR.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav
                minimac4 --compress-reference X_PAR2.vcf.gz > $OUTPUT_DIR/msav/ALL.chrX_PAR2.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.msav

                # Chromosome X has mixed ploidy issues with older minimac4 versions
                print_info "Skipping chromosome X msav conversion (requires newer minimac4 with --refSex support)"
                print_info "You can use a pre-built reference panel or upgrade minimac4"
            else
                minimac4 --compress-reference "$INPUT_VCF" \
                         --output "$OUTPUT_MSAV" \
                         --threads "$THREADS"
            fi
        else
            print_info "msav for chr$chr already exists, skipping..."
        fi
    done

    print_success "msav conversion completed"
else
    print_info "Skipping msav conversion (--skip-msav)"
fi

################################################################################
# Step 4: Create a configuration file for the pipeline
################################################################################

CONFIG_FILE="$OUTPUT_DIR/reference_panel.config"

cat > "$CONFIG_FILE" << EOF
# Reference Panel Configuration
# Generated on $(date)

# Base directory
REF_BASE=$OUTPUT_DIR

# VCF files (for Genotype Harmonizer - Step 2 and Ancestry Analysis - Step 3)
REF_VCF_DIR=$OUTPUT_DIR/vcf

# BCF files (for Phasing - Step 5)
REF_BCF_DIR=$OUTPUT_DIR/bcf

# MSAV files (for Imputation - Step 6)
REF_MSAV_DIR=$OUTPUT_DIR/msav

# Usage in imputation pipeline:
# Set these environment variables before running the pipeline:
#   export REF_PATH=$OUTPUT_DIR/vcf
#   export REF_BCF_BASE=$OUTPUT_DIR/bcf
#   export REF_MSAV_BASE=$OUTPUT_DIR/msav
EOF

print_success "Configuration file created: $CONFIG_FILE"

################################################################################
# Summary
################################################################################

echo
echo "================================================================"
echo "Reference Panel Preparation Complete!"
echo "================================================================"
echo
print_info "Directory structure:"
print_info "  $OUTPUT_DIR/vcf/     - VCF files (for steps 2-3)"
print_info "  $OUTPUT_DIR/bcf/     - BCF files (for step 5)"
print_info "  $OUTPUT_DIR/msav/    - MSAV files (for step 6)"
echo
print_info "Configuration saved to: $CONFIG_FILE"
echo
print_info "To use with the imputation pipeline, update the script or set:"
echo "  export REF_PATH='$OUTPUT_DIR/vcf'"
echo "  export REF_BCF_BASE='$OUTPUT_DIR/bcf'"
echo "  export REF_MSAV_BASE='$OUTPUT_DIR/msav'"
echo
