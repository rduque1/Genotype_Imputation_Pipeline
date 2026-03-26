#!/usr/bin/env bash

################################################################################
# Genotype Harmonizer QC1 Step
#
# Runs Genotype Harmonizer to align genotypes to a reference panel,
# applies call rate filtering, and fixes reference alleles.
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "${SCRIPT_DIR}/utils.sh"

################################################################################
# Usage
################################################################################

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Arguments:
  --input, -i PATH        Input PLINK file prefix (without .bed/.bim/.fam extension)
  --out, -o PATH          Output directory for harmonized files
  --ref-path PATH         Reference panel path (VCF files for Genotype Harmonizer)

Optional Arguments:
  --prefix, -p PREFIX     Output file prefix (default: derived from input)
  --lifted-code CODE      Lifted code suffix (default: GRCh37)
  --chr, -c RANGE         Chromosomes to process (default: 1-22 and X)
  --threads, -t NUM       Number of threads (default: 8)
  --fasta-ref PATH        Reference FASTA file for fixref (optional)
  --help                  Show this help message

Example:
  $0 \\
    --input /path/to/lifted/sample.GRCh37 \\
    --out /path/to/output \\
    --ref-path /path/to/reference/vcf

EOF
    exit 1
}

################################################################################
# Parse Arguments
################################################################################

INPUT_PREFIX=""
OUTROOT=""
REF_PATH=""
PREFIX=""
LIFTED_CODE="GRCh37"
CHR_RANGE=""
THREADS=8
FASTA_REF="/reference_prepared/fasta/human_g1k_v37.fasta.gz"

while [[ $# -gt 0 ]]; do
    case $1 in
        --input|-i)
            INPUT_PREFIX="$2"
            shift 2
            ;;
        --out|-o)
            OUTROOT="$2"
            shift 2
            ;;
        --ref-path)
            REF_PATH="$2"
            shift 2
            ;;
        --prefix|-p)
            PREFIX="$2"
            shift 2
            ;;
        --lifted-code)
            LIFTED_CODE="$2"
            shift 2
            ;;
        --chr|-c)
            CHR_RANGE="$2"
            shift 2
            ;;
        --threads|-t)
            THREADS="$2"
            shift 2
            ;;
        --fasta-ref)
            FASTA_REF="$2"
            shift 2
            ;;
        --help|-h)
            usage
            ;;
        *)
            print_error "Unknown option: $1"
            usage
            ;;
    esac
done

################################################################################
# Validate Arguments
################################################################################

if [[ -z "$OUTROOT" ]]; then
    print_error "Output directory is required (--out)"
    exit 1
fi

if [[ -z "$INPUT_PREFIX" ]]; then
    print_error "Input PLINK file prefix is required (--input)"
    exit 1
fi

if [[ -z "$REF_PATH" ]]; then
    print_error "Reference panel path is required (--ref-path)"
    exit 1
fi

# Derive prefix from input if not provided
if [[ -z "$PREFIX" ]]; then
    PREFIX=$(basename "$INPUT_PREFIX" | sed -e "s/\.${LIFTED_CODE}$//" -e 's/\.GRCh37$//' -e 's/\.lifted_.*$//')
fi

# Set up tool paths
export GH="${SCRIPT_DIR}/required_tools/GenotypeHarmonizer/GenotypeHarmonizer.jar"

################################################################################
# Helper Functions
################################################################################

# Detect which chromosomes are present in the input files
detect_chromosomes() {
    local input_prefix="$1"
    local detected=""

    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
        if [[ -f "${input_prefix}.chr${chr}.bed" ]]; then
            detected="$detected $chr"
        fi
    done

    echo "$detected" | sed 's/^ //' # Trim leading space
}

get_chr_range() {
    if [[ -n "$CHR_RANGE" ]]; then
        echo "$CHR_RANGE"
    else
        detect_chromosomes "$INPUT_PREFIX"
    fi
}

################################################################################
# Main: Genotype Harmonizer QC1
################################################################################

print_step "Genotype Harmonizer and QC1"

mkdir -p "$OUTROOT"
cd "$OUTROOT"

print_info "Input prefix: $INPUT_PREFIX"
print_info "Output directory: $OUTROOT"
print_info "Reference path: $REF_PATH"
print_info "Prefix: $PREFIX"
print_info "Lifted code: $LIFTED_CODE"

# Check FASTA reference
if [[ -n "$FASTA_REF" ]] && [[ ! -f "$FASTA_REF" ]]; then
    print_warning "Reference FASTA file not found at $FASTA_REF. Skipping fixref step."
    FASTA_REF=""
fi

DETECTED_CHRS=$(get_chr_range)
if [[ -z "$DETECTED_CHRS" ]]; then
    print_error "No chromosomes detected. Provide --chr or ensure expected BED files are present."
    exit 1
fi
print_info "Processing chromosomes: $DETECTED_CHRS"

################################################################################
# Function to run Genotype Harmonizer for one chromosome
################################################################################

run_gh() {
    local chr=$1
    print_info "Processing chromosome ${chr}..."

    if [[ "$chr" == "X" || "$chr" == "23" ]]; then
        REF="${REF_PATH}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
        java -Xmx16g -jar "$GH" \
            --keep \
            --input "${INPUT_PREFIX}.chr${chr}" \
            --ref "$REF" \
            --inputType PLINK_BED \
            --refType VCF \
            --update-reference-allele \
            --update-id \
            --callRateFilter 0.90 \
            --output "./${PREFIX}.${LIFTED_CODE}.GH.chr${chr}" \
            --debug
    else
        REF="${REF_PATH}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        java -Xmx16g -jar "$GH" \
            --keep \
            --input "${INPUT_PREFIX}.chr${chr}" \
            --ref "$REF" \
            --inputType PLINK_BED \
            --refType VCF \
            --update-reference-allele \
            --update-id \
            --callRateFilter 0.90 \
            --output "./${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
    fi
}

# Run Genotype Harmonizer for each chromosome
print_info "Running Genotype Harmonizer sequentially..."
for chr in $DETECTED_CHRS; do
    outname="${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
    if [[ -f "${outname}.vcf.gz" ]]; then
        print_info "Chromosome ${chr} already completed (${outname}.vcf.gz), skipping"
        continue
    fi
    if [[ -f "${outname}.bed" && -f "${outname}.bim" && -f "${outname}.fam" ]]; then
        print_info "Chromosome ${chr} already harmonized (${outname}.bed/.bim/.fam), skipping GH"
        continue
    fi
    print_info "Running Genotype Harmonizer for chromosome $chr..."
    run_gh "$chr"
done

################################################################################
# Convert BED to VCF and run fixref
################################################################################

print_info "Converting harmonized BED files to VCF and fixing reference alleles..."

fixref_and_convert() {
    local chr=$1
    local outname="${PREFIX}.${LIFTED_CODE}.GH.chr${chr}"
    local VCF_OUT_PREFIX="${outname}"

    if [[ -f "${outname}.vcf.gz" ]]; then
        print_info "Final VCF already exists for chromosome $chr, skipping: ${outname}.vcf.gz"
        return
    fi

    if [[ -f "${outname}.bed" ]]; then
        print_info "Converting and fixing chromosome $chr..."

        # --set-missing-var-ids @:#:\$1:\$2 : Use chr:pos for variant IDs to avoid issues with duplicates
        # --set-all-var-ids @:#:\$r:\$a : Use ref:alt for variant IDs to avoid issues with duplicates
        # --new-id-max-allele-len 500 : Allow long alleles for indels, default is 50 in plink2 which is too short for some variants in 1000G
        # --output-chr M : ensures chromosomes are named 1-22, X, Y, MT (no "chr").
        plink2 \
            --chr $chr \
            --output-chr M \
            --bfile "${outname}" \
            --export vcf-4.2 bgz \
            --set-all-var-ids @:#:\$r:\$a \
            --new-id-max-allele-len 500 \
            --out "${outname}.0"
        tabix -p vcf "${outname}.0.vcf.gz"

        if [[ -n "$FASTA_REF" ]]; then
            print_info "Running fixref for chr $chr"
            # Fix ref allele
            bcftools +fixref "$VCF_OUT_PREFIX.0.vcf.gz" -O z -o "$VCF_OUT_PREFIX.1.vcf.gz" -- -f "$FASTA_REF" -m ref-alt
            bcftools index -t "$VCF_OUT_PREFIX.1.vcf.gz"
            rm "$VCF_OUT_PREFIX.0.vcf.gz"*

            # Rename the fixed file to the final name
            mv "$VCF_OUT_PREFIX.1.vcf.gz" "$VCF_OUT_PREFIX.vcf.gz"
            mv "$VCF_OUT_PREFIX.1.vcf.gz.tbi" "$VCF_OUT_PREFIX.vcf.gz.tbi"
        else
            # If no reference fasta is provided, just rename the file
            mv "${outname}.0.vcf.gz" "${outname}.vcf.gz"
            mv "${outname}.0.vcf.gz.tbi" "${outname}.vcf.gz.tbi"
        fi
    else
        print_warning "BED file not found for chromosome $chr: ${outname}.bed"
    fi
}

# Run fixref for each chromosome
for chr in $DETECTED_CHRS; do
    fixref_and_convert "$chr"
done

################################################################################
# Complete
################################################################################

print_success "Genotype Harmonizer and fixref completed"
print_info "Output files in: $OUTROOT"
print_info "Output pattern: ${PREFIX}.${LIFTED_CODE}.GH.chr*.vcf.gz"
